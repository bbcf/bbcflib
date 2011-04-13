"""
===============
bbcflib.chipseq
===============

This module provides functions to run a basic ChIP-seq analysis from reads mapped on
a reference genome.
The most important steps are the binding of ``macs`` via the 'add_macs_results' function, 
the peak deconvolution algorithm with the 'run_deconv' function and
the 'parallel_density_sql' function to create quantitative sql tracks from bam files. Below is the script used by the frontend::
    hts_key = 'test_key'
    M = MiniLIMS( '/path/to/chipseq/minilims' )
    gl = use_pickle( M, "global variables" )
    htss = frontend.Frontend( url=gl["hts_url"] )
    job = htss.job( hts_key )
    g_rep = genrep.GenRep( gl["genrep_url"], gl["bwt_root"] )
    with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
        M_ms = MiniLIMS( '/path/to/mapseq/minilims' )
        gl_ms = use_pickle( M_ms, "global variables" )
        mapseq_url = gl_ms['hts_url']
        (ms_files, ms_job) = import_mapseq_results( job.options['mapseq_key'], M_ms, ex.working_directory, gl_ms["hts_url"] )
        g_rep_assembly = g_rep.assembly( ms_job.assembly_id )
        job.groups = ms_job.groups
        files = workflow_groups( ex, job, ms_files, g_rep_assembly.chromosomes, gl['script_path'] )

"""

import sqlite3
from bein import *
from bein.util import *
from bbcflib import frontend, mapseq
from bbcflib.common import *

############ Peaks and annotation ############
@program
def macs( read_length, genome_size, bamfile, ctrlbam=None, shift=80, args=[] ):
    """Binding for the ``macs`` peak caller.

    takes one (optionally two) bam file(s) and
    the 'read_length', 'genome_size' and 'shift' parameters passed to ``macs``. 
    
    Returns the file prefix ('-n' option of ``macs``)
    """
    outname = unique_filename_in()
    macs_args = ["macs14","-t",bamfile]
    if ctrlbam != None:
        macs_args += ["-c",ctrlbam]
    macs_args += ["-n",outname,"-f","BAM",
                  "-g",str(genome_size),"-s",str(read_length)]
    if shift>0:
        macs_args += ["--shiftsize="+str(shift)]
    return {"arguments": macs_args+args, "return_value": outname}

def add_macs_results( ex, read_length, genome_size, bamfile,
                      ctrlbam=[None], name=None, shift=80, alias=None,
                      macs_args=[], via='lsf' ):
    """Calls the ``macs`` function on each possible pair 
    of test and control bam files and adds 
    the respective outputs to the execution repository.
    Returns the set of file prefixes.
    """
    if not(isinstance(bamfile,list)):
        bamfile = [bamfile]
    if not(isinstance(ctrlbam,list)):
        ctrlbam = [ctrlbam]
    futures = {}
    rl = read_length
    gs = genome_size
    for i,bam in enumerate(bamfile):
        n = name['tests'][i]
        if isinstance(read_length,list):
            rl = read_length[i]
        if isinstance(genome_size,list):
            gs = read_length[i]
        for j,cam in enumerate(ctrlbam):
            m = name['controls'][j]
            if m == None:
                nm = (n,)
            else:
                nm = (n,m)
            futures[nm] = macs.nonblocking( ex, rl, gs, bam, cam, shift, 
                                            args=macs_args, via=via )
    prefixes = dict((n,f.wait()) for n,f in futures.iteritems())
    for n,p in prefixes.iteritems():
        description = "_vs_".join(n)
        touch( ex, p )
        ex.add( p, description="none:"+description, alias=alias )
        ex.add( p+"_peaks.xls", description="macs:"+description+"_peaks.xls",
                associate_to_filename=p, template='%s_peaks.xls' )
        ex.add( p+"_peaks.bed", description="macs:"+description+"_peaks.bed",
                associate_to_filename=p, template='%s_peaks.bed' )
        ex.add( p+"_summits.bed", description="macs:"+description+"_summits.bed",
                associate_to_filename=p, template='%s_summits.bed' )
        if len(n)>1:
            ex.add( p+"_negative_peaks.xls", description="macs:"+description+"_negative_peaks.xls",
                    associate_to_filename=p, template='%s_negative_peaks.xls' )
    return prefixes

@program
def sql_prepare_deconv(sql_dict,peaks_bed,chr_name,chr_length,cutoff,read_length):
    """Prepares files for the deconvolution algorithm.
    Calls the stand-alone ``sql_prepare_deconv.py`` script which needs 
    a directory of sql files (quantitative tracks for forward and reverse strand) 
    and a bed file (*_peaks.bed file from ``macs``) of wnriched regions to consider.
    Returns the name of an 'Rdata' file to be passed to the *R* deconvolution script.
    """
    out = unique_filename_in()
    args = [sql_dict['fwd'],sql_dict['rev'],peaks_bed,out,chr_name,
            str(chr_length),str(cutoff),str(read_length)]
    return {"arguments": ["sql_prepare_deconv.py"]+args, 
            "return_value": out}

@program
def run_deconv_r( counts_file, read_length, chr_name, script_path ):
    """Runs the 'deconv.R' script (found under 'script_path') on the 'counts_file' 
    input with parameters 'chr_name' (name of chromosome to process) and 'read_length'.
    Returns the pdf file and the 'Rdata' output. 
    """
    pdf_file = unique_filename_in()
    output_file = unique_filename_in()
    rargs = [counts_file, pdf_file, str(read_length), 
             chr_name, output_file, script_path]
    return {'arguments': ["R","--vanilla","--slave","-f",
                          os.path.join(script_path,"deconv.R"),
                          "--args"] + rargs,
            'return_value': {'pdf': pdf_file, 'rdata': output_file}}

@program
def sql_finish_deconv(sqlout,rdata):
    """Binds the ``sql_finish_deconv.py`` scripts which creates an sqlite file from
    'run_deconv''s output. 
    """
    return {"arguments": ["sql_finish_deconv.py",rdata,sqlout], 
            "return_value": sqlout}

def run_deconv(ex,sql,peaks,chromosomes,read_length,script_path, via='lsf'):
    """Runs the complete deconvolution process for a set of sql files and a bed file,
    parallelized over a set of chromosomes (a dictionary with keys 'chromosome_id' 
    and values a dictionary with chromosome names and lengths).
   
    Returns a dictionary of file outputs with keys file types 
    ('pdf', 'sql' and 'bed') and values the file names.
    """
    prep_futures = dict((c['name'],
                         sql_prepare_deconv.nonblocking( ex, sql, peaks, 
                                                         c['name'], c['length'],
                                                         1500, read_length, 
                                                         via=via ))
                        for c in chromosomes.values())
    rdeconv_futures = dict((c,
                            run_deconv_r.nonblocking( ex, f.wait(), read_length,
                                                      c, script_path, via=via ))
                           for c,f in prep_futures.iteritems())
    rdeconv_out = dict((c, f.wait()) for c,f in rdeconv_futures.iteritems())
    if len(rdeconv_out)>0:
        pdf_future = join_pdf.nonblocking( ex,
                                           [x['pdf'] for x in rdeconv_out.values()
                                            if x != None],
                                           via=via )
    sqlout = create_sql_track( unique_filename_in(), chromosomes )
    for fout in rdeconv_out.values():
        if fout != None:
            f = sql_finish_deconv.nonblocking( ex, sqlout, fout['rdata'], via=via ).wait()
    outfiles = {}
    outfiles['sql'] = sqlout
    outfiles['bed'] = sqlout+'_deconv.bed'
    if len(rdeconv_out)>0:
        outfiles['pdf'] = pdf_future.wait()
    else:
        outfiles['pdf'] = rdeconv_out.values()[0]['pdf']
    return outfiles
    
###################### Workflow ####################

def workflow_groups( ex, job_or_dict, mapseq_files, chromosomes, script_path='', 
                     via='lsf' ):
    """Runs a chipseq workflow over bam files obtained by mapseq. Will optionally run ``macs`` and 'run_deconv'.
    
    Arguments are:

    * ``'ex'``: a 'bein' execution environment to run jobs in,
    
    * ``'job_or_dict'``: a 'Frontend' 'job' object, or a dictionary with key 'groups' and 'options' if applicable,
    
    * ``'chromosomes'``: a dictionary with keys 'chromosome_id' and values a dictionary with chromosome names and lengths,

    * ``'script_path'``: only needed if 'run_deconv' is in the job options, must point to the location of the R scripts.

    Defaults ``macs`` parameters (overriden by job_or_dict['options']['macs_args']) are set as follows:

    * ``'--nomodel'``

    * ``'p'``: .001 (p-value threshold)

    * ``'bw'``: 200 ('bandwith')

    * ``'m'``: 5,60 ('minimum and maximum enrichments relative to background or control)

    Returns a tuple of a dictionary with keys *group_id* from the job groups, *macs* and *deconv* if applicable and values file description dictionaries and a dictionary of *group_ids* to *names* used in file descriptions.
"""
    options = {}
    if isinstance(job_or_dict,frontend.Job):
        options = job_or_dict.options
        groups = job_or_dict.groups
    elif isinstance(job_or_dict,dict) and 'groups' in job_or_dict:
        if 'options' in job_or_dict:
            options = job_or_dict['options']
        groups = job_or_dict['groups']
        for gid in groups.keys():
            if not 'name' in groups[gid]:
                groups[gid]['name'] = gid
    else:
        raise TypeError("job_or_dict must be a frontend.Job object or a dictionary with key 'groups'.")
    merge_strands = -1
    suffixes = ["fwd","rev"]
    if options.get('merge_strands')>=0:
        merge_strands = options['merge_strands']
        suffixes = ["merged"]
    peak_deconvolution = options.get('peak_deconvolution') or False
    macs_args = options.get('macs_args') or ["--nomodel","-m","5,60","--bw=200","-p",".001"]
    b2w_args = options.get('b2w_args') or []
    if not isinstance(mapseq_files,dict):
        raise TypeError("Mapseq_files must be a dictionary.")
    tests = []
    controls = []
    names = {'tests': [], 'controls': []}
    read_length = []
    genome_size = []
    for gid,mapped in mapseq_files.iteritems():
        group_name = groups[gid]['name']
        if not isinstance(mapped,dict):
            raise TypeError("Mapseq_files values must be dictionaries with keys *run_ids* or 'bam'.")
        if 'bam' in mapped:
            mapped = {'_': mapped}
        futures = {}
        for k in mapped.keys():
            if not 'libname' in mapped[k]:
                mapped[k]['libname'] = group_name+"_"+str(k)
            if not 'stats' in mapped[k]:
                futures[k] = mapseq.bamstats.nonblocking( ex, mapped[k]["bam"], via=via )
        for k in futures.keys():
            mapped[k]['stats'] = f.wait()
        if len(mapped)>1:
            bamfile = merge_bam(ex, [m['bam'] for m in mapped.values()])
        else:
            bamfile = mapped.values()[0]['bam']
        if groups[gid]['control']:
            controls.append(bamfile)
            names['controls'].append(group_name)
        else:
            tests.append(bamfile)
            names['tests'].append(group_name)
            read_length.append(mapped[0]['stats']['read_length'])
            genome_size.append(mapped[0]['stats']['genome_size'])
    if len(controls)<1:
        controls = [None]
        names['controls'] = [None]
    processed['macs'] = add_macs_results( ex, read_length, genome_size,
                                          tests, ctrlbam=controls, name=names, 
                                          shift=merge_strands,
                                          macs_args=macs_args, via=via )
    if peak_deconvolution:
        processed['deconv'] = {}
        wigs = {}
        for gid,mapped in mapseq_files.iteritems():
            if groups[gid]['control']:
                continue
            group_name = groups[gid]['name']
            if merge_strands >= 0:
                wig = [parallel_density_sql( ex, m["bam"], chromosomes, 
                                             nreads=m["stats"]["total"], 
                                             merge=-1,
                                             convert=False, 
                                             description=m['libname'], 
                                             b2w_args=b2w_args, via=via ) 
                       for m in mapped.values()]
            else:
                wig = [m['wig'] for m in mapped]
            if len(wig) > 1:
                merged_wig[group_name] = dict((s, 
                                               merge_sql(ex, [x[s] for x in wig],
                                                         [m['libname'] for m in mapped.values()] ,
                                                         description="sql:"+group_name+"_"+s+".sql")) 
                                              for s in suffixes)
            else:
                merged_wig[group_name] = wig[0]
        for i,name in enumerate(names['tests']):
            if names['controls']==[None]:
                macsbed = processed['macs'][(name,)]+"_peaks.bed"
            else:
                macsbed = merge_many_bed(ex,[processed['macs'][(name,x)]+"_peaks.bed" 
                                             for x in names['controls']],via=via)
            deconv = run_deconv( ex, merged_wig[name], macsbed, chromosomes,
                                 read_length[i], script_path, via=via )
            [ex.add(v, description=k+':'+name+'_deconv.'+k)
             for k,v in deconv.iteritems()]
            processed['deconv'][name] = deconv
    return processed

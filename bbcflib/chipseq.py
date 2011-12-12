"""
=======================
Module: bbcflib.chipseq
=======================

This module provides functions to run a basic ChIP-seq analysis from reads mapped on
a reference genome.
The most important steps are the binding of ``macs`` via the 'add_macs_results' function and
the peak deconvolution algorithm with the 'run_deconv' function.

The whole workflow can be run via the function ``workflow_groups`` with appropriate options in 'job.options'.

Below is the script used by the frontend::

    from bbcflib import daflims, genrep, frontend, email, gdv, common
    from bbcflib.mapseq import *
    from bbcflib.chipseq import *
    M = MiniLIMS( '/path/to/chipseq/minilims' )
    ms_limspath = '/path/to/mapseq/minilims'
    working_dir = '/path/to/scratch/on/cluster'
    hts_key = 'test_key'
    assembly_id = 'mm9'
    gl = { 'hts_chipseq': {'url': 'http://htsstation.vital-it.ch/chipseq/'},
           'hts_mapseq': {'url': 'http://htsstation.vital-it.ch/mapseq/'},
           'genrep_url': 'http://bbcftools.vital-it.ch/genrep/',
           'script_path': '/srv/chipseq/lib' }
    htss = frontend.Frontend( url=gl['hts_chipseq']['url'] )
    g_rep = genrep.GenRep( gl['genrep_url'], "" )
    job = htss.job( hts_key )
    g_rep = genrep.GenRep( gl['genrep_url'], "" )
    assembly = g_rep.assembly( assembly_id )
    with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
        (ms_files, job) = get_bam_wig_files( ex, job, ms_limspath, gl['hts_mapseq']['url'], gl['script_path'], via=via )
        files = workflow_groups( ex, job, ms_files, assembly.chromosomes, gl['script_path'] )
    print ex.id
    allfiles = common.get_files( ex.id, M )
    print allfiles
"""

# Built-in modules #
import shutil, pickle, urllib, re, os

# Internal modules #
from bbcflib import frontend, mapseq
from bbcflib.common import merge_sql, merge_many_bed, join_pdf, cat, set_file_descr, unique_filename_in

# Other modules #
from bein import program
from bein.util import touch

################################################################################
# Peaks and annotation #

@program
def macs( read_length, genome_size, bamfile, ctrlbam = None, args = None ):
    """Binding for the ``macs`` peak caller.

    takes one (optionally two) bam file(s) and
    the 'read_length' and 'genome_size' parameters passed to ``macs``.

    Returns the file prefix ('-n' option of ``macs``)
    """
    if args is None:
        args = []
    outname = unique_filename_in()
    macs_args = ["macs14","-t",bamfile]
    if ctrlbam != None:
        macs_args += ["-c",ctrlbam]
    macs_args += ["-n",outname,"-f","BAM","-g",str(genome_size),"-s",str(read_length)]
    return {"arguments": macs_args+args, "return_value": outname}

def add_macs_results( ex, read_length, genome_size, bamfile,
                      ctrlbam=None, name=None, poisson_threshold=None,
                      alias=None, macs_args=None, via='lsf' ):
    """Calls the ``macs`` function on each possible pair
    of test and control bam files and adds
    the respective outputs to the execution repository.

    ``macs`` options can be controlled with `macs_args`.
    If a dictionary of Poisson thresholds for each sample is given, then the enrichment bounds ('-m' option)
    are computed from them otherwise the default is '-m 10,100'.

    Returns the set of file prefixes.
    """
    if not(isinstance(bamfile,list)):
        bamfile = [bamfile]
    if not(isinstance(ctrlbam,list)):
        ctrlbam = [ctrlbam]
    if poisson_threshold is None:
        poisson_threshold = {}
    if macs_args is None:
        macs_args = []
    futures = {}
    rl = read_length
    for i,bam in enumerate(bamfile):
        n = name['tests'][i]
        if poisson_threshold.get(n)>0:
            low = (poisson_threshold.get(n)+1)*5
            enrich_bounds = str(min(30,low))+","+str(10*low)
        else:
            enrich_bounds = "10,100"
        if isinstance(read_length,list):
            rl = read_length[i]
        for j,cam in enumerate(ctrlbam):
            m = name['controls'][j]
            if m == None:
                nm = (n,)
            else:
                nm = (n,m)
            futures[nm] = macs.nonblocking( ex, rl, genome_size, bam, cam,
                                            args=macs_args+["-m",enrich_bounds],
                                            via=via )
    prefixes = dict((n,f.wait()) for n,f in futures.iteritems())
    for n,p in prefixes.iteritems():
        macs_descr0 = {'step':'macs','type':'none','view':'admin'}
        macs_descr1 = {'step':'macs','type':'xls','group':n[0]}
        macs_descr2 = {'step':'macs','type':'bed','group':n[0]}
        filename = "_vs_".join(n)
        touch( ex, p )
        ex.add( p, description=set_file_descr(filename,**macs_descr0), alias=alias )
        ex.add( p+"_peaks.xls",
                description=set_file_descr(filename+"_peaks.xls",**macs_descr1),
                associate_to_filename=p, template='%s_peaks.xls' )
        ex.add( p+"_peaks.bed",
                description=set_file_descr(filename+"_peaks.bed",**macs_descr2),
                associate_to_filename=p, template='%s_peaks.bed' )
        ex.add( p+"_summits.bed",
                description=set_file_descr(filename+"_summits.bed",**macs_descr2),
                associate_to_filename=p, template='%s_summits.bed' )
        if len(n)>1:
            ex.add( p+"_negative_peaks.xls",
                    description=set_file_descr(filename+"_negative_peaks.xls",**macs_descr1),
                    associate_to_filename=p, template='%s_negative_peaks.xls' )
    return prefixes

@program
def sql_prepare_deconv(sql_dict, peaks_bed, chr_name, chr_length, cutoff, read_extension):
    """Prepares files for the deconvolution algorithm.
    Calls the stand-alone ``sql_prepare_deconv.py`` script which needs
    a dictionary of sql files (quantitative tracks for forward and reverse strand)
    and a bed file (*_peaks.bed file from ``macs``) of enriched regions to consider.
    Returns the name of an 'Rdata' file to be passed to the *R* deconvolution script.
    """
    out = unique_filename_in()
    args = [sql_dict['fwd'],sql_dict['rev'],peaks_bed,out,chr_name,
            str(chr_length),str(cutoff),str(read_extension)]
    return {"arguments": ["sql_prepare_deconv.py"]+args,
            "return_value": out}

@program
def run_deconv_r( counts_file, read_extension, chr_name, script_path ):
    """Runs the 'deconv.R' script (found under 'script_path') on the 'counts_file'
    input with parameters 'chr_name' (name of chromosome to process) and 'read_extension'.
    Returns the pdf file and the 'Rdata' output.
    """
    pdf_file = unique_filename_in()
    output_file = unique_filename_in()
    rargs = [counts_file, pdf_file, str(read_extension),
             chr_name, output_file, script_path]
    return {'arguments': ["R","--vanilla","--slave","-f",
                          os.path.join(script_path,"deconv.R"),
                          "--args"] + rargs,
            'return_value': {'pdf': pdf_file, 'rdata': output_file}}

@program
def sql_finish_deconv(sqlout,rdata,chrom):
    """Binds the ``sql_finish_deconv.py`` scripts which creates an sqlite file from
    'run_deconv''s output.
    """
    return {"arguments": ["sql_finish_deconv.py",rdata,sqlout,chrom],
            "return_value": sqlout}

def run_deconv(ex, sql, peaks, chromosomes, read_extension, script_path, via = 'lsf'):
    """Runs the complete deconvolution process for a set of sql files and a bed file,
    parallelized over a set of chromosomes (a dictionary with keys 'chromosome_id'
    and values a dictionary with chromosome names and lengths).

    Returns a dictionary of file outputs with keys file types
    ('pdf', 'sql' and 'bed') and values the file names.
    """
    from .track import new
    rdeconv_futures = {}
    for c in chromosomes.values():
        f = sql_prepare_deconv.nonblocking( ex, sql, peaks, c['name'], c['length'],
                                            1500, read_extension, via='local' ).wait()
        rdeconv_futures[c['name']] =  run_deconv_r.nonblocking( ex, f, read_extension,
                                                                c['name'], script_path, via=via )
    rdeconv_out = dict((c, f.wait()) for c,f in rdeconv_futures.iteritems())
    if len(rdeconv_out)>0:
        pdf_future = join_pdf.nonblocking( ex,
                                           [x['pdf'] for x in rdeconv_out.values()
                                            if x != None],
                                           via=via )
    chrlist = dict((v['name'], {'length': v['length']}) for v in chromosomes.values())
    output = unique_filename_in()
    with new(output, 'sql', datatype="quantitative", chrmeta=chrlist) as track:
        pass
    for c,fout in rdeconv_out.iteritems():
        if fout != None:
            f = sql_finish_deconv.nonblocking( ex, output, fout['rdata'], c, via='local' ).wait()
    outfiles = {}
    outfiles['sql'] = output
    outfiles['bed'] = output+'_deconv.bed'
    if len(rdeconv_out)>0:
        outfiles['pdf'] = pdf_future.wait()
    else:
        outfiles['pdf'] = rdeconv_out.values()[0]['pdf']
    return outfiles

def filter_deconv( bedfile, pval=0.5 ):
    """Filters a bedfile created by deconvolution to select peaks with p-value smaller than 'pval'.
    Returns the filtered file name."""
    outname = unique_filename_in()+".bed"
    with open( bedfile, "r" ) as fin:
        with open( outname, "w" ) as fout:
            for row in fin:
                rowpval = float(re.search(r';FERR=([\d\.]+)\s',row).groups()[0])
                if rowpval <= pval:
                    fout.write(row)
    return outname

def filter_macs( bedfile, ntags=5 ):
    """Filters MACS' summits file to select peaks with a number of tags greater than 'ntags'.
    Returns the filtered file name."""
    outname = unique_filename_in()+".bed"
    with open( bedfile, "r" ) as fin:
        with open( outname, "w" ) as fout:
            for row in fin:
                splitrow = row.split("\t")
                if len(splitrow)>4 and float(splitrow[4])>ntags:
                    fout.write(row)
    return outname

################################################################################
# Workflow #

def workflow_groups( ex, job_or_dict, mapseq_files, chromosomes, script_path='',
                     genrep=None, logfile=None, via='lsf' ):
    """Runs a chipseq workflow over bam files obtained by mapseq. Will optionally run ``macs`` and 'run_deconv'.

    :param ex: a 'bein' execution environment to run jobs in,

    :param job_or_dict: a 'Frontend' 'job' object, or a dictionary with key 'groups' and 'options' if applicable,

    :param chromosomes: a dictionary with keys 'chromosome_id' and values a dictionary with chromosome names and lengths,

    :param script_path: only needed if 'run_deconv' is in the job options, must point to the location of the R scripts.

    Defaults ``macs`` parameters (overriden by ``job_or_dict['options']['macs_args']``) are set as follows:

    * ``'-bw'``: 200 ('bandwith')

    * ``'-m'``: 10,100 ('minimum and maximum enrichments relative to background or control')

    The enrichment bounds will be computed from a Poisson threshold *T*, if available, as *(min(30,5*(T+1)),50*(T+1))*.

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
            if not('name' in groups[gid]):
                groups[gid]['name'] = gid
    else:
        raise TypeError("job_or_dict must be a frontend.Job object or a dictionary with key 'groups'.")
    merge_strands = -1
    suffixes = ["fwd","rev"]
    if options.get('merge_strands') and int(options['merge_strands'])>=0:
        merge_strands = int(options['merge_strands'])
    peak_deconvolution = options.get('peak_deconvolution') or False
    if isinstance(peak_deconvolution,str):
        peak_deconvolution = peak_deconvolution.lower() in ['1','true','t']
    run_meme = options.get('run_meme') or False
    if isinstance(run_meme,str):
        run_meme = run_meme.lower() in ['1','true','t']
    macs_args = options.get('macs_args') or ["--bw=200"]
    b2w_args = options.get('b2w_args') or []
    if not(isinstance(mapseq_files,dict)):
        raise TypeError("Mapseq_files must be a dictionary.")
    tests = []
    controls = []
    names = {'tests': [], 'controls': []}
    read_length = []
    p_thresh = {}
    for gid,mapped in mapseq_files.iteritems():
        group_name = groups[gid]['name']
        if not(isinstance(mapped,dict)):
            raise TypeError("Mapseq_files values must be dictionaries with keys *run_ids* or 'bam'.")
        if 'bam' in mapped:
            mapped = {'_': mapped}
        futures = {}
        ptruns = []
        for k in mapped.keys():
            if not 'libname' in mapped[k]:
                mapped[k]['libname'] = group_name+"_"+str(k)
            if not 'stats' in mapped[k]:
                futures[k] = mapseq.bamstats.nonblocking( ex, mapped[k]["bam"], via=via )
            if mapped[k].get('poisson_threshold')>0:
                ptruns.append(mapped[k]['poisson_threshold'])
        if len(ptruns)>0:
            p_thresh['group_name'] = sum(ptruns)/len(ptruns)
        for k in futures.keys():
            mapped[k]['stats'] = f.wait()
        if len(mapped)>1:
            bamfile = mapseq.merge_bam(ex, [m['bam'] for m in mapped.values()])
        else:
            bamfile = mapped.values()[0]['bam']
        if groups[gid]['control']:
            controls.append(bamfile)
            names['controls'].append(group_name)
        else:
            tests.append(bamfile)
            names['tests'].append(group_name)
            read_length.append(mapped.values()[0]['stats']['read_length'])
    genome_size = mapped.values()[0]['stats']['genome_size']
    if len(controls)<1:
        controls = [None]
        names['controls'] = [None]
    if logfile:
        logfile.write("Starting MACS.\n");logfile.flush()
    processed = {'macs': add_macs_results( ex, read_length, genome_size,
                                           tests, ctrlbam=controls, name=names,
                                           poisson_threshold=p_thresh,
                                           macs_args=macs_args, via=via ) }
    if logfile:
        logfile.write("Done MACS.\n");logfile.flush()
    peak_list = {}
    if run_meme:
        import gMiner
        chrlist = dict((v['name'], {'length': v['length']}) for v in chromosomes.values())
        for i,name in enumerate(names['tests']):
            if names['controls']==[None]:
                macsbed = processed['macs'][(name,)]+"_summits.bed"
            else:
                macsbed = cat([processed['macs'][(name,x)]+"_summits.bed" for x in names['controls']])
            macsbed = filter_macs(macsbed)
            outdir = unique_filename_in()
            os.mkdir(outdir)
            gm_out = gMiner.run(
                track1 = macsbed,
                track1_name = group_name+' peak summits',
                track1_chrs = chrlist,
                operation_type = 'genomic_manip',
                manipulation = 'neighborhood',
                output_location = outdir,
                before_start = -150,
                after_end = 150 )[0]
            outdir = unique_filename_in()
            os.mkdir(outdir)
            gm_out = gMiner.run(
                track1 = gm_out,
                track1_name = group_name+' peak summits +- 150',
                operation_type = 'genomic_manip',
                manipulation = 'merge',
                output_location = outdir )[0]
            peak_list[name] = gm_out
    if peak_deconvolution:
        processed['deconv'] = {}
        merged_wig = {}
        if not('read_extensions' in options and int(options['read_extension'])>0):
            options['read_extension'] = read_length[0]
        if not('-q' in b2w_args):
            b2w_args += ["-q",str(options['read_extension'])]
        for gid,mapped in mapseq_files.iteritems():
            if groups[gid]['control']:
                continue
            group_name = groups[gid]['name']
            wig = []
            for m in mapped.values():
                if merge_strands >= 0 or not('wig' in m) or len(m['wig'])<2:
                    output = mapseq.parallel_density_sql( ex, m["bam"], chromosomes,
                                                          nreads=m["stats"]["total"],
                                                          merge=-1,
                                                          convert=False,
                                                          b2w_args=b2w_args, via=via )
                    wig.append(dict((s,output+s+'.sql') for s in suffixes))
                else:
                    wig.append(m['wig'])
            if len(wig) > 1:
                merged_wig[group_name] = dict((s,merge_sql(ex, [x[s] for x in wig],
                                                           [m['libname'] for m in mapped.values()],
                                                           via='local'))
                                              for s in suffixes)
            else:
                merged_wig[group_name] = wig[0]
        for name in names['tests']:
            if logfile:
                logfile.write(name+" deconvolution.\n");logfile.flush()
            if names['controls']==[None]:
                macsbed = processed['macs'][(name,)]+"_peaks.bed"
            else:
                macsbed = merge_many_bed(ex,[processed['macs'][(name,x)]+"_peaks.bed"
                                             for x in names['controls']],via=via)
            deconv = run_deconv( ex, merged_wig[name], macsbed, chromosomes,
                                 options['read_extension'], script_path, via=via )
            [ex.add(v, description=set_file_descr(name+'_deconv.'+k,type=k,step='deconvolution',group=name))
             for k,v in deconv.iteritems()]
            processed['deconv'][name] = deconv
            peak_list[name] = filter_deconv( deconv['bed'], pval=0.65 )
    if run_meme and not(genrep == None):
        from .motif import parallel_meme
        if logfile:
            logfile.write("Starting MEME.\n");logfile.flush()
        processed['meme'] = parallel_meme( ex, genrep, chromosomes,
                                           peak_list.values(), name=peak_list.keys(),
                                           meme_args=['-nmotifs','4','-revcomp'], via=via )
    return processed

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

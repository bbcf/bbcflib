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
import shutil, pickle, urllib

# Internal modules #
from . import frontend, mapseq
from .common import *

# Other modules #
from bein import *
from bein.util import *

################################################################################
# Peaks and annotation #

@program
def macs( read_length, genome_size, bamfile, ctrlbam=None, args=[] ):
    """Binding for the ``macs`` peak caller.

    takes one (optionally two) bam file(s) and
    the 'read_length' and 'genome_size' parameters passed to ``macs``.

    Returns the file prefix ('-n' option of ``macs``)
    """
    outname = unique_filename_in()
    macs_args = ["macs14","-t",bamfile]
    if ctrlbam != None:
        macs_args += ["-c",ctrlbam]
    macs_args += ["-n",outname,"-f","BAM",
                  "-g",str(genome_size),"-s",str(read_length)]
    return {"arguments": macs_args+args, "return_value": outname}

def add_macs_results( ex, read_length, genome_size, bamfile,
                      ctrlbam=[None], name=None, poisson_threshold={},
                      alias=None, macs_args=[], via='lsf' ):
    """Calls the ``macs`` function on each possible pair
    of test and control bam files and adds
    the respective outputs to the execution repository.

    ``macs`` options can be controlled with `macs_args`.
    If a dictionary of Poisson thresholds for each sample is given, then the enrichment bounds ('-m' option)
    are computed from them otherwise the default is '-m 5,60'.

    Returns the set of file prefixes.
    """
    if not(isinstance(bamfile,list)):
        bamfile = [bamfile]
    if not(isinstance(ctrlbam,list)):
        ctrlbam = [ctrlbam]
    futures = {}
    rl = read_length
    for i,bam in enumerate(bamfile):
        n = name['tests'][i]
        if poisson_threshold.get(n)>0:
            low = poisson_threshold.get(n)
            enrich_bounds = str(min(30,low+1))+","+str((low+1)*30)
        else:
            enrich_bounds = "5,60"
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
def sql_prepare_deconv(sql_dict,peaks_bed,chr_name,chr_length,cutoff,read_extension):
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
def sql_finish_deconv(sqlout,rdata):
    """Binds the ``sql_finish_deconv.py`` scripts which creates an sqlite file from
    'run_deconv''s output.
    """
    return {"arguments": ["sql_finish_deconv.py",rdata,sqlout],
            "return_value": sqlout}

def run_deconv(ex,sql,peaks,chromosomes,read_extension,script_path, via='lsf'):
    """Runs the complete deconvolution process for a set of sql files and a bed file,
    parallelized over a set of chromosomes (a dictionary with keys 'chromosome_id'
    and values a dictionary with chromosome names and lengths).

    Returns a dictionary of file outputs with keys file types
    ('pdf', 'sql' and 'bed') and values the file names.
    """
    prep_futures = dict((c['name'],
                         sql_prepare_deconv.nonblocking( ex, sql, peaks,
                                                         c['name'], c['length'],
                                                         1500, read_extension,
                                                         via=via ))
                        for c in chromosomes.values())
    rdeconv_futures = dict((c,
                            run_deconv_r.nonblocking( ex, f.wait(), read_extension,
                                                      c, script_path, via=via ))
                           for c,f in prep_futures.iteritems())
    rdeconv_out = dict((c, f.wait()) for c,f in rdeconv_futures.iteritems())
    if len(rdeconv_out)>0:
        pdf_future = join_pdf.nonblocking( ex,
                                           [x['pdf'] for x in rdeconv_out.values()
                                            if x != None],
                                           via=via )
    sqlout = create_sql_track( unique_filename_in(), chromosomes.values() )
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

################################################################################
# Workflow #

def get_bam_wig_files( ex, job, minilims=None, hts_url=None,
                       script_path = './', via='lsf' ):
    """
    Will replace file references by actual file paths in the 'job' object.
    These references are either 'mapseq' keys or urls.
    """
    mapped_files = {}
    read_exts = {}
    suffix = ['fwd','rev']
    for gid,group in job.groups.iteritems():
        mapped_files[gid] = {}
        if 'name' in group:
            group_name = re.sub(r'\s+','_',group['name'])
        else:
            group_name = str(gid)
        job.groups[gid]['name'] = group_name
        for rid,run in group['runs'].iteritems():
            file_loc = str(run['url'])
            bamfile = unique_filename_in()
            wig = {}
            name = group_name
            s = None
            p_thresh = None
            if len(group['runs'])>1:
                if all([x in run for x in ['machine','run','lane']]):
                    name += "_".join(['',run['machine'],str(run['run']),str(run['lane'])])
                else:
                    name += "_"+file_loc.split("/")[-1]
            if file_loc.startswith("http://") or file_loc.startswith("https://"):
                urllib.urlretrieve( file_loc, bamfile )
                urllib.urlretrieve( file_loc+".bai", bamfile+".bai" )
            elif os.path.exists(file_loc):
                shutil.copy( file_loc, bamfile )
                shutil.copy( file_loc+".bai", bamfile+".bai" )
            elif os.path.exists(minilims) and os.path.exists(os.path.join(minilims+".files",file_loc)):
                MMS = MiniLIMS(minilims)
                file_loc = os.path.join(minilims+".files",file_loc)
                shutil.copy( file_loc, bamfile )
                shutil.copy( file_loc+".bai", bamfile+".bai" )
                exid = max(MMS.search_executions(with_text=run['key']))
                allfiles = dict((MMS.fetch_file(x)['description'],x)
                                for x in MMS.search_files(source=('execution',exid)))
                with open(MMS.path_to_file(allfiles['py:file_names'])) as q:
                    file_names = pickle.load(q)
#                ms_name = file_names[gid][rid]
                stats_id = allfiles.get("py:"+name+"_filter_bamstat") or allfiles.get("py:"+name+"_full_bamstat")
                with open(MMS.path_to_file(stats_id)) as q:
                    s = pickle.load(q)
                p_thresh = -1
                if "py:"+name+"_Poisson_threshold" in allfiles:
                    pickle_thresh = allfiles["py:"+name+"_Poisson_threshold"]
                    with open(MMS.path_to_file(pickle_thresh)) as q:
                        p_thresh = pickle.load(q)
                if 'py:gdv_json' in allfiles:
                    with open(MMS.path_to_file(allfiles['py:gdv_json'])) as q:
                        job.options['gdv_project'] = pickle.load(q)
                htss = frontend.Frontend( url=hts_url )
                ms_job = htss.job( run['key'] )
                if ms_job.options.get('read_extension')>0 and ms_job.options.get('read_extension')<80:
                    read_exts[rid] = ms_job.options['read_extension']
                else:
                    read_exts[rid] = s['read_length']
                if ms_job.options.get('compute_densities') and ms_job.options.get('merge_strands')<0:
                    wigfile = unique_filename_in()
                    wig_ids = dict(((allfiles['sql:'+name+'_'+s+'.sql'],s),
                                    wigfile+'_'+s+'.sql') for s in suffix)
                    [MMS.export_file(x[0],s) for x,s in wig_ids.iteritems()]
                    wig = dict((x[1],s) for x,s in wig_ids.iteritems())
            else:
                raise ValueError("Couldn't find this bam file anywhere: %s." %file_loc)
            mapped_files[gid][rid] = {'bam': bamfile,
                                      'stats': s or mapseq.bamstats.nonblocking( ex, bamfile, via=via ),
                                      'poisson_threshold': p_thresh,
                                      'libname': name,
                                      'wig': wig}
    if len(read_exts)>0 and not('read_extension' in job.options):
        c = dict((x,0) for x in read_exts.values())
        for x in read_exts.values():
            c[x]+=1
        job.options['read_extension'] = [k for k,v in c.iteritems() if v==max(c.values())][0]
    for gid, group in job.groups.iteritems():
        for rid,run in group['runs'].iteritems():
            if read_exts.get(rid) != job.options['read_extension']:
                mapped_files[gid][rid]['wig'] = []
            if not(isinstance(mapped_files[gid][rid]['stats'],dict)):
                stats = mapped_files[gid][rid]['stats'].wait()
                mapped_files[gid][rid]['stats'] = stats
                pdf = maspeq.add_pdf_stats( ex, {gid:{rid:{'stats':stats}}},
                                            {gid: mapped_files[gid][rid]['libname']},
                                            script_path )
                mapped_files[gid][rid]['p_thresh'] = mapseq.poisson_threshold( 50*stats["actual_coverage"] )
    return (mapped_files,job)


def workflow_groups( ex, job_or_dict, mapseq_files, chromosomes, script_path='',
                     via='lsf' ):
    """Runs a chipseq workflow over bam files obtained by mapseq. Will optionally run ``macs`` and 'run_deconv'.

    Arguments are:

    * ``'ex'``: a 'bein' execution environment to run jobs in,

    * ``'job_or_dict'``: a 'Frontend' 'job' object, or a dictionary with key 'groups' and 'options' if applicable,

    * ``'chromosomes'``: a dictionary with keys 'chromosome_id' and values a dictionary with chromosome names and lengths,

    * ``'script_path'``: only needed if 'run_deconv' is in the job options, must point to the location of the R scripts.

    Defaults ``macs`` parameters (overriden by job_or_dict['options']['macs_args']) are set as follows:

    * ``'-p'``: .001 (p-value threshold)

    * ``'-bw'``: 200 ('bandwith')

    * ``'-m'``: 5,60 ('minimum and maximum enrichments relative to background or control')

    The enrichment bounds will be computed from a Poisson threshold *T*, if available, as *(min(30,T+1),30(T+1))*.

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
    if options.get('merge_strands')>=0:
        merge_strands = options['merge_strands']
    peak_deconvolution = options.get('peak_deconvolution') or False
    macs_args = options.get('macs_args') or ["--bw=200","-p",".001"]
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
            bamfile = merge_bam(ex, [m['bam'] for m in mapped.values()])
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
    processed = {'macs': add_macs_results( ex, read_length, genome_size,
                                           tests, ctrlbam=controls, name=names,
                                           poisson_threshold=p_thresh,
                                           macs_args=macs_args, via=via ) }
    if peak_deconvolution:
        processed['deconv'] = {}
        merged_wig = {}
        if not(options.get('read_extension')>0):
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
                    output = unique_filename_in()
                    touch(ex,output)
                    [create_sql_track( output+s+'.sql', chromosomes.values(),
                                       name=m['libname'] ) 
                     for s in suffixes]
                    mapseq.parallel_density_sql( ex, m["bam"],
                                                 output, chromosomes,
                                                 nreads=m["stats"]["total"],
                                                 merge=-1,
                                                 convert=False,
                                                 b2w_args=b2w_args, via=via )
                    wig.append(dict((s,output+s+'.sql') for s in suffixes))
                else:
                    wig.append(m['wig'])
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
                                 options['read_extension'], script_path, via=via )
            [ex.add(v, description=k+':'+name+'_deconv.'+k)
             for k,v in deconv.iteritems()]
            processed['deconv'][name] = deconv
    return processed

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

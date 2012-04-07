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

    from bbcflib import daflims, genrep, frontend, email, common
    from bbcflib.mapseq import *
    from bbcflib.chipseq import *
    M = MiniLIMS( '/path/to/chipseq/minilims' )
    ms_limspath = '/path/to/mapseq/minilims'
    working_dir = '/path/to/scratch/on/cluster'
    hts_key = 'test_key'
    assembly_id = 'mm9'
    gl = { 'hts_chipseq': {'url': 'http://htsstation.vital-it.ch/chipseq/'},
           'hts_mapseq': {'url': 'http://htsstation.vital-it.ch/mapseq/'},
           'script_path': '/srv/chipseq/lib' }
    htss = frontend.Frontend( url=gl['hts_chipseq']['url'] )
    job = htss.job( hts_key )
    assembly = genrep.Assembly( assembly=assembly_id )
    with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
        (ms_files, job) = get_bam_wig_files( ex, job, ms_limspath, gl['hts_mapseq']['url'], gl['script_path'], via=via )
        files = workflow_groups( ex, job, ms_files, assembly, gl['script_path'] )
    print ex.id
    allfiles = common.get_files( ex.id, M )
    print allfiles
"""

# Built-in modules #
import re, os, gzip, sys

# Internal modules #
from bbcflib import frontend, mapseq, common
from bbcflib import btrack as track
from bbcflib import bFlatMajor as gMiner

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
    outname = common.unique_filename_in()
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
            nm = (n,m)
            futures[nm] = macs.nonblocking( ex, rl, genome_size, bam, cam,
                                            args=macs_args+["-m",enrich_bounds],
                                            via=via )
    prefixes = dict((n,f.wait()) for n,f in futures.iteritems())
    for n,p in prefixes.iteritems():
        macs_descr0 = {'step':'macs','type':'none','view':'admin'}
        macs_descr1 = {'step':'macs','type':'xls','group':n[0]}
        macs_descr2 = {'step':'macs','type':'bed','group':n[0],'ucsc':'1'}
        filename = "_vs_".join([x for x in n if not(x is None)])
        touch( ex, p )
        ex.add( p, description=common.set_file_descr(filename,**macs_descr0), 
                alias=alias )
        ex.add( p+"_peaks.xls",
                description=common.set_file_descr(filename+"_peaks.xls",**macs_descr1),
                associate_to_filename=p, template='%s_peaks.xls' )
        bedzip = gzip.open(p+"_peaks.bed.gz",'wb')
        bedzip.write("track name='"+filename+"_macs_peaks'\n")
        with open(p+"_peaks.bed") as bedinf:
            [bedzip.write(l) for l in bedinf]
        bedzip.close()
        ex.add( p+"_peaks.bed.gz",
                description=common.set_file_descr(filename+"_peaks.bed.gz",**macs_descr2),
                associate_to_filename=p, template='%s_peaks.bed.gz' )
        bedzip = gzip.open(p+"_summits.bed.gz",'wb')
        bedzip.write("track name='"+filename+"_macs_summits'\n")
        with open(p+"_summits.bed") as bedinf:
            [bedzip.write(l) for l in bedinf]
        bedzip.close()
        ex.add( p+"_summits.bed.gz",
                description=common.set_file_descr(filename+"_summits.bed",**macs_descr2),
                associate_to_filename=p, template='%s_summits.bed.gz' )
        if not(n[1] is None):
            ex.add( p+"_negative_peaks.xls",
                    description=common.set_file_descr(filename+"_negative_peaks.xls",**macs_descr1),
                    associate_to_filename=p, template='%s_negative_peaks.xls' )
    return prefixes

@program
def camelPeaks( scores_fwd, scores_rev, peaks, chromosome_name, chromosome_length,
                read_extension, script_path ):
    """Runs the 'camelPeaks.py' wrapper script on the
    'scores_fwd', 'scores_rev' and 'peaks'
    input with parameters 'chromosome_name' (name of chromosome to process)
    and 'read_extension', using functions from 'script_path'/deconv_fcts.R.
    Returns a pdf file and several data tracks.
    """
    output = common.unique_filename_in()
    args = ["-p",peaks,"-f",scores_fwd,"-r",scores_rev,"-o",output,"-c",chromosome_name,
            "-l",str(chromosome_length),"-e",str(read_extension),"-z",script_path,"-s","1500"]
    return {'arguments': ["camelPeaks.py"]+args,
            'return_value': None}

def run_deconv(ex, sql, peaks, chromosomes, read_extension, script_path, via = 'lsf'):
    """Runs the complete deconvolution process for a set of sql files and a bed file,
    parallelized over a set of chromosomes (a dictionary with keys 'chromosome_id'
    and values a dictionary with chromosome names and lengths).

    Returns a dictionary of file outputs with keys file types
    ('pdf', 'sql' and 'bed') and values the file names.
    """
    deconv_futures = {}
    stdout_files = {}
    for cid,chr in chromosomes.iteritems():
        stdout_files[cid] = common.unique_filename_in()
        deconv_futures[cid] = camelPeaks.nonblocking( ex, sql['fwd'], sql['rev'], peaks,
                                                      chr['name'], chr['length'],
                                                      read_extension, script_path,
                                                      via=via, stdout=stdout_files[cid] )
    deconv_out = {}
    for c,f in deconv_futures.iteritems():
        f.wait()
        scan = False
        deconv_out[c] = []
        with open(stdout_files[c]) as sout:
            for row in sout:
                if scan:
                    outf = row.strip()
                    if os.path.exists(outf):
                        deconv_out[c].append(outf)
                if re.search(r'\*\*\*\*\*\*\*\*\*\*\*\*OUTPUT FILES\*\*\*\*\*\*\*\*\*\*',row):
                    scan = True
    if len(deconv_out)>0:
        pdf_future = common.join_pdf.nonblocking( ex,
                                                  [x[0] for x in deconv_out.values() if len(x)>0],
                                                  via=via )
    chrlist = dict((v['name'], {'length': v['length']}) for v in chromosomes.values())
    output = common.unique_filename_in()
    outfiles = {}
    outfiles['bed'] = output+"_peaks.sql"
    outfiles['sql'] = output+"_deconv.sql"
    outbed = track.track(outfiles['bed'],
                         fields=["start","end","score","name"],
                         chrmeta=chrlist,
                         info={'datatype':'qualitative'})
    outwig = track.track(outfiles['sql'],
                         fields=["start","end","score"],
                         chrmeta=chrlist,
                         info={'datatype':'quantitative'})
    for c,fout in deconv_out.iteritems():
        if len(fout) < 3: continue
        outbed.write(track.track(fout[1]).read())
        outwig.write(track.track(fout[2]).read())
    outbed.close()
    outwig.close()
    if len(deconv_out)>0:
        outfiles['pdf'] = pdf_future.wait()
    else:
        outfiles['pdf'] = deconv_out.values()[0][0]
    return outfiles

################################################################################
# Workflow #

def workflow_groups( ex, job_or_dict, mapseq_files, assembly, script_path='',
                     logfile=None, via='lsf' ):
    """Runs a chipseq workflow over bam files obtained by mapseq. Will optionally run ``macs`` and 'run_deconv'.

    :param ex: a 'bein' execution environment to run jobs in,

    :param job_or_dict: a 'Frontend' 'job' object, or a dictionary with key 'groups' and 'options' if applicable,

    :param mapseq_files: a dictionary of files fomr a 'mapseq' execution, imported via 'mapseq.get_bam_wig_files',

    :param assembly: a genrep.Assembly object,

    :param script_path: only needed if 'run_deconv' is in the job options, must point to the location of the R scripts.

    Defaults ``macs`` parameters (overriden by ``job_or_dict['options']['macs_args']``) are set as follows:

    * ``'-bw'``: 200 ('bandwith')

    * ``'-m'``: 10,100 ('minimum and maximum enrichments relative to background or control')

    The enrichment bounds will be computed from a Poisson threshold *T*, if available, as *(min(30,5*(T+1)),50*(T+1))*.

    Returns a tuple of a dictionary with keys *group_id* from the job groups, *macs* and *deconv* if applicable and values file description dictionaries and a dictionary of *group_ids* to *names* used in file descriptions.
"""
    options = {}
    if logfile is None: logfile = sys.stdout
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
    merge_strands = int(options.get('merge_strands',-1))
    suffixes = ["fwd","rev"]
    peak_deconvolution = options.get('peak_deconvolution',False)
    if isinstance(peak_deconvolution,basestring):
        peak_deconvolution = peak_deconvolution.lower() in ['1','true','t']
    run_meme = options.get('run_meme',False)
    if isinstance(run_meme,basestring):
        run_meme = run_meme.lower() in ['1','true','t']
    macs_args = options.get('macs_args',["--bw=200"])
    b2w_args = options.get('b2w_args',[])
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
            if mapped[k].get('poisson_threshold',-1)>0:
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
    logfile.write("Starting MACS.\n");logfile.flush()
    processed = {'macs': add_macs_results( ex, read_length, genome_size,
                                           tests, ctrlbam=controls, name=names,
                                           poisson_threshold=p_thresh,
                                           macs_args=macs_args, via=via ) }
    logfile.write("Done MACS.\n");logfile.flush()
    peak_list = {}
    if run_meme:
        chrlist = assembly.chrmeta
        _select = {'score':(6,sys.maxint)}
        _fields = ['chr','start','end','name','score']
        for i,name in enumerate(names['tests']):
            if len(names['controls']) < 2:
                ctrl = (name,names['controls'][0])
                macsbed = track.track(processed['macs'][ctrl]+"_summits.bed",
                                      chrmeta=chrlist, fields=_fields).read(selection=_select)
            else:
                macsbed = gMiner.stream.concatenate(
                    [track.track(processed['macs'][(name,x)]+"_summits.bed",
                                 chrmeta=chrlist, fields=_fields).read(selection=_select)
                     for x in names['controls']])
            ##############################
            macs_neighb = gMiner.stream.neighborhood(macsbed, before_start=150, after_end=150 )
            peak_list[name] = common.unique_filename_in()+".sql"
            macs_final = track.track( peak_list[name], chrmeta=chrlist,
                                      info={'datatype':'qualitative'},
                                      fields=['start','end','name','score'] )
            macs_final.write(gMiner.common.fusion(macs_neighb))
            macs_final.close()
            ##############################
    if peak_deconvolution:
        processed['deconv'] = {}
        merged_wig = {}
        if int(options.get('read_extension',-1))<=0:
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
                    output = mapseq.parallel_density_sql( ex, m["bam"], assembly.chromosomes,
                                                          nreads=m["stats"]["total"],
                                                          merge=-1,
                                                          convert=False,
                                                          b2w_args=b2w_args, via=via )
                    wig.append(dict((s,output+s+'.sql') for s in suffixes))
                else:
                    wig.append(m['wig'])
            if len(wig) > 1:
                merged_wig[group_name] = dict((s,common.merge_sql(ex, [x[s] for x in wig], via=via))
                                              for s in suffixes)
            else:
                merged_wig[group_name] = wig[0]
        for name in names['tests']:
            logfile.write(name+" deconvolution.\n");logfile.flush()
            if len(names['controls']) < 2:
                ctrl = (name,names['controls'][0])
                macsbed = processed['macs'][ctrl]+"_peaks.bed"
            else:
                macsbed = common.intersect_many_bed( ex, [processed['macs'][(name,x)]+"_peaks.bed"
                                                          for x in names['controls']], via=via )
            
            deconv = run_deconv( ex, merged_wig[name], macsbed, assembly.chromosomes,
                                 options['read_extension'], script_path, via=via )
            ##############################
            def _filter_deconv( stream, pval ):
                for row in stream:
                    rowpatt = re.search(r';FERR=([\d\.]+)\s',row[4])
                    if rowpatt and float(rowpatt.groups()[0]) <= pval: yield row
            ##############################
            peak_list[name] = common.unique_filename_in()+".bed"
            bedfile = track.track(peak_list[name], chrmeta=chrlist,
                                  fields=["chr","start","end","name","score"])
            trbed = track.track(deconv['bed'])
            trfields = ['chr']+trbed.fields
            bedfile.write(track.FeatureStream(
                _filter_deconv(trbed.read(fields=trfields),0.65),fields=trfields))
            bedfile.close()
            ex.add(deconv.pop('bed'), description=common.set_file_descr(name+'_peaks.sql',type='sql',step='deconvolution',group=name))
            [ex.add(v, description=common.set_file_descr(name+'_deconv.'+k,type=k,step='deconvolution',group=name))
             for k,v in deconv.iteritems()]
            processed['deconv'][name] = deconv
    for name, plist in peak_list.iteritems():
        ptrack = track.track(plist)
        annotations = track.track(assembly.sqlite_path())
        peakfile = common.unique_filename_in()
        touch(ex,peakfile)
        peakout = track.track(peakfile, format='txt', 
                              fields=['chr','start','end','name','strand',
                                      'gene','location_type','distance'])
        for chrom in assembly.chrnames:
            peakout.write(gMiner.stream.getNearestFeature(
                    ptrack.read(selection=chrom), 
                    annotations.read(selection=chrom)),mode='append')
        peakout.close()
        common.gzipfile(ex,peakfile)
        ex.add(peakfile+".gz", 
               description=common.set_file_descr(name+'_annotated_peaks.txt.gz',type='text',
                                                 step='annotation',group=name))
    if run_meme:
        from bbcflib.motif import parallel_meme
        logfile.write("Starting MEME.\n");logfile.flush()
        processed['meme'] = parallel_meme( ex, assembly,
                                           peak_list.values(), name=peak_list.keys(),
                                           meme_args=['-nmotifs','4','-revcomp'], via=via )
    return processed

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

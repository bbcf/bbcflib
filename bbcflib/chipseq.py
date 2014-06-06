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
    gl = { 'hts_chipseq': {'url': 'http://htsstation.epfl.ch/chipseq/'},
           'hts_mapseq': {'url': 'http://htsstation.epfl.ch/mapseq/'},
           'script_path': '/srv/chipseq/lib' }
    htss = frontend.Frontend( url=gl['hts_chipseq']['url'] )
    job = htss.job( hts_key )
    assembly = genrep.Assembly( assembly=assembly_id )
    with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
        job = get_bam_wig_files( ex, job, ms_limspath, gl['hts_mapseq']['url'], gl['script_path'], via=via )
        files = chipseq_workflow( ex, job, ms_files, assembly, gl['script_path'] )
    print ex.id
    allfiles = get_files( ex.id, M )
    print allfiles
"""

# Built-in modules #
import re, os, gzip, sys, time

# Internal modules #
from bbcflib import frontend, mapseq
from bbcflib.common import unique_filename_in, set_file_descr, join_pdf, merge_sql, intersect_many_bed, gzipfile
from bbcflib.track import track, FeatureStream, convert
from bbcflib.gfminer.stream import concatenate, neighborhood, getNearestFeature, score_by_feature
from bbcflib.gfminer.common import fusion, apply

# Other modules #
from bein import program
from bein.util import touch

################################################################################
# Peaks and annotation #

@program
def macs( read_length, genome_size, bamfile, ctrlbam=None, args=None ):
    """Binding for the ``macs`` peak caller.
    Takes one (optionally two) bam file(s) and the 'read_length' and 'genome_size' parameters passed to ``macs``.
    Returns the file prefix ('-n' option of ``macs``)
    """
    macs_args = ["macs14","-t",bamfile,"-f","BAM","-g",str(genome_size)]
    if isinstance(args,list): macs_args += args
    if not(ctrlbam is None): macs_args += ["-c",ctrlbam]
    if "-n" in macs_args:
        outname = macs_args[macs_args.index("-n")+1]
    else:
        outname = unique_filename_in()
        macs_args += ["-n",outname]
    if not("-s" in macs_args): macs_args += ["-s",str(read_length)]
    if not("--verbose" in macs_args): macs_args += ["--verbose","1"]
    if not("--keep-dup" in macs_args): macs_args += ["--keep-dup","all"]
    return {"arguments": macs_args, "return_value": outname}

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
        if not("-m" in macs_args): macs_args += ["-m",enrich_bounds]
        if isinstance(read_length,list): rl = read_length[i]
        for j,cam in enumerate(ctrlbam):
            m = name['controls'][j]
            nm = (n,m)
            futures[nm] = macs.nonblocking( ex, rl, genome_size, bam, cam,
                                            args=macs_args, via=via, memory=12 )
    prefixes = {}
    for n,f in futures.iteritems():
        p = f.wait()
        prefixes[n] = p
        macs_descr0 = {'step':'macs','type':'none','view':'admin','groupId':n[0][0]}
        macs_descr1 = {'step':'macs','type':'xls','groupId':n[0][0]}
        macs_descr2 = {'step':'macs','type':'bed','groupId':n[0][0],'ucsc':'1'}
        filename = "_vs_".join([x[1] for x in n if x[0]])
        touch( ex, p )
        ex.add( p, description=set_file_descr(filename,**macs_descr0),
                alias=alias )
        ex.add( p+"_peaks.xls",
                description=set_file_descr(filename+"_peaks.xls",**macs_descr1),
                associate_to_filename=p, template='%s_peaks.xls' )
        bedzip = gzip.open(p+"_peaks.bed.gz",'wb')
        bedzip.write("track name='"+filename+"_macs_peaks'\n")
        with open(p+"_peaks.bed") as bedinf:
            [bedzip.write(l) for l in bedinf]
        bedzip.close()
        ex.add( p+"_peaks.bed.gz",
                description=set_file_descr(filename+"_peaks.bed.gz",**macs_descr2),
                associate_to_filename=p, template='%s_peaks.bed.gz' )
        bedzip = gzip.open(p+"_summits.bed.gz",'wb')
        bedzip.write("track name='"+filename+"_macs_summits'\n")
        with open(p+"_summits.bed") as bedinf:
            [bedzip.write(l) for l in bedinf]
        bedzip.close()
        ex.add( p+"_summits.bed.gz",
                description=set_file_descr(filename+"_summits.bed.gz",**macs_descr2),
                associate_to_filename=p, template='%s_summits.bed.gz' )
        if n[1][0]:
            ex.add( p+"_negative_peaks.xls",
                    description=set_file_descr(filename+"_negative_peaks.xls",**macs_descr0),
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
    output = unique_filename_in()
    args = ["-p",peaks,"-f",scores_fwd,"-r",scores_rev,"-o",output,"-c",chromosome_name,
            "-l",str(chromosome_length),"-e",str(read_extension),"-z",script_path,"-s","1500"]
    return {'arguments': ["camelPeaks.py"]+args, 'return_value': None}

def run_deconv(ex, sql, peaks, chromosomes, read_extension, script_path, via = 'lsf'):
    """Runs the complete deconvolution process for a set of sql files and a bed file,
    parallelized over a set of chromosomes (a dictionary with keys 'chromosome_name'
    and values a dictionary with chromosome lengths).

    Returns a dictionary of file outputs with keys file types
    ('pdf', 'sql' and 'bed') and values the file names.
    """
    deconv_futures = {}
    stdout_files = {}
    for clen,chrom in sorted([(v['length'],k) for k,v in chromosomes.iteritems()],reverse=True):
        stdout_files[chrom] = unique_filename_in()
        deconv_futures[chrom] = camelPeaks.nonblocking( ex, sql['fwd'], sql['rev'], peaks,
                                                        chrom, clen,
                                                        read_extension, script_path,
                                                        queue="long",
                                                        via=via, stdout=stdout_files[chrom] )
        time.sleep(150) ##avoid too many processes reading same sql
    deconv_out = {}
    for c,f in deconv_futures.iteritems():
        f.wait()
        scan = False
        deconv_out[c] = []
        with open(stdout_files[c]) as sout:
            for row in sout:
                if scan:
                    if re.search(r'\*\*\*\*\*\*\*\*\*\*\*\*PARAMETERS\*\*\*\*\*\*\*\*\*\*',row):
                        row = sout.next().strip()
                        deconv_out[c].append(dict(x.split("=") for x in row.split("|")))
                        break
                    outf = row.strip()
                    if os.path.exists(outf):
                        deconv_out[c].append(outf)
                if re.search(r'\*\*\*\*\*\*\*\*\*\*\*\*OUTPUT FILES\*\*\*\*\*\*\*\*\*\*',row):
                    scan = True
    if len(deconv_out)>0:
        pdf_future = join_pdf.nonblocking( ex, [x[0] for x in deconv_out.values() if len(x)>0], via=via )
    output = unique_filename_in()
    outfiles = {}
    outfiles['peaks'] = output+"_peaks.sql"
    outfiles['profile'] = output+"_deconv.sql"
    outbed = track(outfiles['peaks'], chrmeta=chromosomes,
                   fields=["start","end","name","score"])
    info = {'datatype':'qualitative','len':'','mu':[],'lambda':[]}
    outwig = track(outfiles['profile'], chrmeta=chromosomes,
                   fields=["start","end","score"], info={'datatype':'quantitative'})
    outbed.open()
    outwig.open()
    for c,fout in deconv_out.iteritems():
        if len(fout) < 3: continue
        outbed.write(track(fout[1]).read(),chrom=c)
        outwig.write(track(fout[2]).read(),chrom=c)
        if len(fout) > 3:
            info['len'] = int(fout[3].get('len',0))
            info['mu'].append("%.1f"%float(fout[3].get('mu',0)))
            info['lambda'].append("%.1f"%float(fout[3].get('lambda',0)))
    info['mu'] = ",".join(info['mu'])
    info['lambda'] = ",".join(info['lambda'])
    outbed._fix_attributes(info)
    outbed.close()
    outwig.close()
    if len(deconv_out)>0:
        outfiles['pdf'] = pdf_future.wait()
    else:
        outfiles['pdf'] = deconv_out.values()[0][0]
    return outfiles


def parse_MACS_xls(xlsfiles):
    peak_list = [[] for _ in xlsfiles]
    def _row(x):
        return (x[0],)+tuple(int(y) for y in x[1:6])+tuple(float(y) for y in x[6:])
    for nf, xlsf in enumerate(xlsfiles):
        header = None
        with open(xlsf) as _xf:
            for row in _xf:
                if row[0] == "#": continue
                if header:
                    peak_list[nf].append(_row(row.strip().split()))
                elif row.count("\t"):
                    header = row.strip().split()
    return (header,peak_list)




################################################################################
# Workflow #

def chipseq_workflow( ex, job_or_dict, assembly, script_path='', logfile=sys.stdout, via='lsf' ):
    """Runs a chipseq workflow over bam files obtained by mapseq. Will optionally run ``macs`` and 'run_deconv'.

    :param ex: a 'bein' execution environment to run jobs in,

    :param job_or_dict: a 'Frontend' 'job' object, or a dictionary with key 'groups', 'files' and 'options' if applicable,

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
        mapseq_files = job_or_dict.files
    elif isinstance(job_or_dict,dict) and 'groups' in job_or_dict:
        if 'options' in job_or_dict:
            options = job_or_dict['options']
        groups = job_or_dict['groups']
        for gid in groups.keys():
            if not('name' in groups[gid]):
                groups[gid]['name'] = gid
        mapseq_files = job_or_dict.get('files',{})
    else:
        raise TypeError("job_or_dict must be a frontend. Job object or a dictionary with key 'groups'.")
    merge_strands = int(options.get('merge_strands',-1))
    suffixes = ["fwd","rev"]
    peak_deconvolution = options.get('peak_deconvolution',False)
    if isinstance(peak_deconvolution,basestring):
        peak_deconvolution = peak_deconvolution.lower() in ['1','true','t']
    run_meme = options.get('run_meme',False)
    if isinstance(run_meme,basestring):
        run_meme = run_meme.lower() in ['1','true','t']
    macs_args = options.get('macs_args',["--bw","200"])
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
            names['controls'].append((gid,group_name))
        else:
            tests.append(bamfile)
            names['tests'].append((gid,group_name))
            read_length.append(mapped.values()[0]['stats']['read_length'])
    genome_size = mapped.values()[0]['stats']['genome_size']
    if len(controls)<1:
        controls = [None]
        names['controls'] = [(0,None)]
    logfile.write("Starting MACS.\n");logfile.flush()
    processed = {'macs': add_macs_results( ex, read_length, genome_size,
                                           tests, ctrlbam=controls, name=names,
                                           poisson_threshold=p_thresh,
                                           macs_args=macs_args, via=via ) }
    logfile.write("Done MACS.\n");logfile.flush()
    peak_list = {}
    chrlist = assembly.chrmeta
## select only peaks with p-val <= 1e-0.6 = .25 => score = -10log10(p) >= 6
    _select = {'score':(6,sys.maxint)}
    _fields = ['chr','start','end','name','score']
    for i,name in enumerate(names['tests']):
        if len(names['controls']) < 2:
            ctrl = (name,names['controls'][0])
            macsbed = track(processed['macs'][ctrl]+"_summits.bed",
                            chrmeta=chrlist, fields=_fields).read(selection=_select)
        else:
            macsbed = concatenate([apply(track(processed['macs'][(name,x)]+"_summits.bed",
                                         chrmeta=chrlist, fields=_fields).read(selection=_select),
                                         'name', lambda __n,_n=xn: "%s:%i" %(__n,_n))
                                   for xn,x in enumerate(names['controls'])])
        ##############################
        macs_neighb = neighborhood( macsbed, before_start=150, after_end=150 )
        peak_list[name] = unique_filename_in()+".sql"
        macs_final = track( peak_list[name], chrmeta=chrlist,
                            info={'datatype':'qualitative'},
                            fields=['start','end','name','score'] )
        macs_final.write(fusion(macs_neighb),clip=True)
        macs_final.close()
        ##############################

    merged_wig = {}
    options['read_extension'] = int(options.get('read_extension') or read_length[0])
    if options['read_extension'] < 1: options['read_extension'] = read_length[0]
    make_wigs = merge_strands >= 0 or options['read_extension']>100
    if options['read_extension'] > 100: options['read_extension'] = 50
    for gid,mapped in mapseq_files.iteritems():
#            if groups[gid]['control']: continue
        group_name = groups[gid]['name']
        wig = []
        for m in mapped.values():
            if make_wigs or not('wig' in m) or len(m['wig'])<2:
                output = mapseq.parallel_density_sql( ex, m["bam"], assembly.chrmeta,
                                                      nreads=m["stats"]["total"],
                                                      merge=-1, read_extension=options['read_extension'],
                                                      convert=False,
                                                      b2w_args=b2w_args, via=via )
                wig.append(dict((s,output+s+'.sql') for s in suffixes))
            else:
                wig.append(m['wig'])
        if len(wig) > 1:
            merged_wig[group_name] = dict((s,merge_sql(ex, [x[s] for x in wig], via=via))
                                          for s in suffixes)
        else:
            merged_wig[group_name] = wig[0]

    if peak_deconvolution:
        ##############################
        def _filter_deconv( stream, pval ):
            ferr = re.compile(r';FERR=([\d\.]+)$')
            return FeatureStream( ((x[0],)+((x[2]+x[1])/2-150,(x[2]+x[1])/2+150)+x[3:] 
                                   for x in stream 
                                   if "FERR=" in x[3] and float(ferr.search(x[3]).groups()[0]) <= pval), 
                                  fields=stream.fields )
        ##############################
        processed['deconv'] = {}
        for name in names['tests']:
            logfile.write(name[1]+" deconvolution.\n");logfile.flush()
            if len(names['controls']) < 2:
                ctrl = (name,names['controls'][0])
                macsbed = processed['macs'][ctrl]+"_peaks.bed"
            else:
                macsbed = intersect_many_bed( ex, [processed['macs'][(name,x)]+"_peaks.bed"
                                                   for x in names['controls']], via=via )
            deconv = run_deconv( ex, merged_wig[name[1]], macsbed, assembly.chrmeta,
                                 options['read_extension'], script_path, via=via )
            peak_list[name] = unique_filename_in()+".bed"
            trbed = track(deconv['peaks']).read()
            with track(peak_list[name], chrmeta=chrlist, fields=trbed.fields) as bedfile:
                bedfile.write(fusion(_filter_deconv(trbed,0.65)))
            ex.add(deconv['peaks'],
                   description=set_file_descr(name[1]+'_peaks.sql', type='sql',
                                              step='deconvolution', groupId=name[0]))
            ex.add(deconv['profile'],
                   description=set_file_descr(name[1]+'_deconv.sql', type='sql',
                                              step='deconvolution',  groupId=name[0]))
            bigwig = unique_filename_in()
            convert(deconv['profile'],(bigwig,"bigWig"))
            ex.add(bigwig,
                   description=set_file_descr(name[1]+'_deconv.bw', type='bigWig',
                                              ucsc='1', step='deconvolution',
                                              groupId=name[0]))
            ex.add(deconv['pdf'],
                   description=set_file_descr(name[1]+'_deconv.pdf', type='pdf',
                                              step='deconvolution', groupId=name[0]))
            processed['deconv'][name] = deconv

    ##############################
    def _join_macs( stream, xlsl, _f ):
        def _macs_row(_n):
            if len(xlsl) == 1:
                return xlsl[0][(int(_n.split(";")[0][13:]) if _n[:3] == "ID=" else int(_n[10:]))-1][1:]
            else:
                nb = _n.split(";")[0][13:] if _n[:3] == "ID=" else _n[10:]
                nb = nb.split(":")
                return xlsl[int(nb[1])][int(nb[0])-1][1:]
        return FeatureStream( (_p+_macs_row(_p[3]) for _p in stream), fields=_f )
    ##############################
    peakfile_list = []
    for name, plist in peak_list.iteritems():
        ptrack = track(plist,chrmeta=chrlist,fields=["chr","start","end","name","score"])
        peakfile = unique_filename_in()
        xlsh, xlsl = parse_MACS_xls([processed['macs'][(name,_c)]+"_peaks.xls" for _c in names['controls']])
        try:
###### if assembly doesn't have annotations, we skip the "getNearestFeature" but still go through "_join_macs"
            assembly.gene_track()
            _fields = ['chr','start','end','name','score','gene','location_type','distance']\
                +["MACS_%s"%h for h in xlsh[1:5]]+xlsh[5:]
            peakout = track(peakfile, format='txt', chrmeta=chrlist, fields=_fields)
            peakout.make_header("#"+"\t".join(['chromosome','start','end','info','peak_height','gene(s)','location_type','distance']+_fields[8:]))
        except ValueError:
            _fields = ['chr','start','end','name','score']+["MACS_%s"%h for h in xlsh[1:5]]+xlsh[5:]
            peakout = track(peakfile, format='txt', chrmeta=chrlist, fields=_fields)
            peakout.make_header("#"+"\t".join(['chromosome','start','end','info','peak_height']+_fields[8:]))
        for chrom in assembly.chrnames:
            try:
                _feat = assembly.gene_track(chrom)
                peakout.write(_join_macs(getNearestFeature(ptrack.read(selection=chrom),_feat),
                                         xlsl, _fields), mode='append')
            except ValueError:
                peakout.write(_join_macs(ptrack.read(selection=chrom),xlsl, _fields), mode='append')
        peakout.close()
        gzipfile(ex,peakfile)
        peakfile_list.append(track(peakfile+".gz", format='txt', fields=_fields))
        ex.add(peakfile+".gz",
               description=set_file_descr(name[1]+'_annotated_peaks.txt.gz',type='text',
                                          step='annotation',groupId=name[0]))
    stracks = [track(wig,info={'name':name+"_"+st}) 
               for name,wigdict in merged_wig.iteritems() for st,wig in wigdict.iteritems()]
    tablefile = unique_filename_in()
    with open(tablefile,"w") as _tf:
        _tf.write("\t".join(['chr','start','end','name']+[s.name for s in stracks])+"\n")
#### need to do something about peak origin (split names, write to separate columns?)
    for chrom in assembly.chrnames:
        peak_list = [pt.read(chrom,fields=['chr','start','end','name']) for pt in peakfile_list]
        features = fusion(concatenate(peak_list, fields=['chr','start','end','name'], 
                                      remove_duplicates=True, group_by=['chr','start','end']))
        sread = [sig.read(chrom) for sig in stracks]
        quantifs = score_by_feature(sread, features, method='sum')
        with open(tablefile,"a") as _tf:
            for t in quantifs:
                _tf.write("\t".join(str(tt) for tt in t)+"\n")
    gzipfile(ex,tablefile)
    ex.add(tablefile+".gz",
           description=set_file_descr('Combined_peak_quantifications.txt.gz',type='text',
                                      step='summary'))

    if run_meme:
        from bbcflib.motif import parallel_meme
        logfile.write("Starting MEME.\n");logfile.flush()
        processed['meme'] = parallel_meme( ex, assembly,
                                           peak_list.values(), name=peak_list.keys(),
                                           chip=True, meme_args=['-meme-nmotifs','4','-meme-mod','zoops'],
                                           via=via )
    return processed

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

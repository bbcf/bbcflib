"""
========================
Module: bbcflib.dnaseseq
========================
"""
import os, sys, gzip, tarfile
from bbcflib.gfminer.stream import neighborhood, segment_features, score_by_feature
from bbcflib.gfminer.common import fusion, sorted_stream
from bbcflib.gfminer.figure import lineplot
from bbcflib.common import set_file_descr, unique_filename_in, intersect_many_bed, cat
from bbcflib.chipseq import add_macs_results
from bbcflib.mapseq import merge_bam, index_bam
from bbcflib.track import track, convert, FeatureStream
from bbcflib.genrep import GenRep
from bbcflib.motif import save_motif_profile
from bein import program
from bein.util import touch
from numpy import vstack, zeros

_gnrp = GenRep()
_well_flank = 300
_plot_flank = (50,50)

def macs_bedfiles( ex, chrmeta, tests, controls, names, macs_args, via, logfile ):
    missing_beds = [k for k,t in enumerate(tests) if not t[0]]
    if not missing_beds: return tests
    genome_size = sum([x['length'] for x in chrmeta.values()])
    logfile.write("Running MACS.\n");logfile.flush()
    _tts = [tests[k][1] for k in missing_beds]
    _nms = {'tests': [names['tests'][k] for k in missing_beds],
            'controls': names['controls']}
    macsout = add_macs_results( ex, 0, genome_size, _tts, 
                                ctrlbam=controls, name=_nms,
                                poisson_threshold={},
                                macs_args=macs_args, via=via )
    logfile.write("Done MACS.\n");logfile.flush()
    for nbam,bed_bam in enumerate(tests):
        name = names['tests'][nbam]
        if not bed_bam[0]:
            if len(names['controls']) < 2:
                bed_bam = (macsout[(name,names['controls'][0])]+"_summits.bed", bed_bam[1])
            else:
                _beds = [macsout[(name,x)]+"_summits.bed" for x in names['controls']]
                bed_bam = (intersect_many_bed( ex, _beds, via=via ), bed_bam[1])
        tests[nbam] = bed_bam
    return tests

@program
def wellington( bed, bam, output=None, options=[] ):
    if output is None: output = unique_filename_in()
    outdir = unique_filename_in()
    os.mkdir(outdir)
    args = ["wellington_footprints.py","-o",output]
    return {'arguments': args+options+[bed,bam,outdir], 'return_value': (outdir,output)}

def save_wellington( ex, wellout ):
    bedlist = {}
    for name, wlist in wellout.iteritems():
        wellall = unique_filename_in()
#### Dummy file
        touch( ex, wellall )
        ex.add(wellall,
               description=set_file_descr(name[1]+'_wellington_files', type='none', view='admin',
                                          step='footprints', groupId=name[0]))
#### BED at FDR 1%
        bedlist[name[0]] = wellall+"_WellingtonFootprintsFDR01.bed.gz"
        bedzip = gzip.open(bedlist[name[0]],'wb')
        bedzip.write("track name='"+name[1]+"_WellingtonFootprints_FDR_0.01'\n")
        for x in wlist:
            with open(os.path.join(*x)+".WellingtonFootprints.FDR.0.01.bed") as _bed:
                [bedzip.write(l) for l in _bed]
        bedzip.close()
        ex.add(wellall+"_WellingtonFootprintsFDR01.bed.gz",
               description=set_file_descr(name[1]+'_WellingtonFootprintsFDR01.bed.gz', 
                                          type='bed', ucsc='1', step='footprints', groupId=name[0]),
               associate_to_filename=wellall, template='%s_WellingtonFootprintsFDR01.bed.gz')
#### BED at p-values [...]
        bedzip = gzip.open(wellall+"_WellingtonFootprintsPvalCutoffs.bed.gz",'wb')
        for bfile in os.listdir(os.path.join(wlist[0][0],"p_value_cutoffs")):
            cut = os.path.splitext(bfile[:-4])[1][1:] #between . ([1:]) and .bed ([:-4])
            bedzip.write("track name='"+name[1]+"_WellingtonFootprints_Pval_%s'\n" %cut)
            for wdir,wpref in wlist:
                _bedpath = os.path.join(wdir,"p_value_cutoffs",wpref+".WellingtonFootprints."+cut+".bed")
                with open(_bedpath) as _bed:
                    [bedzip.write(l) for l in _bed]
        bedzip.close()
        ex.add(wellall+"_WellingtonFootprintsPvalCutoffs.bed.gz",
               description=set_file_descr(name[1]+'_WellingtonFootprintsPvalCutoffs.bed.gz', 
                                          type='bed', ucsc='1', step='footprints', groupId=name[0]),
               associate_to_filename=wellall, template='%s_WellingtonFootprintsPvalCutoffs.bed.gz')
#### WIG
        cat([os.path.join(*x)+".WellingtonFootprints.wig" for x in wlist], 
            wellall+"_WellingtonFootprints.wig")
        convert(wellall+"_WellingtonFootprints.wig", wellall+"_WellingtonFootprints.bw", 
                chrmeta=assembly.chrmeta)
        ex.add(wellall+"_WellingtonFootprints.bw",
               description=set_file_descr(name[1]+'_WellingtonFootprints.bw', 
                                          type='bigWig', ucsc='1', step='footprints', groupId=name[0]),
               associate_to_filename=wellall, template='%s_WellingtonFootprints.bw')
    return bedlist


def run_wellington( ex, tests, chrnames, via, logfile ):
    futures = {}
    logfile.write("Running Wellington:\n");logfile.flush()
    wellout = {}
    for nbam,bed_bam in enumerate(tests):
        name = names['tests'][nbam]
        wellout[name] = []
        tbed = track(bed_bam[0])
        for chrom in chrnames:
            _chrombed = unique_filename_in()
            with track(_chrombed,format="bed",fields=tbed.fields) as _tt:
                _neighb = neighborhood( tbed.read(chrom), before_start=_well_flank, after_end=_well_flank )
                _tt.write(fusion(_neighb),clip=True)
            if os.path.getsize(_chrombed) > 0:
                futures[(chrom,name)] = wellington.nonblocking(ex, _chrombed, bed_bam[1], via=via, memory=8)
    for chro_name, _fut in futures.iteritems():
        chrom, name = chro_name
        logfile.write(name[1]+" "+chrom+", ");logfile.flush()
        wellout[name].append(_fut.wait())
    logfile.write("\n");logfile.flush()
    bedlist = save_wellington(ex, wellout)
    return bedlist

def motif_scan( ex, bedlist, assembly, groups, via, logfile ):
    logfile.write("Scanning motifs\n");logfile.flush()
    motifbeds = {}
    for gid,bedfile in bedlist.iteritems():
        logfile.write("\n%i: "%gid);logfile.flush()
        group = groups[gid]
        motifs = {}
        for mot in group.get('motif',[]):
            if os.path.exists(mot):
                mname = os.path.basename(os.path.splitext(mot)[0])
                motifs[mname] = mot
            elif os.path.exists(os.path.join(supdir,mot)):
                mname = os.path.basename(os.path.splitext(mot)[0])
                motifs[mname] = os.path.join(supdir,mot)
            else:
                _gnid, mname = mot.split(' ')
                motifs[mname] = _gnrp.get_motif_PWM(int(_gnid), mname, output=unique_filename_in())
            logfile.write(mname+", ");logfile.flush()
        _descr = set_file_descr(group['name']+'_motifs.bed',
                                type='bed', ucsc='1', step='motifs', groupId=gid)
        _out = unique_filename_in()
        motifbeds[gid] = save_motif_profile( ex, motifs, assembly, bedfile,
                                             keep_max_only=True, output=_out,
                                             description=_descr, via=via )
    return motifbeds

def plot_footprint_profile( ex, bedlist, siglist, chrnames, groups, logfile ):
    files = dict((gid,{'pdf':"",'mat':[]}) for gid in bedlist.keys())
    logfile.write("Plotting footprints:\n");logfile.flush()
    for gid, motifbed in bedlist.iteritems():
        signals = [track(sig) for sig in siglist[gid]]
        snames = [sig.name for sig in signals]
        tmotif = track(motifbed)
        data = {}
        numregs = {}
        for chrom in chrnames:
            fread = {}
            for r in tmotif.read(chrom):
                r2 = r[3].split(":")
                key = (r2[0],len(r2[1]))
                if key in fread: fread[key].append(r[1:3])
                else: fread[key] = [r[1:3]]
            for motif, regs in fread.iteritems():
                if motif not in data:
                    data[motif] = zeros(shape=(motif[1]+2*_plot_flank[1],len(signals)))
                    numregs[motif] = 0
                numregs[motif] += len(regs)
                tFeat = sorted_stream(segment_features(FeatureStream(regs,fields=['start','end']),
                                                       nbins=motif[1],upstream=_plot_flank,downstream=_plot_flank))
                for t in score_by_feature([s.read(chrom) for s in signals], tFeat): 
                    data[motif][t[2]] += t[3:]
        files[gid]['pdf'] = unique_filename_in()
        new = True
        last = len(data)
        for motif, dat in data.iteritems():
            last -= 1
            mname, nbins = motif
            dat /= float(numregs[motif])
            X = ['']*dat.shape[0]
            X[0] = "-%i"%_plot_flank[1]
            X[-1] = "%i"%_plot_flank[1]
            for k in range(nbins): X[k+_plot_flank[1]] = str(k+1)
####### Could do a heatmap (sort by intensity)...
            lineplot(X, [dat[:, n] for n in range(dat.shape[-1])],
                     output=files[gid]['pdf'], new=new, last=(last>0), 
                     legend=snames, main=mname)
            new = False
            _datf = unique_filename_in()
            with open(_datf,"w") as dff:
                dff.write("\t".join([""]+[str(x) for x in X])+"\n")
                for n,sn in enumerate(snames):
                    dff.write("\t".join([sn]+[str(x) for x in dat[:, n]])+"\n")
            files[gid]['mat'].append((mname,_datf))
    return files


###############################################################
def dnaseseq_workflow( ex, job, assembly, logfile=sys.stdout, via='lsf' ):
    """
    Main function
    """
    tests = []
    controls = []
    names = {'tests': [], 'controls': []}
    supdir = os.path.split(ex.remote_working_directory)[0]
    for gid,mapped in job.files.iteritems():
        group_name = job.groups[gid]['name']
        if not isinstance(mapped,dict):
            raise TypeError("Mapseq_files values must be dictionaries with keys *run_ids* or 'bam'.")
        if 'bam' in mapped: mapped = {'_': mapped}
        if len(mapped)>1:
            bamfile = merge_bam(ex, [m['bam'] for m in mapped.values()])
            index = index_bam(ex, bamfile)
        else:
            bamfile = mapped.values()[0]['bam']
        if job.groups[gid]['control']:
            controls.append(bamfile)
            names['controls'].append((gid,group_name))
        else:
            if os.path.exists(job.groups[gid].get('bedfile','null')):
                bedfile = job.groups[gid]['bedfile']
            elif os.path.exists(os.path.join(supdir,job.groups[gid].get('bedfile','null'))):
                bedfile = os.path.join(supdir,job.groups[gid]['bedfile'])
            else:
                bedfile = None
            tests.append((bedfile,bamfile))
            names['tests'].append((gid,group_name))
    if len(controls)<1:
        controls = [None]
        names['controls'] = [(0,None)]
    tests = macs_bedfiles( ex, assembly.chrmeta, tests, controls, names, 
                           job.options.get('macs_args',["--keep-dup","10"]), via, logfile )
    bedlist = run_wellington(ex, tests, assembly.chrnames, via, logfile)

######################### Motif scanning / plotting
    if any([gr.get('motif') or None for gr in jobs.groups]):
        motifbeds = motif_scan( ex, bedlist, job.groups, via, logfile )
        siglist = dict((gid,[]) for gid,_ in names['tests'])
        for gid,mapped in job.files.iteritems():
            wig = []
            suffixes = ["fwd","rev"]
            merge_strands = int(job.options.get('merge_strands',-1))
            read_extension = int(options.get('read_extension') or 50)
            make_wigs = merge_strands >= 0 or read_extension != 1
            for m in mapped.values():
                if make_wigs or not('wig' in m) or len(m['wig'])<2:
                    output = mapseq.parallel_density_sql( ex, m["bam"], assembly.chrmeta,
                                                          nreads=m["stats"]["total"],
                                                          merge=-1, read_extension=1,
                                                          convert=False,
                                                          b2w_args=[], via=via )
                    wig.append(dict((s,output+s+'.sql') for s in suffixes))
                else:
                    wig.append(m['wig'])
            if len(wig) > 1:
                wig[0] = dict((s,merge_sql(ex, [x[s] for x in wig], via=via)) for s in suffixes)
            if job.groups[gid]['control']:
                for _g in siglist.keys():
                    siglist[_g].extend(wig[0].values())
            else:
                siglist[gid].extend(wig[0].values())
        plot_files = plot_footprint_profile( ex, motifbeds, siglist, assembly.chrnames, jobs.groups, logfile )
        for gid, flist in plot_files.iteritems():
            gname = job.groups[gid]['name']
            plotall = unique_filename_in()
            touch( ex, plotall )
            ex.add(plotall,
                   description=set_file_descr(gname+'_footprints_plots', type='none', view='admin',
                                              step='motifs', groupId=gid))
            ex.add(flist['pdf'],
                   description=set_file_descr(gname+'_footprints_plots.pdf', 
                                              type='pdf', step='motifs', groupId=gid),
                   associate_to_filename=plotall, template='%s.pdf')
            tarname = unique_filename_in()
            tarfh = tarfile.open(tarname, "w:gz")
            for mname,matf in flist['mat']:
                tarfh.add(matf, arcname="%s_%s.txt" % (gname,mname))
            tarfh.close()
            ex.add( tarname,
                   description=set_file_descr(gname+'_footprints_plots.tar.gz',
                                              type='tar', step='motifs', groupId=gid),
                   associate_to_filename=plotall, template='%s.tar.gz')
    logfile.write("\nDone.\n ");logfile.flush()
    return 0

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

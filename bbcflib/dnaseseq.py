"""
========================
Module: bbcflib.dnaseseq
========================
"""
import os, sys, gzip
from bbcflib.gfminer.stream import neighborhood
from bbcflib.gfminer.common import fusion
from bbcflib.common import set_file_descr, unique_filename_in, intersect_many_bed, cat
from bbcflib.chipseq import add_macs_results
from bbcflib.mapseq import merge_bam, index_bam
from bbcflib.track import track, convert
from bbcflib.genrep import GenRep
from bbcflib.motif import save_motif_profile
from bein import program
from bein.util import touch

_gnrp = GenRep()

@program
def wellington( bed, bam, output=None, options=[] ):
    if output is None: output = unique_filename_in()
    outdir = unique_filename_in()
    os.mkdir(outdir)
    args = ["wellington_footprints.py","-o",output]
    return {'arguments': args+options+[bed,bam,outdir],
            'return_value': (outdir,output)}


###############################################################
def dnaseseq_workflow( ex, job, assembly, logfile=sys.stdout, via='lsf' ):
    """
    Main function
    """
    macs_args = job.options.get('macs_args',["--keep-dup","10"])
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
    missing_beds = [k for k,t in enumerate(tests) if not t[0]]
    if missing_beds:
        genome_size = sum([x['length'] for x in assembly.chrmeta.values()])
        logfile.write("Running MACS.\n");logfile.flush()
        _tts = [tests[k][1] for k in missing_beds]
        _nms = {'tests': [names['tests'][k] for k in missing_beds],
                'controls': names['controls']}
        macsout = add_macs_results( ex, 0, genome_size, _tts, 
                                    ctrlbam=controls, name=_nms,
                                    poisson_threshold={},
                                    macs_args=macs_args, via=via )
        logfile.write("Done MACS.\n");logfile.flush()
    futures = {}
    logfile.write("Running Wellington:\n");logfile.flush()
    wellout = {}
    for nbam,bed_bam in enumerate(tests):
        name = names['tests'][nbam]
        wellout[name] = []
        if not bed_bam[0]:
            if len(names['controls']) < 2:
                bed_bam = (macsout[(name,names['controls'][0])]+"_summits.bed", bed_bam[1])
            else:
                _beds = [macsout[(name,x)]+"_summits.bed" for x in names['controls']]
                bed_bam = (intersect_many_bed( ex, _beds, via=via ), bed_bam[1])
        tbed = track(bed_bam[0])
        for chrom in assembly.chrnames:
            _chrombed = unique_filename_in()
            with track(_chrombed,format="bed",fields=tbed.fields) as _tt:
                _neighb = neighborhood( tbed.read(chrom), before_start=300, after_end=300 )
                _tt.write(fusion(_neighb),clip=True)
            if os.path.getsize(_chrombed) > 0:
                futures[(chrom,name)] = wellington.nonblocking(ex, _chrombed, bed_bam[1], via=via, memory=8)

    for chro_name, _fut in futures.iteritems():
        chrom, name = chro_name
        logfile.write(name[1]+" "+chrom+", ");logfile.flush()
        wellout[name].append(_fut.wait())
    logfile.write("\n");logfile.flush()

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
######################### Motif scanning
    if job.options.get('scan_motifs',False):
        logfile.write("Scanning motifs\n");logfile.flush()
        motifbeds = {}
        for gid,bedfile in bedlist.iteritems():
            logfile.write("\n%i: "%gid);logfile.flush()
            group = job.groups[gid]
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
            _descr = set_file_descr(group['name']+'motifs.bed.gz',
                                    type='bed', ucsc='1', step='motifs', groupId=gid)
            _out = unique_filename_in()
            motifbeds[gid] = save_motif_profile( ex, motifs, assembly, bedfile,
                                                 keep_max_only=True, output=_out,
                                                 description=_descr, via=via )
    logfile.write("\nDone.\n ");logfile.flush()
    return 0

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

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
from bein import program
from bein.util import touch

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
    macs_args = job.options.get('macs_args',[])
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
            if os.path.exists(jobs.groups[gid].get('bedfile')):
                bedfile = jobs.groups[gid]['bedfile']
            elif os.path.exists(os.path.join(supdir,jobs.groups[gid].get('bedfile'))):
                bedfile = os.path.join(supdir,jobs.groups[gid]['bedfile'])
            else:
                bedfile = None
            tests.append((bedfile,bamfile))
            names['tests'].append((gid,group_name))
    if len(controls)<1:
        controls = [None]
        names['controls'] = [(0,None)]
    if any([not t[0] for t in tests]):
        genome_size = sum([x['length'] for x in assembly.chrmeta.values()])
        logfile.write("Running MACS.\n");logfile.flush()
        macsout = add_macs_results( ex, 0, genome_size, [t[1] for t in tests if not t[0]], 
                                    ctrlbam=controls,
                                    name=[n for k,n in enumerate(names) if not tests[k][0]],
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
                bed_bam[0] = macsout[(name,names['controls'][0])]+"_summits.bed"
            else:
                _beds = [macsout[(name,x)]+"_summits.bed" for x in names['controls']]
                bed_bam[0] = intersect_many_bed( ex, _beds, via=via )
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

    for name, wlist in wellout.iteritems():
        wellall = unique_filename_in()
#### Dummy file
        touch( ex, wellall )
        ex.add(wellall,
               description=set_file_descr(name[1]+'_wellington_files', type='none', view='admin',
                                          step='footprints', groupId=name[0]))
#### BED at FDR 1%
        bedzip = gzip.open(wellall+"_WellingtonFootprintsFDR01.bed.gz",'wb')
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
        
    logfile.write("\nDone.\n ");logfile.flush()
    return 0

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

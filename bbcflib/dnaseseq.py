"""
========================
Module: bbcflib.dnaseseq
========================
"""
import os, sys
from bbcflib.common import set_file_descr, unique_filename_in, intersect_many_bed
from bbcflib.chipseq import add_macs_results
from bbcflib.mapseq import merge_bam
from bbcflib.track import track, convert
from bein import program
from bein.util import touch

@program
def wellignton( bed, bam, output=None, options=[] ):
    if output is None: output = unique_filename_in()
    outdir = unique_filename_in()
    os.mkdir(outdir)
    args = ["wellington_footprints.py","-o",output]
    options += [bed,bam,outdir]
    return {'arguments': args+options, 'return_value': os.path.join(outdir,output) }


###############################################################
def dnaseseq_workflow( ex, job, assembly, bedfile=None, logfile=sys.stdout, via='lsf' ):
    """
    Main function
    """
    if bedfile is None:
        macs_args = job.options.get('macs_args',[])
        tests = []
        controls = []
        names = {'tests': [], 'controls': []}
        for gid,mapped in job.files.iteritems():
            group_name = job.groups[gid]['name']
            if not(isinstance(mapped,dict)):
                raise TypeError("Mapseq_files values must be dictionaries with keys *run_ids* or 'bam'.")
            if 'bam' in mapped: mapped = {'_': mapped}
            if len(mapped)>1:
                bamfile = merge_bam(ex, [m['bam'] for m in mapped.values()])
            else:
                bamfile = mapped.values()[0]['bam']
            if job.groups[gid]['control']:
                controls.append(bamfile)
                names['controls'].append((gid,group_name))
            else:
                tests.append(bamfile)
                names['tests'].append((gid,group_name))
        if len(controls)<1:
            controls = [None]
            names['controls'] = [(0,None)]
        genome_size = sum([x['length'] for x in assembly.chrmeta.values()])
        logfile.write("Running MACS.\n");logfile.flush()
        macsout = add_macs_results( ex, 0, genome_size,
                                    tests, ctrlbam=controls, name=names,
                                    poisson_threshold={},
                                    macs_args=macs_args, via=via )
        logfile.write("Done MACS.\n");logfile.flush()
    futures = {}
    logfile.write("Running Wellington:\n");logfile.flush()
    for nbam,bam in enumerate(tests):
        name = names['tests'][nbam]
        if len(names['controls']) < 2:
            macsbed = macsout[(name,names['controls'][0])]+"_peaks.bed"
        else:
            _beds = [macsout[(name,x)]+"_peaks.bed" for x in names['controls']]
            macsbed = intersect_many_bed( ex, _beds, via=via )
        futures[name] = wellington.nonblocking(ex, macsbed, bam, via=via)

    for name, _fut in futures.iteritems():
        logfile.write(name[1]+", ");logfile.flush()
        wellout = _fut.wait()
        touch( ex, wellout )
        ex.add(wellout,
               description=set_file_descr(name[1]+'_wellington_files', type='none', view='admin',
                                          step='footprints', groupId=name[0]))
        ex.add(wellout+".WellingtonFootprints.FDR.0.01.bed",
               description=set_file_descr(name[1]+'.WellingtonFootprints.FDR.0.01.bed', 
                                          type='bed', ucsc='1', step='footprints', groupId=name[0]),
               associate_to_filename=wellout, template='%s.WellingtonFootprints.FDR.0.01.bed')
        twig = track(wellout+".WellingtonFootprints.wig", chrmeta=assembly.chrmeta)
        tbw = track(wellout+".WellingtonFootprints.bigWig", chrmeta=assembly.chrmeta)
        convert(twig,tbw)
        ex.add(wellout+".WellingtonFootprints.bigWig",
               description=set_file_descr(name[1]+'.WellingtonFootprints.bigWig', 
                                          type='bigWig', ucsc='1', step='footprints', groupId=name[0]),
               associate_to_filename=wellout, template='%s.WellingtonFootprints.bigWig')
        
    logfile.write("\nDone.\n ");logfile.flush()
    return 0

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

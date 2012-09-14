"""
=========================
Module: bbcflib.junctions
=========================

Methods of the bbcflib's `junctions` worflow. The main function is ``junctions_workflow()``,
and is usually called by bbcfutils' ``run_junctions.py``, e.g. command-line:

``python run_junctions.py -v local -c gapkowt.config -d junctions``

Uses SOAPsplice ...
"""

# Built-in modules #
import os, pysam

# Internal modules #
from bbcflib.common import set_file_descr, unique_filename_in
from bbcflib import mapseq, genrep
from bbcflib.bFlatMajor.common import cobble, sorted_stream
from bbcflib.btrack import track, FeatureStream, convert
from bein import program


@program
def soapsplice(unmapped_R1, unmapped_R2, index, output=None, **options):
    """Binds 'soapsplice'. Returns a text file containing the list of junctions.

    :param unmapped_R1: (str) path to the fastq file containing the 'left' reads.
    :param unmapped_R2: (str) path to the fastq file containing the 'right' reads.
    :param index: (str) path to the SOAPsplice index.
    :param options: (dict) SOAPsplice options, given as {opt: value}.
    :rtype: str
    """
    args = ['soapsplice','-d',index,'-1',unmapped_R1,'-2',unmapped_R2]
    opts = []
    for k,v in options.iteritems(): opts.extend([str(k),v])
    if not output: output = unique_filename_in()
    return {"arguments": args, "return_value": output}
    #"SOAPsplice-v1.9/bin/soapsplice -d $index -1 $path1$r1 -2 $path2$r2 -I $insertsize -o $outpath$out -f 2 "


def junctions_workflow(ex, job, bam_files, ref_genome, index, via="lsf"):
    """
    Main function of the workflow.

    :rtype: None
    :param ex: the bein's execution Id.
    :param job: a Job object (or a dictionary of the same form) as returned from HTSStation's frontend.
    :param bam_files: a complicated dictionary such as returned by mapseq.get_bam_wig_files.
    :param via: 'local' or 'lsf'.
    """
    assembly = genrep.Assembly(assembly=job.assembly_id,intype=2)
    ref_genome = untar_genome_fasta(assembly, ref_genome, convert=True) # cf snp.py

    group_names={}
    group_ids={}
    for gid,group in job.groups.iteritems():
        gname = str(group['name'])
        group_names[gid] = gname
        group_ids[gname] = gid

    unmapped_fastq = {}
    for gid, group in job.groups.iteritems():
        for rid, run in group['runs'].iteritems():
            # Define pairs of fastq files ...
            unmapped_fastq[gid] = unmapped_fastq.setdefault(gid,[]).append(bam_files[gid][rid].get('unmapped_fastq'))
        if unmapped_fastq[gid]:
            R1 = None
            R2 = None
        junc_file = soapsplice(R1,R2,index)


    return know_junctions, new_junctions


#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#

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
import os, itertools

# Internal modules #
from bbcflib.common import set_file_descr, unique_filename_in, cat
from bbcflib import genrep
from bbcflib.bFlatMajor.common import map_chromosomes, duplicate
from bbcflib.btrack import track, convert
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

def convert_junc_file(filename, assembly):
    """Convert a .junc SOAPsplice output file to sql format,
    and return an iterator on the junction locations.

    :param filename: (str) name of the .junc file to convert.
    :param assembly: genrep.Assembly object.
    :rtype: str
    """
    t = track(filename, fields=['chr','start','end','strand','score'], chrmeta=assembly.chrmeta)
    stream = t.read()
    # Translate chromosome names
    s1 = map_chromosomes(stream, assembly.chromosomes)
    # Add junction IDs
    s2 = duplicate(s1,'strand','name')
    C = itertools.count()
    s3 = apply(s2,'name', lambda x: 'junction'+str(C.next()))
    # Convert to sql and bed formats
    outfile = unique_filename_in()
    sql = outfile + '.sql'
    bed = outfile + '.bed'
    out = track(sql, fields=t.fields, chrmeta=assembly.chrmeta)
    out.write(s1)
    convert(sql,bed)
    return sql, bed

def discovery(juncfile, assembly):
    """Take the intersection between splicing regions found by SOAPsplice, and known exons annotation."""
    known_junctions = {}
    new_junctions = {}
    junc_track = track(juncfile, fields=['chr','start','end','strand','score'], chrmeta=assembly.chrmeta)
    exon_track = assembly.exon_track()
    e_chr,e_start,e_end,exon,e_strand,phase = exon_track.next()
    for chrom in assembly.chrmeta:
        for j in junc_track.read(chrom):
            j_chr,j_start,j_end,name,score,j_strand = j
            while e_end < j_start:
                e_chr,e_start,e_end,exon,e_strand,phase = exon_track.next()
            if e_end == j_start:
                known_junctions[name] = [exon,None]
            if e_start == j_end and name in known_junctions:
                known_junctions[name][1] = exon

    # exon_mapping is a dictionary ``{exon_id: ([transcript_id's],gene_id,gene_name,start,end,chromosome)}``
    #exon_mapping = assembly.get_exon_mapping()
    #for e,v in exon_mapping.iteritems():
    #    start,end = (v[2],v[3])

    return known_junctions, new_junctions


def junctions_workflow(ex, job, bam_files, index, via="lsf"):
    """
    Main function of the workflow.

    :rtype: None
    :param ex: the bein's execution Id.
    :param job: a Job object (or a dictionary of the same form) as returned from HTSStation's frontend.
    :param bam_files: a complicated dictionary such as returned by mapseq.get_bam_wig_files.
    :param via: 'local' or 'lsf'.
    """
    assembly = genrep.Assembly(assembly=job.assembly_id,intype=2)

    unmapped_fastq = {}
    for gid, group in job.groups.iteritems():
        unmapped_fastq[gid] = []
        for rid, run in group['runs'].iteritems():
            # Define pairs of fastq files ...
            unmapped = bam_files[gid][rid].get('unmapped_fastq')
            assert unmapped, "No unmapped reads found."
            assert os.path.exists(unmapped+'_1.gz'), "Pair-end reads required."
            unmapped_fastq[gid].append([unmapped+'_1.gz', unmapped+'_2.gz'])
        R1 = cat(zip(*unmapped_fastq[gid])[0])
        R2 = cat(zip(*unmapped_fastq[gid])[1])
        junc_file = soapsplice(R1,R2,index)
        # Translate to more usual formats
        sql,bed = convert_junc_file(junc_file,assembly)
        sql_descr = set_file_descr('junctions_%s.sql' % group['name'], \
                                   group=group['name'],type='sql',step='1',gdv=1)
        bed_descr = set_file_descr('junctions_%s.bed' % group['name'], \
                                   group=group['name'],type='bed',step='1',ucsc=1)
        ex.add(sql, description=sql_descr)
        ex.add(bed, description=bed_descr)



#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#

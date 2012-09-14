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


def junctions_workflow(ex, job, bam_files, via="lsf"):
    """
    Main function of the workflow.

    :rtype: None
    :param ex: the bein's execution Id.
    :param job: a Job object (or a dictionary of the same form) as returned from HTSStation's frontend.
    :param bam_files: a complicated dictionary such as returned by mapseq.get_bam_wig_files.
    :param via: 'local' or 'lsf'.
    """
    group_names={}; group_ids={}; conditions=[]
    assembly = genrep.Assembly(assembly=job.assembly_id,intype=2)
    groups = job.groups
    for gid,group in groups.iteritems():
        gname = str(group['name'])
        group_names[gid] = gname
        group_ids[gname] = gid

    """ Define conditions """
    for gid,files in bam_files.iteritems():
        k = 0
        for rid,f in files.iteritems():
            k+=1
            cond = group_names[gid]+'.'+str(k)
            conditions.append(cond)
    ncond = len(conditions)

    #Exon labels: ('exonID|geneID|start|end|strand|type', length)
    exons = fetch_labels(bam_files[groups.keys()[0]][groups.values()[0]['runs'].keys()[0]]['bam'])

    """ Extract information from bam headers """
    exonsID=[]; genesID=[]; genesName=[]; starts=[]; ends=[]; strands=[]
    exon_lengths={}; exon_to_gene={}; badexons=[]
    for e in exons:
        length = e[1]
        E = e[0].split('|')
        exon = E[0]; gene = E[1]; start = int(E[2]); end = int(E[3]); strand = int(E[4])
        if end-start>1:
            starts.append(start)
            ends.append(end)
            strands.append(strand)
            exon_lengths[exon] = float(length)
            exonsID.append(exon)
            genesID.append(gene)
            exon_to_gene[exon] = gene
        else: badexons.append(e)
    [exons.remove(e) for e in badexons]

    print "Load mappings"
    """ [0] gene_mapping is a dict ``{gene ID: (gene name,start,end,length,strand,chromosome)}``
        [1] transcript_mapping is a dictionary ``{transcript ID: (gene ID,start,end,length,strand,chromosome)}``
        [2] exon_mapping is a dictionary ``{exon ID: ([transcript IDs],gene ID,start,end,strand,chromosome)}``
        [3] trans_in_gene is a dict ``{gene ID: [IDs of the transcripts it contains]}``
        [4] exons_in_trans is a dict ``{transcript ID: [IDs of the exons it contains]}`` """
    (gene_mapping, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans) = fetch_mappings(assembly)

    """ Map remaining reads to transcriptome """
    additionals = {}; unmapped_bam = {}; unmapped_fastq = {}
    refseq_path = assembly.index_path
    for gid, group in job.groups.iteritems():
        k = 0
        for rid, run in group['runs'].iteritems():
            k +=1
            cond = group_names[gid]+'.'+str(k)
            unmapped_fastq[cond] = bam_files[gid][rid].get('unmapped_fastq')
            if unmapped_fastq[cond] and os.path.exists(refseq_path+".1.ebwt"):
                unmapped_bam[cond] = mapseq.map_reads(ex, unmapped_fastq[cond], {}, refseq_path, \
                      remove_pcr_duplicates=False, bwt_args=[], via=via)['bam']
                if unmapped_bam[cond]:
                    sam = pysam.Samfile(unmapped_bam[cond])
                else: sam = []
                additional = {}
                for read in sam:
                    t_id = sam.getrname(read.tid).split('|')[0]
                    if transcript_mapping.get(t_id) and exons_in_trans.get(t_id):
                        lag = 0
                        r_start = read.pos
                        r_end = r_start + read.rlen
                        E = exons_in_trans[t_id]
                        for e in E:
                            e_start, e_end = exon_mapping[e][2:4]
                            e_len = e_end-e_start
                            if lag <= r_start <= lag+e_len:
                                additional[e] = additional.get(e,0) + 0.5
                            if lag <= r_end <= lag+e_len:
                                additional[e] = additional.get(e,0) + 0.5
                            lag += e_len
                additionals[cond] = additional
                sam.close()

    """ Build exon pileups from bam files """
    print "Build pileups"
    exon_pileups={}; nreads={}
    for gid,files in bam_files.iteritems():
        k = 0
        for rid,f in files.iteritems():
            k+=1
            cond = group_names[gid]+'.'+str(k)
            exon_pileup = build_pileup(f['bam'], exons)
            if unmapped_fastq[cond] and cond in additionals:
                for a,x in additionals[cond].iteritems():
                    exon_pileup[a] = exon_pileup.get(a,0) + x
            exon_pileups[cond] = [exon_pileup[e] for e in exonsID] # {cond1.run1: {pileup}, cond1.run2: {pileup}...}
            nreads[cond] = nreads.get(cond,0) + sum(exon_pileup.values()) # total number of reads
            print "....Pileup", cond, "done"

    """ Print counts for exons """
    hconds = ["counts."+c for c in conditions] + ["rpkm."+c for c in conditions]
    exons_data = [exonsID]+list(counts)+list(rpkm)+[starts,ends,genesID,genesName,strands,echr]
    exons_file = None; genes_file = None; trans_file = None

    print "Get scores of exons"
    exons_data = zip(*exons_data)
    header = ["ExonID"] + hconds + ["Start","End","GeneID","GeneName","Strand","Chromosome"]
    exons_data = zip(*exons_data)
    exons_file = save_results(ex, exons_data, conditions, group_ids, assembly, header=header, feature_type="EXONS")

    return {"exons":exons_file, "genes":genes_file, "transcripts":trans_file}


#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#

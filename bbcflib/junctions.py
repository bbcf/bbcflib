"""
======================
Module: bbcflib.junctions
======================

Methods of the bbcflib's `junctions` worflow. The main function is ``junctions_workflow()``,
and is usually called by bbcfutils' ``run_junctions.py``, e.g. command-line:

``python run_junctions.py -v local -c gapkowt.config -d junctions``

Uses SOAPsplice ...
"""

# Built-in modules #
import os, pysam

# Internal modules #
from bbcflib.common import writecols, set_file_descr, unique_filename_in
from bbcflib import mapseq, genrep
from bbcflib.bFlatMajor.common import cobble, sorted_stream
from bbcflib.btrack import track, FeatureStream, convert
from bein import program


def fetch_mappings(assembly):
    """Given an assembly ID, returns a tuple
    ``(gene_mapping, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans)``

    * [0] gene_mapping is a dict ``{gene ID: (gene name,start,end,length,chromosome)}``
    * [1] transcript_mapping is a dictionary ``{transcript ID: (gene ID,start,end,length,chromosome)}``
    * [2] exon_mapping is a dictionary ``{exon ID: ([transcript IDs],gene ID,start,end,chromosome)}``
    * [3] trans_in_gene is a dict ``{gene ID: [IDs of the transcripts it contains]}``
    * [4] exons_in_trans is a dict ``{transcript ID: [IDs of the exons it contains]}``

    'Lenghts' are always the sum of the lengths of the exons in the gene/transcript.

    :param path_or_assembly_id: can be a numeric or nominal ID for GenRep
    (e.g. 11, 76 or 'hg19' for H.Sapiens), or a path to a file containing a
    pickle object which is read to get the mapping.
    """
    gene_mapping = assembly.get_gene_mapping()
    transcript_mapping = assembly.get_transcript_mapping()
    exon_mapping = assembly.get_exon_mapping()
    exons_in_trans = assembly.get_exons_in_trans()
    trans_in_gene = assembly.get_trans_in_gene()
    mapping = (gene_mapping, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans)
    return mapping

def fetch_labels(bamfile):
    """Returns a list of the exons/transcripts labels in the header of *bamfile*."""
    try: sam = pysam.Samfile(bamfile, 'rb') #bam input
    except ValueError: sam = pysam.Samfile(bamfile,'r') #sam input
    #labels = zip(sam.references,sam.lengths)
    labels = [(t['SN'],t['LN']) for t in sam.header['SQ']]
    sam.close()
    return labels

def build_pileup(bamfile, labels):
    """From a BAM file, returns a dictionary of the form {feature ID: number of reads that mapped to it}.

    :param bamfile: name of a BAM file.
    :param labels: index references - as returned by **fetch_labels()**
    :type bamfile: string
    :type exons: list
    """
    class Counter(object):
        def __init__(self):
            self.n = 0
        def __call__(self, alignment):
            self.n += 1

    counts = {}
    try: sam = pysam.Samfile(bamfile, 'rb')
    except ValueError: sam = pysam.Samfile(bamfile,'r')
    c = Counter()
    for l in labels:
        sam.fetch(l[0], 0, l[1], callback=c) #(name,0,length,Counter())
        #The callback (c.n += 1) is executed for each alignment in a region
        counts[l[0].split('|')[0]] = c.n
        c.n = 0
    sam.close()
    return counts

def save_results(ex, cols, conditions, group_ids, assembly, header=[], feature_type='features'):
    """Save results in a tab-delimited file, one line per feature, one column per run.

    :param ex: bein's execution.
    :param cols: list of iterables, each element being a column to write in the output.
    :param conditions: list of strings corresponding to descriptions of the different samples.
    :param group_ids: dictionary {group name: group ID}.
    :param assembly: a GenRep Assembly object.
    :param header: list of strings, the column headers of the output file.
    :param feature_type: (str) the kind of feature of which you measure the expression.
    """
    conditions = tuple(conditions)
    # Tab-delimited output with all information
    output_tab = unique_filename_in()
    writecols(output_tab,cols,header=header, sep="\t")
    description = set_file_descr(feature_type.lower()+"_expression.tab", step="pileup", type="txt")
    ex.add(output_tab, description=description)
    # Create one track for each group
    if feature_type in ['GENES','EXONS']:
        ncond = len(conditions)
        groups = [c.split('.')[0] for c in conditions]
        start = cols[2*ncond+1]
        end = cols[2*ncond+2]
        chromosomes = cols[-1]
        rpkm = {}; output_sql = {}
        for i in range(ncond):
            group = conditions[i].split('.')[0]
            nruns = groups.count(group)
            # mean of all replicates in the group
            rpkm[group] = asarray(rpkm.get(group,zeros(len(start)))) + asarray(cols[i+ncond+1]) / nruns
            output_sql[group] = output_sql.get(group,unique_filename_in())
        for group,filename in output_sql.iteritems():
            # SQL track
            tr = track(filename+'.sql', fields=['chr','start','end','score'], chrmeta=assembly.chrmeta)
            lines = {}
            for n,c in enumerate(chromosomes):
                if not(c in tr.chrmeta): continue
                if c in lines: lines[c].append((int(start[n]),int(end[n]),rpkm[group][n]))
                else: lines[c] = []
            for chrom, feats in lines.iteritems():
                tr.write(cobble(sorted_stream(FeatureStream(feats, fields=['start','end','score']))),chrom=chrom,clip=True)
            description = set_file_descr(feature_type.lower()+"_"+group+".sql", step="pileup", type="sql", \
                                         groupId=group_ids[group], gdv='1')
            ex.add(filename+'.sql', description=description)
            # UCSC-BED track
            convert(filename+'.sql',filename+'.bedGraph')
            description = set_file_descr(feature_type.lower()+"_"+group+".bedGraph", step="pileup", type="bedGraph", \
                                         groupId=group_ids[group], ucsc='1')
            ex.add(filename+'.bedGraph', description=description)
    print feature_type+": Done successfully."
    return output_tab


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

    """ Treat data """
    print "Process data"
    starts = asarray(starts, dtype=numpy.float_)
    ends = asarray(ends, dtype=numpy.float_)
    nreads = asarray([nreads[cond] for cond in conditions], dtype=numpy.float_)
    counts = asarray([exon_pileups[cond] for cond in conditions], dtype=numpy.float_)

    hconds = ["counts."+c for c in conditions] + ["rpkm."+c for c in conditions]
    genesName, echr = zip(*[(x[0],x[-1]) for x in [gene_mapping.get(g,("NA",)*6) for g in genesID]])
    exons_data = [exonsID]+list(counts)+list(rpkm)+[starts,ends,genesID,genesName,strands,echr]
    exons_file = None; genes_file = None; trans_file = None

    """ Print counts for exons """
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

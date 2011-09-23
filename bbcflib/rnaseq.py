"""
======================
Module: bbcflib.rnaseq
======================

Methods of the bbcflib's RNA-seq worfow
"""

# Built-in modules #
import os, sys, pickle, json, pysam, urllib, math, time
from itertools import combinations

# Internal modules #
from bbcflib.mapseq import map_groups
from bbcflib.genrep import GenRep
from bbcflib.common import timer, results_to_json
from bein import program
from bein.util import unique_filename_in

# Other modules #
import numpy
from rpy2 import robjects
import rpy2.robjects.packages as rpackages
import rpy2.robjects.numpy2ri
import rpy2.rlike.container as rlc
import csv

################################################################################

def fetch_mappings(path_or_assembly_id):
    """Given an assembly ID, return a tuple
    (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans)
    - gene_ids is a list of gene ID's
    - gene_names is a dict {gene ID: gene name}
    - transcript_mapping is a dictionary {transcript ID: gene ID}
    - exon_mapping is a dictionary {exon ID: (transcript ID, gene ID)}
    - trans_in_gene is a dict {gene ID: IDs of the transcripts it contains}
    - exons_in_trans is a dict {transcript ID: IDs of the exons it contains}

    *path_or_assembly_id* can be a numeric or nominal ID for GenRep
    (e.g. 11, 76 or 'hg19' for H.Sapiens), or a path to a file containing a
    pickle object which is read to get the mapping.
    """
    if os.path.exists(str(path_or_assembly_id)):
        with open(path_or_assembly_id, 'rb') as pickle_file:
            mapping = pickle.load(pickle_file)
        print "Mapping found in", os.path.abspath(path_or_assembly_id)
        return mapping
    else:
        grep_root = '/db/genrep'
        grep = GenRep(url='http://bbcftools.vital-it.ch/genrep/',root=grep_root)
        assembly = grep.assembly(path_or_assembly_id)
        mdfive = assembly.md5
        mappings_path = os.path.join(grep_root,'nr_assemblies/exons_pickle/')+mdfive+".pickle"
        with open(mappings_path, 'rb') as pickle_file:
            mapping = pickle.load(pickle_file)
        print "Mapping for assembly", mdfive, "found in /db/"
        return mapping

def map_runs(fun, runs):
    """Parallelization of fun(run) executions"""
    futures = {}
    for group_id, run_list in runs.iteritems():
        futures[group_id] = [fun(run) for run in run_list]
    results = {}
    for group_id, future_list in futures.iteritems():
        results[group_id] = [f.wait() for f in future_list]
    return results

def exons_labels(bamfile):
    """List of the exons labels (SN: 'reference sequence name'; LN: 'sequence length') in *bamfile*."""
    sam = pysam.Samfile(bamfile, 'rb')
    labels = [(t['SN'],t['LN']) for t in sam.header['SQ']]
    sam.close()
    return labels

def pileup_file(bamfile, exons):
    """Return a dictionary {exon ID: count}, the pileup of *exons* from *bamfile*."""
    counts = {}
    sam = pysam.Samfile(bamfile, 'rb')

    class Counter(object):
        def __init__(self):
            self.n = 0
        def __call__(self, alignment):
            self.n += 1

    c = Counter()
    for i,exon in enumerate(exons):
        sam.fetch(exon[0], 0, exon[1], callback=c) #(exon_name,0,exon_length,Counter())
        #The callback (c.n += 1) is executed for each alignment in a region
        counts[exon[0]] = c.n
        c.n = 0
    sam.close()
    return counts

def pairs_to_test(controls):
    """*controls* is a dictionary of group_ids to True/False.

    If all the values are True or all the values are False, then it
    returns all unique pairs of group IDs.  Otherwise it returns all
    combinations of one group ID with value True and one with value False.
    """
    if all(controls.values()) or not(any(controls.values())):
        return list(combinations(controls.keys(), 2))
    else:
        return [(x,y) for x in controls.keys() for y in controls.keys()
                if controls[x] and not(controls[y])]

def translate_gene_ids(fc_ids, dictionary):
    """Replace (unique) gene IDs by (not unique) gene names.
    *fc_ids* is a dict {gene_id: whatever}
    *dictionary* is a dict {gene_id: gene_name} """
    names = []
    for s in fc_ids.keys():
        start = s.find("ENSG")
        if start != -1:
            end = s.split("ENSG")[1].find("|")
            gene_id = "ENSG" + s.split("ENSG")[1].split("|")[0]
            names.append(s.replace(gene_id, dictionary.get(gene_id,gene_id)))
        else: names.append(s)
    fc_names = dict(zip(names,fc_ids.values()))
    return fc_names

def estimate_size_factors(counts):
    """
    The total number of reads may be different between conditions or replicates.
    This treatment makes different count sets being comparable. If rows are features
    and columns are different conditions/replicates, each column is divided by the
    geometric mean of the rows.
    The median of these ratios is used as the size factor for this column. Size
    factors may be used for further variance sabilization.
    
    * *counts* is an array of counts, each line representing a transcript, each
    column a different run.
    """
    counts = numpy.array(counts)
    geo_means = numpy.exp(numpy.mean(numpy.log(counts), axis=0))
    mean = counts/geo_means
    size_factors = numpy.median(mean[:,geo_means>0], axis=1)
    res = counts.transpose()/size_factors
    res = res.transpose()
    return res, size_factors

def save_results(data, conditions=["C1","C2"], filename=None):
    """Save results in a CSV file, one line per transcript, one column per run """
    if not filename:
        filename = unique_filename_in()
    with open(filename,"wb") as f:
        c = csv.writer(f, delimiter='\t')
        c.writerow(["id"]+conditions)
        for k,v in data.iteritems():
            c.writerow([k]+list(v))
    return filename

@timer
def genes_expression(gene_ids, exon_mapping, dexons):
    z = numpy.array([[0]*len(dexons)]*len(gene_ids))
    dgenes = dict(zip(gene_ids,z))
    print dexons.iteritems().next()[0]
    print dexons.iteritems().next()[1]
    for e,c in dexons.iteritems():
        e = e.split('|')[0]
        if exon_mapping.get(e):
            dgenes[exon_mapping[e][1]] += c
    return dgenes

@timer
def transcripts_expression(gene_ids, transcript_mapping, trans_in_gene, exons_in_trans, dexons):
    z = numpy.array([[0]*len(dexons)]*len(transcript_mapping.keys()))
    dtrans = dict(zip(transcript_mapping.keys(),z))
    for g in gene_ids:
        if trans_in_gene.get(g):
            tg = trans_in_gene[g]
            eg = []
            for t in tg:
                if exons_in_trans.get(t):
                    eg.extend(exons_in_trans[t])
            M = numpy.zeros((len(eg),len(tg)))
            exons_1 = numpy.zeros(len(eg))
            exons_2 = numpy.zeros(len(eg))
            for i,e in enumerate(eg):
                for j,t in enumerate(tg):
                    if exons_in_trans.get(t):
                        ebt = exons_in_trans[t]
                    if e in ebt:
                        M[i,j] = 1
                if dexons.get(e):
                    exons_1[i] += dexons[e][0]
                    exons_2[i] += dexons[e][1]
            transcripts_1 = numpy.dot(numpy.linalg.pinv(M),exons_1)
            transcripts_2 = numpy.dot(numpy.linalg.pinv(M),exons_2)
            for k,t in enumerate(tg):
                if dtrans.get(t) is not None:
                    dtrans[t][0] = transcripts_1[k]
                    dtrans[t][1] = transcripts_2[k]
    return dtrans

@timer
def get_expressions(pileups, assembly_id, target=["exons","genes","transcripts"]):
    """
    Get transcripts expression from exons expression values.
    Size factors are computed to make pileups comparable is one had much more
    read counts than the other one.
    
    * *assembly_id* can be a numeric or nominal ID for GenRep
    (e.g. 76 or 'hg19' for H.Sapiens), or a path to a file containing a
    pickle file which is read to get the mapping.
    * *target* is a string or array of strings indicating the features you want to compare.
    Targets can be 'genes', 'transcripts', or 'exons'. E.g. ['genes','transcripts'], or 'genes'.
    * *pileups* is a dictionary of the form {cond1: [{pileup bam1},{pileup bam2},...], cond2:...},
    where each *pileup* is a dictionary {transcript id: number of reads}.
    """

    print "Loading mappings..."
    assembly_id = "../temp/nice_features/nice_mappings" # testing code
    mappings = fetch_mappings(assembly_id)
    """ - gene_ids is a list of gene ID's
        - gene_names is a dict {gene ID: gene name}
        - transcript_mapping is a dictionary {transcript ID: gene ID}
        - exon_mapping is a dictionary {exon ID: (transcript ID, gene ID)}
        - trans_in_gene is a dict {gene ID: IDs of the transcripts it contains}
        - exons_in_trans is a dict {transcript ID: IDs of the exons it contains} """
    print "Loaded."
    (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans) = mappings

    conditions = pileups.keys()
    res = numpy.array([p.values() for p in pileups.values()], dtype=numpy.float32)
    for i in range(len(res.ravel())):
        if res.flat[i]==0: res.flat[i] += 1.0
    res, size_factors = estimate_size_factors(res)

    exon_ids = pileups[conditions[0]].keys()
    gene_ids = [e.split('|')[1] for e in exon_ids]
    coords = [(e.split('|')[2],e.split('|')[3]) for e in exon_ids]
    exon_ids = [e.split('|')[0] for e in exon_ids]
    
    dexons = dict(zip(exon_ids,zip(*res)))
    exons_filename = save_results(translate_gene_ids(dexons, gene_names), conditions=conditions)

    # Get fold change for genes
    genes_filename=None
    if "genes" in target:
        dgenes = genes_expression(gene_ids, exon_mapping, dexons)
        genes_filename = save_results(translate_gene_ids(dgenes, gene_names), conditions=conditions)

    ### Get fold change for the transcripts using pseudo-inverse
    trans_filename=None
    if "transcripts" in target:
        dtrans = transcripts_expression(gene_ids, transcript_mapping, trans_in_gene, exons_in_trans, dexons)
        trans_filename = save_results(translate_gene_ids(dtrans, gene_names), conditions=conditions)
        
    result = {"exons":exons_filename, "genes":genes_filename, "transcripts":trans_filename,
              "conditions": conditions, "size factors": size_factors}
    return result


def rnaseq_workflow(ex, job, assembly, bam_files, target=["genes"], via="lsf", output=None):
    """
    Main function of the workflow. 
    
    * *output*: alternative name for output file. Otherwise it is random.
    * *target*: a string or array of strings indicating the features you want to compare.
    Targets can be 'genes', 'transcripts', or 'exons'. E.g. ['genes','transcripts'], or 'genes'.
    (This part of the workflow may fail if there are too few reads for DESeq to estimate
    variances amongst exons.)
    * *via*: 'local' or 'lsf'
    * *job* is a Job object (or a dictionary of the same form) as returned from
    HTSStation's frontend.
    """
    group_names={}; run_names={}; runs={}; controls={};
    groups = job.groups
    assembly_id = job.assembly_id
    for gid,group in groups.iteritems():
        group_names[gid] = str(group['name']) # group_names = {gid: name}
        for rid, run in group['runs'].iteritems():
            runs[str(gid)+'.'+str(rid)] = (gid, run["url"])    # runs = {gid: [{run info}, {run info}]}
    if isinstance(target,str): target=[target]
    
    # All the bam_files were created against the same index, so
    # they all have the same header in the same order.  I can take
    # the list of exons from just the first one and use it for all of them.
    # Format: ('exonID|geneID|start|end|strand', length)
    exons = exons_labels(bam_files[groups.keys()[0]][groups.values()[0]['runs'].keys()[0]]['bam'])
    
    exon_pileups = {}
    for gid,files in bam_files.iteritems():
        for rid,f in files.iteritems():
            name = group_names[gid]+'.'+str(rid)
            exon_pileups[name] = []
        for rid,f in files.iteritems():
            name = group_names[gid]+'.'+str(rid)
            exon_pileup = pileup_file(f['bam'], exons) #{exon_id: count}
            exon_pileups[name] = exon_pileup #{cond1.run1: {pileup}, cond1.run2: {pileup}...}

    print "Get expressions for genes and transcripts..."
    result = get_expressions(exon_pileups, assembly_id, target)
    
    conditions = tuple(result["conditions"])
    conditions_s = '%s, '*(len(conditions)-1)+' and %s'
    if "exons" in target:
        exons_file = result.get("exons")
        if exons_file:
            ex.add(exons_file,
                   description="csv:Comparison of EXONS in conditions "+conditions_s % conditions)
            print "EXONS: Done successfully."
    if "genes" in target:
        genes_file = result.get("genes")
        if genes_file:
            ex.add(genes_file,
                   description="csv:Comparison of GENES in conditions "+conditions_s % conditions)
            print "GENES: Done successfully."
    if "transcripts" in target:
        trans_file = result.get("transcripts");
        if trans_file:
            ex.add(trans_file,
                   description="csv:Comparison of TRANSCRIPTS in conditions "+conditions_s % conditions)
            print "TRANSCRIPTS: Done successfully."
    print "Done."

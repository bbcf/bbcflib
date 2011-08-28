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
from bbcflib.common import timer
from bein import program, execution, MiniLIMS
from bein.util import unique_filename_in

# Other modules #
import numpy
from scipy import stats
from scipy.interpolate import UnivariateSpline
from rpy2 import robjects
import rpy2.robjects.packages as rpackages
import rpy2.robjects.numpy2ri
import rpy2.rlike.container as rlc
import cogent.db.ensembl as ensembl
import csv

import pdb

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
    (e.g. 76 or 'hg19' for H.Sapiens), or a path to a file containing a
    pickle object which is read to get the mapping.
    """
    assembly_id = path = path_or_assembly_id
    nr_assemblies = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies.json").read()
    nr_assemblies = json.loads(nr_assemblies)
    for a in nr_assemblies:
        if a['nr_assembly']['id'] == assembly_id or a['nr_assembly']['name'] == assembly_id:
            mdfive = a['nr_assembly']['md5']; break

    if os.path.exists(str(path)):
        with open(str(path), 'rb') as pickle_file:
            mapping = pickle.load(pickle_file)
        print "Mapping found in", os.path.dirname(str(path))
        return mapping
    elif os.path.exists("/db/genrep/nr_assemblies/exons_pickle/"+str(mdfive)+".pickle"):
        with open("/db/genrep/nr_assemblies/exons_pickle/"+mdfive+".pickle", 'rb') as pickle_file:
            mapping = pickle.load(pickle_file)
        print "Mapping found in /db/"
        return mapping
    else:
        genrep = GenRep('http://bbcftools.vital-it.ch/genrep/','/db/genrep/nr_assemblies/exons_pickle')
        assembly = genrep.assembly(mdfive)
        with open(assembly.index_path + '.pickle', 'rb') as pickle_file:
            mapping = pickle.load(pickle_file)
        print "Mapping found on GenRep"
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
    '''Replace (unique) gene IDs by (not unique) gene names.
    *fc_ids* is a dict {gene_id: whatever}
    *dictionary* is a dict {gene_id: gene_name} '''
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

def save_results(data, filename=None):
    '''Save results in a CSV file. Data must be of the form {id:(mean,fold_change)} '''
    if not filename:
        filename = unique_filename_in()
    with open(filename,"wb") as f:
        c = csv.writer(f, delimiter='\t')
        c.writerow(["id", "baseMean", "foldChange"])
        for k,v in data.iteritems():
            c.writerow([k,v[0],v[1]])
    return filename

@program
def external_maplot(csv_filename, mode='normal', deg=4, bins=30, assembly_id=None, output=None):
    if not output: output = unique_filename_in()
    call = ["maplot.py", csv_filename, str(mode), str(deg), str(bins), str(assembly_id), str(output)]
    return {"arguments": call, "return_value": output}

@program
def external_deseq(cond1_label, cond1, cond2_label, cond2, assembly_id,
                   target=["exons","genes","transcripts"], method="normal", maplot=None):
    output = unique_filename_in()
    c1 = unique_filename_in()
    with open(c1,'wb') as f:
        pickle.dump(cond1,f,pickle.HIGHEST_PROTOCOL)
    c2 = unique_filename_in()
    with open(c2,'wb') as f:
        pickle.dump(cond2,f,pickle.HIGHEST_PROTOCOL)
    targ = unique_filename_in()
    with open(targ,'wb') as f:
        pickle.dump(target,f,pickle.HIGHEST_PROTOCOL)
    call = ["run_deseq.py", c1, c2, cond1_label, cond2_label, str(assembly_id), 
            targ, method, str(maplot), output]
    return {"arguments": call, "return_value": output}

@timer
def inference(cond1_label, cond1, cond2_label, cond2, assembly_id,
              target=["exons","genes","transcripts"], method="normal", maplot=None):
    """ Writes a CSV file for each selected type of feature,
    each row being of the form: Feature_ID // Mean_expression // Fold_change.
    It calls an R session to use DESeq for size factors and variance stabilization.

    * *cond1* and *cond2* are dictionaries - pileups - of the form
    {feature ID: number of reads mapping to it}.
    * *cond1_label* and *cond2_label* are string which will be used
    to identify the two conditions in R.
    * *assembly_id* can be a numeric or nominal ID for GenRep
    (e.g. 76 or 'hg19' for H.Sapiens), or a path to a file containing a
    pickle file which is read to get the mapping.
    * *target* is a string or array of strings indicating the features you want to compare.
    Targets can be 'genes', 'transcripts', or 'exons'. E.g. ['genes','transcripts'], or 'genes'.
    * *method* can be 'normal' or 'blind', the method used for DESeq variances estimation.
    * *maplot*: MA-plot of data.
    - If 'interactive', one can click on a point (gene or exon) to display its name;
    - if 'normal', name of genes over 99.9%/under 0.1% quantiles are displayed;
    - if None, no figure is produced.
    """

    """ - gene_ids is a list of gene ID's
        - gene_names is a dict {gene ID: gene name}
        - transcript_mapping is a dictionary {transcript ID: gene ID}
        - exon_mapping is a dictionary {exon ID: (transcript ID, gene ID)}
        - trans_in_gene is a dict {gene ID: IDs of the transcripts it contains}
        - exons_in_trans is a dict {transcript ID: IDs of the exons it contains} """
    mappings = fetch_mappings(assembly_id)
    (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans) = mappings

    # Pass the data into R as a data frame
    data_frame_contents = rlc.OrdDict([(cond1_label+'-'+str(i), robjects.IntVector(c.values()))
                                       for i,c in enumerate(cond1)] +
                                      [(cond2_label+'-'+str(i), robjects.IntVector(c.values()))
                                       for i,c in enumerate(cond2)])
    data_frame = robjects.DataFrame(data_frame_contents)
    data_frame.rownames = cond1[0].keys()

    conds = robjects.StrVector([cond1_label for x in cond1] + [cond2_label for x in cond2]).factor()
    
    deseq = rpackages.importr('DESeq')
    cds = deseq.newCountDataSet(data_frame, conds)
    cds = deseq.estimateSizeFactors(cds)
    try: cds = deseq.estimateVarianceFunctions(cds,method=method)
    except : raise rpy2.rinterface.RRuntimeError("Too few reads to estimate variances with DESeq")
    res = deseq.getVarianceStabilizedData(cds)
    
    exon_ids = list(list(res.names)[0]) # 'res.names[0]' kills the python session.
    cond1 = numpy.abs(numpy.array(res.rx(True,1), dtype='f')) # replicates?
    cond2 = numpy.abs(numpy.array(res.rx(True,2), dtype='f'))
          #abs: a zero count may become slightly negative after normalization - see DESeq
    names=[]; means=[]; ratios=[]
    means = numpy.sqrt(cond1*cond2)
    ratios = cond1/cond2
    fc_exons = dict(zip(exon_ids,zip(means,ratios)))
    if "exons" in target:
        exons_filename = save_results(translate_gene_ids(fc_exons, gene_names))
        
    fc_exons = dict(zip([e.split('|')[0] for e in fc_exons.keys()], fc_exons.values()))
        
    # Get fold change for genes
    if "genes" in target:
        fc_genes = dict(zip(gene_ids,
                            numpy.array(zip(numpy.zeros(len(gene_ids)),
                                            numpy.zeros(len(gene_ids))))  ))
        for e,c in fc_exons.iteritems():
            e = e.split('|')[0]
            if exon_mapping.get(e):
                fc_genes[exon_mapping[e][1]] += c
        genes_filename = save_results(translate_gene_ids(fc_genes, gene_names))
            
    ### Get fold change for the transcripts using pseudo-inverse
    if "transcripts" in target:
        fc_trans = dict(zip(transcript_mapping.keys(),
                            numpy.array(zip(numpy.zeros(len(transcript_mapping.keys())),
                                            numpy.zeros(len(transcript_mapping.keys()))))  ))
        for g in gene_ids:
            tg = trans_in_gene[g]
            if len(tg) == 1:
                t = tg[0]
                eg = exons_in_trans[t]
                for e in eg:
                    if fc_exons.get(e):
                        fc_trans[t] += fc_exons[e]
            eg = []
            for t in tg:
                if exons_in_trans.get(t):
                    eg.extend(exons_in_trans[t])
            M = numpy.zeros((len(eg),len(tg)))
            exons_fold = numpy.zeros(len(eg))
            exons_mean = numpy.zeros(len(eg))
            for i,e in enumerate(eg):
                for j,t in enumerate(tg):
                    if exons_in_trans.get(t):
                        ebt = exons_in_trans[t]
                    if e in ebt:
                        M[i,j] = 1
                if fc_exons.get(e):
                    exons_fold[i] += fc_exons[e][1]
                    exons_mean[i] += fc_exons[e][0]
            transcripts_fold = numpy.dot(numpy.linalg.pinv(M),exons_fold)
            transcripts_mean = numpy.dot(numpy.linalg.pinv(M),exons_mean)
            for k,t in enumerate(tg):
                if fc_trans.get(t) != None:
                    fc_trans[t][1] = transcripts_fold[k]
                    fc_trans[t][0] = transcripts_mean[k]
        trans_filename = save_results(translate_gene_ids(fc_trans, gene_names))
        
    result = {"exons":exons_filename, "genes":genes_filename, "transcripts":trans_filename}
    return result


def rnaseq_workflow(ex, job, assembly, target=["genes"], bam_files=False, via="lsf", output=None, maplot="normal"):
    """Run RNASeq inference according to *job_info*.

    * *output*: alternative name for output file. Otherwise it is random.
    * *maplot*: MA-plot of data.
    If 'interactive', one can click on a point (gene or exon) to display its name;
    if 'normal', name of genes over 99.9%/under 0.1% quantiles are displayed;
    if None, no figure is produced.
    * *target*: a string or array of strings indicating the features you want to compare.
    Targets can be 'genes', 'transcripts', or 'exons'. E.g. ['genes','transcripts'], or 'genes'.
    (This part of the workflow may fail if there are too few reads for DESeq to estimate
    variances amongst exons.)
    * *via*: 'local' or 'lsf'
    * *job* is a Job object (or a dictionary of the same form) as returned from
    HTSStation's frontend.
    """
    names = {}; runs = {}; controls = {};
    groups = job.groups
    assembly_id = job.assembly_id
    for i,group in groups.iteritems():
        names[i] = str(group['name'])
        runs[i] = group['runs'].values()
        controls[i] = group['control']
    if isinstance(target,str): target=[target]

    if not bam_files: #shouldn't be there
        print "Alignment..."
        fastq_root = os.path.abspath(ex.working_directory)
        bam_files = map_groups(ex, job, fastq_root, assembly_or_dict = assembly)
        print "Reads aligned."
    else: #testing
        print "Loading bam files..."
        with open("../temp/bam_files","rb") as f:
            bam_files = pickle.load(f)
            
        id1 = ex.lims.import_file("../temp/"+"100k1.bam")
        id2 = ex.lims.import_file("../temp/"+"100k2.bam")
        id3 = ex.lims.import_file("../temp/"+"100k1.bam.bai")
        id4 = ex.lims.import_file("../temp/"+"100k2.bam.bai")
        ex.lims.associate_file(id3,id1,template="%s.bai")
        ex.lims.associate_file(id4,id2,template="%s.bai")
        f1 = ex.use(id1)
        f2 = ex.use(id2)
        bam_files.values()[0].values()[0]['bam'] = f1
        bam_files.values()[1].values()[0]['bam'] = f2
        print "Loaded."

        #assembly_id = "../temp/nice_mappings"

        ## path = '../archive/mapseq_full2_lims.files/'
        ## id1 = ex.lims.import_file(os.path.join(path,'eyFS5XeDLCrzIhVuYItV'))
        ## id2 = ex.lims.import_file(os.path.join(path,'cpCqC1tubV6zoBWAmoqy'))
        ## id3 = ex.lims.import_file(os.path.join(path,'eyFS5XeDLCrzIhVuYItV.bai'))
        ## id4 = ex.lims.import_file(os.path.join(path,'cpCqC1tubV6zoBWAmoqy.bai'))
        ## ex.lims.associate_file(id3,id1,template="%s.bai")
        ## ex.lims.associate_file(id4,id2,template="%s.bai")
        ## f1 = ex.use(id1)
        ## f2 = ex.use(id2)
        ## bam_files.values()[0].values()[0]['bam'] = f1
        ## bam_files.values()[1].values()[0]['bam'] = f2
    
    # All the bam_files were created against the same index, so
    # they all have the same header in the same order.  I can take
    # the list of exons from just the first one and use it for all of them.
    # Format: ('exonID|geneID|start|end|strand', length)
    exons = exons_labels(bam_files.values()[0][1]['bam'])
    
    exon_pileups = {}
    for condition,files in bam_files.iteritems():
        exon_pileups[condition] = []
        for f in files.values():
            exon_pileup = pileup_file(f['bam'], exons) #{exon_id: count}
            exon_pileups[condition].append(exon_pileup) #{cond1: [{pileup bam1},{pileup bam2},...], cond2:...}

    futures = {}
    for (c1,c2) in pairs_to_test(controls):
        if len(runs[c1]) + len(runs[c2]) > 2: method = "normal"
        else: method = "blind"
        print "External DESeq..."
        if 1:
            futures[(c1,c2)] = external_deseq.nonblocking(ex,
                                   names[c1], exon_pileups[c1], names[c2], exon_pileups[c2],
                                   assembly_id, target, method, maplot, via=via)
            
            for c,f in futures.iteritems():
                result_filename = f.wait()
                with open(result_filename,"rb") as f:
                    result = pickle.load(f)

                conditions_desc = (names[c[0]], names[c[1]])
                if "exons" in target:
                    exons_file = result.get("exons")
                    if exons_file:
                        ex.add(exons_file, description="csv:Comparison of EXONS in conditions \
                        '%s' and '%s' " % conditions_desc)
                        print "EXONS: Done successfully."
                    else: print >>sys.stderr, "Exons: Failed during inference, \
                                               probably because of too few reads for DESeq stats."
                if "genes" in target:
                    genes_file = result.get("genes")
                    if genes_file:
                        ex.add(genes_file, description="csv:Comparison of GENES in conditions \
                        '%s' and '%s' " % conditions_desc)
                        print "GENES: Done successfully."
                if "transcripts" in target:
                    trans_file = result.get("transcripts");
                    if trans_file:
                        ex.add(trans_file, description="csv:Comparison of TRANSCRIPTS in conditions \
                        '%s' and '%s' " % conditions_desc)
                    print "TRANSCRIPTS: Done successfully."

        if 0: #testing
            futures[(c1,c2)] = inference(names[c1], exon_pileups[c1], names[c2], exon_pileups[c2],
                                   assembly_id, target, method, maplot)
            for c,f in futures.iteritems():
                result = f

                conditions_desc = (names[c[0]], names[c[1]])
                if "exons" in target:
                    exons_file = result.get("exons")
                    if exons_file:
                        ex.add(exons_file, description="csv:Comparison of EXONS in conditions \
                        '%s' and '%s' " % conditions_desc)
                        print "EXONS: Done successfully."
                    else: print >>sys.stderr, "Exons: Failed during inference, \
                                               probably because of too few reads for DESeq stats."
                if "genes" in target:
                    genes_file = result.get("genes")
                    if genes_file:
                        ex.add(genes_file, description="csv:Comparison of GENES in conditions \
                        '%s' and '%s' " % conditions_desc)
                        print "GENES: Done successfully."
                if "transcripts" in target:
                    trans_file = result.get("transcripts");
                    if trans_file:
                        ex.add(trans_file, description="csv:Comparison of TRANSCRIPTS in conditions \
                        '%s' and '%s' " % conditions_desc)
                    print "TRANSCRIPTS: Done successfully."
            
        print "Done."

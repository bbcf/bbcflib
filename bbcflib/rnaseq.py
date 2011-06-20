#!/bin/env python

"""
======================
Module: bbcflib.rnaseq
======================

"""

import pickle
import json
import pysam
import numpy
import urllib
import time
from itertools import combinations
from bein.util import *
from bbcflib import *
from bbcflib.mapseq import plot_stats, bamstats
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.robjects.numpy2ri
import rpy2.rlike.container as rlc
from maplot import MAplot

def path_to_bowtie_index(ex, assembly_id):
    """Return full path to usable bowtie index.

    If *assembly_id* is a string, it is interpreted as a path to a
    bowtie index.  This makes debugging and testing the workflow
    simple.  If *assembly_id* is an integer, it is taken to be a
    reference to GenRep.
    """
    nr_assemblies = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies.json").read()
    nr_assemblies = json.loads(nr_assemblies)
    for a in nr_assemblies:
        if a['nr_assembly']['id'] == assembly_id:
            mdfive = a['nr_assembly']['md5']; break
    #mdfive = "bb27b89826b88823423282438077cdb836e1e6e5"
            
    if isinstance(assembly_id, str):
        if os.path.exists(assembly_id + ".1.ebwt"):
            print "Index found in", os.path.dirname(assembly_id)
            return assembly_id
        else:
            raise ValueError("No bowtie index at %s" % assembly_id)
    elif os.path.exists("/db/genrep/nr_assemblies/exons_bowtie/"+mdfive+".1.ebwt"):
        print "Index found in /db/"
        return "/db/genrep/nr_assemblies/exons_bowtie/"+mdfive
    else:
        genrep = GenRep('http://bbcftools.vital-it.ch/genrep/',
                        '/db/genrep/nr_assemblies/exons_bowtie')
        assembly = genrep.assembly(mdfive)
        print "Index found on GenRep"
        return assembly.index_path

def fetch_read_files(ex, runs):
    """Return a dictionary of full paths to FASTQ files of reads.

    The input *runs* has type {Int : [String]} or {Int : [Dict]}.  If
    String, then fetch_read_files is being called for local use or
    testing, and the strings are paths to FASTQ files.  If Dict, then
    each dict specifies a run in the DAF LIMS with keys ``facility``,
    ``machine``, ``run``, and ``lane``.  This is exactly what is
    returned as a run in a Job object in bbcflib.frontend.

    The return type is a dictionary {Int : [String]}, where the
    strings are again paths to FASTQ files.
    """
    files = {}
    for group_id, run_list in runs.iteritems():
        files[group_id] = []
        if all(isinstance(r, str) for r in run_list):
            for file_path in run_list:
                if not(os.path.exists(file_path)):
                    raise ValueError("File %s in group %d does not exist" % (file_path,group_id))
                else:
                    files[group_id].append({"path":file_path})
        elif all(isinstance(x, dict) for x in run_list):
            daflims = DAFLIMS(username='jrougemont', password='cREThu6u')
            for run in run_list:
                filename = unique_filename_in()
                files[group_id].append(daflims.fetch_fastq(str(run['facility']), str(run['machine']), 
                                                    run['run'], run['lane'], filename))
    return files


def fetch_transcript_mapping(ex, assembly_id):
    """Given an assembly ID, return a dictionary giving an exon to orf mapping.

    The *assembly_id* can be a string or integer.  If it is an
    integer, the mapping is looked up from GenRep.  Otherwise it is
    assumed to be a file containing a JSON object which is read to get
    the mapping.
    """
    nr_assemblies = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies.json").read()
    nr_assemblies = json.loads(nr_assemblies)
    for a in nr_assemblies:
        if a['nr_assembly']['id'] == assembly_id:
            mdfive = a['nr_assembly']['md5']
    #mdfive = "bb27b89826b88823423282438077cdb836e1e6e5"
    
    if isinstance(assembly_id, str):
        pickle_path = assembly_id + '.pickle'
        if os.path.exists(pickle_path):
            with open(pickle_path, 'rb') as pickle_file:
                mapping = pickle.load(pickle_file)
                return mapping
            print "Exon to ORF mapping found in", os.path.dirname(pickle_path)
        else:
            raise ValueError("No pickle at %s" % pickle_path)
    elif os.path.exists("/db/genrep/nr_assemblies/exons_pickle/"+mdfive+".pickle"):
        print "Exon to ORF mapping found in /db/"
        with open("/db/genrep/nr_assemblies/exons_pickle/"+mdfive+".pickle") as pickle_file:
            mapping = pickle.load(pickle_file)
        return mapping
    else:
        genrep = GenRep('http://bbcftools.vital-it.ch/genrep/',
                        '/db/genrep/nr_assemblies/exons_pickle')
        assembly = genrep.assembly(mdfive)
        with open(assembly.index_path + '.pickle', 'rb') as pickle_file:
            mapping = pickle.load(pickle_file)
        print "Exon to ORF mapping found on GenRep"
        return mapping


def map_runs(fun, runs):
    """Parallelization of fun(run) executions"""
    futures = {}
    for group_id, run_list in runs.iteritems():
        futures[group_id] = [fun(run) for run in run_list]
    results = {}
    for group_id, future_list in futures.iteritems():
        results[group_id] = [f.wait() for f in future_list]
    print "Results:", results.values()
    return results

def align_reads(ex, index, read_files, via="lsf"):
    """Align each member of *read_files* to *index*.

    Returns a dictionary of {integer:['filename','filename',...]},
    each filename pointing to an indexed BAM file of the alignments,
    in the same order as *read_files*.
    """
    print "Aligning reads..."
    def _bowtie(filename):
        return bowtie.nonblocking(ex, index, filename['path'], args=['--best','--strata','-Sqam','100'], via=via)
    def _add_nh(filename):
        return external_add_nh_flag.nonblocking(ex, filename, via=via)
    def _sort(filename):
        return sort_bam.nonblocking(ex, filename, via=via)
    def _index(filename):
        return index_bam.nonblocking(ex, filename, via=via)

    sam_files = map_runs(_bowtie, read_files)
    print "Sam files obtained."
    bam_files = map_runs(_add_nh, sam_files)
    print "Bam files obtained."
    sorted_files = map_runs(_sort, bam_files)
    print "Sorted."
    indexes = map_runs(_index, sorted_files)
    print "Indexed."

    ## try:
    ##     ex.add(sorted_files.values()[0][0], alias="bam1_bb27")
    ##     ex.add(sorted_files.values()[1][0], alias="bam2_bb27")
    ##     ex.add(indexes.values()[0][0], alias="index1_bb27")
    ##     ex.add(indexes.values()[1][0], alias="index2_bb27")
    ## except: print "Failed to add bam files"
    
    return sorted_files

def exons_labels(bamfile):
    """List of the exons labels in *bamfile*."""
    sam = pysam.Samfile(bamfile, 'rb')
    labels = [(t['SN'],t['LN']) for t in sam.header['SQ']] #SN: "reference sequence name"; LN: "sequence length", for each genome sequence SQ.
    sam.close()
    return labels

def pileup_file(bamfile, exons):
    """Return a numpy array of the pileup of *bamfile* in order *exons*."""
    counts = numpy.zeros(len(exons))
    sam = pysam.Samfile(bamfile, 'rb')

    class Counter(object):
        def __init__(self):
            self.n = 0
        def __call__(self, alignment):
            self.n += 1

    c = Counter()
    for i,exon in enumerate(exons):
        sam.fetch(exon[0], 0, exon[1], callback=c)
        counts[i] = c.n
        c.n = 0
    sam.close()
    return counts

def pairs_to_test(controls):
    """*controls* is a dictionary of group_ids to True/False.

    If all the values are True or all the values are False, then it
    returns all unique pairs of group IDs.  Otherwise it returns all
    combinations of one group ID with value True and one with value
    False.
    """
    if all(controls.values()) or not(any(controls.values())):
        return list(combinations(controls.keys(), 2))
    else:
        return [(x,y) for x in controls.keys() for y in controls.keys()
                if controls[x] and not(controls[y])]

@program
def external_deseq(cond1_label, cond1, cond2_label, cond2, transcript_names, method="normal", assembly_id=None, output=None, maplot=None):
    if output:
        result_filename = output
    else: 
        result_filename = unique_filename_in()
    c1 = unique_filename_in()
    with open(c1,'wb') as f:
        pickle.dump(cond1,f,pickle.HIGHEST_PROTOCOL)
    c2 = unique_filename_in()
    with open(c2,'wb') as f:
        pickle.dump(cond2,f,pickle.HIGHEST_PROTOCOL)
    tn = unique_filename_in()
    with open(tn,'wb') as f:
        pickle.dump(transcript_names,f,pickle.HIGHEST_PROTOCOL)
    call = ["run_deseq", c1, c2, tn, cond1_label, cond2_label, method, str(assembly_id), result_filename, str(maplot)]
    return {"arguments": call, "return_value": result_filename}

def results_to_json(lims, exid):
    """Create a JSON string describing the results of execution *exid*.

    The execution is sought in *lims*, and all its output files and
    their descriptions are written to the string.
    """
    produced_file_ids = lims.search_files(source=('execution',exid))
    d = dict([(lims.fetch_file(i)['description'], lims.path_to_file(i))
              for i in produced_file_ids])
    j = json.dumps(d)
    return j

def inference(cond1_label, cond1, cond2_label, cond2, transcript_names, method="normal", assembly_id=None, output=None, maplot=None):
    """Runs DESeq comparing the counts in *cond1* and *cond2*.

    Arguments *cond1* and *cond2* are lists of numpy arrays. Each
    array lists the number of reads mapping to a particular
    transcript. *cond1_label* and *cond2_label* are string which will
    be used to identify the two conditions in R.  *transcript_names*
    is a list of strings, in the same order as the numpy arrays, used
    to name the transcripts in R.

    ``deseq_inference`` writes a tab delimited text file of the
    conditions to R, the first column being the transcript name,
    followed by one column for each numpy array in *cond1*, then one
    column for each numpy array in *cond2*.

    Then it calls DESeq in R, writes out the results in a new,
    randomly named file, and returns that filename.
    """
    ### FAILS IF TOO FEW READS

    # Pass the data into R as a data frame
    data_frame_contents = rlc.OrdDict([(cond1_label+'-'+str(i), robjects.IntVector(c))
                                       for i,c in enumerate(cond1)] +
                                      [(cond2_label+'-'+str(i), robjects.IntVector(c))
                                       for i,c in enumerate(cond2)])
    data_frame = robjects.DataFrame(data_frame_contents)
    data_frame.rownames = transcript_names    
    conds = robjects.StrVector([cond1_label for x in cond1] + [cond2_label for x in cond2]).factor()

    ## DESeq full analysis
    deseq = rpackages.importr('DESeq')
    cds = deseq.newCountDataSet(data_frame, conds)
    cds = deseq.estimateSizeFactors(cds)
    cds = deseq.estimateVarianceFunctions(cds,method=method)
    res = deseq.nbinomTest(cds, cond1_label, cond2_label)

    ## Replace (unique) gene IDs by (not unique) gene names
    try:
        if assembly_id:
            print "Translate gene IDs to gene names..."
            res_ids = list(res[0])
            gene_names = {}; gid = None; l = None
            #assem = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies/" + str(assembly_id) + ".gtf")
            (assem,headers) = urllib.urlretrieve("http://bbcftools.vital-it.ch/genrep/nr_assemblies/" + str(76) + ".gtf")
            with open(assem) as assembly:
                for l in [line for line in assembly.readlines() if (line!='' and line!='\n')]:
                    if l.find("gene_id")!=-1:
                        gid = l.split("gene_id")[1].split(";")[0].strip(" \"")
                    if gid in res_ids:
                        idx = res_ids.index(gid)
                        gene_names[idx] = l.split("gene_name")[1].split(";")[0].strip(" \"")
            urllib.urlcleanup()
            for i in range(len(res_ids)):
                if i not in gene_names.keys(): gene_names[i] = res_ids[i]
            data_frame = res.rx(2)
            for i in range(3,len(res)+1):
                data_frame = data_frame.cbind(res.rx(i))
            data_frame.rownames = [gene_names[i] for i in range(len(gene_names))]
            res = data_frame
    except: print "Failed while changing names"

    ## MA-plot
    try:
        if maplot:
            print "MAplot"
            MAplot(res, mode="normal")
            ex.add("MAplot.png")
    except: print "Failed in MA-plot"

    if output:
        result_filename = output
    else:
        result_filename = unique_filename_in()

    res.to_csvfile(result_filename)
    return result_filename

def rnaseq_workflow(job, lims_path="rnaseq", via="lsf", job_or_dict="job", maplot=None):
    """Run RNASeq inference according to *job_info*.

    Whatever script calls this function should have looked up the job
    description from HTSStation.  The job description from HTSStation
    is returned as JSON, but should have been changed into a
    bbcflib.frontend.Job object (or a dictionary of the same form),
    which is what ``workflow`` expects.
    The object must have the fields:

        * ``assembly_id``, which ``path_to_bowtie_index`` uses to get
          a path from GenRep.

        * ``groups``, a dictionary of the conditions, and the
          information to fetch the files corresponding to those
          conditions from the DAF LIMS.

    If ``assembly_id`` is a string, not an integer, than it is
    interpreted as a path to the bowtie index, and no lookup to GenRep
    is necessary.

    Likewise, if the values of ``groups`` are lists of
    strings, then they are interpreted as a list of filenames and no
    lookup is done on the DAF LIMS.

    Whatever script calls workflow needs to pass back the JSON string
    it returns in some sensible way.  For the usual HTSStation
    frontend, this just means printing it to stdout.
    """
    
    """ Groups as given by the frontend is not that useful. Pull it apart
    into more useful pieces. """
    names = {}
    runs = {}
    controls = {}
    
    if job_or_dict == "job":
        groups = job.groups
        assembly_id = job.assembly_id
    elif job_or_dict == "dict":
        groups = job["groups"]
        assembly_id = job["assembly_id"]

    for i,v in groups.iteritems():
        names[i] = str(v['name'])
        runs[i] = v['runs'].values()
        controls[i] = v['control']

    M = MiniLIMS(lims_path)
    with execution(M) as ex:    #,remote_working_directory="/scratch/cluster/monthly/jdelafon/rnaseq"
        print "Current working directory:", os.getcwd()
        bowtie_index = path_to_bowtie_index(ex, assembly_id)
        print "Bowtie index:", bowtie_index
        
        """ gene_labels is a list whose ith entry is a string giving
        the name of the gene assigned id i. The ids are arbitrary.
        exon_mapping is a list whose ith entry is the integer id in
        gene_labels exon i maps to."""
        (gene_labels,exon_mapping) = fetch_transcript_mapping(ex, assembly_id)
        exon_mapping = numpy.array(exon_mapping)
        
        fastq_files = fetch_read_files(ex, runs)
        print "FastQ files:", [os.path.basename(f[0]['path']) for f in fastq_files.values()]
        bam_files = align_reads(ex, bowtie_index, fastq_files, via=via)
        print "Reads aligned."

        #bam1 = ex.use("bam1"); bam2 = ex.use("bam2")
        #index1 = ex.use("index1"); index2 = ex.use("index2")
        #bam_files = {1:[bam1], 2:[bam2]}
        #os.rename(index1, bam1+".bai"); os.rename(index2, bam2+".bai")
        
        print "Stats..."
        stats = {}
        for group_id,files in bam_files.iteritems():
            for i,f in enumerate(files):
                run_info = runs[group_id][i]
                if isinstance(run_info, str):
                    run_description = 'Group %d, run %d' % (group_id, i)
                else:
                    run_description = 'groupd %d, facility %s, machine %s, run %d, lane %d' % \
                        (group_id, run_info['facility'], run_info['machine'],
                         run_info['run'], run_info['lane'])
                stats[run_description] = bamstats(ex, f)
        stats_plots = plot_stats(ex, stats, '/archive/epfl/bbcf/share/')
        ex.add(stats_plots,
               description="Plot of alignment statistics (PDF)")

        # All the bam_files were created against the same index, so
        # they all have the same header in the same order.  I can take
        # the list of exons from just the first one and use it for all
        # of them.

        exons = exons_labels(bam_files[bam_files.keys()[0]][0])

        exon_pileups = {}
        gene_pileups = {}
        for condition,files in bam_files.iteritems():
            gene_pileups[condition] = []
            exon_pileups[condition] = []
            for f in files:
                exon_pileup = pileup_file(f, exons)
                exon_pileups[condition].append(exon_pileup)
                gene_pileup = numpy.zeros(len(gene_labels))
                for i,c in enumerate(exon_pileup):
                    gene_pileup[exon_mapping[i]] += c
                gene_pileups[condition].append(gene_pileup)
                
        futures = {}
        for (c1,c2) in pairs_to_test(controls):
            if len(runs[c1]) + len(runs[c2]) > 2:
                method = "normal";
            else:
                method = "blind";
            print "Inference..."
            output = None
            
            ## csvfile = inference(names[c1], gene_pileups[c1], 
            ##                               names[c2], gene_pileups[c2],
            ##                               gene_labels, 
            ##                               method, assembly_id,output, maplot)
            ## print csvfile
            ## a = ex.add(csvfile)
            
            try:
                csvfile = external_deseq.nonblocking(ex,
                                          names[c1], gene_pileups[c1], 
                                          names[c2], gene_pileups[c2],
                                          gene_labels, 
                                          method, assembly_id, output, maplot,
                                          via=via)
                print "Wait for results..."
                a = ex.add(csvfile.wait())
                print a
            except: print "Failed in Inference"
            
        ##     futures[(c1,c2)] = [external_deseq.nonblocking(ex,
        ##                                   names[c1], exon_pileups[c1], 
        ##                                   names[c2], exon_pileups[c2],
        ##                                   [x[0].split("|")[0]+"|"+x[0].split("|")[1] for x in exons],
        ##                                   method, assembly_id, output, maplot,
        ##                                   via=via)
        ##                         ,external_deseq.nonblocking(ex,
        ##                                   names[c1], gene_pileups[c1], 
        ##                                   names[c2], gene_pileups[c2],
        ##                                   gene_labels, 
        ##                                   method, assembly_id, output, maplot,
        ##                                   via=via)
        ##                         ]
        ## print "Results....."
        ## #pudb.set_trace()
        ## for c,f in futures.iteritems():
        ##     ex.add(f[0].wait(), description="Comparison of exons in conditions '%s' and '%s' (CSV)" % (names[c[0]], names[c[1]]))
        ##     ex.add(f[1].wait(), description="Comparison of genes in conditions '%s' and '%s' (CSV)" % (names[c[0]], names[c[1]]))
            
    return results_to_json(M, ex.id)

"""
======================
Module: bbcflib.rnaseq
======================

Methods of the bbcflib's RNA-seq worflow. The main function is rnaseq_workflow().
"""

# Built-in modules #
import os, sys, cPickle, json, pysam, urllib, math, time, csv

# Internal modules #
from bbcflib.genrep import GenRep
from bbcflib.common import timer, rstring, writecols

# Other modules #
import numpy

################################################################################

def lsqnonneg(C, d, x0=None, tol=None, itmax_factor=3):
    """Linear least squares with nonnegativity constraints (NNLS), based on MATLAB's lsqnonneg function.

    (x,resnorm,res) = lsqnonneg(C,d) returns
    * the vector x that minimizes norm(d-C*x) subject to x >= 0, C and d must be real
    * the norm of residuals *resnorm*
    * the residuals *res*

    References: Lawson, C.L. and R.J. Hanson, Solving Least-Squares Problems, Prentice-Hall, Chapter 23, p. 161, 1974.
    http://code.google.com/p/diffusion-mri/source/browse/trunk/Python/lsqnonneg.py?spec=svn17&r=17
    """
    eps = 2.22e-16    # from matlab
    def norm1(x):
        return abs(x).sum().max()

    def msize(x, dim):
        s = x.shape
        if dim >= len(s): return 1
        else: return s[dim]

    if tol is None: tol = 10*eps*norm1(C)*(max(C.shape)+1)
    C = numpy.asarray(C)
    (m,n) = C.shape
    P = numpy.zeros(n)
    Z = numpy.arange(1, n+1)
    if x0 is None: x=P
    else:
        if any(x0 < 0): x=P
        else: x=x0
    ZZ=Z
    resid = d - numpy.dot(C, x)
    w = numpy.dot(C.T, resid)
    outeriter=0; it=0
    itmax=itmax_factor*n
    exitflag=1

    # outer loop to put variables into set to hold positive coefficients
    while numpy.any(Z) and numpy.any(w[ZZ-1] > tol):
        outeriter += 1
        t = w[ZZ-1].argmax()
        t = ZZ[t]
        P[t-1]=t
        Z[t-1]=0
        PP = numpy.where(P != 0)[0]+1
        ZZ = numpy.where(Z != 0)[0]+1
        CP = numpy.zeros(C.shape)
        CP[:, PP-1] = C[:, PP-1]
        CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))
        z=numpy.dot(numpy.linalg.pinv(CP), d)
        z[ZZ-1] = numpy.zeros((msize(ZZ,1), msize(ZZ,0)))

        # inner loop to remove elements from the positve set which no longer belong
        while numpy.any(z[PP-1] <= tol):
            it += 1
            if it > itmax:
                max_error = z[PP-1].max()
                raise Exception('Exiting: Iteration count (=%d) exceeded\n Try raising the \
                                 tolerance tol. (max_error=%d)' % (it, max_error))
            QQ = numpy.where((z <= tol) & (P != 0))[0]
            alpha = min(x[QQ]/(x[QQ] - z[QQ]))
            x = x + alpha*(z-x)
            ij = numpy.where((abs(x) < tol) & (P <> 0))[0]+1
            Z[ij-1] = ij
            P[ij-1] = numpy.zeros(max(ij.shape))
            PP = numpy.where(P != 0)[0]+1
            ZZ = numpy.where(Z != 0)[0]+1
            CP[:, PP-1] = C[:, PP-1]
            CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))
            z=numpy.dot(numpy.linalg.pinv(CP), d)
            z[ZZ-1] = numpy.zeros((msize(ZZ,1), msize(ZZ,0)))
        x = z
        resid = d - numpy.dot(C, x)
        w = numpy.dot(C.T, resid)
    return (x, sum(resid*resid), resid)

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
            mapping = cPickle.load(pickle_file)
        print "Mapping found in", os.path.abspath(path_or_assembly_id)
        return mapping
    else:
        grep_root = '/db/genrep'
        grep = GenRep(url='http://bbcftools.vital-it.ch/genrep/',root=grep_root)
        assembly = grep.assembly(path_or_assembly_id)
        mdfive = assembly.md5
        mappings_path = os.path.join(grep_root,'nr_assemblies/exons_pickle/')+mdfive+".pickle"
        with open(mappings_path, 'rb') as pickle_file:
            mapping = cPickle.load(pickle_file)
        print "Mapping for assembly", mdfive, "found in /db/"
        return mapping

def exons_labels(bamfile):
    """List of the exons labels (SN: 'reference sequence name'; LN: 'sequence length') in *bamfile*."""
    sam = pysam.Samfile(bamfile, 'rb')
    labels = [(t['SN'],t['LN']) for t in sam.header['SQ']]
    sam.close()
    return labels

def pileup_file(bamfile, exons):
    """Return a dictionary {exon ID: count}, the pileup of *exons* from *bamfile*."""
    # exons = sam.references #tuple
    counts = []
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
        counts.append(c.n)
        c.n = 0
    sam.close()
    return counts

def translate_gene_ids(fc_ids, dictionary):
    """Replace (unique) gene IDs by (not unique) gene names.
    *fc_ids* is a dict {gene_id: whatever}
    *dictionary* is a dict {gene_id: gene_name} """
    # Maybe better use regexp here
    names = []
    for s in fc_ids.keys():
        start = s.find("ENS")
        if start != -1:
            end = s.split("ENS")[1].find("|")
            gene_id = "ENS" + s.split("ENS")[1].split("|")[0]
            names.append(s.replace(gene_id, dictionary.get(gene_id,gene_id)))
        else: names.append(s)
    fc_names = dict(zip(names,fc_ids.values()))
    return fc_names

def estimate_size_factors(counts):
    """The total number of reads may differ between conditions or replicates.
    This treatment makes different count sets being comparable. If rows are features
    and columns are different conditions/replicates, each column is divided by the
    geometric mean of the rows.
    The median of these ratios is used as the size factor for this column. Each
    column is divided by its size factor in order to normalize the data. Size
    factors may be used for further variance sabilization.

    * *counts* is an array of counts, each line representing a transcript, each
    column a different run.
    """
    numpy.seterr(divide='ignore') # Divisions by zero counts
    counts = numpy.array(counts)
    geo_means = numpy.exp(numpy.mean(numpy.log(counts), axis=0))
    mean = counts/geo_means
    size_factors = numpy.median(mean[:,geo_means>0], axis=1)
    res = counts.transpose()/size_factors
    res = res.transpose()
    return res, size_factors

def save_results(ex, data, conditions=[], name='counts'):
    """Save results in a CSV file, one line per feature, one column per run.
    *data* is a dictionary {ID: [counts in each condition]}
    """
    conditions_s = '%s, '*(len(conditions)-1)+'%s.'
    output = rstring()
    with open(output,"wb") as f:
        c = csv.writer(f, delimiter='\t')
        c.writerow(["id"]+conditions)
        for k,v in data.iteritems():
            c.writerow([k]+list(v))
    conditions = tuple(conditions)
    ex.add(output, description="csv:Expression level of "+name+" in samples "+conditions_s % conditions)
    print name+": Done successfully."

@timer
def genes_expression(gene_ids, exon_mapping, dexons):
    """Get gene counts from exons counts."""
    ncond = len(dexons.iteritems().next()[1])
    z = numpy.array([[0.]*ncond]*len(gene_ids))
    dgenes = dict(zip(gene_ids,z))
    for e,c in dexons.iteritems():
        e = e.split('|')[0]
        if exon_mapping.get(e):
            dgenes[exon_mapping[e][1]] += c
    return dgenes

@timer
def transcripts_expression(gene_ids, transcript_mapping, trans_in_gene, exons_in_trans, dexons):
    """Get transcript counts from exon counts."""
    ncond = len(dexons.iteritems().next()[1])
    z = numpy.array([[0.]*ncond]*len(transcript_mapping.keys()))
    dtrans = dict(zip(transcript_mapping.keys(),z))
    for g in gene_ids:
        if trans_in_gene.get(g):
            tg = trans_in_gene[g]
            # Get all exons in the gene
            eg = set()
            for t in tg:
                if exons_in_trans.get(t):
                    eg = eg.union(set(exons_in_trans[t]))
            # Create the correspondance matrix
            M = numpy.zeros((len(eg),len(tg)))
            cexons = []; ctrans = []
            for c in range(ncond):
                cexons.append( numpy.zeros(len(eg)) )
            for i,e in enumerate(eg):
                for j,t in enumerate(tg):
                    if exons_in_trans.get(t) and e in exons_in_trans[t]:
                        M[i,j] = 1
                # Retrieve exon counts
                if dexons.get(e) is not None:
                    for c in range(ncond):
                        cexons[c][i] += dexons[e][c]
            cexons = numpy.array(cexons)
            # Compute transcript counts
            for c in range(ncond):
                x, resnorm, res = lsqnonneg(M,cexons[c])
                ctrans.append(x)
            # Store results in a dict
            for k,t in enumerate(tg):
                if dtrans.get(t) is not None:
                    for c in range(ncond):
                        dtrans[t][c] = ctrans[c][k]
    return dtrans

@timer
def rnaseq_workflow(ex, job, assembly, bam_files, target=["genes"], via="lsf", output=None):
    """
    Main function of the workflow.

    * *ex*: the bein's execution Id.
    * *job*: a Job object (or a dictionary of the same form) as returned from HTSStation's frontend.
    * *assembly*: the assembly Id of the species, string or int (e.g. 'hg19' or 76).
    * *bam_files*: a complicated dictionary as returned by mapseq.get_bam_wig_files.
    * *target*: a string or array of strings indicating the features you want to compare.
    Targets can be 'genes', 'transcripts', or 'exons'. E.g. ['genes','transcripts'], or 'genes'.
    (This part of the workflow may fail if there are too few reads for DESeq to estimate
    variances amongst exons.)
    * *via*: 'local' or 'lsf'
    * *output*: alternative name for output file. Otherwise it is random.

    To do: -use lsf -pass rpkm, start, end and gene name to output -control with known refseqs -sort by ratios
    """
    group_names={}
    assembly_id = job.assembly_id
    groups = job.groups
    for gid,group in groups.iteritems():
        group_names[gid] = str(group['name']) # group_names = {gid: name}
    if isinstance(target,str): target=[target]

    # All the bam_files were created against the same index, so
    # they all have the same header in the same order.  I can take
    # the list of exons from just the first one and use it for all of them.
    # Format: ('exonID|geneID|start|end|strand', length)
    exons = exons_labels(bam_files[groups.keys()[0]][groups.values()[0]['runs'].keys()[0]]['bam'])

    """ Build pileups from bam files """
    print "Build pileups"
    exon_pileups = {}; conditions = []
    for gid,files in bam_files.iteritems():
        for rid,f in files.iteritems():
            cond = group_names[gid]+'.'+str(rid)
            exon_pileups[cond] = []
            conditions.append(cond)
        for rid,f in files.iteritems():
            cond = group_names[gid]+'.'+str(rid)
            exon_pileup = pileup_file(f['bam'], exons)
            exon_pileups[cond] = exon_pileup #{cond1.run1: {pileup}, cond1.run2: {pileup}...}

    #assembly_id = "../temp/nice_features/nice_mappings" # testing code
    print "Load mappings"
    mappings = fetch_mappings(assembly_id)
    """ [0] gene_ids is a list of gene ID's
        [1] gene_names is a dict {gene ID: gene name}
        [2] transcript_mapping is a dictionary {transcript ID: gene ID}
        [3] exon_mapping is a dictionary {exon ID: ([transcript IDs], gene ID)}
        [4] trans_in_gene is a dict {gene ID: [IDs of the transcripts it contains]}
        [5] exons_in_trans is a dict {transcript ID: [IDs of the exons it contains]} """
    (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans) = mappings

    """ Get counts for exons from bam files """
    print "Normalize"
    res = numpy.array([exon_pileups[cond] for cond in conditions], dtype=numpy.float32)
    res, size_factors = estimate_size_factors(res)
    print "Size factors:", size_factors
    for i in range(len(res.ravel())):
        # if zero counts, add 1 for further comparisons
        if res.flat[i]==0: res.flat[i] += 1.0

    """ Extract information from bam headers """
    exon_lengths=[]; exon_ids=[]; genes=[]; starts=[]; ends=[]
    for e in exons:
        exon, gene, start, end, strand = e[0].split('|')
        exon_lengths.append(e[1])
        exon_ids.append(exon)
        starts.append(start)
        ends.append(end)
        genes.append(gene)

    dexons = dict(zip(exon_ids,zip(*res)))
    if "exons" in target:
        save_results(ex, translate_gene_ids(dexons, gene_names), conditions=conditions, name="EXONS")

    """ Get counts for genes from exons """
    if "genes" in target:
        dgenes = genes_expression(gene_ids, exon_mapping, dexons)
        save_results(ex, translate_gene_ids(dgenes, gene_names), conditions=conditions, name="GENES")

    """ Get counts for the transcripts from exons, using pseudo-inverse """
    if "transcripts" in target:
        dtrans = transcripts_expression(gene_ids, transcript_mapping, trans_in_gene, exons_in_trans, dexons)
        save_results(ex, translate_gene_ids(dtrans, gene_names), conditions=conditions, name="TRANSCRIPTS")

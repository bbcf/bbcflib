"""
======================
Module: bbcflib.rnaseq
======================

Methods of the bbcflib's RNA-seq worflow. The main function is rnaseq_workflow().

From a BAM file produced by an alignement on the *exonome*, gets counts of reads
on the exons, add them to get counts on genes, and uses least-squares to infer
counts on transcripts, avoiding to map on either genome or transcriptome.
Note that the result counts on transcripts are approximate.
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

    ``(x,resnorm,res) = lsqnonneg(C,d)`` returns
    * the vector *x* that minimizes norm(d-C*x) subject to x >= 0
    * the norm of residuals *resnorm*
    * the residuals *res*

    :param x0: Initial point for x.
    :param tol: Tolerance to determine when the error is small enough.
    :param itmax_factor: Maximum number of iterations.

    :rtype x: numpy array
    :rtype resnorm: float
    :rtype res: numpy array

    Reference: Lawson, C.L. and R.J. Hanson, Solving Least-Squares Problems, Prentice-Hall, Chapter 23, p. 161, 1974.
    http://diffusion-mri.googlecode.com/svn/trunk/Python/lsqnonneg.py
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
    Z = ZZ = numpy.arange(1, n+1)
    if x0 is None: x=P
    else:
        if any(x0 < 0): x=P
        else: x=x0
    resid = d - numpy.dot(C, x)
    w = numpy.dot(C.T, resid)
    outeriter=0; it=0
    itmax=itmax_factor*n

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
    #return (x, sum(resid*resid), resid) # norm2
    return (x, sum(numpy.absolute(resid)), resid) # norm1

def fetch_mappings(path_or_assembly_id):
    """Given an assembly ID, return a tuple
    (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans)
    [0] gene_ids is a list of gene ID's
    [1] gene_names is a dict {gene ID: gene name}
    [2] transcript_mapping is a dictionary {transcript ID: gene ID}
    [3] exon_mapping is a dictionary {exon ID: ([transcript IDs], gene ID)}
    [4] trans_in_gene is a dict {gene ID: [IDs of the transcripts it contains]}
    [5] exons_in_trans is a dict {transcript ID: [IDs of the exons it contains]}

    :param path_or_assembly_id: can be a numeric or nominal ID for GenRep
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
    try: sam = pysam.Samfile(bamfile, 'rb')
    except ValueError: sam = pysam.Samfile(bamfile,'r')
    labels = [(t['SN'],t['LN']) for t in sam.header['SQ']]
    sam.close()
    return labels

def pileup_file(bamfile, exons):
    """Return a dictionary {exon ID: count}, the pileup of *exons* from *bamfile*."""
    # exons = sam.references #tuple
    # put it together with exons_labels? It opens the file twice.
    counts = []
    try: sam = pysam.Samfile(bamfile, 'rb')
    except ValueError: sam = pysam.Samfile(bamfile,'r')

    class Counter(object):
        def __init__(self):
            self.n = 0
        def __call__(self, alignment):
            self.n += 1

    c = Counter()
    for exon in exons:
        sam.fetch(exon[0], 0, exon[1], callback=c) #(exon_name,0,exon_length,Counter())
        #The callback (c.n += 1) is executed for each alignment in a region
        counts.append(c.n)
        c.n = 0
    sam.close()
    return counts

def save_results(ex, cols, conditions, header=[], desc='counts'):
    """Save results in a CSV file, one line per feature, one column per run.
    :param data: a dictionary {ID: [counts in each condition]}
    """
    conditions_s = '%s, '*(len(conditions)-1)+'%s.'
    output = rstring()
    writecols(output,cols,header=header, sep="\t")
    conditions = tuple(conditions)
    ex.add(output, description="csv:Expression level of "+desc+" in sample(s) "+conditions_s % conditions)
    print desc+": Done successfully."

@timer
def genes_expression(exons_data, exon_to_gene, ncond):
    """Get gene counts from exons counts."""
    genes = list(set(exons_data[1]))
    z  = [numpy.zeros(ncond,dtype=numpy.float_)]*len(genes)
    zz = [numpy.zeros(ncond,dtype=numpy.float_)]*len(genes)
    gcounts = dict(zip(genes,z))
    grpkms = dict(zip(genes,zz))
    for e,c in zip(exons_data[0],zip(*exons_data[2:ncond+ncond+2])):
        gcounts[exon_to_gene[e]] += c[:ncond]
        grpkms[exon_to_gene[e]] += c[ncond:]
    return gcounts,grpkms

@timer
def transcripts_expression(exons_data, trans_in_gene, exons_in_trans, ncond, method="nnls"):
    """Get transcript counts from exon counts.
    :param method: "nnls" or "pinv", respectively non-negative least-squares and pseudoinverse.
    """
    print "Method:", method
    genes = list(set(exons_data[1]))
    transcripts = []
    for g in genes:
        transcripts.extend(trans_in_gene.get(g,[]))
    transcripts = list(set(transcripts))
    z = numpy.zeros((len(transcripts),ncond))
    tcounts = dict(zip(transcripts,z))
    exons_counts = dict(zip( exons_data[0], zip(*exons_data[2:ncond+2])) )
    totalerror = 0; totalcount = 0; unknown = 0
    pinv = numpy.linalg.pinv; norm = numpy.linalg.norm; zeros = numpy.zeros;
    identity = numpy.identity; sum = numpy.sum; dot = numpy.dot; shape = numpy.shape
    for g in genes:
        if trans_in_gene.get(g): # if the gene is still in the Ensembl database
            # Get all transcripts in the gene
            tg = trans_in_gene[g]
            # Get all exons in the gene
            eg = set()
            for t in tg:
                if exons_in_trans.get(t):
                    eg = eg.union(set(exons_in_trans[t]))
            # Create the correspondance matrix
            M = zeros((len(eg),len(tg)))
            ecnts = zeros((ncond,len(eg)))
            tcnts = []
            for i,e in enumerate(eg):
                for j,t in enumerate(tg):
                    if exons_in_trans.get(t) and e in exons_in_trans[t]:
                        M[i,j] = 1
                # Retrieve exon counts
                if exons_counts.get(e) is not None:
                    for c in range(ncond):
                        ecnts[c][i] += exons_counts[e][c]
            # Compute transcript counts
            for c in range(ncond):
                if method == "pinv":
                    #-----------------------#
                    # Pseudo-inverse method #
                    #-----------------------#
                    x = dot(pinv(M),ecnts[c])
                    tcnts.append(x)
                    totalerror += norm(dot(M,pinv(M))-identity(shape(M)[0]) ,1) # norm 1
                    totalcount += sum(ecnts[c])
                if method == "nnls":
                    #-----------------------------------#
                    # Non-negative least squares method #
                    #-----------------------------------#
                    x, resnorm, res = lsqnonneg(M,ecnts[c],itmax_factor=5)
                    tcnts.append(x)
                    totalerror += resnorm
                    totalcount += sum(ecnts[c])
            # Store results in a dict *tcounts*
            for k,t in enumerate(tg):
                if tcounts.get(t) is not None:
                    for c in range(ncond):
                        tcounts[t][c] = tcnts[c][k]
        else:
            unknown += 1
    print "Evaluation of error for transcript counts:"
    print "\tUnknown transcripts for %d of %d genes (%.2f %%)" \
                       % (unknown, len(genes), 100*float(unknown)/float(len(genes)) )
    print "\tError in number of reads:",int(totalerror)
    print "\tTotal number of reads:",int(totalcount)
    print "\tRatio error/total:",float(totalerror)/totalcount
    return tcounts, totalerror

@timer
def rnaseq_workflow(ex, job, assembly, bam_files, pileup_level=["genes"], via="lsf", output=None):
    """
    Main function of the workflow.

    :param ex: the bein's execution Id.
    :param job: a Job object (or a dictionary of the same form) as returned from HTSStation's frontend.
    :param assembly: the assembly Id of the species, string or int (e.g. 'hg19' or 76).
    :param bam_files: a complicated dictionary as returned by mapseq.get_bam_wig_files.
    :param pileup_level: a string or array of strings indicating the features you want to compare.
    Targets can be 'genes', 'transcripts', or 'exons'.
    :param via: 'local' or 'lsf'
    :param output: alternative name for output file. Otherwise it is random.
    :rtype: None
    """
    group_names={}
    assembly_id = job.assembly_id
    groups = job.groups
    for gid,group in groups.iteritems():
        group_names[gid] = str(group['name']) # group_names = {gid: name}
    if isinstance(pileup_level,str): pileup_level=[pileup_level]

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
            exon_pileups[cond] = exon_pileup # {cond1.run1: {pileup}, cond1.run2: {pileup}...}

    print "Load mappings"
    #assembly_id = "../temp/nice_features/nice_mappings" # testing code
    mappings = fetch_mappings(assembly_id)
    """ [0] gene_ids is a list of gene ID's
        [1] gene_names is a dict {gene ID: gene name}
        [2] transcript_mapping is a dictionary {transcript ID: gene ID}
        [3] exon_mapping is a dictionary {exon ID: ([transcript IDs], gene ID)}
        [4] trans_in_gene is a dict {gene ID: [IDs of the transcripts it contains]}
        [5] exons_in_trans is a dict {transcript ID: [IDs of the exons it contains]} """
    (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans) = mappings

    """ Extract information from bam headers """
    exonsID=[]; genesID=[]; genesName=[]; starts=[]; ends=[]; exon_lengths=[]; exon_to_gene={};
    for e in exons:
        (exon, gene, start, end, strand) = e[0].split('|')
        starts.append(start)
        ends.append(end)
        exon_lengths.append(e[1])
        exonsID.append(exon)
        genesID.append(gene)
        exon_to_gene[exon] = gene
        genesName.append(gene_names.get(gene,gene))

    """ Treat data """
    print "Process data"
    starts = numpy.asarray(starts, dtype=numpy.int_)
    ends = numpy.asarray(ends, dtype=numpy.int_)
    counts = numpy.asarray([exon_pileups[cond] for cond in conditions], dtype=numpy.float_)
    for i in range(len(counts.ravel())):
        if counts.flat[i]==0: counts.flat[i] += 1.0 # if zero counts, add 1 for further comparisons
    rpkms = counts/(starts-ends)

    print "Get counts"
    exons_data = [exonsID,genesID]+list(counts)+list(rpkms)+[starts,ends,genesName]

    """ Print counts for exons """
    if "exons" in pileup_level:
        header = ["ExonID","GeneID"] + conditions*2 + ["Start","End","GeneName"]
        save_results(ex, exons_data, conditions, header=header, desc="EXONS")

    """ Get counts for genes from exons """
    if "genes" in pileup_level:
        header = ["GeneID"] + conditions*2 #+ ["Start","End","GeneName"]
        (gcounts, grpkms) = genes_expression(exons_data, exon_to_gene, len(conditions))
        gnames = gcounts.keys()
        genes_data = [gnames]+list(zip(*gcounts.values()))+list(zip(*grpkms.values()))
        save_results(ex, genes_data, conditions, header=header, desc="GENES")

    """ Get counts for the transcripts from exons, using pseudo-inverse """
    if "transcripts" in pileup_level:
        header = ["TranscriptID","GeneID"] + conditions #*2 + ["Start","End","GeneName"]
        (tcounts, error) = transcripts_expression(exons_data,
                                 trans_in_gene,exons_in_trans,len(conditions))
        tnames = tcounts.keys()
        gnames = [transcript_mapping[t] for t in tnames]
        (gnames, tnames) = zip(*sorted(zip(gnames,tnames))) # sort w.r.t. gene names
        trans_data = [tnames,gnames]+list(zip(*tcounts.values()))
        save_results(ex, trans_data, conditions, header=header, desc="TRANSCRIPTS")

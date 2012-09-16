"""
======================
Module: bbcflib.rnaseq
======================

Methods of the bbcflib's RNA-seq worflow. The main function is ``rnaseq_workflow()``,
and is usually called by bbcfutils' ``run_rnaseq.py``, e.g. command-line:

``python run_rnaseq.py -v local -c gapkowt.txt -d rnaseq -p transcripts,genes``

From a BAM file produced by an alignement on the **exonome**, gets counts of reads
on the exons, add them to get counts on genes, and uses least-squares to infer
counts on transcripts, avoiding to map on either genome or transcriptome.
Note that the resulting counts on transcripts are approximate.
The annotation of the bowtie index has to be consistent to that of the database (same assembly version).
"""

# Built-in modules #
import os, pysam, math

# Internal modules #
from bbcflib.common import writecols, set_file_descr, unique_filename_in
from bbcflib import mapseq, genrep
from bbcflib.bFlatMajor.common import cobble, sorted_stream
from bbcflib.btrack import track, FeatureStream, convert
from bein import program

# Other modules #
import numpy
from numpy import zeros, asarray, nonzero

numpy.set_printoptions(precision=3,suppress=True)
numpy.seterr(divide='ignore')


def positive(x):
    """Set to zero all negative components of an array."""
    for i in range(len(x)):
        if x[i] < 0 : x[i] = 0
    return x

def lsqnonneg(C, d, x0=None, tol=None, itmax_factor=3):
    """Linear least squares with nonnegativity constraints (NNLS), based on MATLAB's lsqnonneg function.

    ``(x,resnorm,res) = lsqnonneg(C,d)`` returns

    * the vector *x* that minimizes norm(d-Cx) subject to x >= 0
    * the norm of residuals *resnorm* = norm(d-Cx)^2
    * the residuals *res* = d-Cx

    :param x0: Initial point for x.
    :param tol: Tolerance to determine what is considered as close enough to zero.
    :param itmax_factor: Maximum number of iterations.

    :type C: nxm numpy array
    :type d: nx1 numpy array
    :type x0: mx1 numpy array
    :type tol: float
    :type itmax_factor: int
    :rtype: *x*: numpy array, *resnorm*: float, *res*: numpy array

    Reference: Lawson, C.L. and R.J. Hanson, Solving Least-Squares Problems, Prentice-Hall, Chapter 23, p. 161, 1974.
    http://diffusion-mri.googlecode.com/svn/trunk/Python/lsqnonneg.py
    """
    eps = 2.22e-16 # from Matlab
    def norm1(x):
        return abs(x).sum().max()

    def msize(x, dim):
        s = x.shape
        if dim >= len(s): return 1
        else: return s[dim]

    if tol is None: tol = 10*eps*norm1(C)*(max(C.shape)+1)
    C = asarray(C)
    (m,n) = C.shape
    P = numpy.zeros(n)            # set P of indices, empty for now
    Z = ZZ = numpy.arange(1, n+1) # set Z of indices, will be transferred to P
    if x0 is None or any(x0 < 0): x=P
    else: x=x0
    resid = d - numpy.dot(C, x)
    w = numpy.dot(C.T, resid) # gradient of (1/2)*||d-Cx||^2
    outeriter=0; it=0
    itmax=itmax_factor*n

    # outer loop to put variables into set to hold positive coefficients
    while numpy.any(Z) and numpy.any(w[ZZ-1] > tol): # if Z is empty or w_j<0 for all j, terminate.
        outeriter += 1
        t = w[ZZ-1].argmax() # find an index t s.t. w_t = max(w)
        t = ZZ[t]
        P[t-1]=t # move the index t from set Z to set P
        Z[t-1]=0 # Z becomes [0] if n=1
        PP = numpy.where(P != 0)[0]+1 # index of first non-zero element of P, +1 (-1 later)
        ZZ = numpy.where(Z != 0)[0]+1 # index of first non-zero element of Z, +1 (-1 later)
        CP = numpy.zeros(C.shape)
        CP[:, PP-1] = C[:, PP-1] # CP[:,j] is C[:,j] if j in P, or 0 if j in Z
        CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))
        z=numpy.dot(numpy.linalg.pinv(CP), d) # solution of least-squares min||d-CPx||
        if isinstance(ZZ,numpy.ndarray) and len(ZZ) == 0: # if Z = [0], ZZ = [] and makes it fail
            return (positive(z), sum(resid*resid), resid)
        else:
            z[ZZ-1] = numpy.zeros((msize(ZZ,1), msize(ZZ,0)))

        # inner loop to remove elements from the positve set which no longer belong
        while numpy.any(z[PP-1] <= tol): # if z_j>0 for all j, set x=z and return to outer loop
            it += 1
            if it > itmax:
                max_error = z[PP-1].max()
                raise Exception('Exiting: Iteration count (=%d) exceeded\n \
                      Try raising the tolerance tol. (max_error=%d)' % (it, max_error))
            QQ = numpy.where((z <= tol) & (P != 0))[0] # find an index q in P s.t. x_q/(x_q-z_q) is min with negative z_q.
            alpha = min(x[QQ]/(x[QQ] - z[QQ]))
            x = x + alpha*(z-x)
            ij = numpy.where((abs(x) < tol) & (P <> 0))[0]+1 # move from P to Z all indices j for which x_j=0
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

#@timer
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

def build_pileup(bamfile, exons):
    """From a BAM file, returns a dictionary of the form {feature ID: number of reads that mapped to it}.

    :param bamfile: name of a BAM file.
    :param exons: list of features of the type (chr,start,end,'exon_id|gene_id|gene_name',strand,phase).
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
    for e in exons:
        sam.fetch(e[0],e[1],e[2], callback=c) #(name,start,end,Counter())
        #The callback (c.n += 1) is executed for each alignment in a region
        counts[e[3].split('|')[0]] = c.n
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
            t=track(filename+'.bedGraph',info={'name':feature_type.lower()+"_"+group})
            t.make_header()
            convert(filename+'.sql',filename+'.bedGraph',mode='append')
            description = set_file_descr(feature_type.lower()+"_"+group+".bedGraph",
                                         step="pileup", type="bedGraph",
                                         groupId=group_ids[group], ucsc='1')
            ex.add(filename+'.bedGraph', description=description)
    print feature_type+": Done successfully."
    return output_tab

#@timer
def genes_expression(exons_data, gene_mapping, exon_mapping, ncond, nreads):
    """Get gene counts from exons counts.

    Returns two dictionaries, one for counts and one for rpkm, of the form ``{gene ID: score}``.

    :param exons_data: list of lists ``[exonsIDs,counts,rpkm,starts,ends,geneIDs,geneNames,strands,chrs]``.
    :param gene_mapping: dictionary ``{gene ID: (gene name,start,end,length,strand,chromosome)}``
    :param exon_mapping: dictionary ``{exon ID: ([transcript IDs],gene ID,start,end,strand,chromosome)}``
    :param ncond: number of samples.
    """
    genes = list(set(exons_data[-4]))
    z = numpy.zeros((len(genes),ncond))
    zz = numpy.zeros((len(genes),ncond))
    gene_counts = dict(zip(genes,z))
    gene_rpkm = dict(zip(genes,zz))
    round = numpy.round
    for e,c in zip(exons_data[0],zip(*exons_data[1:2*ncond+1])):
        g = exon_mapping[e][1]
        gene_counts[g] += round(c[:ncond],2)
    for g in genes:
        if gene_mapping.get(g):
            gene_rpkm[g] = to_rpkm(gene_counts[g], gene_mapping[g][3], nreads)
        else: gene_rpkm[g] = [0]*ncond
    return gene_counts, gene_rpkm

#@timer
def transcripts_expression(exons_data, exon_mapping, transcript_mapping, trans_in_gene, exons_in_trans, ncond, nreads):
    """Get transcript rpkms from exon rpkms.

    Returns two dictionaries, one for counts and one for rpkm, of the form ``{transcript ID: score}``.

    :param exons_data: list of lists ``[exonsIDs,counts,rpkm,starts,ends,geneIDs,geneNames,strands,chrs]``
    :param exon_mapping: dictionary ``{exon ID: ([transcript IDs],gene ID,start,end,strand,chromosome)}``
    :param transcript_mapping: dictionary ``{transcript ID: (gene ID,start,end,length,strand,chromosome)}``
    :param trans_in_gene: dictionary ``{gene ID: [transcript IDs it contains]}``
    :param exons_in_trans: dictionary ``{transcript ID: [exon IDs it contains]}``
    :param ncond: number of samples
    """
    genes = list(set(exons_data[-4]))
    transcripts=[]
    for g in genes:
        transcripts.extend(trans_in_gene.get(g,[]))
    transcripts = list(set(transcripts))
    z = numpy.zeros((len(transcripts),ncond))
    zz = numpy.zeros((len(transcripts),ncond))
    trans_counts = dict(zip(transcripts,z))
    trans_rpkm = dict(zip(transcripts,zz))
    exons_counts = dict(zip( exons_data[0], zip(*exons_data[1:ncond+1])) )
    unknown = 0
    pinv = numpy.linalg.pinv
    norm = numpy.linalg.norm
    for g in genes:
        if trans_in_gene.get(g): # if the gene is (still) in the Ensembl database
            # Get all transcripts in the gene
            tg = trans_in_gene[g]
            # Get all exons in the gene
            eg = set()
            for t in tg:
                if exons_in_trans.get(t):
                    eg = eg.union(set(exons_in_trans[t]))
            # Create the correspondance matrix
            M = zeros((len(eg),len(tg)))
            L = zeros((len(eg),len(tg)))
            ec = zeros((ncond,len(eg)))
            for i,e in enumerate(eg):
                for j,t in enumerate(tg):
                    if exons_in_trans.get(t) and e in exons_in_trans[t]:
                        M[i,j] = 1.
                        if exon_mapping.get(e) and transcript_mapping.get(t):
                            L[i,j] = float((exon_mapping[e][3]-exon_mapping[e][2])) / transcript_mapping[t][3]
                        else:
                            L[i,j] = 1./len(eg)
                # Retrieve exon scores
                if exons_counts.get(e) is not None:
                    for c in range(ncond):
                        ec[c][i] += exons_counts[e][c]
            # If only one transcript, special case for more precision
            if len(tg) == 1:
                for c in range(ncond):
                    trans_counts[t][c] = round(sum(ec[c]),2)
                continue
            # Compute transcript scores
            tc = []
            M = numpy.vstack((M,numpy.ones(M.shape[1]))) # add constraint |E| = |T|
            L = numpy.vstack((L,numpy.ones(M.shape[1])))
            N = M*L
            for c in range(ncond): # - counts
                Ec = ec[c]
                Ec = numpy.hstack((Ec,asarray(sum(Ec))))
                tol = 10*2.22e-16*norm(N,1)*(max(N.shape)+1)
                try: Tc, resnormc, resc = lsqnonneg(N,Ec,tol=100*tol)
                except: Tc = positive(numpy.dot(pinv(N),Ec))
                tc.append(Tc)
            # Store results in a dict *trans_counts*
            for k,t in enumerate(tg):
                if trans_counts.get(t) is not None:
                    for c in range(ncond):
                        trans_counts[t][c] = round(tc[c][k], 2)
        else:
            unknown += 1
    print "Unknown transcripts for %d of %d genes (%.2f %%)" \
           % (unknown, len(genes), 100*float(unknown)/float(len(genes)) )
    # Convert to rpkm
    for t in transcripts:
        trans_rpkm[t] = to_rpkm(trans_counts[t], transcript_mapping[t][3], nreads)
    return trans_counts, trans_rpkm

def estimate_size_factors(counts):
    """
    The total number of reads may be different between conditions or replicates.
    This treatment makes different count sets being comparable. If rows are features
    and columns are different conditions/replicates, each column is divided by the
    geometric mean of the rows.
    The median of these ratios is used as the size factor for this column. Size
    factors may be used for further variance sabilization.

    :param counts: an array of counts, each line representing a transcript, each
                   column a different sample.
    """
    numpy.seterr(divide='ignore')
    counts = numpy.asarray(counts)
    cnts = counts[nonzero(counts[:,0]*counts[:,1])] # none of the counts is zero
    loggeomeans = numpy.mean(numpy.log(cnts), 1)
    size_factors = numpy.exp(numpy.median(numpy.log(cnts).T - loggeomeans, 1))
    res = counts / size_factors
    print "Size factors:",size_factors
    return res, size_factors

def to_rpkm(counts, lengths, nreads):
    if isinstance(counts, numpy.ndarray):
        rpkm = 1000*(1e6*counts.T/nreads).T/lengths
    elif isinstance(counts, float) or isinstance(counts, int):
        rpkm = (1000*1e6*counts)/(nreads*lengths)
    elif isinstance(counts, dict):
        rpkm = {}
        for k,v in counts.iteritems():
            rpkm[k] = 1000*(1e6*v/nreads)/lengths
    return rpkm


#@timer
def rnaseq_workflow(ex, job, bam_files, pileup_level=["exons","genes","transcripts"], via="lsf"):
    """
    Main function of the workflow.

    :rtype: None
    :param ex: the bein's execution Id.
    :param job: a Job object (or a dictionary of the same form) as returned from HTSStation's frontend.
    :param bam_files: a complicated dictionary such as returned by mapseq.get_bam_wig_files.
    :param pileup_level: a string or array of strings indicating the features you want to compare.
                         Targets can be 'genes', 'transcripts', or 'exons'.
    :param via: 'local' or 'lsf'.
    """
    group_names={}; group_ids={}; conditions=[]
    assembly = genrep.Assembly(assembly=job.assembly_id,intype=2)
    groups = job.groups
    for gid,group in groups.iteritems():
        gname = str(group['name'])
        group_names[gid] = gname
        group_ids[gname] = gid
    if isinstance(pileup_level,str): pileup_level=[pileup_level]

    """ Define conditions """
    for gid,files in bam_files.iteritems():
        k = 0
        for rid,f in files.iteritems():
            k+=1
            cond = group_names[gid]+'.'+str(k)
            conditions.append(cond)
    ncond = len(conditions)

    print "Load mappings"
    """ [0] gene_mapping is a dict ``{gene ID: (gene name,start,end,length,strand,chromosome)}``
        [1] transcript_mapping is a dictionary ``{transcript ID: (gene ID,start,end,length,strand,chromosome)}``
        [2] exon_mapping is a dictionary ``{exon ID: ([transcript IDs],gene ID,start,end,strand,chromosome)}``
        [3] trans_in_gene is a dict ``{gene ID: [IDs of the transcripts it contains]}``
        [4] exons_in_trans is a dict ``{transcript ID: [IDs of the exons it contains]}`` """
    (gene_mapping, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans) = fetch_mappings(assembly)

    #Exon track: ('chr1',34553,35174,'ENSE00001727627|ENSG00000237613|FAM138A',-1,-1)
    #Exon map: (exon_id,gene_id,start,end,strand,chr)
    #exons = [(e,v[1:]) for e,v in exon_mapping.iteritems()]
    exons = list(assembly.exon_track())
    #Exon labels: ('exonID|geneID|start|end|strand|type', length)
    #exons = fetch_labels(bam_files[groups.keys()[0]][groups.values()[0]['runs'].keys()[0]]['bam'])

    """ Extract information from bam headers """
    exonsID=[]; genesID=[]; genesName=[]; starts=[]; ends=[]; strands=[]; badexons=[]
    for e in exons:
        chr,start,end,annot,strand,phase = e
        exon_id,gene_id,gene_name = annot.split('|')
        length = end-start
        if length > 1 and len(exon_id) > 0:
            starts.append(start)
            ends.append(end)
            strands.append(strand)
            exonsID.append(exon_id)
            genesID.append(gene_id)
            genesName.append(gene_name)
        else: badexons.append(e)
    [exons.remove(e) for e in badexons]

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
    #counts, sf = estimate_size_factors(counts)
    rpkm = to_rpkm(counts, ends-starts, nreads)

    hconds = ["counts."+c for c in conditions] + ["rpkm."+c for c in conditions]
    genesName, echr = zip(*[(x[0],x[-1]) for x in [gene_mapping.get(g,("NA",)*6) for g in genesID]])
    exons_data = [exonsID]+list(counts)+list(rpkm)+[starts,ends,genesID,genesName,strands,echr]
    exons_file = None; genes_file = None; trans_file = None

    """ Print counts for exons """
    if "exons" in pileup_level:
        print "Get scores of exons"
        exons_data = zip(*exons_data)
        header = ["ExonID"] + hconds + ["Start","End","GeneID","GeneName","Strand","Chromosome"]
        exons_data = zip(*exons_data)
        exons_file = save_results(ex, exons_data, conditions, group_ids, assembly, header=header, feature_type="EXONS")

    """ Get scores of genes from exons """
    if "genes" in pileup_level:
        print "Get scores of genes"
        (gcounts, grpkm) = genes_expression(exons_data, gene_mapping, exon_mapping, ncond, nreads)
        genesID = gcounts.keys()
        genes_data = [[g,gcounts[g],grpkm[g]]+list(gene_mapping.get(g,("NA",)*6)) for g in genesID]
        (genesID,gcounts,grpkm,gname,gstart,gend,glen,gstr,gchr) = zip(*genes_data)
        header = ["GeneID"] + hconds + ["Start","End","GeneName","Strand","Chromosome"]
        genes_data = [genesID]+list(zip(*gcounts))+list(zip(*grpkm))+[gstart,gend,gname,gstr,gchr]
        genes_file = save_results(ex, genes_data, conditions, group_ids, assembly, header=header, feature_type="GENES")

    """ Get scores of transcripts from exons, using non-negative least-squares """
    if "transcripts" in pileup_level:
        print "Get scores of transcripts"
        (tcounts,trpkm) = transcripts_expression(exons_data, exon_mapping,
                   transcript_mapping, trans_in_gene, exons_in_trans, ncond, nreads)
        transID = tcounts.keys()
        trans_data = [[t,tcounts[t],trpkm[t]]+list(transcript_mapping.get(t,("NA",)*6)) for t in transID]
        (transID,tcounts,trpkm,genesID,tstart,tend,tlen,tstr,tchr) = zip(*trans_data)
        genesName = [gene_mapping.get(g,("NA",)*6)[0] for g in genesID]
        header = ["TranscriptID"] + hconds + ["Start","End","GeneID","GeneName","Strand","Chromosome"]
        trans_data = [transID]+list(zip(*tcounts))+list(zip(*trpkm))+[tstart,tend,genesID,genesName,tstr,tchr]
        trans_file = save_results(ex, trans_data, conditions, group_ids, assembly, header=header, feature_type="TRANSCRIPTS")

    return {"exons":exons_file, "genes":genes_file, "transcripts":trans_file}


#-------------------------- DIFFERENTIAL ANALYSIS ----------------------------#

@program
def run_glm(rpath, data_file, options=[]):
    """Run *rpath*/negbin.test.R on *data_file*."""
    output_file = unique_filename_in()
    opts = ["-o",output_file]+options
    script_path = os.path.join(rpath,'negbin.test.R')
    return {'arguments': ["R","--slave","-f",script_path,"--args",data_file]+opts,
            'return_value': output_file}

def clean_before_deseq(filename):
    """Delete all lines of *filename* where counts are 0 in every run."""
    filename_clean = unique_filename_in()
    with open(filename,"rb") as f:
        with open(filename_clean,"wb") as g:
            header = f.readline()
            ncond = sum([h.split('.').count("counts") for h in header.split('\t')])
            header = '\t'.join(header.split('\t')[:1+ncond])+'\n'
            g.write(header)
            for line in f:
                l = line.split('\t')
                scores = l[1:1+ncond]
                if any([float(x) for x in scores]):
                    label = '|'.join([l[0],l[-3],l[-2],l[-1].strip('\r\n')]) # geneID|geneName|strand|chr
                    line = label + '\t' + '\t'.join(scores) + '\n'
                    g.write(line)
    return filename_clean

def clean_deseq_output(filename):
    """Delete all lines of *filename* with NA's everywhere, add 0.5 to zero scores
    before recalculating fold change, and remove rows' ids."""
    filename_clean = unique_filename_in()
    with open(filename,"rb") as f:
        with open(filename_clean,"wb") as g:
            contrast = f.readline().split('-')
            header = f.readline().split('\t')
            header[2] = 'baseMean'+contrast[0].strip()
            header[3] = 'baseMean'+contrast[1].strip()
            g.write('-'.join(contrast))
            g.write('\t'.join(header))
            for line in f:
                line = line.split("\t")[1:]
                if not (line[2]=="0" and line[3]=="0"):
                    meanA = float(line[2]) or 0.5
                    meanB = float(line[3]) or 0.5
                    fold = meanB/meanA
                    log2fold = math.log(fold,2)
                    line[2] = str(meanA); line[3] = str(meanB); line[4] = str(fold); line[5] = str(log2fold)
                    line = '\t'.join(line)
                    g.write(line)
    return filename_clean

def differential_analysis(ex, result, rpath, design=None, contrast=None):
    """For each file in *result*, launches an analysis of differential expression on the count
    values, and saves the output in the MiniLIMS.

    :param ex: the bein's execution.
    :param result: dictionary {type:filename} as returned by rnaseq_workflow.
    :param rpath: path to the R scripts (negbin.test.R).
    :param design: name of the file containing the design matrix.
    :param contrast: name of the file containing the contrasts matrix.
    """
    for type,res_file in result.iteritems():
        if res_file and rpath and os.path.exists(rpath):
            res_file = clean_before_deseq(res_file)
            options = ['-s','tab']
            if design: options += ['-d',design]
            if contrast: options += ['-c', contrast]
            try:
                glmfile = run_glm(ex, rpath, res_file, options)
                output_files = [f for f in os.listdir(ex.working_directory) if glmfile in f]
                for o in output_files:
                    desc = set_file_descr(type+"_differential"+o.split(glmfile)[1]+".txt", step='stats', type='txt')
                    o = clean_deseq_output(o)
                    ex.add(o, description=desc)
            except Exception as exc: print "Skipped differential analysis: %s \n" % exc

#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#

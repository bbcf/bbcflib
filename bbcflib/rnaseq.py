"""
======================
Module: bbcflib.rnaseq
======================

Methods of the bbcflib's RNA-seq worflow. The main function is ``rnaseq_workflow()``.

From a BAM file produced by an alignement on the genome or the exonome, gets counts of reads
on the exons, add them to get counts on genes, and uses least-squares to infer
counts on transcripts, avoiding to map on either genome or transcriptome.
Note that the resulting counts on transcripts are approximate.
"""

# Built-in modules #
import os, sys, pysam, math, itertools
from operator import itemgetter

# Internal modules #
from bbcflib.common import set_file_descr, unique_filename_in, cat, timer
from bbcflib.genrep import Assembly
from bbcflib.mapseq import add_and_index_bam, sam_to_bam, map_reads
from bbcflib.gfminer.common import cobble, sorted_stream, map_chromosomes, duplicate, apply
from bbcflib.track import track, FeatureStream, convert
from bein import program

# Other modules #
import numpy
from numpy import zeros, asarray, nonzero

numpy.set_printoptions(precision=3,suppress=True)
numpy.seterr(invalid='print')
numpy.seterr(divide='ignore')

test = False
#run_htsstation.py rnaseq -c config/gapkowt.txt -p genes --basepath ./

class Counter(object):
    def __init__(self):
        self.n = 0
    def __call__(self, alignment):
        NH = [1.0/t[1] for t in alignment.tags if t[0]=='NH']+[1]
        self.n += NH[0]

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
    `<http://diffusion-mri.googlecode.com/svn/trunk/Python/lsqnonneg.py>`_
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
        # inner loop to remove elements from the positive set which no longer belong
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

def fetch_mappings(assembly):
    """Given an assembly object, returns a tuple
    ``(gene_mapping, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans)``

    * [0] gene_mapping is a dict ``{gene_id: (gene name,start,end,length,chromosome)}``
    * [1] transcript_mapping is a dictionary ``{transcript_id: (gene_id,gene_name,start,end,length,chromosome)}``
    * [2] exon_mapping is a dictionary ``{exon_id: ([transcript_ids],gene_id,gene_name,start,end,chromosome)}``
    * [3] trans_in_gene is a dict ``{gene_id: [IDs of the transcripts it contains]}``
    * [4] exons_in_trans is a dict ``{transcript_id: [IDs of the exons it contains]}``

    :param path_or_assembly_id: can be a numeric or nominal ID for GenRep
    (e.g. 11, 76 or 'hg19' for H.Sapiens), or a path to a file containing a
    json object which is read to get the mapping.
    """
    if test and os.path.exists('../mappings/'):
        import json
        map_path = '../mappings/'
        gene_mapping = json.load(open(os.path.join(map_path,'gene_mapping.json')))
        exon_mapping = json.load(open(os.path.join(map_path,'exon_mapping.json')))
        transcript_mapping = json.load(open(os.path.join(map_path,'transcript_mapping.json')))
        trans_in_gene = json.load(open(os.path.join(map_path,'trans_in_gene.json')))
        exons_in_trans = json.load(open(os.path.join(map_path,'exons_in_trans.json')))
    else:
        gene_mapping = assembly.get_gene_mapping()
        transcript_mapping = assembly.get_transcript_mapping()
        exon_mapping = assembly.get_exon_mapping()
        exons_in_trans = assembly.get_exons_in_trans()
        trans_in_gene = assembly.get_trans_in_gene()
    mapping = (gene_mapping, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans)
    return mapping

@timer
def build_custom_pileup(bamfile, transcript_mapping=None, debugfile=sys.stderr):
    counts = {}
    try: sam = pysam.Samfile(bamfile, 'rb')
    except ValueError: sam = pysam.Samfile(bamfile,'r')
    c = Counter()
    for n,ref in enumerate(sam.references):
        start = 0
        end = transcript_mapping.get(ref.split('|')[0],(sam.lengths[n],)*6)[4]
        sam.fetch(ref, start, end, callback=c)
        counts[ref] = int(.5+c.n)
        c.n = 0
    sam.close()
    return counts

@timer
def build_pileup(bamfile, assembly, gene_mapping, exon_mapping, trans_in_gene, exons_in_trans, debugfile=sys.stderr):
    """From a BAM file, returns a dictionary of the form {feature_id: number of reads that mapped to it}.

    :param bamfile: name of a BAM file.
    :param exons: list of features of the type (exon_id,gene_id,gene_name,start,end,strand,chr).
    :type bamfile: string
    :type exons: list
    """
    class Counter(object):
        def __init__(self):
            self.n = 0 # total number of reads
            self.start = 0 # exon start
            self.counts = [] # vector of counts per non-zero position
        def __call__(self, alignment):
            NH = [1.0/t[1] for t in alignment.tags if t[0]=='NH']+[1]
            self.n += NH[0]
            try:
                self.counts[alignment.pos-self.start] += NH[0]
            except IndexError:
                pass # read overflows but is bigger than the exon, we don't care
        def remove_duplicates(self):
            """Fetches all reads mapped to a transcript, checks if there are mapping positions
            where the number of reads is more than N times the average for this gene.
            If there are, shrink them to the same level as the others."""
            Nsigma = 50 # arbitrary
            argmax = self.counts.index(max(self.counts))
            average = sum(self.counts[:argmax]+self.counts[argmax+1:]) / (len(self.counts)-1)
            #variance = average # Poisson approx
            #limit = round(average+Nsigma*(variance**0.5)+0.5) # avg + N x stdev
            limit = Nsigma*average
            self.counts = [c if c < limit else average for c in self.counts]

    counts = {}
    try: sam = pysam.Samfile(bamfile, 'rb')
    except ValueError: sam = pysam.Samfile(bamfile,'r')
    chromosomes = assembly.chrmeta.keys()
    ref_index = {}
    for ref in sam.references:
        ref_index[ref.split('|')[0]] = ref
    mapped_on = 'genome' if all([ref in chromosomes for ref in sam.references[:100]]) else 'exons'
    c = Counter()
    for g in gene_mapping.iterkeys():
        eg = set()
        for t in trans_in_gene[g]:
            eg.update(exons_in_trans[t])
        eg = sorted([(e,)+exon_mapping[e][3:] for e in eg], key=itemgetter(1,2))
        eg = cobble(FeatureStream(eg,fields=['name','start','end','strand','chr']))
        for e in eg:
            ex = e[0].split('|') # list of cobbled exons spanning the same interval
            origin = ex[-1]      # last representant of this exons set (arbitrary)
            if e[2]-e[1] <= 1: continue
            if mapped_on == 'genome':
                ref = e[-1]; start = e[1]; end = e[2]
            else: # mapped on the exonome
                ostart,oend = exon_mapping[origin][3:5]
                ref = ref_index.get(origin)
                start = e[1]-ostart; end = e[2]-ostart
            if not ref: continue
            try:
                #The callback is executed for each alignment in a region
                c.n = 0
                c.counts = [0]*(2*(end-start)+1) # 2x because of reads overflow (not strict ref interval)
                c.start = start
                sam.fetch(ref,start,end, callback=c)
                c.remove_duplicates()
            except ValueError,ve: # unknown reference
                debugfile.write(str(ve)); debugfile.flush()
            for exon in ex:
                counts[exon] = counts.get(exon,0) + c.n/float(len(ex))
    sam.close()
    return counts

@timer
def save_results(ex, lines, conditions, group_ids, assembly, header, feature_type='features', logfile=sys.stdout):
    """Save results in a tab-delimited file, one line per feature, one column per run.

    :param ex: bein's execution.
    :param lines: list of iterables, each element being a line to write in the output.
    :param conditions: list of sample identifiers, each given as *group.run_id*.
    :param group_ids: dictionary ``{group name: group_id}``.
    :param assembly: a genrep.Assembly object.
    :param header: list of strings, the column headers of the output file.
    :param feature_type: (str) the kind of feature of which you measure the expression.
    """
    # Tab-delimited output with all information
    ncond = len(conditions)
    output_tab = unique_filename_in()
    with open(output_tab,'wb') as f:
        f.write('\t'.join(header)+'\n')
        for l in lines:
            tid = str(l[0])
            counts = ["%d"%x for x in l[1:ncond+1]]
            norm = ["%.2f"%x for x in l[ncond+1:2*ncond+1]]
            rpkm = ["%.2f"%x for x in  l[2*ncond+1:3*ncond+1]]
            rest = [str(x) for x in l[3*ncond+1:]]
            f.write('\t'.join([tid]+counts+norm+rpkm+rest)+'\n')
    description = set_file_descr(feature_type.lower()+"_expression.tab", step="pileup", type="txt")
    ex.add(output_tab, description=description)
    # Create one track for each group
    if feature_type in ['GENES','EXONS']:
        cols = zip(*lines)
        groups = [c.split('.')[0] for c in conditions]
        start = cols[3*ncond+1]
        end = cols[3*ncond+2]
        chromosomes = cols[-1]
        rpkm = {}; output_sql = {}
        for i in range(ncond):
            group = conditions[i].split('.')[0]
            nruns = groups.count(group)
            # Average all replicates in the group
            rpkm[group] = asarray(rpkm.get(group,zeros(len(start)))) + asarray(cols[i+2*ncond+1]) / nruns
            output_sql[group] = output_sql.get(group,unique_filename_in())
        for group,filename in output_sql.iteritems():
            # SQL track - GDV
            tr = track(filename+'.sql', fields=['chr','start','end','score'], chrmeta=assembly.chrmeta)
            towrite = {}
            for n,c in enumerate(chromosomes):
                if c not in tr.chrmeta: continue
                towrite.setdefault(c,[]).append((int(start[n]),int(end[n]),rpkm[group][n]))
            for chrom, feats in towrite.iteritems():
                tr.write(cobble(sorted_stream(FeatureStream(feats, fields=['start','end','score']))),chrom=chrom,clip=True)
            description = set_file_descr(feature_type.lower()+"_"+group+".sql", step="pileup", type="sql",
                                         groupId=group_ids[group], gdv='1')
            ex.add(filename+'.sql', description=description)
            # bigWig track - UCSC
            try: # if the bigWig conversion program fails, the file is not created
                convert(filename+'.sql',filename+'.bw')
                description = set_file_descr(feature_type.lower()+"_"+group+".bw",
                                             step="pileup", type="bw", groupId=group_ids[group], ucsc='1')
                ex.add(filename+'.bw', description=description)
            except IOError: pass
    logfile.write("  %s: Done successfully.\n" % feature_type); logfile.flush()
    return os.path.abspath(output_tab)

@timer
def genes_expression(exon_ids, ecounts_matrix, gene_mapping, exon_mapping, ncond):
    """Get gene counts/rpkm from exon counts (sum).

    Returns a dictionary of the form ``{gene_id: scores}``.

    :param gene_mapping: dictionary ``{gene_id: (gene name,start,end,length,strand,chromosome)}``
    :param exon_mapping: dictionary ``{exon_id: ([transcript_id's],gene_id,gene_name,start,end,strand,chromosome)}``
    :param ncond: (int) number of samples.
    """
    gene_counts={}
    for e,counts in itertools.izip(exon_ids,ecounts_matrix):
        g = exon_mapping[e][1]
        gene_counts[g] = gene_counts.get(g,zeros(ncond)) + counts
    return gene_counts

@timer
def transcripts_expression(exon_ids, ecounts_matrix, exon_mapping, transcript_mapping, trans_in_gene, exons_in_trans,
                           ncond, debugfile=sys.stderr):
    """Get transcript rpkms from exon rpkms.

    Returns a dictionary of the form ``{transcript_id: score}``.

    :param exon_mapping: dictionary ``{exon_id: ([transcript_ids],gene_id,gene_name,start,end,strand,chromosome)}``
    :param transcript_mapping: dictionary ``{transcript_id: (gene_id,gene_name,start,end,length,strand,chromosome)}``
    :param trans_in_gene: dictionary ``{gene_id: [transcript_ids it contains]}``
    :param exons_in_trans: dictionary ``{transcript_id: [exon_ids it contains]}``
    :param ncond: (int) number of samples.
    """
    trans_counts={}
    exons_counts={}; genes=[];
    for e,counts in itertools.izip(exon_ids,ecounts_matrix):
        exons_counts[e] = counts
        genes.append(exon_mapping[e][1])
    del exon_ids, ecounts_matrix
    genes = set(genes)
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
                    trans_counts[t] = zeros(ncond)
            # Create the correspondance matrix
            M = zeros((len(eg),len(tg)))
            L = zeros((len(eg),len(tg)))
            ec = zeros((ncond,len(eg)))
            for i,e in enumerate(eg):
                for j,t in enumerate(tg):
                    if exons_in_trans.get(t) and e in exons_in_trans[t]:
                        M[i,j] = 1.
                        if exon_mapping.get(e) and transcript_mapping.get(t):
                            L[i,j] = float((exon_mapping[e][4]-exon_mapping[e][3])) / transcript_mapping[t][4]
                        else:
                            L[i,j] = 1./len(eg)
                # Retrieve exon scores
                if exons_counts.get(e) is not None:
                    for c in range(ncond):
                        ec[c][i] += exons_counts[e][c]
            # If only one transcript, special case for more precision
            if len(tg) == 1:
                for c in range(ncond):
                    trans_counts[t][c] = sum(ec[c])
                continue
            # Compute transcript scores
            tc = []
            M = numpy.vstack((M,numpy.ones(M.shape[1]))) # add constraint |E| = |T|
            L = numpy.vstack((L,numpy.ones(M.shape[1])))
            N = M*L
            for c in range(ncond):
                Ec = ec[c]
                Ec = numpy.hstack((Ec,asarray(sum(Ec))))
                tol = 10*2.22e-16*norm(N,1)*(max(N.shape)+1)
                try: Tc, resnormc, resc = lsqnonneg(N,Ec,tol=100*tol)
                except: Tc = positive(numpy.dot(pinv(N),Ec))
                tc.append(Tc)
            # Store results in a dict *trans_counts*
            for k,t in enumerate(tg):
                for c in range(ncond):
                    trans_counts[t][c] = tc[c][k]
        else:
            unknown += 1
    if unknown != 0:
        debugfile.write("\tUnknown transcripts for %d of %d genes (%.2f %%)" \
              % (unknown, len(genes), 100*float(unknown)/float(len(genes))) )
    return trans_counts

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
    counts = numpy.asarray(counts)
    cnts = counts[nonzero(numpy.prod(counts,1))] # none of the counts is zero
    logcnt = numpy.log(cnts)
    loggeomeans = numpy.mean(logcnt, 1)
    size_factors = numpy.exp(numpy.median(logcnt.T - loggeomeans, 1))
    res = counts / size_factors
    return res, size_factors

def to_rpkm(counts, lengths):
    """Divides read counts by the transcript length.
    *counts* and *lengths* are numpy arrays."""
    rpkm = 1000.*counts/(lengths[:,numpy.newaxis]) # transpose *lengths* without copying
    return rpkm

def norm_and_format(counts,lengths,map,map_idx,ids=None):
    """Normalize, compute RPKM values and format array lines as they will be printed."""
    if isinstance(counts,dict):
        counts_matrix = asarray(counts.values())
        ids = counts.iterkeys()
    else: counts_matrix = counts
    norm_matrix, sf = estimate_size_factors(counts_matrix)
    rpkm_matrix = to_rpkm(norm_matrix, lengths)
    map_nitems = len(map.iteritems().next())
    data = [(t,) + tuple(counts_matrix[k]) + tuple(norm_matrix[k]) + tuple(rpkm_matrix[k])
                 + itemgetter(*map_idx)(map.get(t,("NA",)*map_nitems))
            for k,t in enumerate(ids)]
    data = sorted(data, key=itemgetter(-1,-map_nitems,-map_nitems+1,-2)) # sort wrt. chr,start,end,strand
    return data

@timer
def rnaseq_workflow(ex, job, assembly=None,
                    pileup_level=["exons","genes","transcripts"], via="lsf",
                    rpath=None, junctions=None, unmapped=None,
                    logfile=sys.stdout, debugfile=sys.stderr):
    """Main function of the workflow.

    :rtype: None
    :param ex: the bein's execution Id.
    :param job: a Frontend.Job object (or a dictionary of the same form).
    :param assembly: a genrep.Assembly object
    :param bam_files: a dictionary such as returned by mapseq.get_bam_wig_files.
    :param rpath: (str) path to the R executable.
    :param junctions: (bool) whether to search for splice junctions using SOAPsplice. [False]
    :param unmapped: (bool) whether to remap to the transcriptome reads that did not map the genome. [False]
    :param via: (str) send job via 'local' or 'lsf'. ["lsf"]
    """
    if test:
        via = 'local'
        logfile = sys.stdout
        debugfile = sys.stderr
        repo_rpath = '/home/jdelafon/repos/bbcfutils/R/'
        if os.path.exists(repo_rpath): rpath = repo_rpath
    group_names={}; group_ids={}; conditions=[]
    if assembly is None:
        assembly = Assembly(assembly=job.assembly_id)
    groups = job.groups
    if len(groups)==0: sys.exit("No groups/runs were given.")
    for gid,group in groups.iteritems():
        gname = str(group['name'])
        group_names[gid] = gname
        group_ids[gname] = gid
    if isinstance(pileup_level,str): pileup_level=[pileup_level]

    """ Define conditions as 'group_name.run_id' """
    for gid,files in job.files.iteritems():
        k = 0
        for rid,f in files.iteritems():
            k+=1
            cond = group_names[gid]+'.'+str(k)
            conditions.append(cond)
    ncond = len(conditions)

    # If the reads were aligned on transcriptome (maybe custom), do that and skip the rest
    if hasattr(assembly,"fasta_origin") or assembly.intype == 2:
        tmap = {}
        if assembly.intype==2: # usual transcriptome
            try:
                tmap = assembly.get_transcript_mapping()
                ftype = "Transcripts"
                header = ["TranscriptID"]
            except:
                pass
        if hasattr(assembly,"fasta_origin"): # build custom transcriptome
            firstbam = job.files.itervalues().next().itervalues().next()['bam']
            firstbamtrack = track(firstbam,format='bam')
            for c,meta in firstbamtrack.chrmeta.iteritems():
                tmap[c] = ('','',0,meta['length'],meta['length'],0,'') #(gene_id,gene_name,start,end,length,strand,chr)
            ftype = "Custom"
            header = ["CustomID"]
        logfile.write("* Build pileups\n"); logfile.flush()
        pileups={}
        for gid,files in job.files.iteritems():
            k = 0
            for rid,f in files.iteritems():
                k+=1
                cond = group_names[gid]+'.'+str(k)
                pileup = build_custom_pileup(f['bam'],tmap,debugfile)
                pileups[cond] = pileup.values()
                logfile.write("  ....Pileup %s done\n" % cond); logfile.flush()
        ids = pileup.keys()
        counts = asarray([pileups[cond] for cond in conditions], dtype=numpy.float_).T
        del pileup, pileups
        tcounts={}
        for k,t in enumerate(ids):
            c = counts[k]
            if sum(c) != 0 and t in tmap:
                tcounts[t] = c
        lengths = asarray([tmap[t][3]-tmap[t][2] for t in tcounts.iterkeys()])
        trans_data = norm_and_format(tcounts,lengths,tmap,(2,3,0,1,5,6))
        header += ["counts."+c for c in conditions] + ["norm."+c for c in conditions] + ["rpkm."+c for c in conditions]
        header += ["Start","End","GeneID","GeneName","Strand","Chromosome"]
        trans_file = save_results(ex,trans_data,conditions,group_ids,assembly,
                                  header=header,feature_type=ftype,logfile=logfile)
        differential_analysis(ex, trans_data, header, rpath, logfile, debugfile, feature_type=ftype.lower(),via=via)
        return 0

    logfile.write("* Load mappings\n"); logfile.flush()
    #    [0] gene_mapping is a dict ``{gene_id: (gene_name,start,end,length,strand,chromosome)}``
    #    [1] transcript_mapping is a dictionary ``{transcript_id: (gene_id,gene_name,start,end,length,strand,chromosome)}``
    #    [2] exon_mapping is a dictionary ``{exon_id: ([transcript_ids],gene_id,gene_name,start,end,strand,chromosome)}``
    #    [3] trans_in_gene is a dict ``{gene_id: [IDs of the transcripts it contains]}``
    #    [4] exons_in_trans is a dict ``{transcript_id: [IDs of the exons it contains]}``
    (gene_mapping, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans) = fetch_mappings(assembly)
    if len(exon_mapping) == 0 or len(gene_mapping) == 0:
        raise ValueError("No genes found for this genome. Abort.")

    # Map remaining reads to transcriptome
    unmapped_fastq = {}
    if unmapped:
        logfile.write("* Align unmapped reads on transcriptome\n"); logfile.flush()
        try:
            unmapped_fastq,additionals = align_unmapped(ex,job,assembly,group_names,
                                                        exon_mapping,transcript_mapping,exons_in_trans,via)
        except Exception, error:
            debugfile.write(str(error)); debugfile.flush()

    # Find splice junctions
    if junctions:
        logfile.write("* Search for splice junctions\n"); logfile.flush()
        try:
            find_junctions(ex,job,assembly,logfile=logfile,debugfile=debugfile,via=via)
        except Exception, error:
            debugfile.write(str(error)); debugfile.flush()

    # Build exon pileups from bam files
    logfile.write("* Build pileups\n"); logfile.flush()
    exon_pileups={}
    for gid,files in job.files.iteritems():
        k = 0
        for rid,f in files.iteritems():
            k+=1
            cond = group_names[gid]+'.'+str(k)
            exon_pileup = build_pileup(f['bam'],assembly,gene_mapping,exon_mapping,trans_in_gene,exons_in_trans,debugfile)
            if unmapped and cond in unmapped_fastq and cond in additionals:
                for a,x in additionals[cond].iteritems():
                    if exon_pileup.get(a):
                        exon_pileup[a] += x
                additionals.pop(cond)
            exon_pileups[cond] = exon_pileup.values()
            logfile.write("  ....Pileup %s done\n" % cond); logfile.flush()
    exon_ids = asarray(exon_pileup.keys()) # same for all conds
    del exon_pileup

    # Arrange exon counts in a matrix
    ecounts_matrix = asarray([exon_pileups[cond] for cond in conditions], dtype=numpy.float_).T
    nonzero_exons = nonzero(numpy.sum(ecounts_matrix,1)) # indices of non-zero lines
    ecounts_matrix = ecounts_matrix[nonzero_exons]
    exon_ids = exon_ids[nonzero_exons]
    del exon_pileups

    hconds = ["counts."+c for c in conditions] + ["norm."+c for c in conditions] + ["rpkm."+c for c in conditions]
    exons_file = None; genes_file = None; trans_file = None

    # Print counts for exons
    if "exons" in pileup_level:
        logfile.write("* Get scores of exons\n"); logfile.flush()
        header = ["ExonID"] + hconds + ["Start","End","GeneID","GeneName","Strand","Chromosome"]
        lengths = asarray([exon_mapping[e][4]-exon_mapping[e][3] for e in exon_ids])
        exons_data = norm_and_format(ecounts_matrix,lengths,exon_mapping,(3,4,1,2,5,6),ids=exon_ids)
        exons_file = save_results(ex, exons_data, conditions, group_ids, assembly, header=header, feature_type="EXONS")
        differential_analysis(ex, exons_data, header, rpath, logfile, debugfile, feature_type='exons',via=via)
        del exons_data

    # Get scores of genes from exons
    if "genes" in pileup_level:
        logfile.write("* Get scores of genes\n"); logfile.flush()
        header = ["GeneID"] + hconds + ["Start","End","GeneName","Strand","Chromosome"]
        gcounts = genes_expression(exon_ids, ecounts_matrix, gene_mapping, exon_mapping, ncond)
        lengths = asarray([gene_mapping[g][3] for g in gcounts.iterkeys()])
        genes_data = norm_and_format(gcounts,lengths,gene_mapping,(1,2,0,4,5))
        genes_file = save_results(ex, genes_data, conditions, group_ids, assembly, header=header, feature_type="GENES")
        differential_analysis(ex, genes_data, header, rpath, logfile, debugfile, feature_type='genes',via=via)
        del genes_data

    # Get scores of transcripts from exons, using non-negative least-squares
    if "transcripts" in pileup_level:
        logfile.write("* Get scores of transcripts\n"); logfile.flush()
        header = ["TranscriptID"] + hconds + ["Start","End","GeneID","GeneName","Strand","Chromosome"]
        tcounts = transcripts_expression(exon_ids, ecounts_matrix, exon_mapping,
                   transcript_mapping, trans_in_gene, exons_in_trans, ncond, debugfile)
        lengths = asarray([transcript_mapping[t][4] for t in tcounts.iterkeys()])
        trans_data = norm_and_format(tcounts,lengths,transcript_mapping,(2,3,0,1,5,6))
        trans_file = save_results(ex, trans_data, conditions, group_ids, assembly, header=header, feature_type="TRANSCRIPTS")
        differential_analysis(ex, trans_data, header, rpath, logfile, debugfile, feature_type='transcripts',via=via)
        del trans_data
    return 0


#-------------------------- DIFFERENTIAL ANALYSIS ----------------------------#

@program
def run_glm(rpath, data_file, options=[]):
    """Run *rpath*/negbin.test.R on *data_file*."""
    output_file = unique_filename_in()
    opts = ["-o",output_file]+options
    script_path = os.path.join(rpath,'negbin.test.R')
    return {'arguments': ["R","--slave","-f",script_path,"--args",data_file]+opts,
            'return_value': output_file}

def clean_before_deseq(data, header, keep=0.6):
    """Delete all lines of *filename* where counts are 0 in every run."""
    norm = 'counts'
    if norm == 'counts': w = 0
    elif norm == 'norm': w = 1
    elif norm == 'rpkm': w = 2
    filename_clean = unique_filename_in()
    ncond = sum([h.split('.').count("counts") for h in header]) # a regexp would be better
    if ncond >1:
        rownames = asarray(['%s|%s|%s|%s' % (x[0],x[-3],x[-2],x[-1]) for x in data])
        M = asarray([x[1+w*ncond:1+(w+1)*ncond] for x in data]) # *norm* columns
        colnames = header[0:1]+header[1:1+ncond] # 'counts' column names, always
        # Remove 40% lowest counts
        sums = numpy.sum(M,1)
        filter = asarray([x[1] for x in sorted(zip(sums,range(len(sums))))])
        filter = filter[:len(filter)/keep]
        M = M[filter]
        rownames = rownames[filter]
        # Create the input tab file
        with open(filename_clean,"wb") as g:
            header = '\t'.join(colnames)+'\n' # ID & *norm* columns
            g.write(header)
            for i,scores in enumerate(M):
                if any(scores):
                    line = rownames[i] + '\t' + '\t'.join([str(x) for x in scores]) + '\n'
                    g.write(line)
    return filename_clean, ncond

def clean_deseq_output(filename):
    """Delete all lines of *filename* with NA's everywhere, add 0.5 to zero scores
    before recalculating the fold change, and remove row ids. Return the new file name."""
    filename_clean = unique_filename_in()
    with open(filename,"rb") as f:
        with open(filename_clean,"wb") as g:
            contrast = f.readline().split('-')
            header = f.readline().split('\t')
            header[2] = 'baseMean'+contrast[0].strip()
            header[3] = 'baseMean'+contrast[1].strip()
            g.write('-'.join(contrast))
            g.write('\t'.join(header[:10]))
            for line in f:
                line = line.split("\t")[1:11] # 1:: remove row ids
                if not (line[2]=="0" and line[3]=="0"):
                    line[1] = '%.2f'%float(line[1])
                    meanA = float(line[2]) or 0.5
                    meanB = float(line[3]) or 0.5
                    fold = meanB/meanA
                    log2fold = math.log(fold,2)
                    line[2]='%.2f'%meanA; line[3]='%.2f'%meanB; line[4]='%.2f'%fold; line[5]='%.2f'%log2fold
                    line = '\t'.join(line)
                    g.write(line)
    return filename_clean

@timer
def differential_analysis(ex, data, header, rpath, logfile, debugfile, feature_type, via='lsf'):
    """For each file in *result*, launch an analysis of differential expression on the count
    values, and saves the output in the MiniLIMS.

    :param rpath: path to the R scripts ("negbin.test.R").
    :param via: (str) send job via 'local' or 'lsf'. ["lsf"]
    """
    if rpath and os.path.exists(rpath):
        res_file, ncond = clean_before_deseq(data, header)
        if ncond < 2:
            logfile.write("  Skipped differential analysis: less than two groups.\n"); logfile.flush()
        else:
            logfile.write("  Differential analysis\n"); logfile.flush()
            debugfile.write("DE: R path: '%s'\n" % rpath); logfile.flush()
        options = ['-s','tab']
        try:
            glmfile = run_glm.nonblocking(ex, rpath, res_file, options, via=via).wait()
        except Exception as exc:
            logfile.write("  Skipped differential analysis: %s \n" % exc); logfile.flush()
            return
        if not glmfile:
            logfile.write("  ....Empty file.\n"); logfile.flush()
            return
        output_files = [f for f in os.listdir(ex.working_directory) if glmfile in f]
        for o in output_files:
            desc = set_file_descr(feature_type+"_differential"+o.split(glmfile)[1]+".txt", step='stats', type='txt')
            o = clean_deseq_output(o)
            ex.add(o, description=desc)
        logfile.write("  ....done.\n"); logfile.flush()


#-------------------------- SPLICE JUNCTIONS SEARCH ----------------------------#

@program
def soapsplice(unmapped_R1, unmapped_R2, index, output=None, path_to_soapsplice=None, options={}):
    """Bind 'soapsplice'. Return a text file containing the list of junctions.

    :param unmapped_R1: (str) path to the fastq file containing the 'left' reads.
    :param unmapped_R2: (str) path to the fastq file containing the 'right' reads.
    :param index: (str) path to the SOAPsplice index.
    :param output: (str) output file name.
    :param path_to_soapsplice: (str) path to the SOAPsplice executable.
        If not specified, the program must be in your $PATH.
    :param options: (dict) SOAPsplice options, given as {opt: value}.
    :rtype: str

    Main options::

    -p: number of threads, <= 20. [1]
    -S: 1: forward strand, 2: reverse strand, 3: both. [3]
    -m: maximum mismatch for one-segment alignment, <= 5. [3]
    -g: maximum indel for one-segment alignment, <= 2. [2]
    -i: length of tail that can be ignored in one-segment alignment. [7]
    -t: longest gap between two segments in two-segment alignment. [500000]
    -a: shortest length of a segment in two-segment alignment. [8]
    -q: input quality type in FASTQ file (0: old Illumina, 1: Sanger). [0]
    -L: maximum distance between paired-end reads. [500000]
    -l: minimum distance between paired-end reads. [50]
    -I: insert length of paired-end reads.

    """
    if not output: output = unique_filename_in()
    path_to_soapsplice = path_to_soapsplice or 'soapsplice'
    args = [path_to_soapsplice,'-d',index,'-1',unmapped_R1,'-2',unmapped_R2,'-o',output,'-f','2']
    opts = []
    for k,v in options.iteritems(): opts.extend([str(k),str(v)])
    return {"arguments": args+opts, "return_value": output}

@timer
def find_junctions(ex,job,assembly,soapsplice_index=None,path_to_soapsplice=None,soapsplice_options={},
                   logfile=sys.stdout,debugfile=sys.stderr,via='lsf'):
    """
    Retrieve unmapped reads from a precedent mapping and runs SOAPsplice on them.
    Return the names of .bed and .sql tracks indicating the junctions positions, as well as
    of a bam file of the alignments attesting the junctions.

    :param soapsplice_index: (str) path to the SOAPsplice index.
    :param path_to_soapsplice: (str) specify the path to the program if it is not in your $PATH.
    :param soapsplice_options: (dict) SOAPsplice options, e.g. {'-m':2}.
    :rtype: str, str, str
    """
    assembly.set_index_path(intype=3)
    soapsplice_index = soapsplice_index or assembly.index_path
    soapsplice_options.update(job.options.get('soapsplice_options',{}))
    soapsplice_options.setdefault('-p',16) # number of threads
    soapsplice_options.setdefault('-q',1)  # Sanger format
    unmapped_fastq = {}
    for gid, group in job.groups.iteritems():
        unmapped_fastq[gid] = []
        for rid, run in group['runs'].iteritems():
            unmapped = job.files[gid][rid].get('unmapped_fastq')
            if not unmapped:
                logfile.write("No unmapped reads found for group %s, run %d. Skip.\n" % (gid,rid)); logfile.flush()
                continue
            elif not isinstance(unmapped,tuple):
                logfile.write("Pair-end reads required. Skip.\n"); logfile.flush()
                continue
            unmapped_fastq[gid].append(unmapped)
        if len(unmapped_fastq[gid]) == 0:
            continue
        R1 = cat(zip(*unmapped_fastq[gid])[0])
        R2 = cat(zip(*unmapped_fastq[gid])[1])
        future = soapsplice.nonblocking(ex,R1,R2,soapsplice_index,
                                        path_to_soapsplice=path_to_soapsplice,
                                        options=soapsplice_options,
                                        via=via, memory=8, threads=soapsplice_options['-p'])
        template = future.wait()
        if not template: return
        junc_file = template+'.junc'
        bed = convert_junc_file(junc_file,assembly)
        bed_descr = set_file_descr('junctions_%s.bed' % group['name'],
                                   groupId=gid,type='bed',step='1',ucsc=1)
        bam_descr = set_file_descr('junctions_%s.bam' % group['name'],
                                   groupId=gid,type='bam',step='1')
        sam = template+'.sam'
        try:
            bam = sam_to_bam(ex,sam,reheader=assembly.name)
            add_and_index_bam(ex, bam, description=bam_descr)
            ex.add(bam, description=bam_descr)
        except Exception, e:
            debugfile.write("%s\n(Qualities may be in the wrong format, try with '-q 0'.)\n" %str(e)); debugfile.flush()
        ex.add(bed, description=bed_descr)

def convert_junc_file(filename, assembly):
    """Convert a .junc SOAPsplice output file to sql and bed formats. Return the two respective file names.

    :param filename: (str) name of the .junc file to convert.
    :param assembly: genrep.Assembly object.
    """
    t = track(filename, format='txt', fields=['chr','start','end','strand','score'], chrmeta=assembly.chrmeta)
    stream = t.read()
    # Translate chromosome names
    s1 = map_chromosomes(stream, assembly.chromosomes)
    # Add junction IDs
    s2 = duplicate(s1,'strand','name')
    C = itertools.count()
    s3 = apply(s2,'name', lambda x: 'junction'+str(C.next()))
    # Convert to bed format
    outfile = unique_filename_in()
    bed = outfile + '.bed'
    out = track(bed, fields=s3.fields, chrmeta=assembly.chrmeta)
    out.write(s3)
    return bed


#-------------------------- UNMAPPED READS ----------------------------#

@timer
def align_unmapped( ex, job, assembly, group_names,
                    exon_mapping, transcript_mapping, exons_in_trans, via ):
    """
    Map reads that did not map to the exons to a collection of annotated transcripts,
    in order to add counts to pairs of exons involved in splicing junctions.

    Return a dictionary ``unmapped_fastq`` of the form ``{sample_name:bam_file}``,
    and a second ``additionals`` of the form ``{sample_name:{exon:additional_counts}}``.

    :param group_names: dict of the form ``{group_id: group_name}``.
    """
    assembly.set_index_path(intype=2)
    additionals = {}
    unmapped_fastq = {}
    refseq_path = assembly.index_path
    bwt2 = job.options.get("bowtie2",True)
    for gid, group in job.groups.iteritems():
        k = 0
        for rid, run in group['runs'].iteritems():
            k += 1
            cond = group_names[gid]+'.'+str(k)
            _fastq = job.files[gid][rid].get('unmapped_fastq')
            if _fastq and os.path.exists(refseq_path+".1.ebwt"):
                try:
                    _bam = map_reads( ex, _fastq, {}, refseq_path, bowtie_2=bwt2,
                                      remove_pcr_duplicates=False, via=via )['bam']
                except:
                    continue
                if _bam:
                    sam = pysam.Samfile(_bam)
                else:
                    continue
                additional = {}
                for read in sam:
                    t_id = sam.getrname(read.tid).split('|')[0]
                    if transcript_mapping.get(t_id) and exons_in_trans.get(t_id):
                        lag = 0
                        r_start = read.pos
                        r_end = r_start + read.rlen
                        E = exons_in_trans[t_id]
                        for e in E:
                            e_start, e_end = exon_mapping[e][3:5]
                            e_len = e_end-e_start
                            if r_start <= lag+e_len and lag <= r_end:
                                additional[e] = additional.get(e,0) + 0.5
                            lag += e_len
                additionals[cond] = additional
                unmapped_fastq[cond] = _fastq
                sam.close()
    return unmapped_fastq,additionals


#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#

"""
======================
Module: bbcflib.rnaseq
======================

Methods of the bbcflib's RNA-seq worflow. The main function is ``rnaseq_workflow()``,
and is usually called by bbcfutils' ``run_rnaseq.py``, e.g. command-line:

``python run_rnaseq.py -v lsf -c config_files/gapkowt.txt -d rnaseq -p transcripts,genes``
``python run_rnaseq.py -v lsf -c config_files/rnaseq2.txt -d rnaseq -p transcripts -m ../mapseq2``

From a BAM file produced by an alignement on the **exonome**, gets counts of reads
on the exons, add them to get counts on genes, and uses least-squares to infer
counts on transcripts, avoiding to map on either genome or transcriptome.
Note that the resulting counts on transcripts are approximate.
The annotation of the bowtie index has to be consistent to that of the database (same assembly version).
"""

# Built-in modules #
import os, pysam, math
from operator import itemgetter

# Internal modules #
from bbcflib.common import writecols, set_file_descr, unique_filename_in
from bbcflib import mapseq, genrep
import track
from bein import program

# Other modules #
import numpy
from numpy import zeros, asarray

numpy.set_printoptions(precision=3,suppress=True)
numpy.seterr(divide='ignore')


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
    P = numpy.zeros(n)
    Z = ZZ = numpy.arange(1, n+1)
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
        Z[t-1]=0
        PP = numpy.where(P != 0)[0]+1 # set P of indices
        ZZ = numpy.where(Z != 0)[0]+1 # set Z of indices
        CP = numpy.zeros(C.shape)
        CP[:, PP-1] = C[:, PP-1] # CP[:,j] is C[:,j] if j in P, or 0 if j in Z
        CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))
        z=numpy.dot(numpy.linalg.pinv(CP), d) # solution of least-squares min||d-CPx||
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

def fusion(X):
    """ Takes a track-like generator with items of the form (chromosome,start,end,score)
    and returns a merged track-like generator with summed scores.
    This is to avoid having overlapping coordinates of features from both DNA strands,
    which some genome browsers cannot handle for quantitative tracks.
    """
    x = X.next()
    c = x[0]
    last = x[1]
    last_was_alone = True
    last_had_same_end = False
    for y in X:
        if y[1] < x[2]:             # y intersects x
            last_had_same_end = False
            if y[1] >= last: # cheating, needed in case 3 features overlap, thanks to the annotation...
                if y[1] != last:    # y does not have same start as x
                    yield (c,last,y[1],x[3])
                if y[2] < x[2]:     # y is embedded in x
                    yield (c,y[1],y[2],x[3]+y[3])
                    last = y[2]
                elif y[2] == x[2]:  # y has same end as x
                    yield (c,y[1],y[2],x[3]+y[3])
                    x = y
                    last_had_same_end = True
                else:               # y exceeds x
                    yield (c,y[1],x[2],x[3]+y[3])
                    last = x[2]
                    x = y
            last_was_alone = False
        else:                       # y is outside of x
            if last_was_alone:
                yield x
            elif last_had_same_end:
                pass
            else:
                yield (c,last,x[2],x[3])
            x = y
            last = x[1]
            last_was_alone = True
            last_had_same_end = False
    if last_was_alone:
        yield x
    elif last_had_same_end:
        pass
    else:
        yield (c,last,x[2],x[3])

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
            rpkm[group] = asarray(rpkm.get(group,zeros(len(start)))) + asarray(cols[i+ncond+1]) / nruns
            output_sql[group] = output_sql.get(group,unique_filename_in())
        for group,filename in output_sql.iteritems():
            lines = zip(*[chromosomes,start,end,rpkm[group]])
            # SQL track
            with track.new(filename+'.sql') as t:
                t.chrmeta = assembly.chrmeta
                for chr in t.chrmeta:
                    goodlines = [l for l in lines if (l[3]!=0.0 and l[0]==chr)]
                    [lines.remove(l) for l in goodlines]
                    if goodlines:
                        goodlines.sort(key=lambda x: itemgetter(1,2)) # sort w.r.t start
                        goodlines = fusion(iter(goodlines))
                        for x in goodlines:
                            t.write(x[0],[(x[1],x[2],x[3])],fields=["start","end","score"])
            description = set_file_descr(feature_type.lower()+"_"+group+".sql", step="pileup", type="sql", \
                                         groupId=group_ids[group], gdv='1')
            ex.add(filename+'.sql', description=description)
            # UCSC-BED track
            track.convert(filename+'.sql',filename+'.bedGraph')
            description = set_file_descr(feature_type.lower()+"_"+group+".bedGraph", step="pileup", type="bedGraph", \
                                         groupId=group_ids[group], ucsc='1')
            ex.add(filename+'.bedGraph', description=description)
    print feature_type+": Done successfully."
    return output_tab

#@timer
def genes_expression(exons_data, exon_lengths, gene_mapping, exon_to_gene, ncond):
    """Get gene counts from exons counts.

    Returns two dictionaries, one for counts and one for rpkm, of the form ``{gene ID: score}``.

    :param exons_data: list of lists ``[exonsID, counts, rpkm, start, end, geneID, geneName]``.
    :param exon_lengths: dictionary ``{exon ID: length}``
    :param gene_mapping: dictionary ``{gene ID: (gene name,start,end,length,chromosome)}``
    :param exon_to_gene: dictionary ``{exon ID: gene ID}``.
    :param ncond: number of samples.
    """
    genes = list(set(exons_data[-3]))
    z = numpy.zeros((len(genes),ncond))
    zz = numpy.zeros((len(genes),ncond))
    gcounts = dict(zip(genes,z))
    grpkm = dict(zip(genes,zz))
    round = numpy.round
    for e,c in zip(exons_data[0],zip(*exons_data[1:2*ncond+1])):
        g = exon_to_gene[e]
        gcounts[g] += round(c[:ncond],2)
        try: # Let's avoid errors of that kind until the article is done
            ratio = exon_lengths[e]/gene_mapping[g][3]
            grpkm[g] += ratio*round(c[ncond:],2)
        except KeyError, ZeroDivisionError: pass
    return gcounts, grpkm

#@timer
def transcripts_expression(exons_data, exon_lengths, transcript_mapping, trans_in_gene, exons_in_trans, ncond):
    """Get transcript rpkms from exon rpkms.

    Returns two dictionaries, one for counts and one for rpkm, of the form ``{transcript ID: score}``.

    :param exons_data: list of lists ``[exonsID, counts, rpkm, start, end, geneID, geneName, chromosome]``
    :param exon_lengths: dictionary ``{exon ID: length}``
    :param transcript_mapping: dictionary ``{transcript ID: (gene ID, start, end, length, chromosome)}``
    :param trans_in_gene: dictionary ``{gene ID: [transcript IDs it contains]}``
    :param exons_in_trans: dictionary ``{transcript ID: [exon IDs it contains]}``
    :param ncond: number of samples
    """
    genes = list(set(exons_data[-3]))
    transcripts=[]
    for g in genes:
        transcripts.extend(trans_in_gene.get(g,[]))
    transcripts = list(set(transcripts))
    z = numpy.zeros((len(transcripts),ncond))
    zz = numpy.zeros((len(transcripts),ncond))
    trans_counts = dict(zip(transcripts,z))
    trans_rpk = dict(zip(transcripts,zz))
    exons_counts = dict(zip( exons_data[0], zip(*exons_data[1:ncond+1])) )
    exons_rpk = dict(zip( exons_data[0], zip(*exons_data[ncond+1:2*ncond+1])) )
    totalerror = 0; unknown = 0; #alltranscount=0; allexonscount=0;
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
            M = zeros((len(eg),len(tg))); L = zeros((len(eg),len(tg)))
            ec = zeros((ncond,len(eg))); er = zeros((ncond,len(eg)))
            for i,e in enumerate(eg):
                for j,t in enumerate(tg):
                    if exons_in_trans.get(t) and e in exons_in_trans[t]:
                        M[i,j] = 1.
                        if exon_lengths.get(e) and transcript_mapping.get(t):
                            L[i,j] = exon_lengths[e]/transcript_mapping[t][3]
                        else:
                            L[i,j] = 1./len(eg)
                # Retrieve exon scores
                if exons_counts.get(e) is not None:
                    for c in range(ncond):
                        ec[c][i] += exons_counts[e][c]
                        er[c][i] += exons_rpk[e][c]
            # Compute transcript scores
            tc = []; tr = []
            for c in range(ncond):  # - rpkm
                Er = er[c]
                tol = 10*2.22e-16*numpy.linalg.norm(M,1)*(max(M.shape)+1)
                try: Tr, resnormr, resr = lsqnonneg(M,Er,tol=100*tol)
                except: Tr = zeros(len(tg)) # Bad idea
                tr.append(Tr)
            #N = M
            M = numpy.vstack((M,numpy.ones(M.shape[1]))) # add constraint |E| = |T|
            L = numpy.vstack((L,numpy.ones(M.shape[1])))
            N = M*L
            for c in range(ncond): # - counts
                Ec = ec[c]
                Ec = numpy.hstack((Ec,asarray(sum(Ec))))
                tol = 10*2.22e-16*numpy.linalg.norm(N,1)*(max(N.shape)+1)
                try: Tc, resnormc, resc = lsqnonneg(N,Ec,tol=100*tol)
                except: Tc = zeros(len(tg)) # Bad idea
                tc.append(Tc)
                totalerror += math.sqrt(resnormc)
            # Store results in a dict *tcounts*/*trpk*
            for k,t in enumerate(tg):
                if trans_rpk.get(t) is not None:
                    for c in range(ncond):
                        trans_counts[t][c] = round(tc[c][k], 2)
                        trans_rpk[t][c] = round(tr[c][k], 2)
            # Testing
            #total_trans_c = sum([sum(trans_counts[t]) for t in tg]) or 0
            #total_exons_c = sum([sum(ec[c]) for c in range(ncond)]) or 0
            #alltranscount += total_trans_c
            #allexonscount += total_exons_c
        else:
            unknown += 1
    print "Unknown transcripts for %d of %d genes (%.2f %%)" \
           % (unknown, len(genes), 100*float(unknown)/float(len(genes)) )
    # Testing
    #try: print "\t Total transcript counts: %.2f, Total exon counts: %.2f, Ratio: %.2f" \
    #       % (alltranscount,allexonscount,alltranscount/allexonscount)
    #except ZeroDivisionError: pass
    #print "\t Total error (sum of resnorms):", totalerror
    return trans_counts, trans_rpk

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
    counts = numpy.array(counts)
    geo_means = numpy.exp(numpy.mean(numpy.log(counts), axis=0))
    mean = counts/geo_means
    size_factors = numpy.median(mean[:,geo_means>0], axis=1)
    res = counts.T/size_factors
    res = res.T
    print "Size factors:",size_factors
    return res, size_factors

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

    #Exon labels: ('exonID|geneID|start|end|strand', length)
    exons = fetch_labels(bam_files[groups.keys()[0]][groups.values()[0]['runs'].keys()[0]]['bam'])

    """ Extract information from bam headers """
    exonsID=[]; genesID=[]; genesName=[]; starts=[]; ends=[];
    exon_lengths={}; exon_to_gene={}; badexons=[]
    for e in exons:
        (exon, gene, start, end, strand) = e[0].split('|')
        start = int(start); end = int(end)
        if end-start>1:
            starts.append(start)
            ends.append(end)
            exon_lengths[exon] = float(e[1])
            exonsID.append(exon)
            genesID.append(gene)
            exon_to_gene[exon] = gene
        else: badexons.append(e)
    [exons.remove(e) for e in badexons]

    print "Load mappings"
    """ [0] gene_mapping is a dict ``{gene ID: (gene name,start,end,length,chromosome)}``
        [1] transcript_mapping is a dictionary ``{transcript ID: (gene ID,start,end,chromosome)}``
        [2] exon_mapping is a dictionary ``{exon ID: ([transcript IDs],gene ID,start,end,chromosome)}``
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
            if unmapped_fastq[cond] and os.path.exists(unmapped_fastq[cond]):
                print "Add splice junction reads for run %s.%s" % (gid,rid)
                assert os.path.exists(refseq_path+".1.ebwt"), "Refseq index not found: %s" % refseq_path+".1.ebwt"
                unmapped_bam[cond] = mapseq.map_reads(ex, unmapped_fastq[cond], {}, refseq_path, \
                      remove_pcr_duplicates=False, bwt_args=[], via=via)['bam']
                if unmapped_bam[cond]:
                    sam = pysam.Samfile(unmapped_bam[cond])
                else: sam = []
                additional = {}
                for read in sam:
                    t_id = sam.getrname(read.tid).split('|')[0]
                    if transcript_mapping.get(t_id):
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
    rpkm = 1000*(1e6*counts.T/nreads).T/(ends-starts)

    hconds = ["counts."+c for c in conditions] + ["rpkm."+c for c in conditions]
    genesName, echr = zip(*[(x[0],x[-1]) for x in [gene_mapping.get(g,("NA",)*5) for g in genesID]])
    exons_data = [exonsID]+list(counts)+list(rpkm)+[starts,ends,genesID,genesName,echr]
    exons_file = None; genes_file = None; trans_file = None

    """ Print counts for exons """
    if "exons" in pileup_level:
        print "Get scores of exons"
        exons_data = zip(*exons_data)
        exons_data = sorted(exons_data, key=itemgetter(7,4)) # sort w.r.t. chromosome, then start
        header = ["ExonID"] + hconds + ["Start","End","GeneID","GeneName","Chromosome"]
        exons_data = zip(*exons_data)
        exons_file = save_results(ex, exons_data, conditions, group_ids, assembly, header=header, feature_type="EXONS")

    """ Get scores of genes from exons """
    if "genes" in pileup_level:
        print "Get scores of genes"
        (gcounts, grpkm) = genes_expression(exons_data, exon_lengths, gene_mapping, exon_to_gene, len(conditions))
        genesID = gcounts.keys()
        genes_data = [[g,gcounts[g],grpkm[g]]+list(gene_mapping.get(g,("NA",)*5)) for g in genesID]
        genes_data = sorted(genes_data, key=itemgetter(7,4)) # sort w.r.t. chromosome, then start
        (genesID,gcounts,grpkm,gname,gstart,gend,glen,gchr) = zip(*genes_data)
        header = ["GeneID"] + hconds + ["Start","End","GeneName","Chromosome"]
        genes_data = [genesID]+list(zip(*gcounts))+list(zip(*grpkm))+[gstart,gend,gname,gchr]
        genes_file = save_results(ex, genes_data, conditions, group_ids, assembly, header=header, feature_type="GENES")

    """ Get scores of transcripts from exons, using non-negative least-squares """
    if "transcripts" in pileup_level:
        print "Get scores of transcripts"
        (tcounts, trpkm) = transcripts_expression(exons_data, exon_lengths,
                   transcript_mapping, trans_in_gene, exons_in_trans,len(conditions))
        transID = tcounts.keys()
        trans_data = [[t,tcounts[t],trpkm[t]]+list(transcript_mapping.get(t,("NA",)*5)) for t in transID]
        trans_data = sorted(trans_data, key=itemgetter(7,4)) # sort w.r.t. chromosome, then start
        (transID,tcounts,trpkm,genesID,tstart,tend,tlen,tchr) = zip(*trans_data)
        genesName = [gene_mapping.get(g,("NA",)*5)[0] for g in genesID]
        header = ["TranscriptID"] + hconds + ["Start","End","GeneID","GeneName","Chromosome"]
        trans_data = [transID]+list(zip(*tcounts))+list(zip(*trpkm))+[tstart,tend,genesID,genesName,tchr]
        trans_file = save_results(ex, trans_data, conditions, group_ids, assembly, header=header, feature_type="TRANSCRIPTS")

    return {"exons":exons_file, "genes":genes_file, "transcripts":trans_file}


#-------------------------- DIFFERENTIAL ANALYSIS ----------------------------#

@program
def run_glm(rpath, data_file, options=[]):
    """Run negbin.test.R"""
    output_file = unique_filename_in()
    options += ["-o",output_file]
    script_path = os.path.join(rpath,'negbin.test.R')
    return {'arguments': ["R","--slave","-f",script_path,"--args",data_file]+options,
            'return_value': output_file}

def clean_before_deseq(filename):
    """Delete all lines of *filename* where counts are 0 in every run."""
    filename_clean = unique_filename_in()
    with open(filename,"rb") as f:
        with open(filename_clean,"wb") as g:
            header = f.readline()
            g.write(header)
            ncond = sum([h.split('.').count("counts") for h in header.split('\t')])
            for line in f:
                goodline = line.split("\t")[1:]
                scores = [float(x) for x in goodline[:ncond]]
                if any(scores):
                    g.write(line)
    return filename_clean

def clean_deseq_output(filename):
    """Delete all lines of *filename* with NA's everywhere, and add 0.5 to zero scores
    before recalculating fold change."""
    filename_clean = unique_filename_in()
    with open(filename,"rb") as f:
        with open(filename_clean,"wb") as g:
            contrast = f.readline()
            header = f.readline()
            g.write(contrast+header)
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
            except: print "Skipped differential analysis"

#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#

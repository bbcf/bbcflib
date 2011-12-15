"""
======================
Module: bbcflib.rnaseq

Methods of the bbcflib's RNA-seq worflow. The main function is ``rnaseq_workflow()``,
and is usually called by bbcfutils' ``run_rnaseq.py``, e.g. command-line:

``python run_rnaseq.py -v lsf -c config_files/gapdh.txt -d rnaseq -p transcripts``
``python run_rnaseq.py -v lsf -c config_files/rnaseq.txt -d rnaseq -p genes -u -m /scratch/cluster/monthly/jdelafon/mapseq``

From a BAM file produced by an alignement on the **exonome**, gets counts of reads
on the exons, add them to get counts on genes, and uses least-squares to infer
counts on transcripts, avoiding to map on either genome or transcriptome.
Note that the resulting counts on transcripts are approximate.
"""

# Built-in modules #
import os, sys, pysam, math

# Internal modules #
from bbcflib.genrep import GenRep
from bbcflib.common import timer, writecols, set_file_descr
from bbcflib import mapseq
#from mapseq import bowtie, sam_to_bam, sort_bam, index_bam
import track

# Other modules #
import sqlite3
import numpy
import cPickle
from numpy import zeros, asarray

numpy.set_printoptions(precision=3,suppress=True)


def rstring(len=20):
    """Generate a random string of length *len* (usually for filenames)."""
    import string, random
    return "".join([random.choice(string.letters+string.digits) for x in range(len)])

def lsqnonneg(C, d, x0=None, tol=None, itmax_factor=3):
    """Linear least squares with nonnegativity constraints (NNLS), based on MATLAB's lsqnonneg function.

    ``(x,resnorm,res) = lsqnonneg(C,d)`` returns

    * the vector *x* that minimizes norm(d-Cx) subject to x >= 0
    * the norm of residuals *resnorm* = norm(d-Cx)^2
    * the residuals *res* = d-Cx

    :param x0: Initial point for x.
    :param tol: Tolerance to determine what is considered as zero.
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

def get_md5(assembly_id):
    grep_root = '/db/genrep'
    grep = GenRep(url='http://bbcftools.vital-it.ch/genrep/',root=grep_root)
    assembly = grep.assembly(assembly_id)
    mdfive = assembly.md5
    return mdfive

def get_chromosomes(assembly_id):
    grep_root = '/db/genrep'
    grep = GenRep(url='http://bbcftools.vital-it.ch/genrep/',root=grep_root)
    assembly = grep.assembly(assembly_id)
    chromosomes = assembly.chromosomes
    return chromosomes

@timer
def fetch_mappings(assembly_id, path_to_map=None):
    """Given an assembly ID, returns a tuple
    ``(gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans)``

    * [0] gene_mapping is a dict ``{gene ID: (gene name,start,end,chromosome)}``
    * [1] transcript_mapping is a dictionary ``{transcript ID: (gene ID,start,end,length,chromosome)}``
    * [2] exon_mapping is a dictionary ``{exon ID: ([trancript IDs],gene ID,start,end,chromosome)}``
    * [3] trans_in_gene is a dict ``{gene ID: [IDs of the transcripts it contains]}``
    * [4] exons_in_trans is a dict ``{transcript ID: [IDs of the exons it contains]}``

    :param path_or_assembly_id: can be a numeric or nominal ID for GenRep
    (e.g. 11, 76 or 'hg19' for H.Sapiens), or a path to a file containing a
    pickle object which is read to get the mapping.
    """
    mdfive = get_md5(assembly_id)

    # Connect to GTF database
    dbpath = "/db/genrep/nr_assemblies/annot_tracks/"+mdfive+".sql"
    if os.path.exists(dbpath):
        db = sqlite3.connect(dbpath, check_same_thread=False)
    else: raise IOError("Database not found in %s" % dbpath)

    # Get chromosome names from GenRep
    c = db.cursor()
    chromosomes = c.execute("""SELECT name FROM chrNames ORDER BY name;""").fetchall()
    chromosomes = [c[0] for c in chromosomes]

    def get_gene_mapping(db,chromosomes):
        """Return a dictionary {geneID: (geneName, start, end, chromosome)}"""
        c = db.cursor()
        gene_mapping = {}; sql=''
        for chr in chromosomes:
            sql += '''SELECT DISTINCT '%s',gene_id,gene_name,MIN(start),MAX(end) FROM '%s'
                     WHERE (type LIKE 'exon') GROUP BY gene_id UNION ''' % (chr,chr,)
        sql = sql[:-7]+';'
        sql_result = c.execute(sql)
        for chr,g,name,start,end in sql_result:
            gene_mapping[g] = (name,start,end,chr)
        return gene_mapping

    def get_transcript_mapping(db,chromosomes):
        """Return a dictionary ``{transcript ID: (gene ID,start,end,length)}``"""
        c = db.cursor()
        transcript_mapping = {}; sql=''
        lengths = {}; sql=''
        for chr in chromosomes:
            sql += '''SELECT transcript_id,sum(end-start) FROM '%s'
                      WHERE (type LIKE 'exon') GROUP BY transcript_id UNION ''' % (chr,)
        sql = sql[:-7]+';'
        sql_result = c.execute(sql)
        for t,l in sql_result:
            lengths[t] = l
        sql = ''
        for chr in chromosomes:
            sql += '''SELECT DISTINCT '%s',transcript_id,gene_id,MIN(start),MAX(end) FROM '%s'
                      WHERE (type LIKE 'exon') GROUP BY transcript_id UNION ''' % (chr,chr,)
        sql = sql[:-7]+';'
        sql_result = c.execute(sql)
        for chr,t,g,start,end in sql_result:
            transcript_mapping[t] = (g,start,end,lengths[t],chr)
        return transcript_mapping

    def get_exon_mapping(db,chromosomes):
        """Return a dictionary ``{exon ID: ([transcript IDs],gene ID,start,end)}``"""
        c = db.cursor()
        exon_mapping = {}; sql=''; T={}
        for chr in chromosomes:
            sql += '''SELECT DISTINCT exon_id,transcript_id from '%s'
                      WHERE (type LIKE 'exon') UNION ''' % (chr,)
        sql = sql[:-7]+';'
        sql_result = c.execute(sql)
        for e,t in sql_result:
            T.setdefault(e,[]).append(t)
        sql = ''
        for chr in chromosomes:
            sql += '''SELECT DISTINCT '%s',exon_id,gene_id,start,end FROM '%s'
                      WHERE (type LIKE 'exon') UNION ''' % (chr,chr,)
        sql = sql[:-7]+';'
        sql_result = c.execute(sql)
        for chr,e,g,start,end in sql_result:
            exon_mapping[e] = (T[e],g,start,end,chr)
        return exon_mapping

    def get_exons_in_trans(db,chromosomes):
        """Return a dictionary ``{transcript ID: list of exon IDs it contains}``"""
        c = db.cursor()
        exons_in_trans = {}; sql=''
        for chr in chromosomes:
            sql += '''SELECT DISTINCT transcript_id,exon_id from '%s'
                      WHERE (type LIKE 'exon') UNION ''' % (chr,)
        sql = sql[:-7]+';'
        sql_result = c.execute(sql)
        for t,e in sql_result:
            exons_in_trans.setdefault(t,[]).append(e)
        return exons_in_trans

    def get_trans_in_gene(db,chromosomes):
        """Return a dictionary ``{gene ID: list of transcript IDs it contains}``"""
        c = db.cursor()
        trans_in_gene = {}; sql=''
        for chr in chromosomes:
            sql += '''SELECT DISTINCT gene_id,transcript_id FROM '%s'
                      WHERE (type LIKE 'exon') UNION ''' % (chr,)
        sql = sql[:-7]+';'
        sql_result = c.execute(sql)
        for g,t in sql_result:
            trans_in_gene.setdefault(g,[]).append(t)
        return trans_in_gene

    path_to_map = "/scratch/cluster/monthly/jdelafon/temp/mappings"
    if os.path.exists(path_to_map):
        with open(path_to_map,"rb") as f:
            mapping = cPickle.load(f)
    else:
        gene_mapping = get_gene_mapping(db,chromosomes)
        transcript_mapping = get_transcript_mapping(db,chromosomes)
        exon_mapping = get_exon_mapping(db,chromosomes)
        exons_in_trans = get_exons_in_trans(db,chromosomes)
        trans_in_gene = get_trans_in_gene(db,chromosomes)
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
    """From a BAM file *bamfile*, returns a list containing for each exon in *exons*,
    the number of reads that mapped to it, in the same order.

    :param bamfile: name of a BAM file.
    :param labels: index references - as returned by **fetch_labels()**
    :type bamfile: string
    :type exons: list
    :rtype: list
    """
    counts = []
    try: sam = pysam.Samfile(bamfile, 'rb')
    except ValueError: sam = pysam.Samfile(bamfile,'r')

    class Counter(object):
        def __init__(self):
            self.n = 0
        def __call__(self, alignment):
            self.n += 1

    c = Counter()
    for l in labels:
        sam.fetch(l[0], 0, l[1], callback=c) #(name,0,length,Counter())
        #The callback (c.n += 1) is executed for each alignment in a region
        counts.append(c.n)
        c.n = 0
    sam.close()
    return counts

def save_results(ex, cols, conditions, header=[], feature_type='features'):
    """Save results in a tab-delimited file, one line per feature, one column per run.

    :param ex: bein's execution.
    :param cols: list of iterables, each element being a column to write in the output.
    :param conditions: list of strings corresponding to descriptions of the different samples.
    :param header: list of strings, the column headers of the output file.
    :param feature_type: (str) the kind of feature of which you measure the expression.
    """
    conditions_s = '%s, '*(len(conditions)-1)+'%s.'
    conditions = tuple(conditions)
    ncond = len(conditions)
    # Tab-delimited output with all information
    output_tab = rstring()
    writecols(output_tab,cols,header=header, sep="\t")
    description = "Expression level of "+feature_type+" in sample(s) "+conditions_s % conditions
    description = set_file_descr(feature_type.lower()+"_expression.tab", step="pileup", type="txt", comment=description)
    ex.add(output_tab, description=description)
    # SQL track output
    if feature_type=="EXONS" and 0:
        output_sql = rstring()
        def to_bedgraph(cols,ncond,i):
            lines = zip(*cols)
            print lines[0]
            yield (lines[0],lines[2*ncond+1],lines[2*ncond+2],lines[i])
        for i in range(ncond):
            with track.new(output_sql, format='bedGraph') as t:
                t.write(feature_type,to_bedgraph(cols,ncond,i),(feature_type[:-1]+"_id","start","end","score"))
            c = conditions[i]
            description = "SQL track of exons rpkm for sample %s" % c
            description = set_file_descr("exons_"+c+".sql", step="pileup", type="sql", comment=description)
            ex.add(output_sql, description=description)
    print feature_type+": Done successfully."

@timer
def genes_expression(exons_data, exon_to_gene, ncond):
    """Get gene counts from exons counts.

    Returns two dictionaries, one for counts and one for rpkm, of the form ``{gene ID: float}``.

    :param exons_data: list of lists ``[exonsID, counts, rpkm, start, end, geneID, geneName]``.
    :param exon_to_gene: dictionary ``{exon ID: gene ID}``.
    :param ncond: number of samples.
    """
    genes = list(set(exons_data[-3]))
    z = [numpy.zeros(ncond,dtype=numpy.float_) for g in genes]
    zz = [numpy.zeros(ncond,dtype=numpy.float_) for g in genes]
    gcounts = dict(zip(genes,z))
    grpk = dict(zip(genes,zz))
    round = numpy.round
    for e,c in zip(exons_data[0],zip(*exons_data[1:2*ncond+1])):
        g = exon_to_gene[e]
        gcounts[g] += round(c[:ncond],2)
        grpk[g] += round(c[ncond:],2)
    return gcounts, grpk

@timer
def transcripts_expression(exons_data, exon_lengths, transcript_mapping, trans_in_gene, exons_in_trans, ncond):
    """Get transcript rpkms from exon rpkms.

    Returns a dictionary of the form ``{transcript ID: rpkm}``.

    :param exons_data: list of lists ``[exonsID, counts, rpkm, start, end, geneID, geneName]``
    :param exon_lengths: dictionary ``{exon ID: length}``
    :param transcript_mapping: dictionary ``{transcript ID: (gene ID, start, end, length, chromosome)}``
    :param trans_in_gene: dictionary ``{gene ID: [transcript IDs it contains]}``
    :param exons_in_trans: dictionary ``{transcript ID: [exon IDs it contains]}``
    :param ncond: number of samples
    """
    genes = list(set(exons_data[-3]))
    transcripts = []
    for g in genes:
        transcripts.extend(trans_in_gene.get(g,[]))
    transcripts = list(set(transcripts))
    z = numpy.zeros((len(transcripts),ncond))
    zz = numpy.zeros((len(transcripts),ncond))
    trans_counts = dict(zip(transcripts,z))
    trans_rpk = dict(zip(transcripts,zz))
    exons_counts = dict(zip( exons_data[0], zip(*exons_data[1:ncond+1])) )
    exons_rpk = dict(zip( exons_data[0], zip(*exons_data[ncond+1:2*ncond+1])) )
    totalerror = 0; unknown = 0; alltranscount=0; allexonscount=0;
    #filE = open("../error_stats.numbers","wb")
    #filE.write("gene \t nbExons \t nbTrans \t ratioNbExonsNbTrans \t totExons \t totTrans \t ratioExonsTrans \t lsqError \n")
    for g in genes:
        if trans_in_gene.get(g): # if the gene is (still) in the Ensembl database
            # Get all transcripts in the gene
            tg = trans_in_gene[g]
            # Get all exons in the gene
            eg = set()
            for t in tg:
                if exons_in_trans.get(t):
                    eg = eg.union(set(exons_in_trans[t]))
            print eg
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
            print M
            for c in range(ncond):  # - rpkm
                Er = er[c]
                tol = 10*2.22e-16*numpy.linalg.norm(M,1)*(max(M.shape)+1)
                try: Tr, resnormr, resr = lsqnonneg(M,Er,tol=100*tol)
                except: Tr = zeros(len(tg)) # Bad idea
                tr.append(Tr)
            N = M
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
            total_trans_c = sum([sum(trans_counts[t]) for t in tg]) or 0
            total_exons_c = sum([sum(ec[c]) for c in range(ncond)]) or 0
            alltranscount += total_trans_c
            allexonscount += total_exons_c
            #try: filE.write("%s\t%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n" \
            #        % (g,len(eg),len(tg),1.*len(eg)/len(tg),total_exons_c,total_trans_c,total_exons_c/total_trans_c,resnormc))
            #except ZeroDivisionError: pass
        else:
            unknown += 1

    #filE.close()
    print "\t Unknown transcripts for %d of %d genes (%.2f %%)" \
           % (unknown, len(genes), 100*float(unknown)/float(len(genes)) )
    try: print "\t Total transcript counts: %.2f, Total exon counts: %.2f, Ratio: %.2f" \
           % (alltranscount,allexonscount,alltranscount/allexonscount)
    except ZeroDivisionError: pass
    print "\t Total error (sum of resnorms):", totalerror
    return trans_rpk, trans_counts, totalerror

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

def to_rpkm(counts, starts, ends, nreads):
    counts, starts, ends, nreads = [asarray(x,dtype=numpy.float_) for x in (counts,starts,ends,nreads)]
    return 1000*(1e6*counts.T/nreads).T/(ends-starts)

@timer
def rnaseq_workflow(ex, job, assembly, bam_files, pileup_level=["exons","genes","transcripts"], via="lsf", unmapped=False):
    """
    Main function of the workflow.

    :rtype: None

    :param ex: the bein's execution Id.
    :param job: a Job object (or a dictionary of the same form) as returned from HTSStation's frontend.
    :param assembly: the assembly Id of the species, string or int (e.g. 'hg19' or 76).
    :param bam_files: a complicated dictionary such as returned by mapseq.get_bam_wig_files.
    :param unmapped: the name or path to a fastq file containing the unmapped reads.
    :param pileup_level: a string or array of strings indicating the features you want to compare.
                         Targets can be 'genes', 'transcripts', or 'exons'.
    :param via: 'local' or 'lsf'.
    """
    group_names={}
    assembly_id = job.assembly_id
    groups = job.groups
    for gid,group in groups.iteritems():
        group_names[gid] = str(group['name']) # group_names = {gid: name}
    if isinstance(pileup_level,str): pileup_level=[pileup_level]

    #Exon labels: ('exonID|geneID|start|end|strand', length)
    exons = fetch_labels(bam_files[groups.keys()[0]][groups.values()[0]['runs'].keys()[0]]['bam'])

    """ Extract information from bam headers """
    exonsID=[]; genesID=[]; genesName=[]; starts=[]; ends=[]; exon_lengths={}; exon_to_gene={}; badexons=[]
    for e in exons:
        (exon, gene, start, end, strand) = e[0].split('|')
        start = int(start); end = int(end)
        if start!=end:
            starts.append(start)
            ends.append(end)
            exon_lengths[exon] = float(e[1])
            exonsID.append(exon)
            genesID.append(gene)
            exon_to_gene[exon] = gene
        else: badexons.append(e)
    [exons.remove(e) for e in badexons]

    """ Map remaining reads to transcriptome """
    junction_pileups={}
    if unmapped:
        print "Get unmapped reads"
        conditions=[]; unmapped_bam={}
        mdfive = get_md5(assembly_id)
        refseq_path = os.path.join("/db/genrep/nr_assemblies/cdna_bowtie/", mdfive)
        assert os.path.exists(refseq_path+".1.ebwt"), "Refseq index not found: %s" % refseq_path+".1.ebwt"
        for gid, group in job.groups.iteritems():
            k = 0
            for rid, run in group['runs'].iteritems():
                k+=1
                cond = group_names[gid]+'.'+str(k)
                conditions.append(cond)
                unmapped_fastq = bam_files[gid][rid].get('unmapped_fastq')
                unmapped_bam[cond] = mapseq.map_reads(ex, unmapped_fastq, {}, refseq_path, \
                                                      remove_pcr_duplicates=False)['bam']
        if unmapped_bam.get(conditions[0]):
            junctions = fetch_labels(unmapped_bam[conditions[0]]) #list of (transcript ID, length)
            for cond in conditions:
                junction_pileup = build_pileup(unmapped_bam[cond], junctions)
                junction_pileup = dict(zip([j[0] for j in junctions], junction_pileup)) # {transcript ID: count}
                junction_pileups[cond] = junction_pileup

    """ Build exon pileups from bam files """
    print "Build pileups"
    exon_pileups = {}; nreads = {}; conditions = []
    for gid,files in bam_files.iteritems():
        k = 0
        for rid,f in files.iteritems():
            k+=1
            cond = group_names[gid]+'.'+str(k)
            conditions.append(cond)
            exon_pileup = build_pileup(f['bam'], exons)
            exon_pileups[cond] = exon_pileup # {cond1.run1: [pileup], cond1.run2: [pileup]...}
            nreads[cond] = nreads.get(cond,0) + sum(exon_pileup) # total number of reads    :w
            print "....Pileup", cond, "done"

    print "Load mappings"
    mappings = fetch_mappings(assembly_id)
    """ [0] gene_mapping is a dict ``{gene ID: (gene name,start,end,chromosome)}``
        [1] transcript_mapping is a dictionary ``{transcript ID: (gene ID,start,end,length,chromosome)}``
        [2] exon_mapping is a dictionary ``{exon ID: ([trancript IDs],gene ID,start,end,chromosome)}``
        [3] trans_in_gene is a dict ``{gene ID: [IDs of the transcripts it contains]}``
        [4] exons_in_trans is a dict ``{transcript ID: [IDs of the exons it contains]}`` """
    (gene_mapping, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans) = mappings

    """ Treat data """
    print "Process data"
    starts = asarray(starts, dtype=numpy.float_)
    ends = asarray(ends, dtype=numpy.float_)
    nreads = asarray([nreads[cond] for cond in conditions], dtype=numpy.float_)
    counts = asarray([exon_pileups[cond] for cond in conditions], dtype=numpy.float_)
    counts, sf = estimate_size_factors(counts)
    rpkm = to_rpkm(counts, starts, ends, nreads)
    #for i in range(len(counts.ravel())):
    #    if counts.flat[i]==0: counts.flat[i] += 1.0 # if zero counts, add 1 for further comparisons

    print "Get scores"
    hconds = ["counts."+c for c in conditions] + ["rpkm."+c for c in conditions]
    genesName, echr = zip(*[(x[0],x[-1]) for x in [gene_mapping.get(g,"NA") for g in genesID]])
    exons_data = [exonsID]+list(counts)+list(rpkm)+[starts,ends,genesID,genesName,echr]

    """ Print counts for exons """
    if "exons" in pileup_level:
        header = ["ExonID"] + hconds + ["Start","End","GeneID","GeneName","Chromosome"]
        save_results(ex, exons_data, conditions, header=header, feature_type="EXONS")

    """ Get scores of genes from exons """
    if "genes" in pileup_level:
        header = ["GeneID"] + hconds + ["Start","End","GeneName","Chromosome"]
        (gcounts, grpkm) = genes_expression(exons_data, exon_to_gene, len(conditions))
        genesID = gcounts.keys()
        genes_data = [[g,gcounts[g],grpkm[g]]+list(gene_mapping.get(g,("NA",)*4,)) for g in genesID]
        (genesID,gcounts,grpkm,gname,gstart,gend,gchr) = zip(*genes_data)
        genes_data = [genesID]+list(zip(*gcounts))+list(zip(*grpkm))+[gstart,gend,gname,gchr]
        save_results(ex, genes_data, conditions, header=header, feature_type="GENES")

    """ Get scores of transcripts from exons, using non-negative least-squares """
    if "transcripts" in pileup_level:
        header = ["TranscriptID"] + hconds + ["Start","End","GeneID","GeneName","Chromosome"]
        (tcounts, trpkm, error) = transcripts_expression(exons_data, exon_lengths,
                   transcript_mapping, trans_in_gene, exons_in_trans,len(conditions))
        transID = tcounts.keys()
        trans_data = [[t,tcounts[t],trpkm[t]]+list(transcript_mapping.get(t,("NA",)*5,)) for t in transID]
        trans_data = sorted(trans_data, key=lambda x: x[-5]) # sort w.r.t. gene IDs
        (transID,tcounts,trpkm,genesID,tstart,tend,tlen,tchr) = zip(*trans_data)
        genesName = [gene_mapping.get(g,("NA",)*4)[0] for g in genesID]
        trans_data = [transID]+list(zip(*tcounts))+list(zip(*trpkm))+[tstart,tend,genesID,genesName,tchr]

        a = trans_data
        print a[1]
        print to_rpkm(a[1],a[3],a[4],nreads[0])
        print a[2]
        # results should be the same for ENSMUST00000073605, ENSMUST00000117757, ENSMUST00000118875,
        # which contain the same exons.

        save_results(ex, trans_data, conditions, header=header, feature_type="TRANSCRIPTS")

    # TEST
    # "ENSMUSG00000057666" "ENSG00000111640": TEST Gapdh
    # for g in list(set(genesID)):
    #     print sum(gcounts[g]), sum([sum(tcounts[t]) for t in trans_in_gene[g]])

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#



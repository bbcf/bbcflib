"""
======================
Module: bbcflib.rnaseq
======================

Methods of the bbcflib's RNA-seq worflow. The main function is ``rnaseq_workflow()``.

From a BAM file produced by an alignement on the genome or the exonome, gets counts of reads
on the exons, add them to get counts on genes, and uses least-squares to infer
counts on transcripts, avoiding to map on the transcriptome.
Note that the resulting counts on transcripts are approximate.
"""

# Built-in modules #
import os, sys, pysam, math, itertools
from operator import itemgetter

# Internal modules #
from bbcflib.common import set_file_descr, unique_filename_in, cat, timer, program_exists
from bbcflib.genrep import Assembly, Transcript
from bbcflib.mapseq import add_and_index_bam, sam_to_bam, map_reads, plot_stats
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


#--------------------------------- MATHS ---------------------------------#


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

    Reference: C.L. Lawson and R.J. Hanson, Solving Least-Squares Problems, Prentice-Hall, Chapter 23, p. 161, 1974.
    `<http://diffusion-mri.googlecode.com/svn/trunk/Python/lsqnonneg.py>`_
    """
    def norm1(x):
        return abs(x).sum().max()

    def msize(x, dim):
        s = x.shape
        if dim >= len(s): return 1
        else: return s[dim]

    eps = sys.float_info.epsilon
    if tol is None: tol = 10*eps*norm1(C)*(max(C.shape)+1)
    C = asarray(C)
    (m,n) = C.shape
    P = numpy.zeros(n)                       # set P of indices where ultimately x_j > 0 & w_j = 0
    Z = ZZ = numpy.arange(1, n+1)            # set Z of indices where ultimately x_j = 0 & w_j <= 0
    if x0 is None or any(x0 < 0): x = P
    else: x = x0
    resid = d - numpy.dot(C, x)
    w = numpy.dot(C.T, resid)                # n-vector C'(d-Cx), "dual" of x, gradient of (1/2)*||d-Cx||^2
    outeriter = it = 0
    itmax = itmax_factor*n
    # Outer loop "A" to hold positive coefficients
    while numpy.any(Z) and numpy.any(w[ZZ-1] > tol): # if Z is empty or w_j<0 for all j, terminate.
        outeriter += 1
        t = w[ZZ-1].argmax()                 # find index t s.t. w_t = max(w), w_t in Z. So w_t > 0.
        t = ZZ[t]
        P[t-1]=t                             # move the index t from set Z to set P
        Z[t-1]=0                             # Z becomes [0] if n=1
        PP = numpy.where(P != 0)[0]+1        # non-zero elements of P for indexing, +1 (-1 later)
        ZZ = numpy.where(Z != 0)[0]+1        # non-zero elements of Z for indexing, +1 (-1 later)
        CP = numpy.zeros(C.shape)
        CP[:, PP-1] = C[:, PP-1]             # CP[:,j] is C[:,j] if j in P, or 0 if j in Z
        CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))
        z=numpy.dot(numpy.linalg.pinv(CP), d)                 # n-vector solution of least-squares min||d-CPx||
        if isinstance(ZZ,numpy.ndarray) and len(ZZ) == 0:     # if Z = [0], ZZ = [] and makes it fail
            return (positive(z), sum(resid*resid), resid)
        else:
            z[ZZ-1] = numpy.zeros((msize(ZZ,1), msize(ZZ,0))) # define z_j := 0 for j in Z
        # Inner loop "B" to remove negative elements from z if necessary
        while numpy.any(z[PP-1] <= tol): # if z_j>0 for all j, set x=z and return to outer loop
            it += 1
            if it > itmax:
                max_error = z[PP-1].max()
                raise Exception('Exiting: Iteration count (=%d) exceeded\n \
                      Try raising the tolerance tol. (max_error=%d)' % (it, max_error))
            QQ = numpy.where((z <= tol) & (P != 0))[0]        # indices j in P s.t. z_j < 0
            alpha = min(x[QQ]/(x[QQ] - z[QQ]))                # step chosen as large as possible s.t. x remains >= 0
            x = x + alpha*(z-x)                               # move x by this step
            ij = numpy.where((abs(x) < tol) & (P <> 0))[0]+1  # indices j in P for which x_j = 0
            Z[ij-1] = ij                                      # Add to Z, remove from P
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


#----------------------------- GLOBAL CLASSES -----------------------------#


class RNAseq(object):
    """Abstract, inherited by different parts of the workflow."""
    def __init__(self,ex,via,job,assembly,conditions,debugfile,logfile,
                 pileup_level,rpath,juliapath,junctions,unmapped,mappings,stranded):
        self.ex = ex                     # bein.Execution object
        self.job = job                   # frontend.Job object
        self.assembly = assembly         # genrep.Assembly object
        self.conditions = conditions     # list of <group_name>.<run_id>
        self.pileup_level = pileup_level # ["genes","exons","transcripts"], or another combination
        self.via = via                   # "local"/"lsf"
        self.rpath = rpath               # Path to R scripts
        self.juliapath = juliapath       # Path to Julia scripts
        self.junctions = junctions       # True/False, want to find or not
        self.unmapped = unmapped         # True/False, want to remap or not
        self.debugfile = debugfile       # debug file name
        self.logfile = logfile           # log file name
        self.M = mappings                # Mappings object
        self.stranded = stranded         # True/False, strand-specific or not
        self.ncond = 2*len(self.conditions) if self.stranded else len(self.conditions)
            # number of columns in the output, = number of runs if not stranded, twice otherwise

    def write_log(self,s):
        self.logfile.write(s+'\n'); self.logfile.flush()

    def write_debug(self,s):
        self.debugfile.write(s+'\n'); self.debugfile.flush()


class Mappings():
    def __init__(self, assembly):
        self.assembly = assembly

    @timer
    def fetch_mappings(self):
        """
        * gene_mapping is a dict ``{gene_id: genrep.Gene}``
        * transcript_mapping is a dictionary ``{transcript_id: genrep.Transcript}``
        * exon_mapping is a dictionary ``{exon_id: genrep.Exon}``
        """
        map_path = '/archive/epfl/bbcf/jdelafon/mappings/test_%s/' % self.assembly.name
        if test and os.path.exists(map_path):
            import pickle
            self.gene_mapping = pickle.load(open(os.path.join(map_path,'gene_mapping.pickle')))
            self.exon_mapping = pickle.load(open(os.path.join(map_path,'exon_mapping.pickle')))
            self.transcript_mapping = pickle.load(open(os.path.join(map_path,'transcript_mapping.pickle')))
        else:
            self.gene_mapping = self.assembly.get_gene_mapping()
            self.transcript_mapping = self.assembly.get_transcript_mapping()
            self.exon_mapping = self.assembly.get_exon_mapping()


class Counter(object):
    def __init__(self, stranded=False):
        self.n = 0 # read count
        self.start = 0 # exon start
        self.counts = [] # vector of counts per non-zero position
        self.strand = 0 # exon strand
        if stranded: self.count_fct = self.count_stranded
        else: self.count_fct = self.count

    def __call__(self, alignment):
        self.count_fct(alignment)

    def count(self, alignment):
        NH = [1.0/t[1] for t in alignment.tags if t[0]=='NH']+[1]
        self.n += NH[0]
        try:
            self.counts[alignment.pos-self.start] += NH[0]
        except IndexError:
            pass # read overflows

    def count_stranded(self, alignment):
        if self.strand == 1 and alignment.is_reverse == False \
        or self.strand == -1 and alignment.is_reverse == True:
            self.count(alignment)

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


#----------------------------- MAIN CALL -----------------------------#


@timer
def rnaseq_workflow(ex, job, assembly=None, pileup_level=["exons","genes","transcripts"], via="lsf",
                    rpath=None, junctions=False, unmapped=False, stranded=False,
                    logfile=sys.stdout, debugfile=sys.stderr):
    """Main function of the workflow.

    :rtype: None
    :param ex: the bein's execution Id.
    :param job: a Frontend.Job object (or a dictionary of the same form).
    :param assembly: a genrep.Assembly object
    :param rpath: (str) path to the R executable.
    :param junctions: (bool) whether to search for splice junctions using SOAPsplice. [False]
    :param unmapped: (bool) whether to remap to the transcriptome reads that did not map the genome. [False]
    :param via: (str) send job via 'local' or 'lsf'. ["lsf"]
    """
    juliapath='/home/jdelafon/repos/bbcfutils/Julia/'

    if test:
        via = 'local'
        logfile = sys.stdout
        debugfile = sys.stdout
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

    # Define conditions as 'group_name.run_id' and store bamfiles
    bamfiles = {}
    for gid,files in job.files.iteritems():
        k = 0
        for rid,f in files.iteritems():
            k+=1
            cond = group_names[gid]+'.'+str(k)
            conditions.append(cond)
            bamfiles[cond] = f['bam']
    ncond = len(conditions)

    M = Mappings(assembly)
    rnaseq_args = (ex,via,job,assembly,conditions,debugfile,logfile,
                   pileup_level,rpath,juliapath,junctions,unmapped,M,stranded)
    PU = Pileups(*rnaseq_args)
    DE = DE_Analysis(*rnaseq_args)
    UN = Unmapped(*rnaseq_args)
    PCA = Pca(*rnaseq_args)
    JN = Junctions(*rnaseq_args)

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
            firstbam = bamfiles.items()[0][1]
            firstbamtrack = track(firstbam,format='bam')
            for c,meta in firstbamtrack.chrmeta.iteritems():
                tmap[c] = Transcript(end=meta['length'], length=meta['length'])
            ftype = "Custom"
            header = ["CustomID"]
        logfile.write("* Build pileups\n"); logfile.flush()
        pileups = {}
        for cond in conditions:
            pileup = PU.build_custom_pileup(bamfiles[cond],tmap)
            pileups[cond] = pileup.values()
            logfile.write("  ....Pileup %s done\n" % cond); logfile.flush()
        ids = pileup.keys()
        counts = asarray([pileups[cond] for cond in conditions], dtype=numpy.float_).T
        del pileup, pileups
        tcounts = {}
        for k,t in enumerate(ids):
            c = counts[k]
            if sum(c) != 0 and t in tmap:
                tcounts[t] = c
        lengths = asarray([tmap[t].end-tmap[t].start for t in tcounts.iterkeys()])
        trans_data = PU.norm_and_format(tcounts,lengths,tmap)
        header += ["counts."+c for c in conditions] + ["norm."+c for c in conditions] + ["rpkm."+c for c in conditions]
        header += ["Start","End","GeneID","GeneName","Strand","Chromosome"]
        PU.save_results(trans_data,group_ids,header=header,feature_type=ftype)
        DE.differential_analysis(trans_data, header, feature_type=ftype.lower())
        return 0

    logfile.write("* Load mappings\n"); logfile.flush()
    M.fetch_mappings()
    if len(M.exon_mapping) == 0 or len(M.gene_mapping) == 0:
        raise ValueError("No genes found for this genome. Abort.")


    # Build exon pileups from bam files
    logfile.write("* Build pileups\n"); logfile.flush()
    exon_pileups = {}
    for i,cond in enumerate(conditions):
        exon_pileup, exon_pileup_rev = PU.build_pileup(bamfiles[cond])
        if stranded:
            for e,count in exon_pileup_rev.iteritems():
                exon_pileups.setdefault(e,zeros(2*ncond))[i] = count
                exon_pileups[e][ncond+i] = exon_pileup_rev[e]
        else:
            for e,count in exon_pileup.iteritems():
                exon_pileups.setdefault(e,zeros(ncond))[i] = count
    del exon_pileup

    if stranded:
        hconds = ["counts."+c for c in conditions] + ["counts_rev."+c for c in conditions] \
               + ["norm."+c for c in conditions] + ["norm_rev."+c for c in conditions] \
               + ["rpkm."+c for c in conditions] + ["rpkm_rev."+c for c in conditions]
    else:
        hconds = ["counts."+c for c in conditions] + ["norm."+c for c in conditions] + ["rpkm."+c for c in conditions]

    # Print counts for exons
    if "exons" in pileup_level:
        logfile.write("* Get scores of exons\n"); logfile.flush()
        header = ["ExonID"] + hconds + ["Start","End","GeneID","GeneName","Strand","Chromosome"]
        lengths = asarray([M.exon_mapping[e].end-M.exon_mapping[e].start for e in exon_pileups])
        exons_data = PU.norm_and_format(exon_pileups,lengths,M.exon_mapping)
        PU.save_results(exons_data, group_ids, header=header, feature_type="EXONS")
        DE.differential_analysis(exons_data, header, feature_type='exons')
        del exons_data

    # Map remaining reads to transcriptome
    unmapped_fastq = {}
    if unmapped:
        logfile.write("* Align unmapped reads on transcriptome\n"); logfile.flush()
        try:
            unmapped_fastq,additionals = UN.align_unmapped(group_names)
        except Exception, error:
            debugfile.write("Remapping failed: %s\n"%str(error)); debugfile.flush()

    # Add junction reads to exon pileups
    for i,cond in enumerate(conditions):
        if unmapped and (cond in unmapped_fastq) and (cond in additionals):
            for e,x in additionals[cond].iteritems():
                if e in exon_pileups:
                    exon_pileups[e][i] += x
                else:
                    exon_pileups.setdefault(e,zeros(ncond))[i] = x
    if unmapped: del additionals

    # Get scores of genes from exons
    if "genes" in pileup_level:
        logfile.write("* Get scores of genes\n"); logfile.flush()
        header = ["GeneID"] + hconds + ["Start","End","GeneName","Strand","Chromosome"]
        gcounts = PU.genes_expression(exon_pileups)
        lengths = asarray([M.gene_mapping[g].length for g in gcounts.iterkeys()])
        genes_data = PU.norm_and_format(gcounts,lengths,M.gene_mapping)
        genes_file = PU.save_results(genes_data, group_ids, header=header, feature_type="GENES")
        DE.differential_analysis(genes_data, header, feature_type='genes')
        del genes_data

        # PCA of groups ~ gene expression
        if ncond >= 2:
            try:
                PCA.pca_rnaseq(genes_file)
            except Exception, error:
                debugfile.write("PCA failed: %s\n"%str(error)); debugfile.flush()

    # Get scores of transcripts from exons, using non-negative least-squares
    if "transcripts" in pileup_level:
        logfile.write("* Get scores of transcripts\n"); logfile.flush()
        header = ["TranscriptID"] + hconds + ["Start","End","GeneID","GeneName","Strand","Chromosome"]
        tcounts = PU.transcripts_expression(exon_pileups)
        lengths = asarray([M.transcript_mapping[t].length for t in tcounts.iterkeys()])
        trans_data = PU.norm_and_format(tcounts,lengths,M.transcript_mapping)
        PU.save_results(trans_data, group_ids, header=header, feature_type="TRANSCRIPTS")
        DE.differential_analysis(trans_data, header, feature_type='transcripts')
        del trans_data

    # Find splice junctions
    if junctions:
        logfile.write("* Search for splice junctions\n"); logfile.flush()
        try:
            JN.find_junctions()
        except Exception, error:
            debugfile.write("Junctions search failed: %s\n"%str(error)); debugfile.flush()

    return 0


#----------------------------- READ COUNTING -----------------------------#


class Pileups(RNAseq):
    def __init__(self,*args):
        RNAseq.__init__(self,*args)

    @timer
    def build_custom_pileup(self, bamfile, transcript_mapping):
        counts = {}
        try: sam = pysam.Samfile(bamfile, 'rb')
        except ValueError: sam = pysam.Samfile(bamfile,'r')
        c = Counter()
        for n,ref in enumerate(sam.references):
            start = 0
            end = transcript_mapping.get(ref.split('|')[0], Transcript(end=0) ).end
            sam.fetch(ref, start, end, callback=c)
            counts[ref] = int(.5+c.n)
            c.n = 0
        sam.close()
        return counts

    @timer
    def build_pileup(self, bamfile):
        """From a BAM file, returns a dictionary of the form {feature_id: number of reads that mapped to it}."""
        counts = {}
        counts_rev = {}
        try: sam = pysam.Samfile(bamfile, 'rb')
        except ValueError: sam = pysam.Samfile(bamfile,'r')
        chromosomes = self.assembly.chrmeta.keys()
        ref_index = {}
        for ref in sam.references:
            ref_index[ref.split('|')[0]] = ref
        mapped_on = 'genome' if all([ref in chromosomes for ref in sam.references[:100]]) else 'exons'
        c = Counter(stranded=self.stranded)

        def _count(etrack):
            for (chrom,start,end,names,strand,phase) in etrack:
                if strand == 0: continue  # skip overlapping exons of different strands to avoid typeI errors
                exon_ids = [e.split('|')[0] for e in names.split('$')]  # all spanning this interval - maybe more
                if mapped_on == 'genome':
                    ref = chrom
                else: # mapped on 'exons': need to subtract exon start coordinate.
                    # All have the same sequence in this interval. Take any of them,
                    # just choose the shift from its start accordingly.
                    oneofthem = exon_ids[0]
                    ostart = self.M.exon_mapping[oneofthem].start
                    shift = start-ostart  # dist from start of this exon to our region start
                    ref = ref_index.get(oneofthem)
                    start = shift; end = shift+(end-start)
                if not ref: continue  # exon in bam header but not in exon_track
                try:
                    c.n = 0
                    c.counts = [0]*(end-start+1)  # do not allow read overflow
                    c.start = start
                    c.strand = strand
                    sam.fetch(ref,start,end, callback=c)  # calls Counter.__call__ for each alignment: updates c.n
                    #c.remove_duplicates()
                except ValueError,ve:  # unknown reference
                    self.write_debug("Unknown reference: %s",str(ve))
                nexons = float(len(exon_ids))
                for exon in exon_ids:
                    if c.n > 0:  # remove this if want to keep all references in output
                        counts[exon] = counts.get(exon,0) + c.n/nexons
            return counts

        for chrom in chromosomes:
            etrack = self.assembly.exon_track(chromlist=[chrom], biotype=None)
            etrack = cobble(etrack, aggregate={'name':lambda n:'$'.join(n)})
            counts = _count(etrack)
            if self.stranded:
                etrack_rev = self.assembly.exon_track(chromlist=[chrom], biotype=None)
                etrack_rev = apply(etrack_rev, 'strand', lambda x:-x)  # reverse strand info
                etrack_rev = cobble(etrack_rev, aggregate={'name':lambda n:'$'.join(n)})
                counts_rev = _count(etrack_rev)
        sam.close()
        return counts, counts_rev

    @timer
    def save_results(self, lines, group_ids, header, feature_type='features'):
        """Save results in a tab-delimited file, one line per feature, one column per run.

        :param lines: list of iterables, each element being a line to write in the output.
        :param group_ids: dictionary ``{group name: group_id}``.
        :param header: list of strings, the column headers of the output file.
        :param feature_type: (str) the kind of feature of which you measure the expression.
        """
        # Tab-delimited output with all information
        output_tab = unique_filename_in()
        n = self.ncond
        with open(output_tab,'wb') as f:
            f.write('\t'.join(header)+'\n')
            for l in lines:
                tid = str(l[0])
                counts = ["%d"%x for x in l[1:n+1]]
                norm = ["%.2f"%x for x in l[n+1:2*n+1]]
                rpkm = ["%.2f"%x for x in  l[2*n+1:3*n+1]]
                rest = [str(x) for x in l[3*n+1:]]
                f.write('\t'.join([tid]+counts+norm+rpkm+rest)+'\n')
        description = set_file_descr(feature_type.lower()+"_expression.tab", step="pileup", type="txt")
        self.ex.add(output_tab, description=description)
        # Create one track for each group
        if feature_type in ['GENES','EXONS']:
            cols = zip(*lines)
            groups = [c.split('.')[0] for c in self.conditions]
            start = cols[3*n+1]
            end = cols[3*n+2]
            chromosomes = cols[-1]
            rpkm = {}; output_sql = {}
            for i,cond in enumerate(self.conditions):
                group = cond.split('.')[0]
                nruns = groups.count(group)
                # Average all replicates in the group
                rpkm[group] = asarray(rpkm.get(group,zeros(len(start)))) + asarray(cols[i+2*n+1]) / nruns
                output_sql[group] = output_sql.get(group,unique_filename_in())
            for group,filename in output_sql.iteritems():
                # SQL track - GDV
                tr = track(filename+'.sql', fields=['chr','start','end','score'], chrmeta=self.assembly.chrmeta)
                towrite = {}
                for n,c in enumerate(chromosomes):
                    if c not in tr.chrmeta: continue
                    towrite.setdefault(c,[]).append((int(start[n]),int(end[n]),rpkm[group][n]))
                for chrom, feats in towrite.iteritems():
                    tr.write(cobble(sorted_stream(FeatureStream(feats, fields=['start','end','score']))),
                             chrom=chrom,clip=True)
                description = set_file_descr(feature_type.lower()+"_"+group+".sql", step="pileup", type="sql",
                                             groupId=group_ids[group])
                self.ex.add(filename+'.sql', description=description)
                # bigWig track - UCSC
                try: # if the bigWig conversion program fails, the file is not created
                    convert(filename+'.sql',filename+'.bw')
                    description = set_file_descr(feature_type.lower()+"_"+group+".bw",
                                                 step="pileup", type="bw", groupId=group_ids[group], ucsc='1')
                    self.ex.add(filename+'.bw', description=description)
                except IOError: pass
        self.write_log("  %s: Done successfully." % feature_type)
        return os.path.abspath(output_tab)

    @timer
    def genes_expression(self, exon_pileups):
        """Get gene counts/rpkm from exon counts (sum).
        Returns a dictionary of the form ``{gene_id: scores}``.
        """
        gene_counts = {}
        for e,counts in exon_pileups.iteritems():
            g = self.M.exon_mapping[e].gene_id
            gene_counts[g] = gene_counts.get(g,zeros(self.ncond)) + counts
        return gene_counts

    @timer
    def transcripts_expression(self, exon_pileups):
        """Get transcript rpkms from exon rpkms using NNLS.
        Returns a dictionary of the form ``{transcript_id: score}``.
        """
        emap = self.M.exon_mapping
        tmap = self.M.transcript_mapping
        gmap = self.M.gene_mapping
        trans_counts={}
        exons_counts={}; genes=[];
        for e,counts in exon_pileups.iteritems():
            exons_counts[e] = counts
            genes.append(emap[e].gene_id)
        genes = set(genes)
        unknown = 0
        pinv = numpy.linalg.pinv
        for g in genes:
            if g in gmap: # if the gene is (still) in the Ensembl database
                # Get all transcripts in the gene
                tg = gmap[g].transcripts
                # Get all exons in the gene
                eg = set()
                for t in tg:
                    if t in tmap:
                        eg = eg.union(set(tmap[t].exons))
                        trans_counts[t] = zeros(self.ncond)
                # Create the correspondance matrix
                M = zeros((len(eg),len(tg)))
                L = zeros((len(eg),len(tg)))
                ec = zeros((self.ncond,len(eg)))
                for i,e in enumerate(eg):
                    for j,t in enumerate(tg):
                        if t in tmap and e in tmap[t].exons:
                            M[i,j] = 1.
                            if e in emap and t in tmap:
                                L[i,j] = float((emap[e].end-emap[e].start)) / tmap[t].length
                            else:
                                L[i,j] = 1./len(eg)
                    # Retrieve exon scores
                    if e in exons_counts:
                        for c in range(self.ncond):
                            ec[c][i] += exons_counts[e][c]
                # If only one transcript, special case for more precision
                if len(tg) == 1:
                    for c in range(self.ncond):
                        trans_counts[t][c] = sum(ec[c])
                    continue
                # Compute transcript scores
                tc = []
                M = numpy.vstack((M,numpy.ones(M.shape[1]))) # add constraint |E| = |T|
                L = numpy.vstack((L,numpy.ones(M.shape[1])))
                N = M*L
                for c in range(self.ncond):
                    Ec = ec[c]
                    Ec = numpy.hstack((Ec,asarray(sum(Ec))))
                    try: Tc, resnormc, resc = lsqnonneg(N,Ec,itmax_factor=100)
                    except: Tc = positive(numpy.dot(pinv(N),Ec))
                    tc.append(Tc)
                # Store results in a dict *trans_counts*
                for k,t in enumerate(tg):
                    for c in range(self.ncond):
                        trans_counts[t][c] = tc[c][k]
            else:
                unknown += 1
        if unknown != 0:
            self.write_debug("\tUnknown transcripts for %d of %d genes (%.2f %%)" \
                  % (unknown, len(genes), 100*float(unknown)/float(len(genes))) )
        return trans_counts

    def norm_and_format(self,counts,lengths,tmap,ids=None):
        """Normalize, compute RPKM values and format array lines as they will be printed."""

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
            cnts = counts[nonzero(numpy.prod(counts,1))]  # none of the counts is zero
            logcnt = numpy.log(cnts)
            loggeomeans = numpy.mean(logcnt, 1)
            size_factors = numpy.exp(numpy.median(logcnt.T - loggeomeans, 1))
            res = counts / size_factors
            return res, size_factors

        if isinstance(counts,dict):
            counts_matrix = asarray(counts.values())
            ids = counts.iterkeys()
        else:
            counts_matrix = counts

        def feature_info(x):
            info = tmap.get(x, Transcript())
            return (info.start,info.end,info.gene_id,info.gene_name,info.strand,info.chrom)

        norm_matrix, sf = estimate_size_factors(counts_matrix)
        rpkm_matrix = 1000.*norm_matrix/(lengths[:,numpy.newaxis])
        data = [(t,) + tuple(counts_matrix[k]) + tuple(norm_matrix[k]) + tuple(rpkm_matrix[k])
                     + feature_info(t)  for k,t in enumerate(ids)]
        data = sorted(data, key=itemgetter(-1,0,1,-2)) # sort wrt. chr,start,end,strand
        return data


#-------------------------- DIFFERENTIAL ANALYSIS ----------------------------#


class DE_Analysis(RNAseq):
    def __init__(self,*args):
        RNAseq.__init__(self,*args)

    def clean_before_deseq(self, data, header, keep=0.6):
        """Delete all lines of *filename* where counts are 0 in every run.

        :param keep: fraction of highest counts to keep. [0.6]
        """
        norm = 'counts'  # the cols we want to send to DESeq
        if norm == 'counts': w = 0
        elif norm == 'norm': w = 1
        elif norm == 'rpkm': w = 2
        filename_clean = unique_filename_in()
        ncond = sum([h.split('.').count("counts") for h in header]) \
              + sum([h.split('.').count("counts_rev") for h in header])
        if ncond >1:
            rownames = asarray(['%s|%s|%s|%s' % (x[0],x[-3],x[-2],x[-1]) for x in data])
            M = asarray([x[1+w*ncond:1+(w+1)*ncond] for x in data])
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

    def clean_deseq_output(self, filename):
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
    def differential_analysis(self, data, header, feature_type):
        """For each file in *data*, launch an analysis of differential expression on the count
        values, and saves the output in the MiniLIMS."""

        @program
        def run_glm(data_file, script_path, options=[]):
            """Run *rpath*/negbin.test.R on *data_file*."""
            output_file = unique_filename_in()
            opts = ["-o",output_file]+options
            return {'arguments': ["R","--slave","-f",script_path,"--args",data_file]+opts,
                    'return_value': output_file}

        res_file, ncond = self.clean_before_deseq(data, header)
        if ncond < 2:
            self.write_log("  Skipped differential analysis: less than two groups.")
        else:
            self.write_log("  Differential analysis")
            options = ['-s','tab']
            script_path = os.path.join(self.rpath,'negbin.test.R')
            if not os.path.exists(script_path):
                self.write_debug("DE: R path not found: '%s'" % self.rpath)
            try:
                glmfile = run_glm.nonblocking(self.ex, res_file, script_path, options, via=self.via).wait()
            except Exception as exc:
                self.write_log("  Skipped differential analysis: %s" % exc)
                return
            if not glmfile:
                self.write_log("  ....Empty file.\n")
                return
            output_files = [f for f in os.listdir(self.ex.working_directory) if glmfile in f]
            for o in output_files:
                desc = set_file_descr(feature_type+"_differential"+o.split(glmfile)[1]+".txt", step='stats', type='txt')
                o = self.clean_deseq_output(o)
                self.ex.add(o, description=desc)
            self.write_log("  ....done.")


#-------------------------- SPLICE JUNCTIONS SEARCH ----------------------------#


class Junctions(RNAseq):
    def __init__(self,*args):
        RNAseq.__init__(self,*args)

    @timer
    def find_junctions(self, soapsplice_index=None, path_to_soapsplice=None, soapsplice_options={}):
        """
        Retrieve unmapped reads from a precedent mapping and runs SOAPsplice on them.
        Return the names of a .bed track indicating the junctions positions, as well as
        of a bam file of the alignments attesting the junctions.

        :param soapsplice_index: (str) path to the SOAPsplice index.
        :param path_to_soapsplice: (str) specify the path to the program if it is not in your $PATH.
        :param soapsplice_options: (dict) SOAPsplice options, e.g. {'-m':2}.
        :rtype: str, str
        """

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

        assert program_exists('soapsplice')
        self.assembly.set_index_path(intype=3)
        soapsplice_index = soapsplice_index or self.assembly.index_path
        soapsplice_options.update(self.job.options.get('soapsplice_options',{}))
        soapsplice_options.setdefault('-p',16) # number of threads
        soapsplice_options.setdefault('-q',1)  # Sanger format
        unmapped_fastq = {}
        for gid, group in self.job.groups.iteritems():
            unmapped_fastq[gid] = []
            for rid, run in group['runs'].iteritems():
                unmapped = self.job.files[gid][rid].get('unmapped_fastq')
                if not unmapped:
                    self.write_log("No unmapped reads found for group %s, run %d. Skip." % (gid,rid))
                    continue
                elif not isinstance(unmapped,tuple):
                    self.write_log("Pair-end reads required. Skip.")
                    continue
                unmapped_fastq[gid].append(unmapped)
            if len(unmapped_fastq[gid]) == 0:
                continue
            R1 = cat(zip(*unmapped_fastq[gid])[0])
            R2 = cat(zip(*unmapped_fastq[gid])[1])
            future = soapsplice.nonblocking(self.ex,R1,R2,soapsplice_index,
                                            path_to_soapsplice=path_to_soapsplice,
                                            options=soapsplice_options,
                                            via=self.via, memory=8, threads=soapsplice_options['-p'])
            template = future.wait()
            if not template: return
            junc_file = template+'.junc'
            bed = self.convert_junc_file(junc_file,self.assembly)
            bed_descr = set_file_descr('junctions_%s.bed' % group['name'],
                                       groupId=gid,type='bed',step='junctions',ucsc=1)
            bam_descr = set_file_descr('junctions_%s.bam' % group['name'],
                                       groupId=gid,type='bam',step='junctions')
            sam = template+'.sam'
            try:
                bam = sam_to_bam(self.ex,sam,reheader=self.assembly.name)
                add_and_index_bam(self.ex, bam, description=bam_descr)
                self.ex.add(bam, description=bam_descr)
            except Exception, e:
                self.write_debug("%s\n(Qualities may be in the wrong format, try with '-q 0'.)" %str(e))
            self.ex.add(bed, description=bed_descr)
        return bed, bam

    def convert_junc_file(self, filename):
        """Convert a .junc SOAPsplice output file to bed format. Return the file name.

        :param filename: (str) name of the .junc file to convert.
        """
        t = track(filename, format='txt', fields=['chr','start','end','strand','score'],
                  chrmeta=self.assembly.chrmeta)
        stream = t.read()
        # Translate chromosome names
        s1 = map_chromosomes(stream, self.assembly.chromosomes)
        # Add junction IDs
        s2 = duplicate(s1,'strand','name')
        C = itertools.count()
        s3 = apply(s2,'name', lambda x: 'junction'+str(C.next()))
        # Convert to bed format
        outfile = unique_filename_in()
        bed = outfile + '.bed'
        out = track(bed, fields=s3.fields, chrmeta=self.assembly.chrmeta)
        out.write(s3)
        return bed


#-------------------------- RE-MAPPING ON TRANSCRIPTOME ----------------------------#


class Unmapped(RNAseq):
    def __init__(self,*args):
        RNAseq.__init__(self,*args)

    @timer
    def align_unmapped(self, group_names):
        """
        Map reads that did not map to the exons to a collection of annotated transcripts,
        in order to add counts to pairs of exons involved in splicing junctions.

        Return a dictionary ``unmapped_fastq`` of the form ``{sample_name:bam_file}``,
        and a second ``additionals`` of the form ``{sample_name:{exon:additional_counts}}``.

        :param group_names: dict of the form ``{group_id: group_name}``.
        """
        self.assembly.set_index_path(intype=2)
        additionals = {}
        unmapped_fastq = {}
        refseq_path = self.assembly.index_path
        bwt2 = self.job.options.get("bowtie2",True)
        for gid, group in self.job.groups.iteritems():
            k = 0
            for rid, run in group['runs'].iteritems():
                k += 1
                cond = group_names[gid]+'.'+str(k)
                _fastq = self.job.files[gid][rid].get('unmapped_fastq')
                if _fastq and os.path.exists(refseq_path+".1.ebwt"):
                    try:
                        _bam = map_reads( self.ex, _fastq, {}, refseq_path, bowtie_2=bwt2,
                                          remove_pcr_duplicates=False, via=self.via )
                    except Exception, error:
                        self.write_debug("Sample %s: map_reads (on transcriptome) failed: \n%s." % (cond,str(error)))
                        continue
                    if not _bam: continue
                    bamfile = _bam['bam']
                    stats = _bam['stats']  # result of mapseq.bamstats
                    #description = set_file_descr(cond+"_transcriptome_mapping.bam",step="pileup",type="bam")
                    #ex.add(bamfile, description=description)
                    try:
                        pdfstats = plot_stats(self.ex, stats, script_path=self.rpath)
                        description = set_file_descr(cond+"_transcriptome_mapping_report.pdf",step="pileup",type="pdf")
                        self.ex.add(pdfstats, description=description)
                    except Exception, error:
                        self.write_debug("Sample %s: bamstats failed: \n%s" % (cond,str(error)))
                    sam = pysam.Samfile(bamfile)
                    additional = {}
                    for read in sam:
                        t_id = sam.getrname(read.tid).split('|')[0]
                        if t_id in self.M.transcript_mapping:
                            E = self.M.trasncript_mapping[t_id].exons
                            if len(E) < 2: continue  # need 2 exons to make a junction
                            lag = 0
                            r_start = read.pos
                            r_end = r_start + read.rlen
                            NH = [1.0/t[1] for t in read.tags if t[0]=='NH']+[1]
                            for e in E:
                                e_start = self.M.exon_mapping[e].start
                                e_end = self.M.exon_mapping[e].end
                                e_len = e_end-e_start
                                if r_start < lag <= r_end <= lag+e_len \
                                or lag <= r_start <= lag+e_len < r_end:
                                    additional[e] = additional.get(e,0) + 0.5*NH[0]
                                lag += e_len
                    additionals[cond] = additional
                    unmapped_fastq[cond] = _fastq
                    sam.close()
        return unmapped_fastq,additionals


#---------------------------------- PCA ----------------------------------#


class Pca(RNAseq):
    def __init__(self,*args):
        RNAseq.__init__(self,*args)

    def pca_rnaseq(self,counts_table_file):
        @program
        def pcajl(counts_table_file):
            assert program_exists('julia')
            outprefix = unique_filename_in()
            args = ['julia', os.path.join(self.juliapath,'pca_rnaseq.jl'), counts_table_file, outprefix]
            return {"arguments": args, "return_value": outprefix}

        #outprefix = pcajl.nonblocking(self.ex, counts_table_file, via=self.via).wait()
        outprefix = pcajl(self.ex, counts_table_file)
        pca_descr_png = set_file_descr('pca_groups.png', type='png', step='pca')
        pca_descr_js = set_file_descr('pca_groups.js', type='txt', step='pca', view='admin')
        pcaeigv_descr_png = set_file_descr('pca_groups_sdev.png', type='png', step='pca')
        pcaeigv_descr_js = set_file_descr('pca_groups_sdev.js', type='txt', step='pca', view='admin')
        self.ex.add(outprefix+'.png', description=pca_descr_png)
        self.ex.add(outprefix+'.js', description=pca_descr_js)
        self.ex.add(outprefix+'_sdev.png', description=pcaeigv_descr_png)
        self.ex.add(outprefix+'_sdev.js', description=pcaeigv_descr_js)


#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#

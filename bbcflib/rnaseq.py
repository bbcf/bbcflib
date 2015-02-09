"""
======================
Module: bbcflib.rnaseq
======================

Methods of the bbcflib's RNA-seq worflow. The main function is ``rnaseq_workflow()``.

From a BAM file produced by an alignement on the genome, gets counts of reads
on the exons, and uses least-squares to infer counts on genes and transcripts.
"""

# Built-in modules #
import os, sys, math, itertools, shutil
from operator import itemgetter

# Internal modules #
from bbcflib.common import set_file_descr, unique_filename_in, cat, timer, program_exists
from bbcflib.mapseq import add_and_index_bam, sam_to_bam
from bbcflib.gfminer.common import map_chromosomes, duplicate, apply
from bbcflib.track import track
from bein import program

# Other modules #
import numpy
from numpy import asarray
import pysam

numpy.set_printoptions(precision=3,suppress=True)
numpy.seterr(invalid='print')
numpy.seterr(divide='ignore')

_TEST_ = 0


#----------------------------- GLOBAL CLASSES -----------------------------#


class RNAseq(object):
    """Abstract, inherited by different parts of the workflow."""
    def __init__(self,ex,via,job,assembly,conditions,debugfile,logfile,
                 pileup_level,junctions,stranded):
        self.ex = ex                     # bein.Execution object
        self.job = job                   # frontend.Job object
        self.assembly = assembly         # genrep.Assembly object
        self.conditions = conditions     # list of <group_name>.<run_id>
        self.pileup_level = pileup_level # ["genes","exons","transcripts"], or another combination
        self.via = via                   # "local"/"lsf"
        self.junctions = junctions       # True/False, want to find or not
        self.debugfile = debugfile       # debug file name
        self.logfile = logfile           # log file name
        self.stranded = stranded         # True/False, strand-specific or not
            # number of columns in the output, = number of runs if not stranded, twice otherwise

    def write_log(self,s):
        self.logfile.write(s+'\n'); self.logfile.flush()

    def write_debug(self,s):
        self.debugfile.write(s+'\n'); self.debugfile.flush()


#--------------------------- GTF CREATION ----------------------------#


def gtf_from_bam_header(bam):
    """In case of alignment on a custom sequence."""
    bamtrack = track(bam,format='bam')
    gtf = unique_filename_in()+'.gtf'
    with open(gtf,"wb") as g:
        for c,meta in bamtrack.chrmeta.iteritems():
            n = c.split("|")[1] if "|" in c else c
            gtfline = '\t'.join([c,'','exon','1',str(meta['length']),'.','.','.',
                    'exon_id "%s"; transcript_id "%s"; gene_id "%s"; gene_name "%s"' % (c,c,c,n)])+'\n'
            g.write(gtfline)
    bamtrack.close()
    return gtf

def transcriptome_gtf_from_genrep(assembly):
    """In case of mapping on the transcriptome - it if still ever happens."""
    tmap = assembly.get_transcript_mapping()
    gtf = unique_filename_in()
    gtflines = []
    smap = {1:'+', -1:'-'}
    with open(gtf,"wb") as g:
        for tid,t in tmap.iteritems():
            gtflines.append([t.chrom,'Ensembl','exon',t.start,t.end,'.',smap.get(t.strand,'.'),'.',
                           'exon_id "%s"; transcript_id "%s"; gene_id "%s"; gene_name "%s"' \
                           % (t.id,t.id,t.gene_id,t.gene_name)])
        del tmap
        gtflines.sort(key=itemgetter(0,3,4))  # chrom,start,end
        for gtfline in gtflines:
            g.write('\t'.join([str(x) for x in gtfline])+'\n')
    return gtf


#----------------------------- MAIN CALL -----------------------------#


@timer
def rnaseq_workflow(ex, job, pileup_level=["genes","transcripts"],
                    via="lsf", junctions=False, stranded=False,
                    logfile=sys.stdout, debugfile=sys.stderr):
    """Main function of the workflow.

    :rtype: None
    :param ex: a bein execution.
    :param job: a Frontend.Job object (or a dictionary of the same form).
    :param assembly: a genrep.Assembly object
    :param junctions: (bool) whether to search for splice junctions using SOAPsplice. [False]
    :param via: (str) send job via 'local' or 'lsf'. ["lsf"]
    """
    group_names={}; conditions=[]
    groups = job.groups
    assembly = job.assembly
    assert len(groups) > 0, "No groups/runs were given."
    for gid,group in groups.iteritems():
        gname = str(group['name'])
        group_names[gid] = gname
    if isinstance(pileup_level,basestring): pileup_level=[pileup_level]

    # Define conditions as 'group_name.run_id' and store bamfiles in the same order
    bamfiles = []
    for gid,files in job.files.iteritems():
        k = 0
        for rid,f in files.iteritems():
            k+=1
            cond = group_names[gid]+'.'+str(k)
            conditions.append(cond)
            bamfiles.append(f['bam'])
    ncond = len(conditions)

    # Get the assembly's GTF
    # ...from fasta origin
    logfile.write("* Prepare GTF\n"); logfile.flush()
    if hasattr(assembly,"fasta_origin"):
        logfile.write("  ... from fasta origin\n"); logfile.flush()
        gtf = gtf_from_bam_header(bamfiles[0])
        pileup_level = ["transcripts"]
        if stranded:
            stranded=False
            logfile.write("  ... Cannot exploit strand information from custom fasta reference.\n"); logfile.flush()
    # ... or from (wrong) mapping on the transcriptome
    elif assembly.intype==2:
        logfile.write("  ... from mapping on the transcriptome\n"); logfile.flush()
        gtf = transcriptome_gtf_from_genrep(assembly)
    # ... or from config file
    else:
        gtf = job.options.get('annot_file')
        if gtf and os.path.exists(os.path.join('..', gtf)):
            gtf = os.path.join('..', gtf)
            logfile.write("  ... from config file: %s\n" % gtf); logfile.flush()
        elif gtf and os.path.exists(gtf):
            gtf = os.path.abspath(gtf)
            logfile.write("  ... from config file: %s\n" % gtf); logfile.flush()
    # ... or from GenRep
        else:
            logfile.write("  ... from GenRep\n"); logfile.flush()
            gtf = assembly.create_exome_gtf()
    descr = set_file_descr(gtf, type='txt', step='pileup', view='admin')
    ex.add(gtf, description=descr)
    #shutil.copy(gtf,"../")

    # Build controllers
    rnaseq_args = (ex,via,job,assembly,conditions,debugfile,logfile,
                   pileup_level,junctions,stranded)
    CNT = Counter(*rnaseq_args)
    DE = DE_Analysis(*rnaseq_args)
    PCA = Pca(*rnaseq_args)
    JN = Junctions(*rnaseq_args)

    # Count reads on genes, transcripts with "rnacounter"
    count_files = CNT.count_reads(bamfiles, gtf)

    def differential_analysis(counts_file, feature_type):
        counts_file = DE.clean_before_deseq(counts_file)
        diff_files = DE.differential_analysis(counts_file)
        if diff_files is not None:
            diff_files = [DE.clean_deseq_output(o) for o in diff_files]
            for diff in diff_files:
                head = open(diff).readline().strip().replace(' ','')
                oname = feature_type + "_differential_"+ head + ".txt"
                desc = set_file_descr(oname, step='stats', type='txt')
                ex.add(diff, description=desc)

    # DE and PCA
    if "genes" in pileup_level:
        # PCA of groups ~ gene expression
        description = set_file_descr("genes_expression.txt", step="pileup", type="txt")
        ex.add(count_files['genes'], description=description)
        differential_analysis(count_files['genes'], "genes")
        if stranded:
            description = set_file_descr("genes_antisense_expression.txt", step="pileup", type="txt")
            ex.add(count_files['genes_anti'], description=description)
            differential_analysis(count_files['genes_anti'], "genes_antisense")
        if ncond > 2:
            PCA.pca_rnaseq(count_files['genes'])

    if "transcripts" in pileup_level:
        description = set_file_descr("transcripts_expression.txt", step="pileup", type="txt")
        ex.add(count_files['transcripts'], description=description)
        differential_analysis(count_files['transcripts'], "transcripts")
        if stranded:
            description = set_file_descr("transcripts_antisense_expression.txt", step="pileup", type="txt")
            ex.add(count_files['transcripts_anti'], description=description)
            differential_analysis(count_files['transcripts_anti'], "transcripts_antisense")

    # Find splice junctions
    if junctions:
        logfile.write("* Search for splice junctions\n"); logfile.flush()
        JN.find_junctions()

    return 0


#--------------------------------- COUNTING ----------------------------------#


@program
def rnacounter(bam,gtf, options=[]):
    """Counts reads from *bam* in every feature of *gtf* and infers
    genes and transcripts expression."""
    opts = [] + options
    args = ["rnacounter"] + opts + [os.path.abspath(bam), os.path.abspath(gtf)]
    return {"arguments": args, "return_value": None}

@program
def rnacounter_join(files):
    """Run `rnacounter join` to merge the output files"""
    args = ["rnacounter", "join"] + files
    return {"arguments": args, "return_value": None}

class Counter(RNAseq):
    def __init__(self,*args):
        RNAseq.__init__(self,*args)

    @timer
    def count_reads(self, bamfiles, gtf):
        self.write_log("* Counting reads")

        # Count reads on genes, transcripts with "rnacounter"
        ncond = len(self.conditions)
        tablenames = [None]*ncond
        futures = [None]*ncond
        max_rlen = 0
        counter_options = []
        for bam in bamfiles:
            sam = pysam.Samfile(bam,'rb')
            max_rlen = max(max_rlen, sam.next().rlen)
        counter_options += ["--exon_cutoff", str(max_rlen)]
        bwt_args = self.job.options.get('map_args',{}).get('bwt_args',[])
        if not "--local" in bwt_args:
            counter_options += ["--nh"]
        if hasattr(self.assembly,"fasta_origin") or self.assembly.intype==2:
            counter_options += ["--type","transcripts", "--method","raw"]
        else:
            counter_options += ["--type","genes,transcripts", "--method","raw,nnls"]
        if self.stranded:
            counter_options += ["--stranded"]
        for i,c in enumerate(self.conditions):
            tablenames[i] = unique_filename_in()
            futures[i] = rnacounter.nonblocking(self.ex, bamfiles[i], gtf, stdout=tablenames[i], via=self.via,
                               options=counter_options)

        # Put samples together
        for i,c in enumerate(self.conditions):
            try:
                futures[i].wait()
            except Exception as err:
                self.write_debug("Counting failed: %s." % str(err))
                raise err
            if futures[i] is None:
                self.write_debug("Counting failed.")
                raise ValueError("Counting failed.")
            # Keep intermediate tables
            #shutil.copy(tablenames[i], "../counts%d.txt"%i)
            descr = set_file_descr(self.conditions[i]+'_'+tablenames[i]+'.txt', type='txt', step='pileup', view='admin')
            self.ex.add(tablenames[i], description=descr)
        joined = unique_filename_in()
        rnacounter_join.nonblocking(self.ex, tablenames, stdout=joined, via=self.via).wait()

        # Split genes and transcripts into separate files
        genes_filename = unique_filename_in()
        trans_filename = unique_filename_in()
        genes_file = open(genes_filename,"wb")
        trans_file = open(trans_filename,"wb")
        if self.stranded:
            genes_anti_filename = unique_filename_in()
            trans_anti_filename = unique_filename_in()
            genes_anti_file = open(genes_anti_filename,"wb")
            trans_anti_file = open(trans_anti_filename,"wb")
        with open(joined) as jfile:
            header = jfile.readline()
            hconds = ["counts."+c for c in self.conditions] + ["rpkm."+c for c in self.conditions]
            hinfo = header.strip().split('\t')[2*ncond+1:]
            header = '\t'.join(["ID"] + hconds + hinfo)+'\n'
            genes_file.write(header)
            trans_file.write(header)
            type_idx = header.split('\t').index("Type")
            if self.stranded:
                genes_anti_file.write(header)
                trans_anti_file.write(header)
                sense_idx = header.split('\t').index("Sense")
                for line in jfile:
                    L = line.split('\t')
                    ftype = L[type_idx].lower()
                    sense = L[sense_idx].lower()
                    if ftype == 'gene':
                        if sense == 'antisense':
                            genes_anti_file.write(line)
                        else:
                            genes_file.write(line)
                    elif ftype == 'transcript':
                        if sense == 'antisense':
                            trans_anti_file.write(line)
                        else:
                            trans_file.write(line)
            else:
                for line in jfile:
                    L = line.split('\t')
                    ftype = L[type_idx].lower()
                    if ftype == 'gene':
                        genes_file.write(line)
                    elif ftype == 'transcript':
                        trans_file.write(line)
        genes_file.close()
        trans_file.close()
        if self.stranded:
            count_files = {'genes':genes_filename, 'transcripts':trans_filename,
                           'genes_anti':genes_anti_filename, 'transcripts_anti':trans_anti_filename}
        else:
            count_files = {'genes':genes_filename, 'transcripts':trans_filename}
        return count_files


#-------------------------- DIFFERENTIAL ANALYSIS ----------------------------#


class DE_Analysis(RNAseq):
    def __init__(self,*args):
        RNAseq.__init__(self,*args)

    def clean_before_deseq(self, filename, keep=0.6):
        """Delete all lines of *filename* where counts are 0 in every run.
        Return only the first column with added information concatenated,
        and the Count columns.

        :param keep: fraction of highest counts to keep. [0.6]
        """
        filename_clean = unique_filename_in()
        ncond = len(self.conditions)
        if ncond >1:
            with open(filename) as f:
                header = f.readline().split('\t')
                data = [x.strip().split('\t') for x in f]
                if not data: return
                rownames = asarray(['%s|%s|%s' % (x[0],x[-1],x[-2]) for x in data])  # id,gene_name,strand
                M = asarray([x[1:1+ncond] for x in data], dtype=numpy.float_)
                colnames = header[0:1]+header[1:1+ncond] # 'counts' column names, always
                # Remove 40% lowest counts
                sums = numpy.sum(M,1)
                filter = asarray([x[1] for x in sorted(zip(sums,range(len(sums))))])
                filter = filter[:len(filter)/keep]
                M = M[filter]
                rownames = rownames[filter]
                # Create the input tab file
                with open(filename_clean,"wb") as g:
                    header = '\t'.join(colnames)+'\n'
                    g.write(header)
                    for i,scores in enumerate(M):
                        if any(scores):
                            line = rownames[i] + '\t' + '\t'.join([str(x) for x in scores]) + '\n'
                            g.write(line)
        return filename_clean

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
    def differential_analysis(self, filename):
        """Launch an analysis of differential expression on the count
        values, and saves the output in the MiniLIMS."""

        @program
        def run_deseq(data_file, options=[]):
            """Run negbin.test.R on *data_file*."""
            output_file = unique_filename_in()
            opts = ["-o",output_file]+options
            arguments = ["negbin.test.R", data_file]+opts
            return {'arguments': arguments, 'return_value': output_file}

        if not program_exists('negbin.test.R'):
            self.write_debug("Skipped DE analysis: negbin.test.R not found.")
            return
        if filename is None:
            self.write_log("  Skipped differential analysis: empty counts file.")
            return
        ncond = len(self.conditions)
        if ncond < 2:
            self.write_log("  Skipped differential analysis: less than two groups.")
            return
        else:
            self.write_log("* Differential analysis")
            options = ['-s','tab']
            try:
                deseqfile = run_deseq.nonblocking(self.ex, filename, options, via=self.via).wait()
            except Exception as err:
                self.write_debug("DE analysis failed: %s." % str(err))
                return
            if deseqfile is None:
                self.write_debug("DE analysis failed.")
                return
            output_files = [f for f in os.listdir(self.ex.working_directory) if deseqfile in f]
            if (deseqfile is None) or isinstance(deseqfile,Exception) or len(output_files)==0:
                self.write_debug("Skipped differential analysis: %s." % str(deseqfile))
                return
            self.write_log("  ....done.")
            return output_files


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

        if not program_exists('soapsplice'):
            self.write_debug("Skipped junctions search: soapsplice not found.")
            return
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
            try:
                template = future.wait()
            except Exception as err:
                self.write_debug("SOAPsplice failed: %s." % str(err))
                return
            if template is None:
                self.write_debug("SOAPsplice failed.")
                return
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
            except Exception as e:
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


#---------------------------------- PCA ----------------------------------#


class Pca(RNAseq):
    def __init__(self,*args):
        RNAseq.__init__(self,*args)

    def pca_rnaseq(self,counts_table_file):
        @program
        def pca(counts_table_file):
            outprefix = unique_filename_in()
            args = ['pca.R', counts_table_file, outprefix, "rpkm"]
            return {"arguments": args, "return_value": outprefix}

        if not program_exists('pca.R'):
            self.write_debug("Skipped PCA: pca.R not found.")
            return
        try:
            self.write_log("* PCA")
            outprefix = pca.nonblocking(self.ex, counts_table_file, via=self.via).wait()
        except Exception as err:
            self.write_debug("PCA failed: %s." % str(err))
            return
        if outprefix is None:
            self.write_debug("PCA failed.")
            return
        pca_descr_pdf = set_file_descr('pca.pdf', type='pdf', step='pca')
        self.ex.add(outprefix+'.pdf', description=pca_descr_pdf)


#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#

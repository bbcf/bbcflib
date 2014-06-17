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
import json, urllib2

# Internal modules #
from bbcflib.common import set_file_descr, unique_filename_in, cat, timer, program_exists
from bbcflib.common import gunzipfile
from bbcflib.genrep import Assembly
from bbcflib.mapseq import add_and_index_bam, sam_to_bam
from bbcflib.gfminer.common import map_chromosomes, duplicate, apply
from bbcflib.track import track
from bein import program

# Other modules #
import numpy
from numpy import asarray

numpy.set_printoptions(precision=3,suppress=True)
numpy.seterr(invalid='print')
numpy.seterr(divide='ignore')

test = False


#----------------------------- GLOBAL CLASSES -----------------------------#


class RNAseq(object):
    """Abstract, inherited by different parts of the workflow."""
    def __init__(self,ex,via,job,assembly,conditions,debugfile,logfile,
                 pileup_level,rpath,juliapath,junctions,stranded):
        self.ex = ex                     # bein.Execution object
        self.job = job                   # frontend.Job object
        self.assembly = assembly         # genrep.Assembly object
        self.conditions = conditions     # list of <group_name>.<run_id>
        self.pileup_level = pileup_level # ["genes","exons","transcripts"], or another combination
        self.via = via                   # "local"/"lsf"
        self.rpath = rpath               # Path to R scripts
        self.juliapath = juliapath       # Path to Julia scripts
        self.junctions = junctions       # True/False, want to find or not
        self.debugfile = debugfile       # debug file name
        self.logfile = logfile           # log file name
        self.stranded = stranded         # True/False, strand-specific or not
        self.ncond = 2*len(self.conditions) if self.stranded else len(self.conditions)
            # number of columns in the output, = number of runs if not stranded, twice otherwise

    def write_log(self,s):
        self.logfile.write(s+'\n'); self.logfile.flush()

    def write_debug(self,s):
        self.debugfile.write(s+'\n'); self.debugfile.flush()


#----------------------------- MAIN CALL -----------------------------#


@timer
def rnaseq_workflow(ex, job, assembly=None, pileup_level=["genes","transcripts"], via="lsf",
                    rpath=None, juliapath=None, junctions=False, stranded=False,
                    logfile=sys.stdout, debugfile=sys.stderr):
    """Main function of the workflow.

    :rtype: None
    :param ex: the bein's execution Id.
    :param job: a Frontend.Job object (or a dictionary of the same form).
    :param assembly: a genrep.Assembly object
    :param rpath: (str) path to the R executable.
    :param junctions: (bool) whether to search for splice junctions using SOAPsplice. [False]
    :param via: (str) send job via 'local' or 'lsf'. ["lsf"]
    """

    # While environment not properly set on V_IT
    if juliapath is None:
        juliapath='/home/jdelafon/repos/bbcfutils/Julia/'

    if test:
        via = 'local'
        logfile = sys.stdout
        debugfile = sys.stdout
        repo_rpath = '/home/jdelafon/repos/bbcfutils/R/'
        juliapath = '/home/jdelafon/repos/bbcfutils/Julia/'
        if os.path.exists(repo_rpath): rpath = repo_rpath

    if assembly is None:
        assembly = Assembly(assembly=job.assembly_id)

    group_names={}; group_ids={}; conditions=[]
    groups = job.groups
    assert len(groups) > 0, "No groups/runs were given."
    for gid,group in groups.iteritems():
        gname = str(group['name'])
        group_names[gid] = gname
        group_ids[gname] = gid
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
        logfile.write("* ... from fasta origin\n"); logfile.flush()
        firstbam = bamfiles[0]
        firstbamtrack = track(firstbam,format='bam')
        gtf = unique_filename_in()+'.gtf'
        with open(gtf,"wb") as g:
            for c,meta in firstbamtrack.chrmeta.iteritems():
                g.write('\t'.join([c,'','exon','0',str(meta['length']),'.','.','.',
                        'exon_id %s; transcript_id %s; gene_id %s; gene_name %s' % (c,)*4])+'\n')
        firstbamtrack.close()
    # ... or from config file
    else:
        gtf = job.files.itervalues().next().itervalues().next().get('gtf')  # from config file
        if gtf: logfile.write("* ... from config file: %s\n" % gtf); logfile.flush()
    # ... or from GenRep
    if gtf is None:
        genrep_root = assembly.genrep.root
        gtf = os.path.join(genrep_root,"nr_assemblies/gtf/%s_%i.gtf.gz"%(assembly.md5,0))
        logfile.write("* ... from GenRep: %s\n" % gtf); logfile.flush()
        for n in range(1,100):
            gtf_n = os.path.join(genrep_root,"nr_assemblies/gtf/%s_%i.gtf.gz"%(assembly.md5,n))
            if not os.path.exists(gtf_n): break
            gtf = gtf_n
        shutil.copy(gtf,ex.working_directory)
        gtf = os.path.join(ex.working_directory, os.path.basename(gtf))
        gtf = gunzipfile(ex, gtf)[0]

        # Resolve chromosome names between assembly (BAM) and GTF
        chromosomes = ','.join(assembly.chrmeta.keys())
        url = os.path.join(assembly.genrep.url,
                  "chromosomes/get_mapping.json?identifier=%s" % chromosomes \
                + "&assembly_name=%s&source_name=Ensembl" % assembly.name)
        request = urllib2.Request(url)
        chrom_json = json.load(urllib2.urlopen(request))
        chrom_mapping = dict((v['chr_name'][0],k) for k,v in chrom_json.iteritems())
        gtf_renamed = os.path.splitext(gtf)[0]+'_renamed.gtf'
        with open(gtf) as f:
            with open(gtf_renamed,"wb") as g:
                for line in f:
                    L = line.strip().split('\t')
                    chrom = chrom_mapping.get(L[0])
                    if chrom is None: continue  # chrom name not found in GenRep
                    g.write('\t'.join([chrom]+L[1:])+'\n')
        gtf = gtf_renamed
    assert os.path.exists(gtf), "GTF not found: %s" %gtf

    # Build controllers
    rnaseq_args = (ex,via,job,assembly,conditions,debugfile,logfile,
                   pileup_level,rpath,juliapath,junctions,stranded)
    CNT = Counter(*rnaseq_args)
    DE = DE_Analysis(*rnaseq_args)
    PCA = Pca(*rnaseq_args)
    JN = Junctions(*rnaseq_args)

    # Count reads on genes, transcripts with "rnacounter"
    genes_filename,trans_filename = CNT.count_reads(bamfiles, gtf)

    # DE and PCA
    if "genes" in pileup_level:
        DE.differential_analysis(genes_filename, "genes")
        # PCA of groups ~ gene expression
        if ncond > 2:
            try:
                PCA.pca_rnaseq(genes_filename)
            except Exception, error:
                debugfile.write("PCA failed: %s\n"%str(error)); debugfile.flush()
        description = set_file_descr("genes_expression.tab", step="pileup", type="txt")
        ex.add(genes_filename, description=description)

    if "transcripts" in pileup_level:
        DE.differential_analysis(trans_filename, "transcripts")
        description = set_file_descr("transcripts_expression.tab", step="pileup", type="txt")
        ex.add(trans_filename, description=description)

    # Find splice junctions
    if junctions:
        logfile.write("* Search for splice junctions\n"); logfile.flush()
        try:
            JN.find_junctions()
        except Exception, error:
            debugfile.write("Junctions search failed: %s\n"%str(error)); debugfile.flush()

    return 0


#--------------------------------- COUNTING ----------------------------------#


@program
def rnacounter(bam,gtf, options=[]):
    """Counts reads from *bam* in every feature of *gtf* and infers
    genes and transcripts expression."""
    opts = [] + options
    args = ["rnacounter"] + opts + [os.path.abspath(bam), os.path.abspath(gtf)]
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
        counter_options = ["--type","genes,transcripts", "--method","raw,nnls", "--multiple"]
        if self.stranded: counter_options += ["--stranded"]
        for i,c in enumerate(self.conditions):
            tablenames[i] = unique_filename_in()
            futures[i] = rnacounter.nonblocking(self.ex, bamfiles[i], gtf, stdout=tablenames[i], via=self.via,
                               options=counter_options)

        # Construct file headers
        if self.stranded:
            hconds = ["counts."+c for c in self.conditions] + ["counts_rev."+c for c in self.conditions] \
                   + ["rpkm."+c for c in self.conditions] + ["rpkm_rev."+c for c in self.conditions]
        else:
            hconds = ["counts."+c for c in self.conditions] + ["rpkm."+c for c in self.conditions]
        hinfo = ["Chromosome","Start","End","Strand","GeneName"]
        header = ["ID"] + hconds + hinfo

        # Put samples together, split genes and transcripts
        tracks = [None]*ncond
        for i,c in enumerate(self.conditions):
            futures[i].wait()
            tracks[i] = track(tablenames[i], format="txt", fields=["ID","Count","RPKM"]+hinfo+['type']).read()
        genes_filename = unique_filename_in()
        trans_filename = unique_filename_in()
        genes_file = open(genes_filename,"wb")
        trans_file = open(trans_filename,"wb")
        genes_file.write('\t'.join(header)+'\n')
        trans_file.write('\t'.join(header)+'\n')
        while 1:
            info = tuple(t.next() for t in tracks)
            if not info: break
            fid = info[0][0]
            rest = info[0][3:-1]
            ftype = info[0][-1]
            counts = tuple(x[1] for x in info)
            rpkm = tuple(x[2] for x in info)
            towrite = '\t'.join((fid,)+counts+rpkm+rest)+'\n'
            if ftype == 'gene':
                genes_file.write(towrite)
            elif ftype=='transcript':
                trans_file.write(towrite)
        genes_file.close()
        trans_file.close()
        return genes_filename, trans_filename


#-------------------------- DIFFERENTIAL ANALYSIS ----------------------------#


class DE_Analysis(RNAseq):
    def __init__(self,*args):
        RNAseq.__init__(self,*args)

    def clean_before_deseq(self, filename, keep=0.6):
        """Delete all lines of *filename* where counts are 0 in every run.

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
    def differential_analysis(self, filename, feature_type):
        """Launch an analysis of differential expression on the count
        values, and saves the output in the MiniLIMS."""

        @program
        def run_deseq(data_file, script_path, options=[]):
            """Run *rpath*/negbin.test.R on *data_file*."""
            output_file = unique_filename_in()
            opts = ["-o",output_file]+options
            arguments = ["R","--slave","-f",script_path,"--args",data_file]+opts
            return {'arguments': arguments, 'return_value': output_file}

        ncond = len(self.conditions)
        res_file = self.clean_before_deseq(filename)
        if res_file is None:
            self.write_log("  Skipped differential analysis: empty counts file for %s." % feature_type)
        elif ncond < 2:
            self.write_log("  Skipped differential analysis: less than two groups.")
        else:
            self.write_log("* Differential analysis")
            options = ['-s','tab']
            script_path = os.path.join(self.rpath,'negbin.test.R')
            if not os.path.exists(script_path):
                self.write_debug("DE: R path not found: '%s'" % self.rpath)
            try:
                deseqfile = run_deseq.nonblocking(self.ex, res_file, script_path, options, via=self.via).wait()
            except Exception as exc:
                self.write_log("  Skipped differential analysis: %s" % exc)
                return
            output_files = [f for f in os.listdir(self.ex.working_directory) if deseqfile in f]
            if (not deseqfile) or len(output_files)==0:
                self.write_log("  ....Empty file.")
                return
            for o in output_files:
                oname = feature_type+"_differential"+o.split(deseqfile)[1]+".txt"
                desc = set_file_descr(oname, step='stats', type='txt')
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

        outprefix = pcajl.nonblocking(self.ex, counts_table_file, via=self.via).wait()
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

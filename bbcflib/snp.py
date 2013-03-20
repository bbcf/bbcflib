"""
===================
Module: bbcflib.snp
===================

From a set of BAM files produced by an alignement on the genome, calls snps and annotates them 
with respect to a set of coding genes on the same genome.
"""
# Built-in modules #
import re, os, sys, tarfile
from itertools import product

# Internal modules #
from bbcflib.common import unique_filename_in, set_file_descr, sam_faidx
from bbcflib import mapseq
from bbcflib.btrack import FeatureStream, track
from bbcflib.bFlatMajor.common import concat_fields
from bbcflib.bFlatMajor import stream as gm_stream
from bein import program

_iupac = {'M':'AC', 'Y':'CT', 'R':'AG', 'S':'CG', 'W': 'AT', 'K':'GT', 
          'B':'CGT', 'D':'AGT', 'H':'ACT', 'V':'ACG',
          'N': 'ACGT'}

_translate = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L/START",
              "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
              "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M & START",
              "GTT": "V", "GTA": "V", "GTC": "V", "GTG": "V/START",
              "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
              "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
              "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
              "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
              "TAT": "Y", "TAC": "Y", "TAA": "STOP", "TAG": "STOP",
              "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
              "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
              "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
              "TGT": "C", "TGC": "C", "TGA": "STOP", "TGG": "W",
              "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
              "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
              "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G" }

def _ploidy(assembly):
    if str(assembly.name) in ["EB1_e_coli_k12","MLeprae_TN","mycoSmeg_MC2_155",
                              "mycoTube_H37RV","NA1000","vibrChol1","TB40-BAC4"]:
        ploidy = 1 # procaryote
    elif str(assembly.name) in ["araTha1"]:
        ploidy = 3
    else:
        ploidy = 2 # eucaryote
    return ploidy

@program
def sam_pileup(assembly,bamfile,refGenome):
    """Binds 'samtools pileup'.

    :param assembly: Genrep.Assembly object.
    :param bamfile: path to the BAM file.
    :param refGenome: path to the reference genome fasta file.
    """
    ploidy = _ploidy(assembly)
    return {"arguments": ["samtools","pileup","-Bcvsf",refGenome,"-N",str(ploidy),bamfile],
            "return_value": None}

def find_snp(info,mincov,minsnp,assembly):
    """Parse the output of samtools pileup, provided as a list of already split *info*::

        info = ['chrV', '91668', 'G', 'R', '4', '4', '60', '4', 'a,.^~.', 'LLLL', '~~~~\n']
        info = [chr, pos, ref, consensus, cons_qual, snp_qual, max_map_qual, nreads, 'a,.^~.', 'LLLL', '~~~~\n']
                 0    1    2       3          4         5           6           7       8         9       10
        cons_qual : consensus quality is the Phred-scaled probability that the consensus is wrong.
        snp_qual: SNP quality is the Phred-scaled probability that the consensus is identical to the reference.

        In case of indels:
        info = ['chrV', '1128659', '*', '*/+CT', '5', '5', '37', '11', '*', '+CT', '10', '1', '0', '0', '0']
        info = [chr, pos, *, c1/c2, cons_qual, snp_qual, max_map_qual, nreads, c1, c2, nindel, nref, n3rd]
                 0    1   2    3        4         5           6           7    8   9     10     11    12
        c1/c2: consensus on strand1/strand2, of the form \[+-][0-9]+[ACGTNacgtn]+ .
        nindel/nref/n3rd: number of reads supporting the indel/reference/3rd allele.
        `http://samtools.sourceforge.net/cns0.shtml`_
    """
    ref = info[2].upper()  # reference base
    cons = info[3].upper() # consensus base
    nreads = int(info[7])  # total coverage at this position
    denom = 100/float(nreads)
    ploidy = _ploidy(assembly)
    #snp_qual = 10**(-float(info[5])/10)
    consensus = []
    if ref == "*": # indel
        nindel = int(info[10])
        nref = int(info[11])
        if nindel >= mincov and 100*nindel*ploidy >= minsnp*nreads:
            consensus.append("%s/%s (%.2f%% of %s)" %(info[8],info[9],denom*nindel,nreads))
    else:
#        nref = info[8].count(".") + info[8].count(",") # number of reads supporting wild type (.fwd and ,rev)
        for char in _iupac.get(cons,cons):
            if char == ref: continue
            nvar = info[8].count(char)+info[8].count(char.lower())
            if nvar >= mincov and 100*nvar*ploidy >= minsnp*nreads:
                consensus.append("%s (%.2f%% of %s)" %(char,denom*nvar,nreads))
    consensus = ", ".join(consensus) or ref
    return consensus

def all_snps(ex,chrom,dictPileup,outall,assembly,sample_names,mincov,minsnp):
    """For a given chromosome, returns a summary file containing all SNPs identified
    in at least one of the samples.
    Each row contains: chromosome id, SNP position, reference base, SNP base (with proportions)

    :param chrom: (str) chromosome name.
    :param dictPileup: (dict) dictionary of the form {filename: (bein.Future, sample_name, bam_filename)}.
    :param outall: (str) name of the file that will contain the list of all SNPs.
    :param sample_names: (list of str) list of sample names.
    :param mincov: (int) Minimum number of reads supporting an SNP at a position for it to be considered. [5]
    :param minsnp: (int) Minimum percentage of reads supporting the SNP for it to be returned.
        N.B.: Effectively, half of it on each strand for diploids. [40]
    """
    snames = []
    pileup_files = []
    bam_tracks = []
    allsnps = []
    tarname = unique_filename_in()
    tarfh = tarfile.open(tarname, "w:gz")
    for pileup_filename,trio in dictPileup.iteritems():
        trio[0].wait() #file p is created from samtools pileup
        snames.append(trio[1])
        pileup_files.append(open(pileup_filename))
        tarfh.add(pileup_filename,arcname=trio[1]+"_"+chrom+".pileup")
        bam_tracks.append(track(trio[2],format='bam'))
    tarfh.close()
    ex.add( tarname, description=set_file_descr("pileups_%s.tgz"%chrom,step="pileup",type="tar",view='admin') )
    nsamples = len(snames)
    sorder = [snames.index(x) for x in sample_names]
    current = [pileup_files[i].readline().split('\t') for i in sorder]
    lastpos = 0
    while any([len(x)>1 for x in current]):
        current_pos = [int(x[1]) if len(x)>1 else sys.maxint for x in current]
        pos = min(current_pos)
        current_snp_idx = set(i for i in range(nsamples) if current_pos[i]==pos)
        current_snps = [None]*nsamples
        for i in current_snp_idx:
            info = current[i]
            chrbam = info[0]
            ref = info[2].upper()
            current_snps[i] = find_snp(info,mincov,minsnp,assembly)
        for i in set(range(nsamples))-current_snp_idx:
            try: coverage = bam_tracks[sorder[i]].coverage((chrbam,pos-1,pos)).next()[-1]
            except StopIteration: coverage = 0
            current_snps[i] = ref if coverage else "0"
        if not all([s in ["0",ref] for s in current_snps]):
            if pos != lastpos: # indel can be located at the same position as an SNP
                allsnps.append((chrom,pos-1,pos,ref)+tuple(current_snps))
            lastpos = pos
        for i in current_snp_idx:
            current[i] = pileup_files[sorder[i]].readline().split('\t')
    for f in pileup_files: f.close()
    for b in bam_tracks: b.close()

    snp_read = FeatureStream(allsnps, fields=['chr','start','end','name']+sample_names)
    annotation = assembly.gene_track(chrom)
    annotated_stream = gm_stream.getNearestFeature(snp_read,annotation,
                                                   thresholdPromot=3000,thresholdInter=3000,thresholdUTR=10)
    with open(outall,"a") as fout:
        for snp in annotated_stream:
            # snp: ('chrV',154529, 154530, 'T', 'A', '* A', 'YER002W|NOP16_YER001W|MNN1', 'Upstream_Included', '2271_1011')
            fout.write('\t'.join([str(x) for x in (snp[0],)+snp[2:]])+'\n')
            # chrV  1606    T   43.48% C / 56.52% T YEL077W-A|YEL077W-A_YEL077C|YEL077C 3UTR_Included   494_-2491
    return allsnps

def exon_snps(chrom,outexons,allsnps,assembly,sample_names,genomeRef={}):
    """Annotates SNPs described in `filedict` (a dictionary of the form {chromosome: filename}
    where `filename` is an output of parse_pileupFile).
    Adds columns 'gene', 'location_type' and 'distance' to the output of parse_pileupFile.
    Returns two files: the first contains all SNPs annotated with their position respective to genes in
    the specified assembly, and the second contains only SNPs found within CDS regions.

    :param chrom: (str) chromosome name.
    :param outexons: (str) name of the file containing the list of SNPs on exons.
    :param allsnps: list of tuples (chr,start,end,ref) as returned by all_snps().
    :param assembly: genrep.Assembly object
    :param sample_names: list of sample names.
    :param genomeRef: dict of the form {'chr1': filename}, where filename is the name of a fasta file
        containing the reference sequence for the chromosome.
    """
    def _revcomp(seq):
        cmpl = dict((('A','T'),('C','G'),('T','A'),('G','C'),
                     ('a','t'),('c','g'),('t','a'),('g','c'),
                     ('M','K'),('K','M'),('Y','R'),('R','Y'),('S','S'),('W','W'),
                     ('m','k'),('k','m'),('y','r'),('r','y'),('s','s'),('w','w'),
                     ('B','V'),('D','H'),('H','D'),('V','B'),
                     ('b','v'),('d','h'),('h','d'),('v','b'),
                     ('N','N'),('n','n')))
        return "".join(reversed([cmpl.get(x,x) for x in seq]))

    def _write_buffer(_buffer, outex):
        new_codon = None
        for chr,pos,refbase,variants,cds,strand,ref_codon,shift in _buffer:
            if refbase == "*":
                break
            varbase = list(variants)
            if new_codon is None:
                new_codon = [[ref_codon] for _ in range(len(varbase))]
            variants = []
            for variant in varbase:
                if variant == '0':
                    variants.append([refbase])
                else: # 'C (43%),G (12%)' : heterozygous double snp
                    variants.append([v[0] for v in variant.split(',')])
            for k,v in enumerate(variants):
                cnumb = len(new_codon[k])
                newc = new_codon[k]*len(v)
                for i,vari in enumerate(v):
                    for j in range(cnumb):
                        newc[i*cnumb+j] = newc[i*cnumb+j][:shift]  +vari +newc[i*cnumb+j][shift+1:]
                        assert ref_codon[shift] == refbase, "bug with shift within codon"
                new_codon[k] = newc
        if new_codon is None: return
        if strand == -1:
            ref_codon = _revcomp(ref_codon)
            new_codon = [[_revcomp(s) for s in c] for c in new_codon]
        for chr,pos,refbase,variants,cds,strand,dummy,shift in _buffer:
            refc = [[y for y in _iupac.get(x,x)] for x in ref_codon]
            ref_codon = [''.join(x) for x in product(*refc)]
            newc = [[[[y for y in _iupac.get(x,x)] for x in variant] for variant in sample]
                    for sample in new_codon]
            new_codon = [[''.join(x) for codon in sample for x in product(*codon)] for sample in newc]
            if refbase == "*":
                result = [chr, pos+1, refbase] + list(variants) + [cds, strand] \
                         + [','.join([_translate.get(refc,'?') for refc in ref_codon])] + ["indel"]
            else:
                result = [chr, pos+1, refbase] + list(variants) + [cds, strand] \
                         + [','.join([_translate.get(refc,'?') for refc in ref_codon])] \
                         + [','.join([_translate.get(s,'?') for s in newc]) for newc in new_codon]
            outex.write("\t".join([str(r) for r in result])+"\n")

    #############################################################
    outex = open(outexons,"a")
    outex.write('#'+'\t'.join(['chromosome','position','reference']+sample_names+['exon','strand','ref_aa'] \
                                  + ['new_aa_'+s for s in sample_names])+'\n')
    snp_stream = FeatureStream(allsnps, fields=['chr','start','end','ref']+sample_names)
    inclstream = concat_fields(snp_stream, infields=snp_stream.fields[3:], as_tuple=True)
    snp_stream = FeatureStream(allsnps, fields=['chr','start','end','ref']+sample_names)
    inclstream = concat_fields(snp_stream, infields=snp_stream.fields[3:], as_tuple=True)

    annotstream = concat_fields(assembly.annot_track('CDS',chrom),
                                infields=['name','strand','frame'], as_tuple=True)
    annotstream = FeatureStream((x[:3]+(x[1:3]+x[3],) for x in annotstream),fields=annotstream.fields)
    _buffer = {1:[], -1:[]}
    last_start = {1:-1, -1:-1}
    for x in gm_stream.intersect([inclstream, annotstream]):
        # x = ('chrV',1606,1607, ('T','C (43%)', 1612,1724,'YEL077C|YEL077C',-1,0, 1712,1723,'YEL077W-A|YEL077W-A',1,0))
        nsamples = len(sample_names)
        chr = x[0]; pos = x[1]; rest = x[3]
        refbase = rest[0]
        annot = [rest[5*i+nsamples+1 : 5*i+5+nsamples+1]
                 for i in range(len(rest[nsamples+1:])/5)] # list of (start,end,cds,strand,phase)
        for es,ee,cds,strand,phase in annot:
            if strand == 1:
                shift = (pos-es-phase) % 3
            elif strand == -1:
                shift = (pos-ee+phase) % 3
            else:
                print "***",es,ee,cds,strand,phase
                continue
            codon_start = pos-shift
            ref_codon = assembly.fasta_from_regions({chr: [[codon_start,codon_start+3]]}, out={},
                                                    path_to_ref=genomeRef.get(chr))[0][chr][0]
            info = [chr,pos,refbase,list(rest[1:nsamples+1]),cds,strand,ref_codon,shift]
            # Either the codon is the same as the previous one on this strand, or it will never be.
            # Only if one codon is passed, can write its snps to a file.
            if codon_start == last_start[strand]:
                _buffer[strand].append(info)
            else:
                _write_buffer(_buffer[strand],outex)
                _buffer[strand] = [info]
                last_start[strand] = codon_start
    for strand in [1,-1]: _write_buffer(_buffer[strand],outex)
    outex.close()
    return outexons

def create_tracks(ex, outall, sample_names, assembly):
    """Write BED tracks showing SNPs found in each sample."""
    ###### TODO: tracks for coverage + quality + heterozyg. 
    infields = ['chromosome','position','reference']+sample_names+['gene','location_type','distance']
    intrack = track(outall, format='text', fields=infields, chrmeta=assembly.chrmeta,
                    intypes={'position':int})
    instream = intrack.read(fields=infields[:-3])
    outtracks = {}
    for sample_name in sample_names:
        out = unique_filename_in()+'.bed.gz'
        t = track(out,fields=['name'])
        t.make_header(name=sample_name+"_SNPs")
        outtracks[sample_name] = (t,out)

    def _row_to_annot(x,ref,n):
        if x[3+n][0] == ref: return None
        else: return "%s>%s"%(ref,x[3+n][0])

    for x in instream:
        coord = (x[0],x[1]-1,x[1])
        ref = x[2]
        snp = dict((name, _row_to_annot(x,ref,n)) for n,name in enumerate(sample_names))
        for name, tr in outtracks.iteritems():
            if snp[name]: tr[0].write([coord+(snp[name],)],mode='append')
#       convert(out+'.sql',out+'.bed.gz',mode='append')
    for name, tr in outtracks.iteritems():
        tr[0].close()
        description = set_file_descr(name+"_SNPs.bed.gz",type='bed',step='tracks',gdv='1',ucsc='1')
        ex.add(tr[1], description=description)


def snp_workflow(ex, job, assembly, minsnp=40, mincov=5, path_to_ref='', via='local', logfile=sys.stdout, debugfile=sys.stderr):
    """Main function of the workflow"""
    if path_to_ref is None:
        path_to_ref = os.path.join(assembly.genrep.root,'nr_assemblies/fasta',assembly.md5+'.tar.gz')
    assert os.path.exists(path_to_ref), "Reference sequence not found: %s." % path_to_ref
    genomeRef = assembly.untar_genome_fasta(path_to_ref, convert=True)
    pileup_dict = dict((chrom,{}) for chrom in genomeRef.keys()) # {chr: {}}
    [g.wait() for g in [sam_faidx.nonblocking(ex,f,via=via) \
                            for f in genomeRef.values()]]
    sample_names = []
    bam = {}
    print >> logfile, "* Run samtools pileup"; logfile.flush()
    for gid in sorted(job.files.keys()):
        files = job.files[gid]
        sample_name = job.groups[gid]['name']
        sample_names.append(sample_name)
        # Merge all bams belonging to the same group
        runs = [r['bam'] for r in files.itervalues()]
        if len(runs) > 1:
            bam = mapseq.merge_bam(ex,runs)
            mapseq.index_bam(ex,bam)
        else: bam = runs[0]
        # Samtools pileup
        for chrom,ref in genomeRef.iteritems():
            pileup_filename = unique_filename_in()
            future = sam_pileup.nonblocking(ex,assembly,bam,ref,via=via,stdout=pileup_filename )
            pileup_dict[chrom][pileup_filename] = (future,sample_name,bam)
        print >> logfile, "....Group %s done." % sample_name; logfile.flush()

    print >> logfile, "* Annotate SNPs"; logfile.flush()
    outall = unique_filename_in()
    outexons = unique_filename_in()
    with open(outall,"w") as fout:
        fout.write('#'+'\t'.join(['chromosome','position','reference']+sample_names+['gene','location_type','distance'])+'\n')
    for chrom in assembly.chrnames:
        dictPileup = pileup_dict[chrom]
        allsnps = all_snps(ex,chrom,dictPileup,outall,assembly,sample_names,mincov,minsnp)
        exon_snps(chrom,outexons,allsnps,assembly,sample_names,genomeRef)
    description = set_file_descr("allSNP.txt",step="SNPs",type="txt")
    ex.add(outall,description=description)
    description = set_file_descr("exonsSNP.txt",step="SNPs",type="txt")
    ex.add(outexons,description=description)

    print >> logfile, "* Create tracks"; logfile.flush()
    create_tracks(ex,outall,sample_names,assembly)
    return 0

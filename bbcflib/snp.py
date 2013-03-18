# Built-in modules #
import re, os, sys, tarfile
from operator import itemgetter
from itertools import product

# Internal modules #
from bbcflib.common import unique_filename_in, set_file_descr, sam_faidx
from bbcflib import mapseq
from bbcflib.btrack import FeatureStream, track, convert
from bbcflib.bFlatMajor.common import concat_fields
from bbcflib.bFlatMajor import stream as gm_stream

# Other modules #
from bein import program

_iupac = {'M':['A','a','C','c'],'Y':['T','t','C','c'],'R':['A','a','G','g'],
          'S':['G','g','C','c'],'W':['A','a','T','t'],'K':['T','t','G','g']}

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
    return {"arguments": ["samtools","pileup","-B","-cvsf",refGenome,"-N",str(ploidy),bamfile],
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
    def _parse_info8(readbase,cons):
        var1_fwd,var1_rev,var2_fwd,var2_rev = _iupac[cons]
        nvar1_fwd = readbase.count(var1_fwd)
        nvar1_rev = readbase.count(var1_rev)
        nvar2_fwd = readbase.count(var2_fwd)
        nvar2_rev = readbase.count(var2_rev)
        return var1_fwd,var2_fwd, nvar1_fwd,nvar1_rev,nvar2_fwd,nvar2_rev

    ref = info[2].upper()  # reference base
    cons = info[3].upper() # consensus base
    nreads = int(info[7])  # total coverage at this position
    #snp_qual = 10**(-float(info[5])/10)
    if ref == "*": # indel
        nindel = int(info[10])
        nref = int(info[11])
        if nreads-nref < mincov:
            return ref
        consensus = "%s/%s (%.2f%% of %s)" % (info[8],info[9],nindel*100./nreads,nreads)
    elif re.search(r'[ACGTN]',cons): # if bases are encoded normally (~100% replacement)
        consensus = "%s (%s)" % (cons,nreads)
    elif re.search(r'[MSYWRK]',cons):
        ploidy = _ploidy(assembly)
        nref = info[8].count(".") + info[8].count(",") # number of reads supporting wild type (.fwd and ,rev)
        if nreads-nref < mincov:
            return ref
        var1,var2, nvar1_fwd,nvar1_rev, nvar2_fwd,nvar2_rev = _parse_info8(info[8],cons)
        nvar1pc = (100/float(nreads)) * (nvar1_fwd + nvar1_rev) # % of consensus 1
        nvar2pc = (100/float(nreads)) * (nvar2_fwd + nvar2_rev) # % of consensus 2
        check1 = nvar1_fwd >= 5 and nvar1_rev >= 5 and nvar1pc >= minsnp/ploidy
        check2 = nvar2_fwd >= 5 and nvar2_rev >= 5 and nvar2pc >= minsnp/ploidy
        check0 = nvar2_fwd >= 3 and nvar2_rev >= 3 and nvar2pc >= 0.5*minsnp/ploidy \
             and nvar1_fwd >= 3 and nvar1_rev >= 3 and nvar1pc >= 0.5*minsnp/ploidy
        if ref.upper() == var1 and check2:   # if ref base == first variant
            consensus = "%s (%.2f%% of %s)" % (var2,nvar2pc,nreads)
        elif ref.upper() == var2 and check1: # if ref base == second variant
            consensus = "%s (%.2f%% of %s)" % (var1,nvar1pc,nreads)
        elif check0:                         # if third allele
            consensus = "%s (%.2f%%),%s (%.4g%% of %s)" % (var1,nvar1pc,var2,nvar2pc,nreads)
        else: consensus = ref
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
    ex.add( tarname, description=set_file_descr(tarname,step="pileup",type="tar",view='admin') )
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
            try: coverage = bam_tracks[i].coverage((chrbam,pos,pos+1)).next()[-1]
            except StopIteration: coverage = 0
            current_snps[i] = ref if coverage else "0"
        if not all([s == ref for s in current_snps]):
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
                        if strand < 0:
                            newc[i*cnumb+j] = newc[i*cnumb+j][:2-shift]+vari +newc[i*cnumb+j][3-shift:]
                            assert ref_codon[2-shift] == refbase, "bug with shift within codon"
                        else:
                            newc[i*cnumb+j] = newc[i*cnumb+j][:shift]  +vari +newc[i*cnumb+j][shift+1:]
                            assert ref_codon[shift] == refbase, "bug with shift within codon"
                new_codon[k] = newc
        if new_codon is None: return
        if strand < 0:
            ref_codon = _revcomp(ref_codon)
            new_codon = [[_revcomp(s) for s in c] for c in new_codon]
        for chr,pos,refbase,variants,cds,strand,ref_codon,shift in _buffer:
            refc = [itemgetter(0,2)(_iupac[x]) if x in _iupac else (x,) for x in ref_codon]
            ref_codon = [''.join(x) for x in product(*refc)]
            newc = [[[itemgetter(0,2)(_iupac[x]) if x in _iupac else (x,) for x in variant]
                    for variant in sample] for sample in new_codon]
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
                shift = (pos - (es + phase)) % 3
                codon_start = pos - shift
            elif strand == -1:
                shift = (ee - phase - pos) % 3
                codon_start = pos + shift - 2
            else:
                print "***",es,ee,cds,strand,phase
                continue
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
    """Write SQL and BED tracks showing all SNPs found."""
    with open(outall,'rb') as f:
        infields = f.readline().lstrip('#').split()
    intrack = track(outall, format='text', fields=infields, chrmeta=assembly.chrmeta)
    def _add_start(stream):
        for x in stream:
            end_idx = stream.fields.index('end')
            end = int(x[end_idx])
            yield x[:end_idx]+(end-1,)+x[end_idx:]
    for sample_name in sample_names:
        instream = intrack.read(fields=['chromosome','position','reference',sample_name])
        try: instream.next() # skip header
        except StopIteration: continue
        instream = concat_fields(instream,['reference',sample_name],'snp',separator=' -> ')
        instream.fields = ['chr','end','snp']
        instream = FeatureStream(_add_start(instream), fields=['chr','start','end','snp'])
        # SQL track
        out = unique_filename_in()
        outtrack = track(out+'.sql', format='sql', fields=['chr','start','end','snp'], chrmeta=assembly.chrmeta)
        outtrack.write(instream)
        description = set_file_descr("allSNP_track_"+sample_name+".sql" ,type='sql',step='tracks',gdv='1')
        ex.add(out+'.sql', description=description)
        # BED track
        t = track(out+'.bed.gz')
        t.make_header(name="allSNP_track_"+sample_name)
        convert(out+'.sql',out+'.bed.gz',mode='append')
        description = set_file_descr("allSNP_track_"+sample_name+".bed.gz",type='bed',step='tracks',ucsc='1')
        ex.add(out+'.bed.gz', description=description)


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
    for gid, files in job.files.iteritems():
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

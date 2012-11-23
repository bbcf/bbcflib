# Built-in modules #
import re, os, sys

# Internal modules #
from bbcflib.common import unique_filename_in, set_file_descr
from bbcflib import mapseq
from bbcflib.btrack import FeatureStream, track, convert
from bbcflib.bFlatMajor.common import concat_fields
from bbcflib.bFlatMajor import stream as gm_stream

# Other modules #
from bein import program

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
    else:
        ploidy = 2 # eucaryote
    return ploidy

@program
def sam_pileup(assembly,bamfile,refGenome,via='lsf'):
    """Binds 'samtools pileup'.

    :param assembly: Genrep.Assembly object.
    :param bamfile: path to the BAM file.
    :param refGenome: path to the reference genome fasta file.
    """
    ploidy = _ploidy(assembly)
    return {"arguments": ["samtools","pileup","-B","-cvsf",refGenome,"-N",str(ploidy),bamfile],
            "return_value": None}

def all_snps(chrom,outall,dictPileup,sample_names,allSNPpos,assembly):
    """For a given chromosome, returns a summary file containing all SNPs identified
    in at least one of the samples.
    Each row contains: chromosome id, SNP position, reference base, SNP base (with proportions)

    :param dictPileup: (dict) dictionary of the form {filename: (bein.Future, sample_name, bam_filename)}.
    :param sample_names: (list of str) list of sample names.
    :param allSNPpos: dict fo the type {3021: 'A'} as returned by posAllUniqSNP(...)[0].
    :param chrom: (str) chromosome name.
    """
    # Note: sample_names is redundant with dictPileup, can find a way to get rid of it
    allSamples = {}

    def _parse_info8(readbase,cons):
        iupac = {'M':['A','a','C','c'],'Y':['T','t','C','c'],'R':['A','a','G','g'],
                 'S':['G','g','C','c'],'W':['A','a','T','t'],'K':['T','t','G','g']}
        indels = []
        irb = iter(readbase)
        rb = ''
        for x in irb:
            if x in ['+','-']:
                n_indel = int(irb.next())
                indels.append(''.join([irb.next() for _ in range(n_indel)]))
            elif x == '$': pass
            elif x == '^': irb.next()
            else: rb += x
        readbase = ''.join(rb)
        var1_fwd,var1_rev,var2_fwd,var2_rev = iupac[cons]
        nvar1_fwd = readbase.count(var1_fwd)
        nvar1_rev = readbase.count(var1_rev)
        nvar2_fwd = readbase.count(var2_fwd)
        nvar2_rev = readbase.count(var2_rev)
        return var1_fwd, var2_fwd, nvar1_fwd, nvar1_rev, nvar2_fwd, nvar2_rev, indels

    for pileup_filename,trio in dictPileup.iteritems():
        allpos = sorted(allSNPpos.keys(),reverse=True) # list of positions [int] with an SNP across all groups
        sname = trio[1]
        bamtrack = track(trio[2],format='bam')
        ploidy = _ploidy(assembly)
        pos = -1
        ref = None # In case sample is empty
        allSamples[sname] = {}
        with open(pileup_filename) as sample:
            for line in sample:
                info = line.split("\t")
                chrbam = info[0]
                ref = info[2].upper() # reference base
                cons = info[3].upper() # consensus base
                snp_qual = 10**(-float(info[5])/10)
                nreads = info[7] # coverage at this position
                # info = ['chrV', '91668', 'G', 'R', '4', '4', '60', '4', 'a,.^~.', 'LLLL', '~~~~\n']
                # info = [chr, pos, ref, consensus, cons_qual, snp_qual, max_map_qual, nreads, 'a,.^~.', 'LLLL', '~~~~\n']
                #          0    1    2       3          4         5           6           7       8         9       10
                # cons_qual : consensus quality is the Phred-scaled probability that the consensus is wrong.
                # snp_qual: SNP quality is the Phred-scaled probability that the consensus is identical to the reference.
                while int(info[1]) > pos and len(allpos) > 0: # while pos not found in the given sample
                    pos = allpos.pop()
                    try: coverage = bamtrack.coverage((chrbam,pos,pos+1)).next()[-1]
                    except StopIteration: coverage = 0
                    allSamples[sname][pos] = allSNPpos[pos] if coverage else "0"  # "0" if not covered, ref base otherwise
                if not(int(info[1]) == pos): continue
                # SNP found in allpos, treat:
                elif int(nreads) < 10:
                    star = "* " # add a star *A if the snp is supported by less than 10 reads
                else:
                    star = ""
                if re.search(r'[ACGTN]',cons): # if bases are encoded normally (~100% replacement)
                    star += cons
                    allSamples[sname][pos] = star
                elif re.search(r'[MSYWRK]',cons): # avoid cons='*/*' as it happened once
                    mincov = 40 # half of it for diploids. Used to be set to 80.
                    var1,var2, nvar1_fwd,nvar1_rev,nvar2_fwd,nvar2_rev, indels = _parse_info8(info[8],cons)
                    nvar1 = (100/float(nreads)) * (nvar1_fwd + nvar1_rev) # % of consensus 1
                    nvar2 = (100/float(nreads)) * (nvar2_fwd + nvar2_rev) # % of consensus 2
                    check1 = nvar1_fwd >= 5 and nvar1_rev >= 5 and nvar1 >= mincov/ploidy
                    check2 = nvar2_fwd >= 5 and nvar2_rev >= 5 and nvar2 >= mincov/ploidy
                    check0 = nvar2_fwd >= 3 and nvar2_rev >= 3 and nvar2 >= 0.5*mincov/ploidy \
                         and nvar1_fwd >= 3 and nvar1_rev >= 3 and nvar1 >= 0.5*mincov/ploidy
                    if ref.upper() == var1 and check2:   # if ref base == first variant
                        allSamples[sname][pos] = star + "%s (%.4g%%)" % (var2,nvar2)
                    elif ref.upper() == var2 and check1: # if ref base == second variant
                        allSamples[sname][pos] = star + "%s (%.4g%%)" % (var1,nvar1)
                    elif check0:                         # if third allele
                        allSamples[sname][pos] = star + "%s (%.4g%%),%s (%.4g%%)" % (var1,nvar1,var2,nvar2)
                    else:
                        allSamples[sname][pos] = ref
            while allpos:  # write '0' for all pos after the last one of this sample
                pos = allpos.pop()
                try: coverage = bamtrack.coverage((chrom,pos,pos+1)).next()[-1]
                except StopIteration: coverage = 0
                allSamples[sname][pos] = allSNPpos[pos] if coverage else "0"
        bamtrack.close()

    allpos = sorted(allSNPpos.keys(),reverse=False) # re-init
    allsamples = []
    for pos in allpos:
        nbNoSnp = 0 # Check if at least one sample still has the SNP (after filtering)
        for sname in sample_names:
            if allSamples[sname][pos] in ("0",allSNPpos[pos]): nbNoSnp += 1
        if nbNoSnp != len(allSamples):
            refbase = allSNPpos[pos]
            allsamples.append((chrom,pos-1,pos,refbase) + tuple([allSamples[s][pos] for s in sample_names]))
            # append: chr start end ref_base    sample1 ... sampleN
            # sampleX is one of '0', 'A', 'T (28%)', 'T (28%),G (13%)', with or without star

    snp_read = FeatureStream(allsamples, fields=['chr','start','end','name']+sample_names )
    annotation = assembly.gene_track(chrom)
    annotated_stream = gm_stream.getNearestFeature(snp_read, annotation,
                                                   thresholdPromot=3000, thresholdInter=3000, thresholdUTR=10)
    with open(outall,"a") as fout:
        fout.write('#'+'\t'.join(['chromosome','position','reference'] + sample_names
                             + ['gene','location_type','distance']) + '\n')
        for snp in annotated_stream:
            # snp: ('chrV',154529, 154530, 'T', 'A', '* A', 'YER002W|NOP16_YER001W|MNN1', 'Upstream_Included', '2271_1011')
            fout.write('\t'.join([str(x) for x in (snp[0],)+snp[2:]])+'\n')
            # chrV  1606    T   43.48% C / 56.52% T YEL077W-A|YEL077W-A_YEL077C|YEL077C 3UTR_Included   494_-2491
    return allsamples

def exon_snps(chrom,outexons,allsamples, sample_names, assembly, genomeRef={}):
    """Annotates SNPs described in `filedict` (a dictionary of the form {chromosome: filename}
    where `filename` is an output of parse_pileupFile).
    Adds columns 'gene', 'location_type' and 'distance' to the output of parse_pileupFile.
    Returns two files: the first contains all SNPs annotated with their position respective to genes in
    the specified assembly, and the second contains only SNPs found within CDS regions.

    :param sample_names: list of sample names
    :param assembly: genrep.Assembly object
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
            varbase = [r.strip('* ') for r in variants]
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
            result = [chr, pos+1, refbase] + list(variants) + [cds, strand] \
                     + [_translate[ref_codon]] \
                     + [','.join([_translate[s] for s in c]) for c in new_codon]
            outex.write("\t".join([str(r) for r in result])+"\n")

    outex = open(outexons,"a")
    outex.write('#'+'\t'.join(['chromosome','position','reference'] + sample_names + ['exon','strand','ref_aa'] \
                          + ['new_aa_'+s for s in sample_names])+'\n')

    snp_stream = FeatureStream(allsamples, fields=['chr','start','end','ref']+sample_names)
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
        instream.next() # skip header
        instream = concat_fields(instream,['reference',sample_name],'snp',separator=' -> ')
        instream.fields = ['chr','end','snp']
        instream = FeatureStream(_add_start(instream), fields=['chr','start','end','snp'])
        # SQL track
        out = unique_filename_in()
        outtrack = track(out+'.sql', format='sql', fields=['chr','start','end','snp'], chrmeta=assembly.chrmeta)
        outtrack.write(instream)
        description = set_file_descr("allSNP_track_"+sample_name+".sql" ,type='sql',step='SNPs',gdv='1')
        ex.add(out+'.sql', description=description)
        # BED track
        t = track(out+'.bed.gz')
        t.make_header(name="allSNP_track_"+sample_name)
        convert(out+'.sql',out+'.bed.gz',mode='append')
        description = set_file_descr("allSNP_track_"+sample_name+".bed.gz",type='bed',step='SNPs',ucsc='1')
        ex.add(out+'.bed.gz', description=description)

def posAllUniqSNP(dictPileup):
    """
    Retrieve the results from samtools pileup and return a dict of the type {3021: 'A'}.

    :param dictPileup: (dict) dictionary of the form {filename: (bein.Future,sample_name,bam_file)}
    """
    d={}
    for filename,trio in dictPileup.iteritems():
        trio[0].wait() #file p is created from samtools pileup
        with open(filename,'rb') as f:
            for l in f:
                data = l.split("\t")
                cpt = data[8].count(".") + data[8].count(",") # number of reads supporting wild type (.fwd and ,rev)
                if int(data[7])-cpt >= 5:
                    d[int(data[1])] = data[2].upper()
    return d


def snp_workflow(ex, job, bam_files, assembly, path_to_ref, via):
    """Main function of the workflow"""
    if path_to_ref is None:
        path_to_ref = os.path.join(assembly.genrep.root,'nr_assemblies/fasta',assembly.md5+'.tar.gz')
    assert os.path.exists(path_to_ref), "Reference sequence not found: %s." % path_to_ref
    genomeRef = assembly.untar_genome_fasta(path_to_ref, convert=True)
    pileup_dict = dict((chrom,{}) for chrom in genomeRef.keys()) # {chr: {}}
    sample_names = []
    bam = {}
    for gid, files in bam_files.iteritems():
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
            pileupFilename = unique_filename_in()
            future = sam_pileup.nonblocking(ex,assembly,bam,ref,via,stdout=pileupFilename )
            pileup_dict[chrom][pileupFilename] = (future,sample_name,bam)

    outall = unique_filename_in()
    outexons = unique_filename_in()
    for chrom, dictPileup in pileup_dict.iteritems():
        # Get the results from sam_pileup
        # Write the list of all snps of THIS chromosome, from ALL samples
        allSNPpos = posAllUniqSNP(dictPileup) # {3021:'A'}
        if len(allSNPpos) == 0: continue
        # Store the info from the pileup file for this chromosome in a dictionary
        allsamples = all_snps(chrom,outall,assembly,sample_names,dictPileup,allSNPpos)
        # Add exon & codon information & write in file
        exon_snps(chrom,outexons,allsamples,sample_names,assembly,genomeRef)

    description = set_file_descr("allSNP.txt",step="SNPs",type="txt")
    ex.add(outall,description=description)
    description = set_file_descr("exonsSNP.txt",step="SNPs",type="txt")
    ex.add(outexons,description=description)
    # Create tracks for UCSC and GDV:
    create_tracks(ex,outall,sample_names,assembly)



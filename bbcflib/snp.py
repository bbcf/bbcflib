# Built-in modules #
import re, tarfile, os, sys

# Internal modules #
from bbcflib.common import unique_filename_in, set_file_descr
from bbcflib.btrack import FeatureStream, track, convert
from bbcflib.bFlatMajor.common import select, concat_fields
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

def untar_genome_fasta(assembly, path_to_ref=None, convert=True):
    """Untar and concatenate reference sequence fasta files.
    Returns a dictionary {chr_name: file_name}

    :param assembly: the GenRep.Assembly instance of the species of interest.
    :param convert: (bool) True if chromosome names need conversion from id to chromosome name.
    """
    if convert:
        # Create a dict of the form {'2704_NC_001144.4': chrXII}
        chrlist = dict((str(k[0])+"_"+str(k[1])+"."+str(k[2]),v['name'])
                       for k,v in assembly.chromosomes.iteritems())
    else:
        chrlist = {}
    archive = tarfile.open(path_to_ref)
    genomeRef = {}
    for f in archive.getmembers():
        if f.isdir(): continue
        inf = archive.extractfile(f)
        header = inf.readline()
        headpatt = re.search(r'>(\S+)\s',header)
        if headpatt: chrom = headpatt.groups()[0]
        else: continue
        if convert and chrom in chrlist:
            header = re.sub(chrom,chrlist[chrom],header)
            chrom = chrlist[chrom]
        genomeRef[chrom] = unique_filename_in()
        with open(genomeRef[chrom],"wb") as outf:
            outf.write(header)
            [outf.write(l) for l in inf]
        inf.close()
    archive.close()
    return genomeRef

@program
def sam_pileup(assembly,bamfile,refGenome,via='lsf',minSNP=10,minCoverage=80):
    """Binds 'samtools pileup'.

    :param assembly: Genrep.Assembly object.
    :param bamfile: path to the BAM file.
    :param refGenome: path to the reference genome fasta file.
    """
    if str(assembly.name) in ['MLeprae_TN','MSmeg_MC2_155','MTb_H37Rv','NA1000','TB40-BAC4']:
        ploidy=1
        minSNP=5 #10
        minCoverage=80
    else:
        ploidy=2
        minSNP=10 #20
        minCoverage=40
    return {"arguments": ["samtools","pileup","-B","-cvsf",refGenome,"-N",str(ploidy),bamfile],
            "return_value": [minCoverage,minSNP]}

def write_pileupFile(dictPileup,sample_names,allSNPpos,chrom,minCoverage=80,minSNP=10):
    """For a given chromosome, returns a summary file containing all SNPs identified
    in at least one of the samples.
    Each row contains: chromosome id, SNP position, reference base, SNP base (with proportions)

    :param dictPileup: (dict) dictionary of the form {filename: (bein.Future, sample_name)}.
    :param sample_names: (list of str) list of sample names.
    :param allSNPpos: dict fo the type {3021: 'A'} as returned by posAllUniqSNP(...)[0].
    :param chrom: (str) chromosome name.
    :param minCoverage: (int) the minimal percentage of reads supporting a SNP to reject a sequencing error.
    :param minSNP: (int) the minimal coverage of the SNP position to accept the SNP.
    """
    # Note: sample_names is redundant with dictPileup, can find a way to get rid of it
    formattedPileupFilename = unique_filename_in()
    allSamples = {}
    iupac = {'M':['A','a','C','c'],'Y':['T','t','C','c'],'R':['A','a','G','g'],
             'S':['G','g','C','c'],'W':['A','a','T','t'],'K':['T','t','G','g']}

    for pileup_filename,pair in dictPileup.iteritems():
        allpos = sorted(allSNPpos.keys(),reverse=True) # list of positions [int] with an SNP across all groups
        sname = pair[1]
        with open(pileup_filename) as sample:
            pos = -1
            allSamples[sname] = {}
            for line in sample:
                info = line.split("\t")
                ref = info[2].upper() # reference base
                cons = info[3].upper() # consensus base
                if cons == 'N': continue
                nreads = info[7]; # coverage at this position
                # info = ['chrV', '91668', 'G', 'R', '4', '4', '60', '4', 'a,.^~.', 'LLLL', '~~~~\n']
                # info = [chr, pos, ref, consensus, cons_qual, snp_qual, max_map_qual, nreads, 'a,.^~.', 'LLLL', '~~~~\n']
                #          0    1    2       3          4         5           6           7       8         9       10
                while int(info[1]) > pos and len(allpos) > 0: # write '0' for all pos until the one on this line of the sample
                    pos = allpos.pop()
                    allSamples[sname][pos] = "0"
                if not(int(info[1]) == pos): continue
                # SNP found in allpos, treat
                if int(nreads) == 0:
                    allSamples[sname][pos] = "0"
                elif int(nreads) < minSNP:
                    star = "* " # add a star *A if the snp is supported by less than minSNP reads
                else:
                    star = ""
                if re.search(r'[ACGT]',cons): # if bases are encoded normally (~100% replacement?)
                    star += cons
                    allSamples[sname][pos] = star
                else:
                    snp = info[8].count(iupac[cons][0]) + info[8].count(iupac[cons][1]) # nreads with one variant
                    snp2 = info[8].count(iupac[cons][2]) + info[8].count(iupac[cons][3]) # nreads with the other
                    # If the number of useful reads > minCoverage % of total reads
                    if (snp+snp2)*100 > minCoverage*int(nreads):
                        cov = 100 / float(nreads)
                        if ref == iupac[cons][0]:   # if ref base == first variant
                            allSamples[sname][pos] = star+"%s (%.4g%%)" % (iupac[cons][2], snp2*cov)
                        elif ref == iupac[cons][2]: # if ref base == second variant
                            allSamples[sname][pos] = star+"%s (%.4g%%)" % (iupac[cons][0], snp*cov)
                        else:                       # if third allele
                            allSamples[sname][pos] = star+"%s (%.4g%%),%s (%.4g%%)" \
                                    % (iupac[cons][0], snp*cov, iupac[cons][2], snp2*cov)
            while allpos:  # write '0' for all pos after the last one of this sample
                pos = allpos.pop()
                allSamples[sname][pos] = "0"

    #firstSample = allSamples.itervalues().next()
    allpos = sorted(allSNPpos.keys(),reverse=False) # list of positions [int] with an SNP across all groups
    with open(formattedPileupFilename,'wb') as outfile:
        for pos in allpos:
            nbNoSnp = 0 # Check if at least one sample had the SNP (??)
            for sname in sample_names:
                nbNoSnp += allSamples[sname][pos].count("0")
            if nbNoSnp != len(allSamples):
                refbase = allSNPpos[pos]
                outfile.write("\t".join([chrom,str(pos),refbase] + [allSamples[s][pos] for s in sample_names])+"\n")
                # Write: chr    pos    ref_base    sample1 ... sampleN
                # sampleX is one of '0', 'A', 'T (28%)', 'T (28%),G (13%)', with or without star

    return formattedPileupFilename

def annotate_snps(filedict, sample_names, assembly, genomeRef=None ):
    """Annotates SNPs described in `filedict` (a dictionary of the form {chromosome: filename}
    where `filename` is an output of parse_pileupFile).
    Adds columns 'gene', 'location_type' and 'distance' to the output of parse_pileupFile.
    Returns two files: the first contains all SNPs annotated with their position respective to genes in
    the specified assembly, and the second contains only SNPs found within CDS regions.

    :param filedict: {chr: formattedPileupFilename}
    :param sample_names: list of sample names
    :param assembly: genrep.Assembly object
    :param genomeRef: dict of the form {'chr1': filename}, where filename is the name of a fasta file
        containing the reference sequence for the chromosome.
    """
    def _process_annot(stream, fname):
        """Write (append) features from the *stream* in a file *fname*."""
        with open(fname,'a') as fout:
            for snp in stream:
                fout.write("\t".join([str(x) for x in (snp[2],snp[1])+snp[3:]])+"\n")
                # chrV  1606    T   43.48% C / 56.52% T YEL077W-A|YEL077W-A_YEL077C|YEL077C 3UTR_Included   494_-2491
                if "Included" in snp[-2].split("_"):
                    yield (snp[2],int(snp[0]),int(snp[1]))+snp[3:-3]
                    # ('chrV', 1606, 1607, 'T', 'C (43%)')

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
            variants = []
            if new_codon is None: new_codon = [[ref_codon] for _ in range(len(varbase))]
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
        for chr,pos,refbase,variants,cds,strand,refcodon,shift in _buffer:
            result = [chr, pos+1, refbase] + list(variants) + [cds, strand] \
                     + [_translate[ref_codon]] \
                     + [','.join([_translate[s] for s in c]) for c in new_codon]
            outex.write("\t".join([str(r) for r in result])+"\n")

    if genomeRef is None: genomeRef = {}
    output = unique_filename_in()
    outall = output+"_all_snps.txt"
    outexons = output+"_exon_snps.txt"
    outex = open(outexons,"w")
    with open(outall,"w") as fout:
        fout.write("\t".join(['chromosome','position','reference'] + sample_names
                             + ['gene','location_type','distance']) + "\n")
    outex.write("\t".join(['chromosome','position','reference'] + sample_names + ['exon','strand','ref_aa'] \
                          + ['new_aa_'+s for s in sample_names])+"\n")

    for chrom, filename in filedict.iteritems():
        # For each chromosome, read the result of parse_pileupFile and make a track
        snp_file = track( filename, format='text',
                          fields=['chr','end','name']+sample_names, chrmeta=assembly.chrmeta )
        # Add a 'start' using end-1 and make an iterator from the track
        snp_read = FeatureStream( ((y[0],y[1]-1)+y[1:] for y in snp_file.read(chrom)),
                                        fields=['chr','start','end','name']+sample_names)
        annotation = assembly.gene_track(chrom)
        annotated_stream = gm_stream.getNearestFeature(snp_read, annotation,
                                                       thresholdPromot=3000, thresholdInter=3000, thresholdUTR=10)
        # Write a line in outall at each iteration; yield if the snp is Included in a gene only.
        inclstream = concat_fields(FeatureStream(_process_annot(annotated_stream, outall),
                                                             fields=snp_read.fields),
                                         infields=['name']+sample_names, as_tuple=True)
        annotstream = concat_fields(assembly.annot_track('CDS',chrom),
                                    infields=['name','strand','frame'], as_tuple=True)
        annotstream = FeatureStream((x[:3]+(x[1:3]+x[3],) for x in annotstream),fields=annotstream.fields)
        _buffer = {1:[], -1:[]}
        last_start = {1:-1, -1:-1}
        for x in gm_stream.combine([inclstream, annotstream], gm_stream.intersection):
            # x = (1606,1607,'chrV', ('T','C (43%), 1612,1724,'YEL077C|YEL077C',-1,0, 1712,1723,'YEL077W-A|YEL077W-A',1,0))
            nsamples = len(sample_names)
            pos = x[0]; chr = x[2]; rest = x[3]
            refbase = rest[0]
            annot = [rest[5*i+nsamples+1 : 5*i+5+nsamples+1]
                     for i in range(len(rest[nsamples+1:])/5)] # list of [start,end,cds,strand,phase]
            for es,ee,cds,strand,phase in annot:
                if strand == 1:
                    shift = (pos - (es + phase)) % 3
                    codon_start = pos - shift
                elif strand == -1:
                    shift = (ee - phase - pos) % 3
                    codon_start = pos + shift - 2
                ref_codon = assembly.fasta_from_regions({chr: [[codon_start,codon_start+3]]}, out={},
                                                        path_to_ref=genomeRef.get(chr))[0][chr][0]
                info = [chr, pos, refbase, list(rest[1:1+nsamples]), cds, strand, ref_codon, shift]
                # Either the codon is the same as the previous one on this strand, or it will never be.
                # Only if one codon is passed, can write its snps to a file.
                if codon_start == last_start[strand]:
                    _buffer[strand].append(info)
                else:
                    _write_buffer(_buffer[strand],outex)
                    _buffer[strand] = [info]
                    last_start[strand] = codon_start
        for strand in [1,-1]:
            _write_buffer(_buffer[strand],outex)
    outex.close()
    return (outall, outexons)


def create_tracks(ex, filename, sample_names, assembly):
    with open(filename,'rb') as f:
        infields = f.readline().split()
    intrack = track(filename, format='text', fields=infields, chrmeta=assembly.chrmeta)
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
        # BED track
        out = unique_filename_in()+'.bed'
        outtrack = track(out, format='bed', fields=['chr','start','end','snp'], chrmeta=assembly.chrmeta)
        outtrack.write(instream)
        description = set_file_descr("allSNP_track_"+sample_name+".bed" ,type='bed',step='SNPs')
        ex.add(out, description=description)
        # SQL track
        convert(out,out+'.sql')
        description = set_file_descr("allSNP_track_"+sample_name+".sql" ,type='sql',step='SNPs')
        ex.add(out+'.sql', description=description)


def posAllUniqSNP(PileupDict):
    """
    Retrieve the results from samtools pileup and store them into a couple of the type
    ({3021: 'A'}, (minCoverage, minSNP))

    :param PileupDict: (dict) dictionary of the form {filename: bein.Future}
    """
    d={}
    for filename,pair in PileupDict.iteritems():
        parameters = pair[0].wait() #file p is created and samtools pileup returns its own parameters
        minCoverage = parameters[0]
        minSNP = parameters[1]
        with open(filename,'rb') as f:
            for l in f:
                data = l.split("\t")
                cpt = data[8].count(".") + data[8].count(",")
                if cpt*100 < int(data[7]) * int(minCoverage) and int(data[7]) >= 8:
                    d[int(data[1])] = data[2].upper()
    return (d,parameters)


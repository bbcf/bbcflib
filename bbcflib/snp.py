# Built-in modules #
import re, tarfile, os, sys

# Internal modules #
from bbcflib.common import unique_filename_in
from bbcflib import btrack as track
from bbcflib.bFlatMajor import stream as gm_stream

# Other modules #
from bein import program

translate = { "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L/START",
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

def write_pileupFile(samples,allSNPpos,chrom,minCoverage=80,minSNP=10):
    """For a given chromosome, returns a summary file containing all SNPs identified 
    in at least one of the samples from *samples*.
    Each row contains: chromosome id, SNP position, reference base, SNP base (with proportions)

    :param ex: a bein.Execution instance.
    :param samples: (dict) dictionary of the form {filename: sample_name}
    :param allSNPpos: dict fo the type {3021: 'A'} as returned by posAllUniqSNP(...)[0].
    :param chrom: (str) chromosome name.
    :param minCoverage: (int) the minimal percentage of reads supporting a SNP to reject a sequencing error.
    :param minSNP: (int) the minimal coverage of the SNP position to accept the SNP.
    """
    formattedPileupFilename = unique_filename_in()
    allSamples = {}
    iupac = {'M':['A','a','C','c'],'Y':['T','t','C','c'],'R':['A','a','G','g'],
             'S':['G','g','C','c'],'W':['A','a','T','t'],'K':['T','t','G','g']}

    for pileup_filename,sname in samples.iteritems():
        allpos = sorted(allSNPpos.keys(),reverse=True) # list of positions [int] with an SNP across all groups
        with open(pileup_filename) as sample:
            pos = -1
            allSamples[sname] = {}
            for line in sample:
                info = line.split("\t")
                ref = info[2].upper() # reference base
                cons = info[3].upper() # consensus base
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
            for sname in allSamples: 
                nbNoSnp += allSamples[sname][pos].count("0")
            if nbNoSnp != len(allSamples):
                refbase = allSNPpos[pos]
                outfile.write("\t".join([chrom,str(pos),refbase] + [allSamples[s][pos] for s in allSamples])+"\n")
                # Write: chr    pos    ref_base    sample1 ... sampleN
                # sampleX is one of '0', 'A', 'T (28%)', 'T (28%),G (13%)', with or without star

    return formattedPileupFilename

def annotate_snps( filedict, sample_names, assembly ):
    """Annotates SNPs described in `filedict` (a dictionary of the form {chromosome: filename}
    where `filename` is an output of parse_pileupFile).
    Adds columns 'gene', 'location_type' and 'distance' to the output of parse_pileupFile.
    Returns two files: the first contains all SNPs annotated with their position respective to genes in
    the specified assembly, and the second contains only SNPs found within CDS regions.

    :param filedict: {chr: formattedPileupFilename}
    :sample_names: list of sample names
    :assembly: genrep.Assembly object
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

    trans_in_gene = assembly.get_trans_in_gene() # {gene_id: [transcript_ids]}
    exons_in_trans = assembly.get_exons_in_trans() # {transcript_id: [exon_ids]}
    transcript_mapping = assembly.get_transcript_mapping() # {transcript_id: (gene_id,start,end,length,strand,chromosome)}
    exon_mapping = assembly.get_exon_mapping() # {exon_id: ([transcript_ids],gene_id,start,end,strand,chromosome)}

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
        snp_file = track.track( filename, format='text',
                                fields=['chr','end','name']+sample_names, chrmeta=assembly.chrmeta )
        # Add a 'start' using end-1 and make an iterator from the track
        snp_read = track.FeatureStream( ((y[0],y[1]-1)+y[1:] for y in snp_file.read(chrom)),
                                        fields=['chr','start','end','name']+sample_names)
        annotation = assembly.gene_track(chrom) # gene annotation from genrep, FeatureStream instance
        # Find the nearest gene from each snp
        annotated_stream = gm_stream.getNearestFeature(snp_read, annotation,
                                thresholdPromot=3000, thresholdInter=3000, thresholdUTR=10)
        # Put together in a tuple the reference base and the base of every sample
        # Write a line in outall at each iteration; yield if the snp is in a CDS only
        inclstream = track.concat_fields(track.FeatureStream(_process_annot(annotated_stream, outall),
                                                                 fields=snp_read.fields),
                                         infields=['name']+sample_names, as_tuple=True)
        # import existing CDS annotation from genrep as a track, and join some fields
        annotstream = track.concat_fields(assembly.annot_track('CDS',chrom),
                                          infields=['name','strand','frame'], as_tuple=True)

        buffer = {1:{}, -1:{}}
        for x in gm_stream.combine([inclstream, annotstream], gm_stream.intersection):
            # x = (1606, 1607, 'chrV', ('T', '43.48% C / 56.52% T', 'YEL077C|YEL077C', -1, 0, 'YEL077W-A|YEL077W-A', 1, 0))
            nsamples = len(sample_names)
            pos = x[0]; chr = x[2]; rest = x[3]
            refbase = rest[0]
            varbase = [r.strip('* ') for r in rest[1:1+nsamples]]
            variants = []
            for variant in varbase:
                if variant == '0':
                    variants.append(refbase)
                elif len(variant.split(',')) > 1: # 'C (43%),G (12%)'
                    variants.append(variant.split()[0])
                else:                             # 'C (43%)' or 'C'
                    variants.append(variant.split()[0])
            annot = [rest[3*i+nsamples+1 : 3*i+3+nsamples+1] 
                     for i in range(len(rest[nsamples+1:])/3)] # list of [cds,strand,phase]
            for cds,strand,phase in annot:
                exon_id, gene_id, gene_name = cds.split('|')
                # If no exon_id, let's find one (sometimes)
                transcripts = trans_in_gene[gene_id]
                exons = None
                for t in transcripts: 
                    ts, te = transcript_mapping[t][1:3]
                    if ts <= pos and pos <= te:
                        exons = exons_in_trans[t]
                        break # found one; the phase should be the same for the other transcripts
                if exons:
                    for e in exons:
                        es, ee = exon_mapping[e][2:4]
                        if es <= pos and pos <= ee:
                            exon_id = e
                            break
                else: continue
                if exon_id == '': exon_id = '-'
                if strand == 1:
                    shift = (pos - (es + phase)) % 3
                    codon_start = pos - shift
                elif strand == -1:
                    shift = (ee - phase - pos) % 3
                    codon_start = pos + shift - 2
                ref_codon = assembly.fasta_from_regions({chr: [[codon_start,codon_start+3]]}, out={})[0][chr][0]
                new_codon = [list(ref_codon) for v in variants]
                for k in range(nsamples):
                    if strand == 1:
                        new_codon[k][shift] = variants[k]
                        assert ref_codon[shift] == refbase, "bug with shift within codon"
                    elif strand == -1:
                        new_codon[k][2-shift] = variants[k]
                        assert ref_codon[2-shift] == refbase, "bug with shift within codon"
                    new_codon[k] = "".join(new_codon[k])
                result = [x[2], str(x[0]), refbase] + list(rest[1:1+nsamples]) + [exon_id+'|'+gene_id+'|'+gene_name, strand] \
                          + [translate[ref_codon]] + [translate[n] for n in new_codon]

                #print strand, "buffer in", buffer[strand]
                if codon_start in buffer[strand]:
                    #print "add to existing buffer:", result[:4], strand, pos, codon_start
                    buffer[strand][codon_start].append(result)
                else:
                    for c,x in buffer[strand].iteritems():
                        for res in x:
                            outex.write("\t".join([str(r) for r in res])+"\n")
                    buffer[strand].clear()
                    buffer[strand][codon_start] = [result]
                    #print "add new to buffer:", result[:4], strand, pos, codon_start
                #print strand,"buffer out", buffer[strand],"\n"
                
                #print result,'\n'
        for strand in [1,-1]:
            for start, x in buffer[strand].iteritems():
                for res in x:
                    outex.write("\t".join([str(r) for r in res])+"\n")
        print "\n\n"

    outex.close()
    return (outall, outexons)


def posAllUniqSNP(PileupDict):
    """
    Retrieve the results from samtools pileup and store them into a couple of the type
    ({3021: 'A'}, (minCoverage, minSNP))

    :param PileupDict: (dict) dictionary of the form {filename: bein.Future}
    """
    d={}
    for filename,future in PileupDict.iteritems():
        parameters = future.wait() #file p is created and samtools pileup returns its own parameters
        minCoverage = parameters[0]
        minSNP = parameters[1]
        with open(filename,'rb') as f:
            for l in f:
                data = l.split("\t")
                cpt = data[8].count(".") + data[8].count(",")
                if cpt*100 < int(data[7]) * int(minCoverage) and int(data[7]) >= 8:
                    d[int(data[1])] = data[2].upper()
    return (d,parameters)


# Built-in modules #
import re, tarfile

# Internal modules #
from bbcflib.common import unique_filename_in
from bbcflib import btrack as track
from bbcflib.bFlatMajor import stream as gm_stream

# Other modules #
from bein import program


def untar_genome_fasta(assembly, convert=True):
    """Untar and concatenate reference sequence fasta files.

    :param assembly: the GenRep.Assembly instance of the species of interest.
    :param convert: (bool) True if chromosome names need conversion RefSeq -> Ensembl, 
        False otherwise.
    """
    if convert:
        chrlist = dict((str(k[0])+"_"+str(k[1])+"."+str(k[2]),v['name'])
                       for k,v in assembly.chromosomes.iteritems())
    else:
        chrlist = {}
    archive = tarfile.open(assembly.fasta_path())
    genomeRef = {}
    for f in archive.getmembers():
        if f.isdir(): continue
        inf = archive.extractfile(f)
        header = inf.readline()
        headpatt = re.search(r'>(\S+)\s',header)
        if headpatt: chrom = headpatt.groups()[0]
        else: continue
        if chrom in chrlist:
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
def sam_pileup(assembly,bamfile,refGenome,via='lsf'):
    """Launches 'samtools pileup' on command-line.
    
    :param assembly: Genrep.Assembly object for the species of interest.
    :param bamfile: path to the BAM file to run the samtool pileup command on.
    :param refGenome: path to the species' reference genome (fasta file).
    """
    if str(assembly.name) in ['MLeprae_TN','MSmeg_MC2_155','MTb_H37Rv','NA1000','TB40-BAC4']:
        ploidy=1
        minSNP=10
        minCoverage=80
    else:
        ploidy=2
        minSNP=20
        minCoverage=40
    return {"arguments": ["samtools","pileup","-B","-cvsf",refGenome,"-N",str(ploidy),bamfile],
            "return_value": [minCoverage,minSNP]}

def parse_pileupFile(dictPileupFile,allSNPpos,chrom,minCoverage=80,minSNP=10):
    """ return filename of summary file containing all SNPs significantly found in samples provide in dictPileupFile.
    Each raw contains chromosome id, SNP position, reference base, SNP base (with quantification)

    :param ex: a bein.Execution instance.
    :param dictPileupFile: (dict) dictionary of the form {[]}
    :param allSNPPos: (str) name of a file as returned by posAllUniqSNP (it contains position of all SNPs found in all samples)  
    :param chrom: (str) chromosome name.
    :param minCoverage: (int) the minimal percentage of reads with SNP to considere SNP as true and not as sequencing error
    :param minSNP: (int) the minimal coverage of the SNP position to considere SNP
 
    """
    formatedPileupFilename = unique_filename_in()
    allSample={}
    iupac={'M':['A','a','C','c'],'Y':['T','t','C','c'],'R':['A','a','G','g'],
           'S':['G','g','C','c'],'W':['A','a','T','t'],'K':['T','t','G','g']}

    for p,sname in dictPileupFile.iteritems():
        cpt=0
        allpos = sorted(allSNPpos.keys(),reverse=True)
        position = -1
        allSample[sname]={}
        with open(p) as sample:
            for line in sample:
                info=line.split("\t")
                while int(info[1])>position:
                    if not(allpos): break
                    position = allpos.pop()
                    allSample[sname][position]="-"
                if not(int(info[1]) == position): continue
                if int(info[7])<minSNP:
                    string="* "
                else:
                    string=""
                if re.search(r'[ACGT]',info[3]):
                    string+=info[3]
                    allSample[sname][position]=string
                else:
                    snp=0
                    snp2=0
                    snp=info[8].count(iupac[info[3]][0])+info[8].count(iupac[info[3]][1])
                    snp2=info[8].count(iupac[info[3]][2])+info[8].count(iupac[info[3]][3])
                    if (snp+snp2)*100 > minCoverage*int(info[7]):
                        cov = 100/float(info[7])
                        if info[2] == iupac[info[3]][0]:
                            allSample[sname][position]=string+"%.4g%% %s / %.4g%% %s" \
                                    %(snp2*cov,iupac[info[3]][2],100-snp2*cov,iupac[info[3]][0])
                        elif info[2] == iupac[info[3]][2]:
                            allSample[sname][position]=string+"%.4g%% %s / %.4g%% %s" \
                                    %(snp*cov, iupac[info[3]][0],100-snp*cov, iupac[info[3]][2])
                        else:
                            allSample[sname][position]=string+"%.4g%% %s / %.4g%% %s" \
                                    %(snp*cov, iupac[info[3]][0],    snp2*cov,iupac[info[3]][2])
            while allpos:
                position = allpos.pop()
                allSample[sname][position]="-"

    firstSample = allSample.values()[0]
    with open(formatedPileupFilename,'w') as outfile:
        for p in sorted(firstSample):
            nbNoSnp=0
            for s in allSample:
                nbNoSnp+=allSample[s][p].count("-")
            if nbNoSnp!=len(allSample.keys()):
                outfile.write("\t".join([chrom,str(p),allSNPpos[p]]+[allSample[s][p] for s in allSample])+"\n")

    return formatedPileupFilename

def annotate_snps( filedict, sample_names, assembly ):
    def _process_annot( stream, fname ):
        with open(fname,'a') as fout:
            for snp in stream:
                fout.write("\t".join([str(x) for x in (snp[2],snp[1])+snp[3:]])+"\n")
                if "Included" in snp[-2].split("_"):
                    yield (snp[2],int(snp[0]),int(snp[1]))+snp[3:-3]

    output = unique_filename_in()
    outall = output+"_all_snps.txt"
    outexons = output+"_exon_snps.txt"
    outex = open(outexons,"w")
    with open(outall,"w") as fout:
        newcols = sample_names+['gene','location_type','distance']
        fout.write("chromosome\tposition\treference\t"+"\t".join(newcols)+"\n")
    for chrom, filename in filedict.iteritems():
        snp_file = track.track( filename, format='text',
                                fields=['chr','end','name']+sample_names, chrmeta=assembly.chrmeta )
        snp_read = track.FeatureStream( ((y[0],y[1]-1)+y[1:] for y in snp_file.read(chrom)),
                                        fields=['chr','start','end']+snp_file.fields[2:])
        annotations = assembly.gene_track(chrom)
        fstream = gm_stream.getNearestFeature(snp_read, annotations,3000,3000)
        inclstream = track.concat_fields(track.FeatureStream(_process_annot(fstream, outall),
                                                             fields=snp_read.fields),
                                         infields=['name']+sample_names, as_tuple=True)
        annotstream = track.concat_fields(assembly.annot_track('CDS',chrom),
                                          infields=['name','strand','frame'], as_tuple=True)
        for x in gm_stream.combine([inclstream, annotstream], gm_stream.intersection):
            outex.write("\t".join([str(y) for y in x])+"\n")
    outexons.close()
    return (outall, outexons)


def synonymous(job,allSnp):
    """Writes the first line of the file *allSnp* to the file *allCodon*.

    :param job: a Frontend.Job object.
    :param allSnp: path to the file summarizing the localization of SNPs.
    """
    allCodon = unique_filename_in()
#    infile = track.rtack(allSnp)
#    outtrack = track.track(allCodon)
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

    return allCodon


def posAllUniqSNP(PileupFile,minCoverage=80):
    """ 
    :param PileupFile: (dict) dictionary of the form {filename: [?, bein.Future]} 
    :param minCoverage: (int) 
    """
    d={}
    for p,v in PileupFile.iteritems():
        parameters=v[1].wait()
        PileupFile[p]=v[0]
        with open(p) as f:
            for l in f:
                data=l.split("\t")
                cpt=data[8].count(".")+data[8].count(",")
                if cpt*100 < int(data[7])*int(minCoverage) and int(data[7])>9:
                    d[int(data[1])]=data[2]
    return (d,parameters)


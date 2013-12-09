"""
===================
Module: bbcflib.snp
===================

From a set of BAM files produced by an alignement on the genome, calls snps and annotates them
with respect to a set of coding genes on the same genome.
"""
# Built-in modules #
import os, sys, tarfile
from itertools import product

# Internal modules #
from bbcflib.common import unique_filename_in, set_file_descr, sam_faidx, timer, iupac, translate, unique
from bbcflib.mapseq import merge_bam, index_bam, replace_bam_header
from bbcflib.track import FeatureStream, track
from bbcflib.gfminer.common import concat_fields
from bbcflib.gfminer import stream as gm_stream
from bein import program

import pysam

_DEBUG_ = False


@program
def pileup(bams,path_to_ref,seq_depth=1000):
    """Use 'samtools mpileup' followed by 'bcftools view' and 'vcfutils'
    to call for SNPs. Ref.: http://samtools.sourceforge.net/mpileup.shtml
    !! No indel call, to speed it up, until their output is implemented."""
    if not isinstance(bams,(list,tuple)): bams = [bams]
    bams = ' '.join(bams)
    bcf = unique_filename_in()
    vcf = unique_filename_in()
    script_name = unique_filename_in()
    if _DEBUG_:
        pileup = unique_filename_in()
        vcf_raw = unique_filename_in()
        script = """
samtools mpileup -uDS -f %s %s > %s
bcftools view -bvcg %s > %s
bcftools view %s > %s
vcfutils.pl varFilter -D%d %s > %s;
""" %(path_to_ref, bams, pileup, pileup, bcf, bcf, vcf_raw, seq_depth, vcf_raw, vcf)
        print script
    else:
        with open(script_name,'w') as s:
            s.write("""
samtools mpileup -uDS -I -f %s %s | bcftools view -vcg - | vcfutils.pl varFilter -D%d > %s 
""" % (path_to_ref,bams,seq_depth,vcf))              # ^  the dash is important
    return {'arguments': ["sh",script_name], 'return_value': vcf}

def parse_vcf(vcf_line):
    """cf. ``http://samtools.sourceforge.net/mpileup.shtml``"""
    snp_info = {}
    if not vcf_line.strip(): return ''
    line = vcf_line.strip().split('\t')
    chrbam,pos,id,ref,alt,qual,filter,info,format = line[:9]
    #   qual  : -10 log_10 P(call in ALT is wrong) - the higher, the best.
    #   filter: semicolon-sep list of filters that failed to call snp.
    #   info  : added by bcftools about SNP calling.
    #   format: genotype fields - for parsing the next columns.
    general = (chrbam,int(pos),ref.upper(),alt,float(qual))
    for x in info.split(';'):
        if x != 'INDEL':
            a,b = x.split('=')
            snp_info[a] = b
        else:
            snp_info['INDEL'] = True
    sample_stats = (format,line[9]) # for us, always one sample at a time
    return (general,snp_info,sample_stats)

def filter_snp(general,snp_info,sample_stats,mincov,minsnp,ploidy):
    ref = general[2]
    dp4 = map(int, snp_info.get('DP4','0,0,0,0').split(',')) # fw.ref, rev.ref, fw.alt, rev.alt (filtered)
    total_reads = sum(dp4)
    if dp4[2]+dp4[3] < mincov:  # too few supporting alt
        return '-'  # '/'.join([ref]*ploidy)
    ratio = 100.0*(dp4[2]+dp4[3])/total_reads
    if ratio < minsnp/ploidy:
        return "0"
    alt = general[3]
    alts = [ref]+alt.split(',')
    sample = sample_stats[1]  # ex: "0/1:48,0,53:14:0:50"
    genotype = sample.split(':')[0]  # todo: better according to *format* = GT:PL:DP:SP:GQ
    sep = '/' if '/' in genotype else '|'  # phased if |, unphased if /
    genotype = [alts[int(i)] for i in genotype.split(sep)]
    #   Diploid:  GT: '0/1' -> 'T/A (50% of 8)'  [if ref='T' and alt='A']
    #   Haploid:  GT: '0/1' -> 'A (50% of 8)'
    if ploidy == 1:
        genotype = unique(genotype)
        if ref in genotype: genotype.remove(ref)
    genotype = sep.join(genotype)
    genotype = "%s (%.0f%% of %d)" % (genotype,ratio,total_reads)
    return genotype

@timer
def all_snps(ex,chrom,vcfs,bams, outall,assembly,sample_names, mincov,minsnp,
             logfile=sys.stdout,debugfile=sys.stderr):
    """For a given chromosome, returns a summary file containing all SNPs identified
    in at least one of the samples.
    Each row contains: chromosome id, SNP position, reference base, SNP base (with proportions)

    :param chrom: (str) chromosome name.
    :param outall: (str) name of the file that will contain the list of all SNPs.
    :param sample_names: (list of str) list of sample names.
    :param mincov: (int) Minimum number of reads supporting an SNP at a position for it to be considered. [5]
    :param minsnp: (int) Minimum percentage of reads supporting the SNP for it to be returned.
        N.B.: Effectively, half of it on each strand for diploids. [40]
    """
    ploidy = assembly.ploidy
    allsnps = []
    nsamples = len(sample_names)
    sorder = range(len(sample_names))
    current = [None]*nsamples
    bam_tracks = [track(v,format='bam') for k,v in sorted(bams.items())]
    vcf_handles = [open(v) for k,v in sorted(vcfs.items())]
    for i,vh in enumerate(vcf_handles):
        line = '#'
        while line and line[0]=='#':
            line = vh.readline()
        current[i] = parse_vcf(line)
    lastpos = 0
    pos = -1
    while any(current):
        current_pos = [int(x[0][1]) if x else sys.maxint for x in current]
        pos = min(current_pos)
        if pos == sys.maxint: break
        current_snp_idx = set(i for i in range(nsamples) if current_pos[i]==pos)
        current_snps = ["0"]*nsamples
        for i in current_snp_idx:
            general,snp_info,sample_stats = current[i]
            chrbam = general[0]
            ref = general[2]
            current_snps[i] = filter_snp(general,snp_info,sample_stats, mincov,minsnp,ploidy)
        # If there were still snp called at this position after filtering
        if any(current_snps[i] not in ("-","0") for i in current_snp_idx):
            for i in set(range(nsamples))-current_snp_idx:
                for coverage in bam_tracks[sorder[i]].coverage((chrbam,pos-1,pos)):
                    if coverage[-1] > 0:
                        current_snps[i] = '-'  # '/'.join([ref]*ploidy)
            if pos != lastpos: # indel can be located at the same position as an SNP
                allsnps.append((chrom,pos-1,pos,ref)+tuple(current_snps))
            lastpos = pos
        for i in current_snp_idx:
            current[i] = parse_vcf(vcf_handles[i].readline())
    for f in vcf_handles: f.close()
    for b in bam_tracks: b.close()

    logfile.write("  Annotate all SNPs\n"); logfile.flush()
    snp_read = FeatureStream(allsnps, fields=['chr','start','end','name']+sample_names)
    try:
        annotation = assembly.gene_track(chrom)
        annotated_stream = gm_stream.getNearestFeature(snp_read,annotation,
            thresholdPromot=3000, thresholdInter=3000, thresholdUTR=10)
    except:
        annotated_stream = snp_read
    logfile.write("  Write all SNPs\n"); logfile.flush()
    with open(outall,"a") as fout:
        for snp in annotated_stream:
            # snp: ('chrV',154529, 154530,'T','A (50% of 10)','A (80% of 10)',
            #       'YER002W|NOP16_YER001W|MNN1','Upstream_Included','2271_1011')
            fout.write('\t'.join(str(x) for x in (snp[0],)+snp[2:])+'\n')  # just remove end coord
    return allsnps


@timer
def exon_snps(chrom,outexons,allsnps,assembly,sample_names,genomeRef={},
              logfile=sys.stdout,debugfile=sys.stderr):
    """Annotates SNPs described in `filedict` (a dictionary of the form {chromosome: filename}
    where `filename` is an output of parse_pileupFile).
    Adds columns 'gene', 'location_type' and 'distance' to the output of parse_pileupFile.
    Returns two files: the first contains all SNPs annotated with their position respective to genes in
    the specified assembly, and the second contains only SNPs found within CDS regions.

    :param chrom: (str) chromosome name.
    :param outexons: (str) name of the file containing the list of SNPs on exons.
    :param allsnps: list of tuples (chr,start,end,ref,alt1..altN) as returned by all_snps().
        Ex: [('chr', 3684115, 3684116, 'G', 'G', 'G', 'G', 'T (56% of 167)', 'G'), ...]
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
        # One position at a time
        for chr,pos,refbase,variants,cds,strand,ref_codon,shift in _buffer:
            varbase = list(variants)  # Ex: ['G','G','G','T (56% of 167)','G'],  ['A/A','G/G (100% of 7)']
            if new_codon is None:
                new_codon = [[ref_codon] for _ in range(len(varbase))]
            variants = []  # [[variants sample1], [variants sample2], ...]
            # One sample at a time
            for variant in varbase:
                if variant in ['0','-']:
                    variants.append([refbase])
                else: # Ex: 'C/G (80% of 10)' : heterozygous simple (ref is C) or double snp (ref is not C)
                    v = variant.split()[0]  # Ex: C/G
                    v = unique(v.split('/'))  # Ex: G/G -> 'G'
                    if refbase in v: v.remove(refbase)  # Ex: C/G -> 'G' (if ref is C)
                    variants.append(v)
            # One sample at a time
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
            refc = [iupac.get(x,x) for x in ref_codon]
            ref_codon = [''.join(x) for x in product(*refc)]
            newc = [[[iupac.get(x,x) for x in variant] for variant in sample]
                    for sample in new_codon]
            new_codon = [[''.join(x) for codon in sample for x in product(*codon)] for sample in newc]
            if refbase == "*":
                result = [chr, pos+1, refbase] + list(variants) + [cds, strand] \
                         + [','.join([translate.get(refc,'?') for refc in ref_codon])] + ["indel"]
            else:
                result = [chr, pos+1, refbase] + list(variants) + [cds, strand] \
                         + [','.join([translate.get(refc,'?') for refc in ref_codon])] \
                         + [','.join([translate.get(s,'?') for s in newc]) for newc in new_codon]
            outex.write("\t".join([str(r) for r in result])+"\n")

    #############################################################
    snp_stream = FeatureStream(allsnps, fields=['chr','start','end','ref']+sample_names)
    inclstream = concat_fields(snp_stream, infields=snp_stream.fields[3:], as_tuple=True)
    snp_stream = FeatureStream(allsnps, fields=['chr','start','end','ref']+sample_names)
    inclstream = concat_fields(snp_stream, infields=snp_stream.fields[3:], as_tuple=True)
    try:
        annotstream = concat_fields(assembly.annot_track('CDS',chrom),
                                    infields=['name','strand','frame'], as_tuple=True)
        annotstream = FeatureStream((x[:3]+(x[1:3]+x[3],) for x in annotstream),fields=annotstream.fields)
    except:
        return False
    _buffer = {1:[], -1:[]}
    last_start = {1:-1, -1:-1}
    logfile.write("  Intersection with CDS - codon changes\n"); logfile.flush()
    outex = open(outexons,"a")
    for x in gm_stream.intersect([inclstream, annotstream]):
        # x = (chr,   start,end, (alt1,alt2,   , start1,end1,cds1,strand1,phase1,  start2,end2,cds2,strand2,phase2    ))
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
                continue
            codon_start = pos-shift
            ref_codon = assembly.fasta_from_regions({chr: [[codon_start,codon_start+3]]}, out={},
                                                    path_to_ref=genomeRef.get(chr))[0][chr][0]
            info = [chr,pos,refbase,list(rest[1:nsamples+1]),cds,strand,
                    ref_codon.upper(),shift]
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
    return True

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
    for name, tr in outtracks.iteritems():
        tr[0].close()
        description = set_file_descr(name+"_SNPs.bed.gz",type='bed',step='tracks',gdv='1',ucsc='1')
        ex.add(tr[1], description=description)


@timer
def snp_workflow(ex, job, assembly, minsnp=40., mincov=5, path_to_ref=None, via='local',
                 logfile=sys.stdout, debugfile=sys.stderr):
    """Main function of the workflow"""
    logfile.write('Current working dir: %s'%os.getcwd()); logfile.flush()
    ref_genome = assembly.fasta_by_chrom

    sample_names = []
    for gid in sorted(job.files.keys()):
        sample_name = job.groups[gid]['name']
        sample_names.append(sample_name)

    logfile.write("\n* Generate vcfs for each chrom/group\n"); logfile.flush()
    vcfs = dict((chrom,{}) for chrom in ref_genome.keys()) # {chr: {}}
    bams = dict((chrom,{}) for chrom in ref_genome.keys()) # {chr: {}}
    bam = {}
    # Launch the jobs
    for gid in sorted(job.files.keys()):
        sample_name = job.groups[gid]['name']
        # Merge all bams belonging to the same group
        runs = [r['bam'] for r in job.files[gid].itervalues()]
        bam = pysam.Samfile(runs[0])
        header = bam.header
        headerfile = unique_filename_in()
        for h in header["SQ"]:
            if h["SN"] in assembly.chrmeta:
                h["SN"] = assembly.chrmeta[h["SN"]]["ac"]
        head = pysam.Samfile( headerfile, "wh", header=header )
        head.close()
        if len(runs) > 1:
            bam = merge_bam(ex,runs)
        else:
            bam = replace_bam_header(ex,headerfile,runs[0],stdout=unique_filename_in())
        index_bam(ex,bam)
        # Samtools mpileup + bcftools + vcfutils.pl
        for chrom,ref in ref_genome.iteritems():
            future = pileup.nonblocking(ex, bam, ref, via=via)
            vcfs[chrom][gid] = future
            bams[chrom][gid] = bam
        logfile.write("  ...Group %s running.\n" % sample_name); logfile.flush()
    # Wait for vcfs to finish and store them in *vcfs[chrom][gid]*
    for gid in sorted(job.files.keys()):
        sample_name = job.groups[gid]['name']
        for chrom,ref in ref_genome.iteritems():
            vcf = vcfs[chrom][gid].wait()
            vcfs[chrom][gid] = vcf
        logfile.write("  ...Group %s done.\n" % sample_name); logfile.flush()
    # Targz the pileup files (vcf)
    for chrom,v in vcfs.iteritems():
        tarname = unique_filename_in()
        tarfh = tarfile.open(tarname, "w:gz")
        for gid,vcf in v.iteritems():
            tarfh.add(vcf, arcname="%s_%s.vcf" % (job.groups[gid]['name'],chrom))
        tarfh.close()
        ex.add( tarname, description=set_file_descr("vcfs_%s.tar.gz"%chrom,step="pileup",type="tar",view='admin') )

    logfile.write("\n* Merge info from vcf files\n"); logfile.flush()
    outall = unique_filename_in()
    outexons = unique_filename_in()
    with open(outall,"w") as fout:
        fout.write('#'+'\t'.join(['chromosome','position','reference']+sample_names+ \
                                 ['gene','location_type','distance'])+'\n')
    with open(outexons,"w") as fout:
        fout.write('#'+'\t'.join(['chromosome','position','reference']+sample_names+['exon','strand','ref_aa'] \
                                  + ['new_aa_'+s for s in sample_names])+'\n')
    for chrom,v in vcfs.iteritems():
        logfile.write("  > Chromosome '%s'\n" % chrom); logfile.flush()
    # Put together info from all vcf files
        logfile.write("  - All SNPs\n"); logfile.flush()
        allsnps = all_snps(ex,chrom,vcfs[chrom],bams[chrom],outall,assembly,
                           sample_names,mincov,float(minsnp),logfile,debugfile)
    # Annotate SNPs and check synonymy
        logfile.write("  - Exonic SNPs\n"); logfile.flush()
        exon_snps(chrom,outexons,allsnps,assembly,sample_names,ref_genome,logfile,debugfile)
    description = set_file_descr("allSNP.txt",step="SNPs",type="txt")
    ex.add(outall,description=description)
    description = set_file_descr("exonsSNP.txt",step="SNPs",type="txt")
    ex.add(outexons,description=description)
    # Create UCSC bed tracks
    logfile.write("\n* Create tracks\n"); logfile.flush()
    create_tracks(ex,outall,sample_names,assembly)
    return 0




import pysam, re
import json
import os
import pickle
from bbcflib import frontend
from numpy import *
from scipy.misc import factorial
from bein import *
from bein.util import *

############ Preprocessing ############
@program
def bamstats(bamfile):
    """Runs the Bamstat program and parses the output. Returns a dictionary.
    """
    def extract_pairs(s,head,foot):
        m=re.search(head+r'\n([\d\s]*)\n'+foot,s,
                    flags=re.MULTILINE).groups()[0].splitlines()
        def f(x):
            (a,b) = re.search(r'(\d+)\s(\d+)',x).groups()
            return (int(a),int(b))
        return dict([f(x) for x in m])
    def coverage_stats(p):
        results = {}        
        s=''.join(p.stdout)
        results["read_length"]=int(re.search(r'Read length (\d+)',s).groups()[0])
        results["genome_size"]=int(re.search(r'Genome size (\d+)',s).groups()[0])
        results["nb_positions"]=int(re.search(r'Nb positions (\d+)',s).groups()[0])
        results["multi_hits"]=extract_pairs(s,"Hits Reads","Total")
        results["total"]=int(re.search(r'Total (\d+)',s).groups()[0])
        [total,fwd,rev]=re.search(r'Alignments (\d+)\s*\(fwd:\s+(\d+)/rev:\s+(\d+)\)',
                                  s,flags=re.MULTILINE).groups()
        results["alignments"]={"total": int(total),
                                "fwd": int(fwd),
                                "rev": int(rev)}
        results["unmapped"]=int(re.search(r'Unmapped ([\d.]+)',s).groups()[0])
        results["expected_coverage"]=float(re.search(r'Expected coverage ([\d.]+)',
                                                     s).groups()[0])
        results["actual_coverage"]=float(re.search(r'Actual coverage ([\d.]+)',
                                                   s).groups()[0])
        results["mismatches"]=extract_pairs(s,"Mismatches Reads","")
        return results
    return {"arguments": ["bamstat",bamfile], "return_value": coverage_stats}

@program
def plot_stats(sample_stats,script_path="./"):
    json_string = json.dumps(sample_stats)
    pdf_file = unique_filename_in()
    return {'arguments': ["R","--vanilla","--slave","-f",
                          os.path.join(script_path,"pdfstats.R"),
                          "--args"] + [json_string,pdf_file],
            'return_value': pdf_file}

def poisson_threshold(mu, cutoff=0.95, max_terms=100):
    """Calculate confidence threshold for Poisson distributions.

    Returns the largest integer k such that, for a Poisson
    distribution random value X with mean 'mu', P(X <= k) <= 'cutoff'.
    It will calculate no farther than k = 'max_terms'.  If it reaches
    that point, it raises an exception.
    """
    p = cumsum( exp(-mu) * array([mu**k / float(factorial(k)) for k in range(0,max_terms)] ))
    n = len(p[p <= cutoff])
    if n == max_terms:
        raise ValueError("In poisson_threshold, reached max_terms. Try raising max_terms.")
    else:
        return n
    
def remove_duplicate_reads( bamfile, chromosomes,
                            maxhits=None, pilesize=1, convert=False ):
    """Filters a bam file for multi-hits above 'maxhits' and for duplicate reads beyond 'pilesize'.

    Reads with NH tag > maxhits are discarded, each genomic position 
    will have at most 'pilesize' reads per library and per strand. 
    If the 'convert' flag is True, the reference sequence ids are replaced by 
    their names as provided in 'chromosomes'.
    """
    infile = pysam.Samfile( bamfile, "rb" )
    outname = unique_filename_in()
    header = infile.header
    pilesize = max(1,pilesize)
    if convert:
        for h in header["SQ"]:
            if h["SN"] in chromosomes:
                h["SN"] = chromosomes[h["SN"]]["name"]
    outfile = pysam.Samfile( outname, "wb", header=header )
    count_per_lib = {}
    pos_per_lib = {}
    for read in infile:
        nh = dict(read.tags).get('NH')
        if nh == None:
            nh = 1
        if nh < 1:
            continue
        lname = re.search(r'^(.*?:.*?):',read.qname).groups()[0]
        lib = lname+":"+(read.is_reverse and '1' or '0')
        pos = "%s:%d" % (read.rname, read.pos)
        if pos != pos_per_lib.get(lib):
            pos_per_lib[lib] = pos
            count_per_lib[lib] = 0
        if (maxhits == None or nh <= maxhits) and count_per_lib[lib] < pilesize:
            outfile.write(read)
        count_per_lib[lib] += 1
    outfile.close()
    infile.close()
    return outname
############################################################
 
def map_reads( ex, fastq_file, chromosomes, bowtie_index,
               maxhits=5, antibody_enrichment=50, name='',
               remove_pcr_duplicates=True ):
    if len(name) > 0:
        name += ' '
    bwtarg = ["-Sam",str(max(20,maxhits)),"--best","--strata"]
    if count_lines( ex, fastq_file )>10000000:
        bam = parallel_bowtie_lsf( ex, bowtie_index, fastq_file,
                                   n_lines=8000000,
                                   bowtie_args=bwtarg,
                                   add_nh_flags=True )
    else:
        future = bowtie.lsf( ex, bowtie_index, fastq_file, bwtarg )
        samfile = future.wait()
        bam = add_nh_flag( samfile )
    sorted_bam = add_and_index_bam( ex, bam, name+"complete bam file from bowtie" )
    stats = bamstats( ex, sorted_bam )
    add_pickle( ex, stats, name+"bamstat output" )
    if remove_pcr_duplicates:
        thresh = poisson_threshold( antibody_enrichment*stats["actual_coverage"] )
        add_pickle( ex, thresh, name+"Poisson threshold" )
        bam2 = remove_duplicate_reads( sorted_bam, chromosomes,
                                       maxhits, thresh, convert=True )
        reduced_bam = add_and_index_bam( ex, bam2, name+"filtered bam file" )
    else:
        reduced_bam = sorted_bam
    return {"full": sorted_bam, "reduced": reduced_bam, "stats": stats}

############################################################ 

def map_groups( ex, job, daflims, globals, genrep ):
    processed = {}
    g_rep_assembly = genrep.assembly( job.assembly_id )
    chromosomes = dict([(str(k[0])+"_"+k[1]+"."+str(k[2]),v) 
                        for k,v in g_rep_assembly.chromosomes.iteritems()])
    pcr_dupl = True
    if 'discard_pcr_duplicates' in job.options:
        pcr_dupl = job.options['discard_pcr_duplicates']
    for gid,group in job.groups.iteritems():
        processed[gid] = {}
        for rid,run in group['runs'].iteritems():
            file = daflims.fetch( str(run['facility']), str(run['machine']),
                                  run['run'], run['lane'], to=globals["fastq_root"] )
            m = map_reads( ex, file.values()[0], chromosomes, 
                           g_rep_assembly.index_path, name=str(rid),
                           remove_pcr_duplicates=pcr_dupl )
            processed[gid][file.keys()[0]] = m
    return processed

def add_pdf_stats( ex, processed, group_names, script_path, description = "mapping pdf report" ):
    all_stats = {}
    for gid, p in processed.iteritems():
        for lib_id, mapped in p.iteritems():
            all_stats[group_names[gid]+": "+lib_id] = mapped['stats']
    pdf = plot_stats(ex, all_stats, script_path=script_path)
    ex.add(pdf,description)
    return pdf

def get_mapseq_files( hts_key, minilims, lib_root, url ):
    def path_from_lib( M, id, root ):
        return os.path.relpath(M.path_to_file(id),root)
    bamfiles = {}
    htss = frontend.Frontend( url=url )
    job = htss.job( hts_key )
    try: 
        exid = max(minilims.search_executions(with_text=hts_key))
    except ValueError, v:
        raise ValueError("No execution with key "+hts_key)
    for group in job.groups.values():
        name = group['name']
        for rid,run in group['runs'].iteritems():
            both = minilims.search_files(with_text=str(rid)+" filtered bam file", source=('execution',exid))
            file = "_".join([name,run['machine'],str(run['run']),str(run['lane'])])+".bam"
            bamfiles[file+".bai"] = path_from_lib(minilims,both.pop(),lib_root)
            bamfiles[file] = path_from_lib(minilims,both.pop(),lib_root)
    names = "_".join([g['name'] for g in job.groups.values()])
    pdffile = {names+".pdf": path_from_lib(minilims,minilims.search_files(
        with_text="mapping pdf report",
        source=('execution',exid)).pop(),lib_root)}
    return {"bam": bamfiles, "pdf": pdffile}

def import_mapseq_results( hts_key, minilims, ex_root, url ):
    processed = {}
    htss = frontend.Frontend( url=url )
    job = htss.job( hts_key )
    try: 
        exid = max(minilims.search_executions(with_text=hts_key))
    except ValueError, v:
        raise ValueError("No execution with key "+hts_key)
    for gid, group in job.groups.iteritems():
        name = group['name']
        processed[gid] = {}
        for rid,run in group['runs'].iteritems():
            bams = minilims.search_files(with_text=str(rid)+" filtered bam file", source=('execution',exid))
            rname = "_".join([name,run['machine'],str(run['run']),str(run['lane'])])
            bamfile = os.path.join(ex_root, unique_filename_in(ex_root))
            minilims.export_file(bams.pop(),bamfile+".bai")
            minilims.export_file(bams.pop(),bamfile)
            stats = minilims.search_files(with_text=str(rid)+" bamstat output", source=('execution',exid))[0]
            with open(minilims.path_to_file(stats)) as q:
                s = pickle.load(q)
            processed[gid][rname] = {'reduced': bamfile, 'stats': s}
    return (processed,job)

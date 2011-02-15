import pysam
import re
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
    sorted_bam = add_and_index_bam( ex, bam, "bam:"+name+"complete.bam" )
    stats = bamstats( ex, sorted_bam )
    add_pickle( ex, stats, "py:"+name+"bamstat" )
    if remove_pcr_duplicates:
        thresh = poisson_threshold( antibody_enrichment*stats["actual_coverage"] )
        add_pickle( ex, thresh, "py:"+name+"Poisson_threshold" )
        bam2 = remove_duplicate_reads( sorted_bam, chromosomes,
                                       maxhits, thresh, convert=True )
        reduced_bam = add_and_index_bam( ex, bam2, "bam:"+name+"filtered.bam" )
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
    file_names = {}
    for gid,group in job.groups.iteritems():
        processed[gid] = {}
        file_names[gid] = {}
        group_name = re.sub(r'\s+','_',group['name'])
        for rid,run in group['runs'].iteritems():
            fq_file = daflims[run['facility']].fetch_fastq( str(run['facility']),
                                                            str(run['machine']),
                                                            run['run'], run['lane'],
                                                            to=globals['fastq_root'] )
            if len(group['runs'])>1:
                name = "_".join([run['machine'],str(run['run']),str(run['lane'])])
            else:
                name = group_name
            m = map_reads( ex, fq_file['path'], chromosomes, 
                           g_rep_assembly.index_path, name=name+"_",
                           remove_pcr_duplicates=pcr_dupl )
            file_names[gid][rid] = name
            m.update({'libname': name})
            processed[gid][rid] = m
    add_pickle( ex, file_names, "py:file_names" )
    return processed

def add_pdf_stats( ex, processed, group_names, script_path,
                   description = "pdf:mapping_report.pdf" ):
    all_stats = {}
    for gid, p in processed.iteritems():
        i=1
        for mapped in p.values():
            name = group_names[gid]
            if 'libname' in mapped:
                name = mapped['libname']
            if name in all_stats:
                name += ":"+str(i)
            all_stats[name] = mapped['stats']
            i += 1
    pdf = plot_stats(ex, all_stats, script_path=script_path)
    ex.add(pdf,description)
    return pdf

def import_mapseq_results( hts_key, minilims, ex_root, url ):
    processed = {}
    htss = frontend.Frontend( url=url )
    job = htss.job( hts_key )
    try: 
        exid = max(minilims.search_executions(with_text=hts_key))
    except ValueError, v:
        raise ValueError("No execution with key "+hts_key)
    allfiles = dict((minilims.fetch_file(x)['description'],x)
                    for x in minilims.search_files(source=('execution',exid)))
    with open(minilims.path_to_file(allfiles['py:file_names'])) as q:
        file_names = pickle.load(q)
    for gid, group in job.groups.iteritems():
        group_name = re.sub(r'\s+','_',group['name'])
        processed[gid] = {}
        for rid,run in group['runs'].iteritems():
            bamfile = os.path.join(ex_root, unique_filename_in(ex_root))
            name = file_names[gid][rid] 
            bam_id = allfiles['bam:'+name+'filtered.bam']
            bam_bai_id = allfiles['bam:'+name+'filtered.bam (BAM index)']
            minilims.export_file(bam_id,bamfile)
            minilims.export_file(bam_bai_id,bamfile+".bai")
            stats_id = allfiles["py:"+name+"bamstat"]
            with open(minilims.path_to_file(stats_id)) as q:
                s = pickle.load(q)
            processed[gid][rid] = {'reduced': bamfile, 'stats': s, 
                                   'libname': name}
    return (processed,job)

def get_files( id_or_key, minilims ):
    if isinstance(id_or_key, str):
        try: 
            exid = max(minilims.search_executions(with_text=id_or_key))
        except ValueError, v:
            raise ValueError("No execution with key "+id_or_key)
    else:
        exid = id_or_key
    file_dict = {}
    all = dict((y['repository_name'],y['description']) for y in
               [minilims.fetch_file(x) for x in minilims.search_files(source=('execution',exid))])
    for f,d in all.iteritems():
        cat,name = re.search(r'([^:]+):(.*)$',d).groups()
        name = re.sub(r'\s+\(BAM index\)','.bai',name)
        if cat in file_dict:
            file_dict[cat].update({f: name})
        else:
            file_dict[cat] = {f: name}
    return file_dict

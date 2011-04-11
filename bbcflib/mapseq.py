"""
===============
bbcflib.mapseq
===============

This module provides functions useful to map raw reads using bowtie. 
The most common workflow will use ``map_reads`` which takes the following arguments:

  * ``'ex'``: an execution environment to run jobs in,

  * ``'fastq_file'``: the raw reads as a fastq file, 

  * ``'chromosomes'``: a dictionary with keys 'chromosome_id' as used in the bowtie indexes and values a dictionary with common chromosome names and lengths,

  * ``'bowtie_index'``: the file prefix to the bowtie index,

  * ``'maxhits'``: the maximum number of hits per read to keep in the output (default *5*),

  * ``'antibody_enrichment'``: an approximate enrichement ratio of protein bound sequences over controls (default *50*),

  * ``'name'``: sample name to be used in descriptions and file names,

  * ``'remove_pcr_duplicates'``: whether to remove probable PCR amplification artifacts based on a Poisson confidence threshold (default *True*).

The function ``map_groups`` will take a collection of sample as described in a *job* object from the ``frontend`` module and run fetch fastq files for each of them through using a ``daflims`` object, use ``genrep`` to get the bowtie indexes and run ``map_reads`` for each sample.
Below is the script used by the frontend::
    M = MiniLIMS( limspath )
    gl = use_pickle(M, "global variables")
    htss = frontend.Frontend( url=gl["hts_url"] )
    job = htss.job( hts_key )
    g_rep = genrep.GenRep( gl["genrep_url"], gl["bwt_root"] )
    daflims = dict((loc,daflims.DAFLIMS( username=gl['lims']['user'], 
                                         password=gl['lims']['passwd'][loc] ))
                   for loc in gl['lims']['passwd'].keys())
    with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
        processed = map_groups( ex, job, daflims, ex.working_directory, g_rep )
        pdf = add_pdf_stats( ex, processed,
                             dict((k,v['name']) for k,v in job.groups.iteritems()),
                             gl['script_path'] )

"""

import pysam
import re
import json
import os
import pickle
from bbcflib import frontend, genrep, daflims, common
from numpy import *
from scipy.misc import factorial
from bein import *
from bein.util import *

############ Preprocessing ############
@program
def bamstats(bamfile):
    """Wrapper to the ``bamstat`` program.

    This program computes read mapping statistics on a bam file. The output will 
    be parsed and converted to a dictionary.
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
    """Wrapper to the ``pdfstats.R`` script which generates 
    a pdf report of the mapping statistics.

    The input is the dictionary return by the ``bamstats`` call. 
    This is passed as a json string to the R script. 
    Returns the pdf file created by the script.
    """
    stats_file = unique_filename_in()
    with open( stats_file, 'w' ) as f:
        f.write(json.dumps(sample_stats))
    pdf_file = unique_filename_in()
    return {'arguments': ["R","--vanilla","--slave","-f",
                          os.path.join(script_path,"pdfstats.R"),
                          "--args"] + [stats_file,pdf_file],
            'return_value': pdf_file}

def poisson_threshold(mu, cutoff=0.95, max_terms=100):
    """Calculate confidence threshold for Poisson distributions.

    Returns the largest integer *k* such that, for a Poisson
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
               remove_pcr_duplicates=True, bwt_args=[], via='lsf' ):
    """Runs ``bowtie`` in parallel over lsf for the `fastq_file` input. 
    Returns the full bamfile, its filtered version (see 'remove_duplicate_reads') 
    and the mapping statistics dictionary (see 'bamstats').

    The input file will be split into subfiles if it contains more than 10M lines.
    The 'add_nh_flag' function will be called to add the number of hits per read 
    in the bowtie output.
    If 'remove_pcr_duplicates' is *True*, the 'chromosomes' and 'maxhits' arguments
    are passed to the 'remove_duplicate_reads' 
    function and the 'antibody_enrichment' will be used as input to 
    the 'poisson_threshold' function to compute its 'pilesize' argument. 

    The mapping statistics dictionary is pickled and added to the execution's 
    repository, as well as both the full and filtered bam files.
    """
    bwtarg = ["-Sam",str(max(20,maxhits)),"--best","--strata"]+bwt_args
    if count_lines( ex, fastq_file )>10000000:
        bam = parallel_bowtie( ex, bowtie_index, fastq_file,
                               n_lines=8000000,
                               bowtie_args=bwtarg,
                               add_nh_flags=True, via=via )
    else:
        future = bowtie.nonblocking( ex, bowtie_index, fastq_file, 
                                     bwtarg, via=via )
        samfile = future.wait()
        bam = add_nh_flag( samfile )
    sorted_bam = add_and_index_bam( ex, bam, "bam:"+name+"complete.bam" )
    full_stats = bamstats( ex, sorted_bam )
    add_pickle( ex, full_stats, "py:"+name+"full_bamstat" )
    return_dict = {"fullbam": sorted_bam}
    if remove_pcr_duplicates:
        thresh = poisson_threshold( antibody_enrichment*full_stats["actual_coverage"] )
        add_pickle( ex, thresh, "py:"+name+"Poisson_threshold" )
        bam2 = remove_duplicate_reads( sorted_bam, chromosomes,
                                       maxhits, thresh, convert=True )
        reduced_bam = add_and_index_bam( ex, bam2, "bam:"+name+"filtered.bam" )
        filtered_stats = bamstats( ex, reduced_bam )
        add_pickle( ex, filtered_stats, "py:"+name+"filter_bamstat" )
        return_dict['bam'] = reduced_bam 
        return_dict['fullstats'] = full_stats
        return_dict['stats'] = filtered_stats
    else:
        infile = pysam.Samfile( sorted_bam, "rb" )
        reduced_bam = unique_filename_in()
        header = infile.header
        for h in header["SQ"]:
            if h["SN"] in chromosomes:
                h["SN"] = chromosomes[h["SN"]]["name"]
        outfile = pysam.Samfile( reduced_bam, "wb", header=header )
        for read in infile:
            nh = dict(read.tags).get('NH')
            if nh == None:
                nh = 1
            if nh < 1:
                continue
            outfile.write(read)
        outfile.close()
        infile.close()
        return_dict['bam'] = reduced_bam 
        return_dict['stats'] = full_stats
    return return_dict

############################################################ 

def map_groups( ex, job_or_dict, daflims_or_files, fastq_root, genrep_or_dict, map_args={} ):
    """Fetches fastq files and bowtie indexes, and runs the 'map_reads' function for 
    a collection of samples described in a 'Frontend' 'job'.

    Arguments are:

    * ``'ex'``: a 'bein' execution environment to run jobs in,
    
    * ``'job_or_dict'``: a 'Frontend' 'job' object, or a dictionary with keys 'groups' and 'assembly_id',
    
    * ``'daflims_or_files'``: a dictionary of 'Daflims' objects (keys are the facility names), or a dictionary with run_ids as keys (same as in job.groups) and values file names,

    * ``'fastq_root'``: path where to download raw fastq files,

    * ``'genrep_or_dict'``: a 'GenRep' object, or a dictionary of 'chromosomes' and 'index_path'.

    * ``'map_args'``: a dictionary of arguments passed to map_reads.

    Returns a dictionary with keys *group_id* from the job object and values dictionaries mapping *run_id* to the corresponding return value of the 'map_reads' function.
    """
    processed = {}
    file_names = {}
    options = {}
    if isinstance(job_or_dict,frontend.Job):
        options = job_or_dict.options
        groups = job_or_dict.groups
        assembly_id = job_or_dict.assembly_id
    elif isinstance(job_or_dict,dict) and 'groups' in job_or_dict and 'assembly_id' in job_or_dict:
        if 'options' in job_or_dict:
            options = job_or_dict['options']
        groups = job_or_dict['groups']
        assembly_id = job_or_dict['assembly_id']
    else:
        raise TypeError("job_or_dict must be a frontend.Job object or a dictionary with keys 'groups' and 'assembly_id'.")
    pcr_dupl = True
    if 'discard_pcr_duplicates' in options:
        pcr_dupl = options['discard_pcr_duplicates']
    if isinstance(genrep_or_dict,genrep.GenRep):
        g_rep_assembly = genrep_or_dict.assembly( assembly_id )
        chromosomes = dict([(str(k[0])+"_"+k[1]+"."+str(k[2]),v) 
                            for k,v in g_rep_assembly.chromosomes.iteritems()])
        index_path = g_rep_assembly.index_path 
    elif isinstance(genrep_or_dict,dict) and 'chromosomes' in genrep_or_dict:
        chromosomes = genrep_or_dict['chromosomes']
        index_path = genrep_or_dict['index_path ']
    else:
        raise TypeError("genrep_or_dict must be a genrep.GenRep object or a dictionary with keys 'chromosomes' and 'index_path'.")
    for gid,group in groups.iteritems():
        processed[gid] = {}
        file_names[gid] = {}
        if 'name' in group:
            group_name = re.sub(r'\s+','_',group['name'])
        else:
            group_name = gid
        if not 'runs' in group:
            group = {'runs': group}
        for rid,run in group['runs'].iteritems():
            name = group_name
            if isinstance(daflims_or_files.values()[0],daflims.DAFLIMS):
                dafl = daflims_or_files[run['facility']]
                fq_file = dafl.fetch_fastq( str(run['facility']), str(run['machine']),
                                            run['run'], run['lane'], to=fastq_root )['path']
                if len(group['runs'])>1:
                    name += "_".join(['',run['machine'],str(run['run']),str(run['lane'])])
            else:
                fq_file = os.path.join(fastq_root,daflims_or_files[gid][rid])
                if len(group['runs'])>1:
                    name += "_"+str(rid)
            m = map_reads( ex, fq_file, chromosomes, index_path, name=name+"_",
                           remove_pcr_duplicates=pcr_dupl, **map_args )
            file_names[gid][rid] = name
            m.update({'libname': name})
            processed[gid][rid] = m
    add_pickle( ex, file_names, "py:file_names" )
    return processed

def add_pdf_stats( ex, processed, group_names, script_path,
                   description = "pdf:mapping_report.pdf" ):
    """Runs the 'plot_stats' function and adds its pdf output to the execution's repository.
    
    Arguments are the output of 'map_groups' ('processed'),
    a dictionary of group_id to names used in the display, 
    the path to the script used by 'plot_stats', 
    and the 'description' to use in the repository.

    Returns the name of the pdf file.
    """
    all_stats = {}
    for gid, p in processed.iteritems():
        for i,mapped in enumerate(p.values()):
            name = group_names[gid]
            if 'libname' in mapped:
                name = mapped['libname']
            if name in all_stats:
                name += ":"+str(i+1)
            if 'fullstats' in mapped:
                all_stats[name+":full"] = mapped['fullstats']
                all_stats[name+":filter"] = mapped['stats']
            else:
                all_stats[name] = mapped['stats']
    pdf = plot_stats(ex, all_stats, script_path=script_path)
    ex.add(pdf,description)
    return pdf

############################################################ 
@program
def wigToBigWig( sql ):
    """Binds ``wigToBigWig`` from the UCSC tools.
    """
    chrsizes = unique_filename_in()
    chromosomes = []
    connection = sqlite3.connect( sql )
    cur = connection.cursor()
    with open(chrsizes,'w') as f:
        cur = connection.cursor()
        cur.execute('select * from chrNames')
        connection.commit()
        for sql_row in cur:
            chromosomes.append(sql_row[0])
            f.write(' '.join([str(x) for x in sql_row])+"\n")
        cur.close()
    bedgraph = unique_filename_in()
    with open(bedgraph,'w') as f:
        for c in chromosomes:
            cur.execute('select * from "'+c+'"')
            connection.commit()
            for sql_row in cur:
                f.write("\t".join([c]+[str(x) for x in sql_row])+"\n")
            cur.close()
    bigwig = unique_filename_in()
    return {"arguments": ['wigToBigWig',bedgraph,chrsizes,bigwig], 
            "return_value": bigwig}

@program
def bam_to_density( bamfile, chromosome_accession, chromosome_name, output,
                    nreads=1, merge=-1, read_length=-1, convert=True, sql=False,
                    args=[] ):
    """Runs the ``bam2wig`` program on a bam file and 
    normalizes for the total number of reads
    provided as argument 'nreads'. 

    Returns the name of the output wig or sql file(s) (if 'sql' is True).

    Use 'convert'=False if the bam already uses chromosome names instead of ids.
    """
    b2w_args = ["-w",str(nreads),"-s",bamfile,"-o",output]
    if convert:
        b2w_args += ["-a",chromosome_accession,"-n",chromosome_name]
    else:
        b2w_args += ["-a",chromosome_name]
    if merge>=0:
        b2w_args += ["-p",str(merge)]
    if read_length>0:
        b2w_args += ["-q",str(read_length)]
    if sql:
        b2w_args += ["-d"]
        if merge<0:
            files = [output+"fwd.sql",output+"rev.sql"]
        else:
            files = [output+"merged.sql"]
    else:
        if merge<0:
            b2w_args += ["-6"]
        files = output
    b2w_args += args
    return {"arguments": ["bam2wig"]+b2w_args, "return_value": files}

def compact_chromosome_name(key):
    if isinstance(key,str):
        return key
    elif isinstance(key,tuple) and len(key)>2:
        return str(key[0])+"_"+key[1]+"."+str(key[2])
    else:
        raise ValueError("Can't handle this chromosomes key ",key)

def parallel_density_wig( ex, bamfile, chromosomes, 
                          nreads=1, merge=-1, read_length=-1, convert=True, 
                          description="", alias=None, 
                          b2w_args=[], via='lsf' ):
    """Runs 'bam_to_density' in parallel 
    for every chromosome in the 'chromosomes' list with 'sql' set to False.
    Returns produces a single text wig file.
    """
    futures = [bam_to_density.nonblocking( ex, bamfile, compact_chromosome_name(k),
                                           v['name'], unique_filename_in(), nreads, merge, 
                                           read_length, convert, False, 
                                           args=b2w_args, via=via )
               for k,v in chromosomes.iteritems()]
    results = []
    for f in futures:
        try: 
            results.append(f.wait())
        except ProgramFailed:
            pass
    output = common.cat(results)
    ex.add( output, description=description, alias=alias )
    return output

def parallel_density_sql( ex, bamfile, chromosomes, 
                          nreads=1, merge=-1, read_length=-1, convert=True, 
                          description="", alias=None, b2w_args=[], via='lsf' ):
    """Runs 'bam_to_density' for every chromosome in the 'chromosomes' list.
    
    Returns one or two sqlite files depending 
    if 'merge'>=0 (shift and merge strands into one tracks) 
    or 'merge'<0 (keep seperate tracks for each strand).
    """
    output = unique_filename_in()
    touch(ex,output)
    if merge<0:
        suffix = ['fwd','rev']
    else:
        suffix = ['merged']
    for s in suffix:
        conn1 = sqlite3.connect( output+s+'.sql' )
        conn1.execute('create table chrNames (name text, length integer)')
        conn1.execute('create table attributes (key text, value text)')
        conn1.execute('insert into attributes (key,value) values ("datatype","quantitative")')
        vals = [(v['name'],str(v['length'])) for v in chromosomes.values()]
        conn1.executemany('insert into chrNames (name,length) values (?,?)',vals)
        [conn1.execute('create table "'+v['name']+'" (start integer, end integer, score real)') 
         for v in chromosomes.values()]
        [conn.execute('create index '+v['name']+'_range_idx on "'+v['name']+'" (start, end)') 
         for v in chromosomes.values()]
        conn1.commit()
        conn1.close()
    results = []
    for k,v in chromosomes.iteritems():
        future = bam_to_density.nonblocking( ex, bamfile, compact_chromosome_name(k),
                                             v['name'], output, nreads, merge,
                                             read_length, convert, True, 
                                             args=b2w_args, via=via )
        try: 
            results.append(future.wait())
        except ProgramFailed:
            pass
    ex.add( output, description='none:'+description+'.sql', alias=alias )
    [ex.add( output+s+'.sql', description='sql:'+description+'_'+s+'.sql',
             associate_to_filename=output, template='%s_'+s+'.sql' ) 
     for s in suffix]
    return dict((s,output+s+'.sql') for s in suffix)

############################################################ 
def import_mapseq_results( key_or_id, minilims, ex_root, url_or_dict ):
    """Imports all files created by a previous 'mapseq' workflow into the current execution environement.

    * ``'key_or_id'`` is the previous execution's description in the MiniLIMS or its id,

    * ``'minilims'`` is the MiniLIMS where that execution was saved,

    * ``'ex_root'`` is the current execution's directory (where files will be copied to),

    * ``'url_or_dict'`` is either the 'Frontend' url to fetch the mapseq job's description, or a job description dictionary (see ``map_groups``)

    Returns a tuple with a dictionary similar to the output of 'map_groups' and a 'job' object describing the mapseq runs.
    """
    processed = {}
    if isinstance(url_or_dict, str):
        htss = frontend.Frontend( url=url_or_dict )
        job = htss.job( key_or_id )
        job_groups = job.groups
    else:
        job = url_or_dict
        job_groups = job['groups']
    if isinstance(key_or_id, str):
        try: 
            exid = max(minilims.search_executions(with_text=key_or_id))
        except ValueError, v:
            raise ValueError("No execution with key "+key_or_id)
    else:
        exid = key_or_id
    allfiles = dict((minilims.fetch_file(x)['description'],x)
                    for x in minilims.search_files(source=('execution',exid)))
    with open(minilims.path_to_file(allfiles['py:file_names'])) as q:
        file_names = pickle.load(q)
    for gid, group in job_groups.iteritems():
        if 'name' in group:
            group_name = re.sub(r'\s+','_',group['name'])
        else:
            group_name = gid
        if not 'runs' in group:
            group = {'runs': group}
        processed[gid] = {}
        for rid,run in group['runs'].iteritems():
######## sql files....
            bamfile = os.path.join(ex_root, unique_filename_in(ex_root))
            name = file_names[gid][rid] 
            bam_id = allfiles['bam:'+name+'_filtered.bam']
            bam_bai_id = allfiles['bam:'+name+'_filtered.bam (BAM index)']
            minilims.export_file(bam_id,bamfile)
            minilims.export_file(bam_bai_id,bamfile+".bai")
            stats_id = allfiles["py:"+name+"_filter_bamstat"]
            with open(minilims.path_to_file(stats_id)) as q:
                s = pickle.load(q)
            processed[gid][rid] = {'bam': bamfile, 'stats': s, 'libname': name}
    return (processed,job)


"""
======================
Module: bbcflib.mapseq
======================

This module provides functions useful to map raw reads using bowtie.
The most common workflow will first use ``map_reads`` which takes the following arguments:

  * ``'ex'``: an execution environment to run jobs in,

  * ``'fastq_file'``: the raw reads as a fastq file,

  * ``'chromosomes'``: a dictionary with keys 'chromosome_id' as used in the bowtie indexes and values a dictionary with usual chromosome names and their lengths,

  * ``'bowtie_index'``: the file prefix to the bowtie index,

  * ``'maxhits'``: the maximum number of hits per read to keep in the output (default *5*),

  * ``'antibody_enrichment'``: an approximate enrichement ratio of protein bound sequences over controls (default *50*),

  * ``'name'``: sample name to be used in descriptions and file names,

  * ``'remove_pcr_duplicates'``: whether to remove probable PCR amplification artifacts based on a Poisson confidence threshold (default *True*),

  * ``'bwt_args'``: a list of specific arguments for bowtie (in addition to ["-Sam",str(max(20,maxhits)),"--best","--strata"]).

The function ``map_groups`` will take a collection of sample as described in a *job* object from the ``frontend`` module and run fetch fastq files for each of them through using a ``daflims`` object, use an 'Assembly' object from ``genrep`` to get the bowtie indexes and run ``map_reads`` for each sample.

The second step in the workflow consists in generating genome-wide density maps with ``densities_groups`` which needs a 'bein' execution environment, the same *job* and *assembly* as above and a dictionary 'file_dict' as returned by ``map_groups``. The function then runs ``parallel_density_sql`` on each 'bam' file to obtain a normalized read density profile as a 'sqlite' file.

Below is the script used by the frontend::

    from bbcflib import daflims, genrep, frontend, gdv, common
    from bbcflib.mapseq import *
    M = MiniLIMS( limspath )
    working_dir = '/path/to/scratch/on/cluster'
    hts_key = 'test_key'
    gl = { 'hts_mapseq': {'url': 'http://htsstation.vital-it.ch/mapseq/'},
           'genrep_url': 'http://bbcftools.vital-it.ch/genrep/',
           'bwt_root': '/db/genrep/nr_assemblies',
           'script_path': '/srv/chipseq/lib',
           'lims': {'user': 'alice',
                     'passwd': {'lgtf': 'bob123',
                                'gva': 'joe456'}},
            'gdv': {'url': 'http://svitsrv25.epfl.ch/gdv',
                    'email': 'alice.ecila@somewhere.edu',
                    'key': 'xxxxxxxxxxxxxxxxxxxxxxxx'} }
    assembly_id = 'mm9'
    htss = frontend.Frontend( url=gl['hts_mapseq']['url'] )
    job = htss.job( hts_key )
    g_rep = genrep.GenRep( gl['genrep_url'], gl['bwt_root'] )
    assembly = g_rep.assembly( assembly_id )
    daflims1 = dict((loc,daflims.DAFLIMS( username=gl['lims']['user'],
                                          password=gl['lims']['passwd'][loc] ))
                    for loc in gl['lims']['passwd'].keys())
    job.options['ucsc_bigwig'] = True
    with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
        job = get_fastq_files( job, ex.working_directory, daflims1 )
        mapped_files = map_groups( ex, job, ex.working_directory, assembly )
        pdf = add_pdf_stats( ex, mapped_files,
                             dict((k,v['name']) for k,v in job.groups.iteritems()),
                             gl['script_path'] )
        density_files = densities_groups( ex, job, mapped_files, assembly.chromosomes )
        gdv_project = gdv.create_gdv_project( gl['gdv']['key'], gl['gdv']['email'],
                                              job.description, hts_key,
                                              assembly.nr_assembly_id,
                                              gdv_url=gl['gdv']['url'],
                                              public=True )
        add_pickle( ex, gdv_project, description='py:gdv_json' )
    print ex.id
    allfiles = common.get_files( ex.id, M )
    print allfiles
"""

# Built-in modules #
import os, re, json, shutil, gzip, tarfile, pickle

# Internal modules #
from . import frontend, genrep, daflims, common

# Other modules #
import pysam
from numpy import *
from scipy.misc import factorial
from bein import *
from bein.util import *

################################################################################
# Preprocessing #
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
    This is passed as a json file to the R script.
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

def pprint_bamstats(sample_stats) :
    """Pretty stdout-print for sample_stats.
    
    The input is the dictionary return by the ``bamstats`` call.
    """
    span = 5
    width_left = max([len(x) for x in sample_stats.keys()]) + span
    width_right = max([len(str(x)) for x in sample_stats.values()]) + span
    width_table = width_left + width_right +7
    print "{0:->{twh}}".format("", twh=width_table)
    for k, v in sample_stats.iteritems() :
        print "* {0:{lwh}} | {1:>{rwh}} *".format(k,v, lwh=width_left,rwh=width_right)
    print "{0:->{twh}}".format("", twh=width_table)
    return(0)



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
        future = bowtie.nonblocking( ex, bowtie_index, fastq_file, bwtarg, via=via )
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
        bam2 = unique_filename_in()
        header = infile.header
        for h in header["SQ"]:
            if h["SN"] in chromosomes:
                h["SN"] = chromosomes[h["SN"]]["name"]
        outfile = pysam.Samfile( bam2, "wb", header=header )
        for read in infile:
            nh = dict(read.tags).get('NH')
            if nh == None:
                nh = 1
            if nh < 1:
                continue
            if maxhits == None or nh <= maxhits:
                outfile.write(read)
        outfile.close()
        infile.close()
        reduced_bam = add_and_index_bam( ex, bam2, "bam:"+name+"filtered.bam" )
        return_dict['bam'] = reduced_bam
        return_dict['stats'] = full_stats
    return return_dict

############################################################

def get_fastq_files( job, fastq_root, dafl=None, set_seed_length=True ):
    """
    Will replace file references by actual file paths in the 'job' object.
    These references are either 'dafl' run descriptions or urls.
    Argument 'dafl' is a dictionary of 'Daflims' objects (keys are the facility names).
    If 'set_seed_length' is true, a dictionary job.groups[gid]['seed_lengths']
    of seed lengths is constructed with values corresponding to 70% of the
    read length.
    """
    if  not(dafl == None or isinstance(dafl.values()[0],daflims.DAFLIMS)):
        raise ValueError("Need DAFLIMS objects in get_fastq_files.")
    for gid,group in job.groups.iteritems():
        job.groups[gid]['seed_lengths'] = {}
        job.groups[gid]['run_names'] = {}
        for rid,run in group['runs'].iteritems():
            if len(run['url'])>0:
                run = str(run['url'])
            if isinstance(run,dict) and all([x in run for x in ['facility','machine','run','lane']]):
                dafl1 = dafl[run['facility']]
                daf_data = dafl1.fetch_fastq( str(run['facility']), str(run['machine']),
                                             run['run'], run['lane'], to=fastq_root )
                job.groups[gid]['runs'][rid] = daf_data['path']
                if (set_seed_length):
                    job.groups[gid]['seed_lengths'][rid] = max(28,int(0.7*daf_data['cycle']))
                job.groups[gid]['run_names'][rid] = "_".join(['',run['machine'],str(run['run']),str(run['lane'])])
            elif isinstance(run,str):
                fq_file = unique_filename_in(fastq_root)
                target = os.path.join(fastq_root,fq_file)
                if run.startswith("http://") or run.startswith("https://"):
                    urllib.urlretrieve( run, target )
                else:
                    shutil.copy(run,target)
#                run = re.sub('.seq.gz','_seq.tar',run)
                is_gz = run.endswith(".gz") or run.endswith(".gzip")
                run_strip = run
                if is_gz:
                    run_strip = re.sub('.gz[ip]*','',run)
                is_tar = run_strip.endswith(".tar")
                if is_tar:
                    mode = 'r'
                    if is_gz:
                        mode += '|gz'
                    fq_file2 = unique_filename_in(fastq_root)
                    target2 = os.path.join(fastq_root,fq_file2)
                    with open(target,'rb') as tf:
                        tar = tarfile.open(fileobj=tf, mode=mode)
                        tar_filename = tar.next()
                        input_file = tar.extractfile(tar_filename)
                        with open(target2, 'w') as output_file:
                            while True:
                                chunk = input_file.read(4096)
                                if chunk == '':
                                    break
                                else:
                                    output_file.write(chunk)
                        input_file.close()
                        tar.close()
                    fq_file = fq_file2
                elif is_gz:
                    fq_file2 = unique_filename_in(fastq_root)
                    target2 = os.path.join(fastq_root,fq_file2)
                    with open(target2,'w') as output_file:
                        with gzip.open(target, 'rb') as input_file:
                            output_file.write(input_file.read())
                    fq_file = fq_file2
                job.groups[gid]['runs'][rid] = fq_file
                job.groups[gid]['run_names'][rid] = os.path.basename(run.split("/")[-1])
    return job

############################################################

def map_groups( ex, job_or_dict, fastq_root, assembly_or_dict, map_args={} ):
    """Fetches fastq files and bowtie indexes, and runs the 'map_reads' function for
    a collection of samples described in a 'Frontend' 'job'.

    Arguments are:

    * ``'ex'``: a 'bein' execution environment to run jobs in,

    * ``'job_or_dict'``: a 'Frontend' 'job' object, or a dictionary with keys 'groups',

    * ``'fastq_root'``: path to raw fastq files,

    * ``'assembly_or_dict'``: an 'Assembly' object, or a dictionary of 'chromosomes' and 'index_path'.

    * ``'map_args'``: a dictionary of arguments passed to map_reads.

    Returns a dictionary with keys *group_id* from the job object and values dictionaries
    mapping *run_id* to the corresponding return value of the 'map_reads' function.
    """
    processed = {}
    file_names = {}
    options = {}
    if isinstance(job_or_dict,frontend.Job):
        options = job_or_dict.options
        groups = job_or_dict.groups
    elif isinstance(job_or_dict,dict) and 'groups' in job_or_dict:
        if 'options' in job_or_dict:
            options = job_or_dict['options']
        groups = job_or_dict['groups']
    else:
        raise TypeError("job_or_dict must be a frontend.Job object or a dictionary with keys 'groups'.")
    pcr_dupl = True
    if 'discard_pcr_duplicates' in options:
        pcr_dupl = options['discard_pcr_duplicates']
        if isinstance(pcr_dupl,str):
            pcr_dupl = pcr_dupl.lower() in ['1','true','t']
    if isinstance(assembly_or_dict,genrep.Assembly):
        chromosomes = dict([(str(k[0])+"_"+k[1]+"."+str(k[2]),v)
                            for k,v in assembly_or_dict.chromosomes.iteritems()])
        index_path = assembly_or_dict.index_path
    elif isinstance(assembly_or_dict,dict) and 'chromosomes' in assembly_or_dict:
        chromosomes = assembly_or_dict['chromosomes']
        index_path = assembly_or_dict['index_path ']
    else:
        raise TypeError("assembly_or_dict must be a genrep.Assembly object or a dictionary \
                         with keys 'chromosomes' and 'index_path'.")
    for gid,group in groups.iteritems():
        processed[gid] = {}
        file_names[gid] = {}
        if 'name' in group:
            group_name = re.sub(r'\s+','_',group['name'])
        else:
            group_name = str(gid)
        if not 'runs' in group:
            group = {'runs': group}
        for rid,run in group['runs'].iteritems():
            fq_file = os.path.join(fastq_root,run)
            if 'seed_lengths' in group and group['seed_lengths'].get(rid) > 0:
                seed_len = str(group['seed_lengths'][rid])
                if 'bwt_args' in map_args:
                    if "-l" in map_args['bwt_args']:
                        map_args['bwt_args'][map_args['bwt_args'].index("-l")+1] = seed_len
                    else:
                        map_args['bwt_args'] += ["-l",seed_len]
                else:
                    map_args['bwt_args'] = ["-l",seed_len]
            name = group_name
            if len(group['runs'])>1:
                name += "_"
                name += group['run_names'].get(rid) or str(rid)
            m = map_reads( ex, fq_file, chromosomes, index_path, name=name+"_",
                           remove_pcr_duplicates=pcr_dupl, **map_args )
            file_names[gid][rid] = str(name)
            m.update({'libname': str(name)})
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
            name = group_names.get(gid)
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
def bam_to_density( bamfile, output, chromosome_accession=None, chromosome_name=None,
                    nreads=1, merge=-1, read_extension=-1, convert=True, sql=False,
                    args=[] ):
    """Runs the ``bam2wig`` program on a bam file and
    normalizes for the total number of reads
    provided as argument 'nreads'.

    Returns the name of the output wig or sql file(s) (if 'sql' is True).

    Use 'convert'=False if the bam already uses chromosome names instead of ids.
    """
    b2w_args = ["-w",str(nreads),"-s",bamfile,"-o",output]
    if chromosome_accession != None:
        if convert and chromosome_name != None:
            b2w_args += ["-a",chromosome_accession,"-n",chromosome_name]
        else:
            b2w_args += ["-a",chromosome_name]
    if merge>=0 and not('-p' in b2w_args):
        b2w_args += ["-p",str(merge)]
    if read_extension>0 and not('-q' in b2w_args):
        b2w_args += ["-q",str(read_extension)]
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
    if key == None or isinstance(key,str):
        return key
    elif isinstance(key,tuple) and len(key)>2:
        return str(key[0])+"_"+str(key[1])+"."+str(key[2])
    else:
        raise ValueError("Can't handle this chromosomes key ",key)

def parallel_density_wig( ex, bamfile, chromosomes,
                          nreads=1, merge=-1, read_extension=-1, convert=True,
                          description="", alias=None,
                          b2w_args=[], via='lsf' ):
    """Runs 'bam_to_density' in parallel
    for every chromosome in the 'chromosomes' list with 'sql' set to False.
    Returns a single text wig file.
    """
    futures = [bam_to_density.nonblocking( ex, bamfile, unique_filename_in(),
                                           compact_chromosome_name(k), v['name'],
                                           nreads, merge, read_extension, convert,
                                           False, args=b2w_args, via=via )
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

def parallel_density_sql( ex, bamfile, output, chromosomes,
                          nreads=1, merge=-1, read_extension=-1, convert=True,
                          b2w_args=[], via='lsf' ):
    """Runs 'bam_to_density' for every chromosome in the 'chromosomes' list.

    Returns a dictionary with one or two sqlite files depending
    if 'merge'>=0 (shift and merge strands into one tracks)
    or 'merge'<0 (keep seperate tracks for each strand).
    """
    if chromosomes == None:
        chromosomes = {None: {'name': None}}
    for k,v in chromosomes.iteritems():
        future = bam_to_density.nonblocking( ex, bamfile, output,
                                             compact_chromosome_name(k), v['name'],
                                             nreads, merge, read_extension, convert,
                                             True, args=b2w_args, via=via )
        try:
            _ = future.wait()
        except ProgramFailed:
            pass
    return output

############################################################

def densities_groups( ex, job_or_dict, file_dict, chromosomes, via='lsf' ):
    """
    Arguments are:

    * ``'ex'``: a 'bein' execution environment to run jobs in,

    * ``'job_or_dict'``: a 'Frontend' 'job' object, or a dictionary with keys 'groups',

    * ``'file_dict'``: a dictionary of files,

    * ``'chromosomes'``: a dictionary with keys 'chromosome_id' as used in the bowtie indexes and values a dictionary with usual chromosome names and their lengths,

    Returns a dictionary with keys *group_id* from the job object and values the files fo each group ('bam' and 'wig').
    """
    processed = {}
    options = {}
    if isinstance(job_or_dict,frontend.Job):
        options = job_or_dict.options
        groups = job_or_dict.groups
    elif isinstance(job_or_dict,dict) and 'groups' in job_or_dict:
        if 'options' in job_or_dict:
            options = job_or_dict['options']
        groups = job_or_dict['groups']
    else:
        raise TypeError("job_or_dict must be a frontend.Job object or a dictionary with keys 'groups'.")
    merge_strands = -1
    suffixes = ["fwd","rev"]
    if options.get('merge_strands')>=0:
        merge_strands = options['merge_strands']
        suffixes = ["merged"]
    ucsc_bigwig = False
    if 'ucsc_bigwig' in options:
        ucsc_bigwig = options['ucsc_bigwig']
    b2w_args = []
    if 'b2w_args' in options:
        b2w_args = options['b2w_args']
    if options.get('read_extension')>0 and not('-q' in b2w_args):
        b2w_args += ["-q",str(options['read_extension'])]
    processed = {}
    for gid,group in groups.iteritems():
        if 'name' in group:
            group_name = re.sub(r'\s+','_',group['name'])
        else:
            group_name = gid
        mapped = file_dict[gid]
        if not isinstance(mapped,dict):
            raise TypeError("processed values must be dictionaries with keys *run_ids* or 'bam'.")
        if 'bam' in mapped:
            mapped = {'_': mapped}
        for k in mapped.keys():
            if not 'libname' in mapped[k]:
                mapped[k]['libname'] = group_name+"_"+str(k)
            if not 'stats' in mapped[k]:
                mapped[k]['stats'] = bamstats( ex, mapped[k]["bam"] )
            if not(options.get('read_extension')>0):
                options['read_extension'] = mapped[k]['stats']['read_length']
                if not('-q' in b2w_args):
                    b2w_args += ["-q",str(options['read_extension'])]
        wig = []
        for m in mapped.values():
            output = unique_filename_in()
            touch(ex,output)
            [common.create_sql_track( output+s+'.sql', chromosomes.values(), 
                                      name=m['libname'] ) 
             for s in suffixes]
            wig.append(parallel_density_sql( ex, m["bam"], output, chromosomes,
                                             nreads=m["stats"]["total"],
                                             merge=merge_strands,
                                             convert=False,
                                             b2w_args=b2w_args, via=via ))
            ex.add( output, description='none:'+m['libname']+'.sql' )
            [ex.add( output+s+'.sql', description='sql:'+m['libname']+'_'+s+'.sql',
                     associate_to_filename=output, template='%s_'+s+'.sql' )
             for s in suffixes]
        if len(mapped)>1:
            merged_bam = merge_bam(ex, [m['bam'] for m in mapped.values()])
            ids = [m['libname'] for m in mapped.values()]
            merged_wig = dict((s, common.merge_sql(ex, [x+s+".sql" for x in wig], ids,
                                                   description="sql:"+group_name+"_"+s+".sql"))
                              for s in suffixes)
        else:
            merged_bam = mapped.values()[0]['bam']
            merged_wig =  dict((s, wig[0]+s+".sql") for s in suffixes)
        processed[gid] = {'bam': merged_bam, 'wig': merged_wig,
                          'read_length': mapped.values()[0]['stats']['read_length']}
        if ucsc_bigwig:
            bw_futures = [wigToBigWig.nonblocking( ex, merged_wig[s], via=via )
                          for s in suffixes]
            [ex.add(bw_futures[i].wait(),description='bigwig:'+group_name+'_'+s+'.bw')
             for i,s in enumerate(suffixes)]
    processed.update({'read_extension': options.get('read_extension'),
                      'genome_size': mapped.values()[0]['stats']['genome_size']})
    return processed

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
    merge = -1
    if isinstance(url_or_dict, str):
        htss = frontend.Frontend( url=url_or_dict )
        job = htss.job( key_or_id )
        job_groups = job.groups
        if 'merge_strands' in job.options and job.options['merge_strands']>=0:
            merge = job.options['merge_strands']
    else:
        job = url_or_dict
        job_groups = job['groups']
        if 'merge_strands' in job['options'] and job['options']['merge_strands']>=0:
            merge = job['options']['merge_strands']
    if merge<0:
        suffix = ['fwd','rev']
    else:
        suffix = ['merged']
    if isinstance(key_or_id, str):
        try:
            exid = max(minilims.search_executions(with_text=key_or_id))
        except ValueError, v:
            raise ValueError("No execution with key "+key_or_id)
    else:
        exid = key_or_id
    allfiles = dict((minilims.fetch_file(x)['description'],x)
                    for x in minilims.search_files(source=('execution',exid)))
    if 'py:gdv_json' in allfiles:
        with open(minilims.path_to_file(allfiles['py:gdv_json'])) as q:
            gdv_json = pickle.load(q)
        if isinstance(url_or_dict, str):
            job.options['gdv_project'] = gdv_json
        else:
            job['options']['gdv_project'] = gdv_json
    with open(minilims.path_to_file(allfiles['py:file_names'])) as q:
        file_names = pickle.load(q)
    for gid, group in job_groups.iteritems():
        if not 'runs' in group:
            group = {'runs': group}
        processed[gid] = {}
        for rid,run in group['runs'].iteritems():
            bamfile = os.path.join(ex_root, unique_filename_in(ex_root))
            name = file_names[gid][rid]
            bam_id = allfiles['bam:'+name+'_filtered.bam']
            bam_bai_id = allfiles['bam:'+name+'_filtered.bam (BAM index)']
            minilims.export_file(bam_id,bamfile)
            minilims.export_file(bam_bai_id,bamfile+".bai")
            stats_id = allfiles.get("py:"+name+"_filter_bamstat") or allfiles.get("py:"+name+"_full_bamstat")
            with open(minilims.path_to_file(stats_id)) as q:
                stats = pickle.load(q)
            pickle_thresh = allfiles["py:"+name+"_Poisson_threshold"]
            with open(minilims.path_to_file(pickle_thresh)) as q:
                p_thresh = pickle.load(q)
            wigfile = os.path.join(ex_root, unique_filename_in(ex_root))
            wig_ids = dict(((allfiles['sql:'+name+'_'+s+'.sql'],s),
                            wigfile+'_'+s+'.sql') for s in suffix)
            [minilims.export_file(x[0],s) for x,s in wig_ids.iteritems()]
            processed[gid][rid] = {'bam': bamfile,
                                   'stats': stats,
                                   'poisson_threshold': p_thresh,
                                   'libname': name,
                                   'wig': dict((x[1],s)
                                               for x,s in wig_ids.iteritems())}
    return (processed,job)


#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

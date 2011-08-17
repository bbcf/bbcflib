"""
=====================
Module: bbcflib.motif
=====================

No documentation
"""

# Built-in modules #
import re, os
from operator import add

# Internal modules #
from . import track, common

# Other modules #
from bein import unique_filename_in, program

################################################################################
@program
def meme( fasta, maxsize=10000000, args=None ):
    """Binding for the ``meme`` motif finder.
    """
    if args is None:
        args = []
    outname = unique_filename_in()
    call = ["meme", fasta, "-o", outname, "-dna", "-maxsize", str(maxsize)]+args
    return {"arguments": call, "return_value": outname}

def parse_meme_xml( ex, meme_file, chromosomes ):
    """ Parse meme xml file and convert to track """
    from xml.etree import ElementTree as ET
    tree = ET.parse(meme_file)
    ncol = {}
    allmatrices = {}
    for motif in tree.find('motifs').findall('motif'):
        mid = motif.attrib['id']
        ncol[mid] = 0
        allmatrices[mid] = unique_filename_in()
        with open(allmatrices[mid],'w') as mat_out:
            for parray in motif.find('probabilities')[0].findall('alphabet_array'):
                ncol[mid] += 1
                m = {'letter_A':0,'letter_C':0,'letter_G':0,'letter_T':0}
                for col in parray:
                    m[col.attrib['letter_id']] = float(col.text)
                mat_out.write("1\t%f\t%f\t%f\t%f\n" %(m['letter_A'],m['letter_C'],m['letter_G'],m['letter_T']))
    def _xmltree(chr,tree):
        seq_name = {}
        seq_chr = None
        for it in tree.getiterator():
            if it.tag == 'sequence':
                seq_name[it.attrib['id']] = it.attrib['name']
            if it.tag == 'scanned_sites':
                name = seq_name[it.attrib['sequence_id']]
                name,seq_chr,start,end = re.search(r'(.*)\|(.+):(\d+)-(\d+)',name).groups()
            if it.tag == 'scanned_site' and chr == seq_chr:
                start = int(start)+int(it.attrib['position'])-1
                end = str(start+ncol[it.attrib['motif_id']])
                start = str(start)
                strnd = it.attrib['strand'] == 'plus' and 1 or -1
                score = it.attrib['pvalue']
                yield (start,end,it.attrib['motif_id'],score,strnd)
    outsql = unique_filename_in()
    chrlist = dict((v['name'], {'length': v['length']}) for v in chromosomes.values())
    with track.new(outsql, format='sql', chrmeta=chrlist, datatype='qualitative') as t:
        for chr in chrlist.keys():
            t.write(chr,_xmltree(chr,tree),fields=['start','end','name','score','strand'])
    return {'sql':outsql,'matrices':allmatrices}
                

def parallel_meme( ex, genrep, chromosomes, data, name=None, meme_args=None, via='lsf' ):
    """Fetches sequences, then calls ``meme``
    on them and finally saves the results in the repository.
    """
    if meme_args is None:
        meme_args   = []
    if not(isinstance(data,list)):
        data = [data]
    if not(isinstance(name,list)):
        if name == None:
            name = '_'
        name = [name]
    futures = {}
    for i,n in enumerate(name):
        (fasta, size) = genrep.fasta_from_data( chromosomes, data[i], out=unique_filename_in() )
        futures[n] = meme.nonblocking( ex, fasta, maxsize=size*1.5, args=meme_args, via=via )
    all_res = {}
    for n,f in futures.iteritems():
        meme_out = f.wait()
        archive = common.compress( ex, meme_out, 'gz' )
        meme_res = parse_meme_xml( ex, os.path.join(meme_out, "meme.xml"),
                                   chromosomes )
        ex.add( os.path.join(meme_out, "meme.html"),
                description="html:"+n+"_meme.html" )
        ex.add( meme_res['sql'], description="sql:"+n+"_meme_sites.sql" )
        ex.add( archive, description="tar:"+n+"_meme.tgz" )
        for i,motif in enumerate(meme_res['matrices'].keys()):
            ex.add( meme_res['matrices'][motif], description="txt:"+n+"_meme_"+motif+".txt" )
            ex.add( os.path.join(meme_out, "logo"+str(i+1)+".png"), description="png:"+n+"_meme_"+motif+".png" )
        all_res[n] = meme_res
    return all_res

@program
def motif_scan( fasta, motif, background, threshold=0 ):
    """Binding for the ``S1K`` motif scanner.
    """
    call = ["S1K", motif, background, str(threshold), fasta]
    return {"arguments": call, "return_value": None}

def save_motif_profile( ex, motifs, background, genrep, chromosomes, data_path,
                        keep_max_only=False, description='', threshold=0, via='lsf' ):
    """Scan a set of motifs on a set of regions and saves the results as an sql file.
    The 'motifs' argument is a dictionary with keys motif names and values PWM with 'n' rows like:
    "1 p(A) p(C) p(G) p(T)" where the sum of the 'p's is 1 and the first column allows to skip a position with a '0'.
    """
    fasta, size = genrep.fasta_from_data( chromosomes, data_path )
    sqlout = unique_filename_in()
    if not(isinstance(motifs, dict)):
        raise ValueError("'Motifs' must be a dictionary with keys 'motif_names' and values the PWMs.")
    for name, pwm in motifs.iteritems():
        output = unique_filename_in()
        futures[name] = ( output, motif_scan.nonblocking( ex, fasta, pwm, background, threshold, stdout=output, via=via ) )
    chrlist = dict((v['name'], {'length': v['length']}) for v in chromosomes.values())
##############
    def _parse_s1k( file, chr ):
        if keep_max_only:
            maxscore = {}
        with open(file, 'r') as fin:
            for line in fin:
                row = line.split("\t")
                name,seq_chr,start,end = re.search(r'(.*)\|(.+):(\d+)-(\d+)',row[0]).groups()
                if seq_chr != chr: 
                    continue
                start = int(row[3])+int(start)-1
                name = row[1]
                end = start+len(name)
                score = float(row[2])
                strnd = row[4]=='+' and 1 or -1
                if keep_max_only:
                    if not(row[0] in maxscore) or score>maxscore[row[0]][3]:
                        maxscore[row[0]] = (start,end,name,score,strnd)
                else:
                    yield (start,end,name,score,strnd)
        if keep_max_only:
            [yield x for x in maxscore.values()]
##############
    with track.new( sqlout, format="sql", datatype="qualitative", chrmeta=chrlist ) as track_result:
#        track_result.attributes = {'source': 'S1K'}
        for name, future in futures.iteritems():
            _ = future[1].wait()
            for chr in track_result:
                track_result.write( chr, _parse_s1k(future[0],chr) )
    ex.add( sqlout, description = "sql:"+description+"_motif_scan.sql" )
    return sqlout

def false_discovery_rate_p_value(false_positive_list, true_positive_list, nb_false_positive_hypotesis = 1.0):
    """
    Returns false discovery rate
    """
    score_tp = 0
    score_fp = 0
    if len(true_positive_list) >0:
        score_tp = reduce(add, true_positive_list)* nb_false_positive_hypotesis
    if len(false_positive_list) >0:
        score_fp = reduce(add, false_positive_list)
    return score_fp / float(score_tp + score_fp)

def false_discovery_rate(false_positive, true_positive, alpha = 1, nb_false_positive_hypotesis = 1.0):
    """
    Returns false discovery rate
    """
    p_value     = 0
    keys        = []
    if isinstance(false_positive, dict):
        keys    +=  false_positive.keys()
    else:
        raise TypeError(u"Allowed type for false_positive: tuple, list, dict. Type subited: %s" %type(false_positive))
    if isinstance(true_positive, dict):
        keys    += true_positive.keys()
    else:
        raise TypeError(u"Allowed type for true_positive: tuple, list, dict. Type subited: %s" %type(true_positive))
    keys                = tuple(sorted(set(keys)))
    index               = len(keys) - 1
    false_positive_list = []
    true_positive_list  = []

    is_running = True
    while is_running:
        if index < 0:
            is_running = False
        else:
            for key in  keys[index:]:
                if key in false_positive:
                    false_positive_list.append(false_positive[key])
                if key in true_positive:
                    true_positive_list.append(true_positive[key])
            p_value = false_discovery_rate_p_value(
                                                    false_positive_list,
                                                    true_positive_list,
                                                    nb_false_positive_hypotesis = nb_false_positive_hypotesis
                                                  )
            if p_value == alpha:
                is_running = False
            elif p_value > alpha:
                is_running = False
                index = (index + 1 < len(true_positive)) and index + 1 or index
            else:
                index -= 1
    return keys[index]

def sqlite_to_false_discovery_rate  (
                                        ex, motifs, background, genrep, chromosomes,
                                        true_positive_data_path, false_positive_data_path,
                                        description = '',
                                        threshold = 0,    via = 'lsf', keep_max_only = False,
                                        alpha = 0.05, nb_false_positive_hypotesis = 1.0
                                    ):
    """
    sqls or beds take an array (list or tuple) like:
    sqls = (sql, sql_random)
    beds = (bed, bed_random)
    """
    if true_positive_data_path is None or false_positive_data_path is None:
        raise ValueError(u"Variables true_positive_data_path or false_positive_data_path is set to None!")
    true_positive_data_path     = os.path.normcase(os.path.expanduser(true_positive_data_path))
    false_positive_data_path    = os.path.normcase(os.path.expanduser(false_positive_data_path))
    true_positive_result        = None
    false_positive_result       = None

    # original
    true_positive_result    = save_motif_profile(
                                                    ex, motifs, background, genrep, chromosomes,
                                                    true_positive_data_path,
                                                    description     = description, threshold = threshold,
                                                    via             = via, keep_max_only = keep_max_only
                                                )
    # random
    false_positive_result   = save_motif_profile(
                                                    ex, motifs, background, genrep, chromosomes,
                                                    false_positive_data_path,
                                                    description     = description, threshold = threshold,
                                                    via             = via, keep_max_only = keep_max_only
                                                )

    fp_scores = None
    tp_scores = None
    with track.load(false_positive_result,  format = "sql") as data:
        fp_scores = data.get_scores_frequencies()
    with track.load(true_positive_result,  format = "sql") as data:
        tp_scores = data.get_scores_frequencies()
    fdr = false_discovery_rate  (
                                    fp_scores, tp_scores,
                                    alpha = alpha, nb_false_positive_hypotesis = nb_false_positive_hypotesis
                                )
    return true_positive_result, fdr

def parse_meme_html_output(ex, meme_file, fasta, chromosomes):
    """Legacy signature"""
    meme_file = re.sub(r'\.html','.xml',meme_file)
    return parse_meme_xml( ex, meme_file, chromosomes )
    
def add_meme_files( ex, genrep, chromosomes, description = '',
                    bed = None, sql = None, meme_args = None, via = 'lsf' ):
    """Legacy signature"""
    data = sql
    if sql == None:
        data = bed
    return parallel_meme( ex, genrep, chromosomes, data=data, name=description, 
                          meme_args=meme_args, via=via )
    
#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

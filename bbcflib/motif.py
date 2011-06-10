"""
===============
bbcflib.motif
===============
"""
import os
from operator import add
import sqlite3
from bein import *
from bein.util import *
from bbcflib import common, genrep

@program
def meme( fasta, maxsize=10000000, args=[] ):
    """Binding for the ``meme`` motif finder.
    """
    outname = unique_filename_in()
    call = ["meme",fasta,"-o",outname,"-dna","-maxsize",str(maxsize),"-revcomp"]+args
    return {"arguments": call, "return_value": outname}

def add_meme_files( ex, genrep, chromosomes, description='',
                    bed=None, sql=None, meme_args=[], via='lsf' ):
    """Fetches sequences, then calls ``meme``
    on them and finally saves the results in the repository.
    """
    fasta,size = genrep.fasta_from_bed( chromosomes, out=unique_filename_in(),
                                        bed=bed, sql=sql )
    meme_out = meme.nonblocking( ex, fasta, maxsize=size*1.5, args=meme_args, via=via ).wait()
    html = os.path.join(meme_out,"meme.html")
    ex.add( html, description="html:"+description+"meme.html" )
    return meme_out

@program
def motif_scan( fasta, motif, background, description='', threshold=0 ):
    """Binding for the ``S1K`` motif scanner.
    """
    call = ["S1K",motif,background,str(threshold),fasta]
    return {"arguments": call, "return_value": None}

def save_motif_profile( ex, motifs, background, genrep, chromosomes,
                        description='', sql=None, bed=None, threshold=0, via='lsf' ):
    """Scan a set of motifs on a set of chromosomes and saves the results as an sql file.
    The 'motifs' argument is a dictionary with keys motif names and values PWM with 'n' rows like:
    "1 p(A) p(C) p(G) p(T)" where the sum of the 'p's is 1 and the first column allows to skip a position with a '0'.
    """
    fasta,size = genrep.fasta_from_bed( chromosomes, out=unique_filename_in(), bed=bed, sql=sql )
    futures= {}
    if not(isinstance(motifs,dict)):
        raise ValueError("'Motifs' must be a dictionary with keys 'motif_names' and values the PWMs.")
    for name,pwm in motifs.iteritems():
        output = unique_filename_in()
        futures[name] = (output,
                         motif_scan.nonblocking( ex, fasta, pwm, background, threshold,
                                                 via=via, stdout=output ))
    regions = {}
    if bed != None:
        chr_ids = dict((cn['name'],c[0]) for c,cn in chromosomes.iteritems())
        with open(bed,"r") as f:
            for l in f:
                row = l.rstrip('\n').split('\t')
                name = re.search(r'(\S+)\s*',row[3]).groups()[0]
                regions[name] = (row[0],int(row[1]))
    elif sql != None:
        connection = sqlite3.connect( sql )
        for v in chromosomes.values():
            for row in connection.execute('select start,name from "'+v['name']+'"'):
                name = re.search(r'(\S+)\s*',row[1]).groups()[0]
                regions[name] = (v['name'],int(row[0]))
        connection.close()
    else:
        raise TypeError("save_motif_profile requires either a 'sqlite' or a 'bed' file.")
    sqlout = common.create_sql_track( unique_filename_in(), chromosomes, datatype="qualitative" )
    connection = sqlite3.connect( sqlout )
    for name,f in futures.iteritems():
        vals = []
        _ = f[1].wait()
        cur_chr = ''
        with open(f[0],'r') as f:
            for l in f:
                s = l.rstrip('\n').split('\t')
                reg = regions[s[0]]
                start = reg[1]+int(s[3])
                vals.append((start-1,start+len(s[1]),float(s[2]),s[4],name+":"+s[1]))
                if len(cur_chr)>0 and reg[0] != cur_chr:
                    connection.executemany('insert into "'+cur_chr+'" (start,end,score,strand,name) values (?,?,?,?,?)',vals)
                    connection.commit()
                    vals = []
                cur_chr = reg[0]
        if len(vals)>0:
            connection.executemany('insert into "'+cur_chr+'" (start,end,score,strand,name) values (?,?,?,?,?)',vals)
            connection.commit()
    connection.close()
    ex.add( sqlout, description="sql:"+description+"motif_scan.sql" )
    return sqlout

def false_discovery_rate_p_value(false_positive_list, true_positive_list, index):
    """
    Return false discovery rate
    """
    tp = 0
    fn = 0
    if index < len(true_positive_list):
        tp = reduce(add, true_positive_list[index:])
    if index < len(false_positive_list):
        fp = reduce(add, false_positive_list[index:])
    return fp / float(tp + fp)

def false_discovery_rate(false_positive, true_positive, alpha=1, factor_fp=1.0, factor_tp=1.0):
    """
    Return false discovery rate
    """
    if isinstance(false_positive, tuple) or isinstance(false_positive, list):
        if isinstance(false_positive[0], tuple) or isinstance(false_positive[0], list):
            false_positive_list = [len(i)* (1/float(factor_fp)) for i in false_positive]
        else:
            false_positive_list = [i* (1/float(factor_fp)) for i in false_positive]
    elif isinstance(false_positive, dict):
        false_positive_list = [ false_positive[f] * (1/float(factor_fp)) for f in false_positive ]
    else:
        raise TypeError(u"Allowed type for false_positive: tuple, list, dict. Type subited: %s" %type(false_positive))
    if isinstance(true_positive, tuple) or isinstance(true_positive, list):
        if isinstance(true_positive[0], tuple) or isinstance(true_positive[0], list):
            true_positive_list = [len(i)* (1/float(factor_tp)) for i in true_positive]
        else:
            true_positive_list = [i* (1/float(factor_tp)) for i in true_positive]
    elif isinstance(true_positive, dict):
        true_positive_list = [ true_positive[f]* (1/float(factor_tp)) for f in true_positive ]
    else:
        raise TypeError(u"Allowed type for true_positive: tuple, list, dict. Type subited: %s" %type(true_positive))
    index       = min( len(true_positive_list), len(false_positive_list) )
    index       -=  1
    p_value     = 0
    result      = 0

    isRunning = True
    while isRunning:
        if index < 0:
            isRunning = False
        else:
            p_value = false_discovery_rate_p_value(false_positive_list, true_positive_list, index)
            if 1 - p_value == alpha:
                isRunning = False
            elif 1 - p_value > alpha:
                isRunning = False
                index       = (index + 1 < len(true_positive_list)) and index + 1 or index
            else:
                index -= 1
    if isinstance(true_positive, tuple) or isinstance(true_positive, list):
        if isinstance(true_positive[0], tuple) or isinstance(true_positive[0], list):
            result = true_positive[index]
        else:
            result =  index
    elif isinstance(true_positive, dict):
        result = true_positive.keys()[index]
    return result

def sqlite_to_false_discovery_rate  (
                                        ex, motifs, background, genrep, chromosomes,
                                        description='',
                                        sqls=None,      beds=None,
                                        threshold=0,    via='lsf',
                                        alpha=1, factor_fp=1.0, factor_tp=1.0
                                    ):
    """
    sqls or beds take an array (list or tuple) like:
    sqls = (sql, sql_random)
    beds = (bed, bed_random)
    """
    if sqls is None and beds is None:
        raise ValueError(u"Variables sqls and beds is set to None! Set one of these variables")
    elif sqls is not None and beds is not None:
        raise ValueError(u"Variables sqls and beds is not None! Set only one of these variables, sqls or beds")
    background  = os.path.expanduser(background)
    if not os.path.isabs(background):
        background = os.path.normcase("../"+background)
    true_positive_result    = None
    false_positive_result   = None

    if sqls is not None:
        new_sql0 = os.path.expanduser(sqls[0])
        new_sql1 = os.path.expanduser(sqls[1])
        if not os.path.isabs(new_sql0):
            new_sql0 = os.path.normcase("../"+new_sql0)
        if not os.path.isabs(new_sql1):
            new_sql1 = os.path.normcase("../"+new_sql1)
        # original
        true_positive_result    = save_motif_profile(
                                                        ex, motifs, background, genrep, chromosomes,
                                                        description='',
                                                        sql=new_sql0, bed=None, threshold=threshold, via=via
                                                    )
        # random
        false_positive_result   = save_motif_profile(
                                                        ex, motifs, background, genrep, chromosomes,
                                                        description='',
                                                        sql=new_sql1, bed=None, threshold=threshold, via=via
                                                    )
    else:
        new_bed0 = os.path.expanduser(beds[0])
        new_bed1 = os.path.expanduser(beds[1])
        if not os.path.isabs(new_bed0):
            new_bed0 = os.path.normcase("../"+new_bed0)
        if not os.path.isabs(new_bed1):
            new_bed1 = os.path.normcase("../"+new_bed1)
        # original
        true_positive_result    = save_motif_profile(
                                                        ex, motifs, background, genrep, chromosomes,
                                                        description=description,
                                                        sql=None, bed=new_bed0, threshold=threshold, via=via
                                                    )
        # random
        false_positive_result   = save_motif_profile(
                                                        ex, motifs, background, genrep, chromosomes,
                                                        description=description,
                                                        sql=None, bed=new_bed1, threshold=threshold, via=via
                                                    )

    fp_scores = get_score(ex.working_directory+"/"+false_positive_result)
    tp_scores = get_score(ex.working_directory+"/"+true_positive_result)

    fdr = false_discovery_rate  (
                                    fp_scores, tp_scores,
                                    alpha=alpha, factor_fp=factor_fp, factor_tp=factor_tp
                                )
    return fdr

def get_score( sql ):
    connection          = sqlite3.connect(sql)
    cursor              = connection.cursor()
    cursor.execute(u"SELECT name FROM chrNames;")
    chromosomes_names   = [ chromosome[0] for chromosome in cursor ]
    scores              = {}
    for chromosome_name in chromosomes_names:
        cursor.execute(u"SELECT count (*), score FROM '"+chromosome_name+u"' GROUP BY score;")
        for result in cursor:
            if result[1] not in scores:
                scores[ result[1] ] = result[0]
            else:
                scores[ result[1] ] += result[0]
    connection.close()
    return scores

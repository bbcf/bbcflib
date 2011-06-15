"""
===============
bbcflib.motif
===============
"""
from operator import add
import sqlite3, re, os
from bein import *
from bein.util import *
from bbcflib import common, genrep
from BeautifulSoup import BeautifulSoup

@program
def meme( fasta, maxsize=10000000, args=[] ):
    """Binding for the ``meme`` motif finder.
    """
    outname = unique_filename_in()
    call = ["meme",fasta,"-o",outname,"-dna","-maxsize",str(maxsize),"-revcomp"]+args
    return {"arguments": call, "return_value": outname}

def parse_meme_html_output(ex, file_path):
    soup = None
    with open(file_path, "r") as f:
        soup = BeautifulSoup("".join(f.readlines()))
    html_matrices = soup.findAll('input',  id=re.compile('^pspm\d+'))
    #~ raw_matrices = [item.attrs[3][1].lstrip("\n").rstrip("\n\n") for item in html_matrices]
    for item in html_matrices:
        matrix_file = unique_filename_in()
        raws = item.attrs[3][1].lstrip("\n").rstrip("\n\n").split("\n")
        raws[0] = ">" + raws[0]
        raws[1:] = [ "1 "+raws[i] for i in xrange(1,len(raws))]
        with open(matrix_file, "w") as f:
            f.write(",".join(raws))
        ex.add( matrix_file, description="matrix:"+raws[0] )

def add_meme_files( ex, genrep, chromosomes, description='',
                    bed=None, sql=None, meme_args=[], via='lsf' ):
    """Fetches sequences, then calls ``meme``
    on them and finally saves the results in the repository.
    """
    if bed not None:
        bed = os.path.expanduser(bed)
        if not os.path.isabs(bed):
            bed = os.path.normcase("../"+bed)
    elif sql not None:
        sql = os.path.expanduser(sql)
        if not os.path.isabs(sql):
            sql = os.path.normcase("../"+sql)
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
                        description='', sql=None, bed=None, threshold=0, via='lsf', keep_max_only=False ):
    """Scan a set of motifs on a set of chromosomes and saves the results as an sql file.
    The 'motifs' argument is a dictionary with keys motif names and values PWM with 'n' rows like:
    "1 p(A) p(C) p(G) p(T)" where the sum of the 'p's is 1 and the first column allows to skip a position with a '0'.
    """
    if bed not None:
        bed = os.path.expanduser(bed)
        if not os.path.isabs(bed):
            bed = os.path.normcase("../"+bed)
    elif sql not None:
        sql = os.path.expanduser(sql)
        if not os.path.isabs(sql):
            sql = os.path.normcase("../"+sql)
    if motifs not None:
        motifs = os.path.expanduser(motifs)
        if not os.path.isabs(motifs):
            motifs = os.path.normcase("../"+motifs)
    if background not None:
        background = os.path.expanduser(background)
        if not os.path.isabs(background):
            background = os.path.normcase("../"+background)
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
        vals        = []
        _           = f[1].wait()
        cur_chr     = ''
        index       = 0
        previous_feature=""
        with open(f[0],'r') as f:
            for l in f:
                s       = l.rstrip('\n').split('\t')
                reg     = regions[s[0]]
                start   = reg[1]+int(s[3])
                if cur_chr != '' and reg[0] != cur_chr:
                    connection.executemany('insert into "'+cur_chr+'" (start,end,score,strand,name) values (?,?,?,?,?)',vals)
                    connection.commit()
                    index           = 0
                    vals            = [(start-1,start+len(s[1]),float(s[2]),s[4],name+":"+s[1])]
                    previous_feature= s[0]
                cur_chr = reg[0]
                if not keep_max_only:
                    vals.append((start-1,start+len(s[1]),float(s[2]),s[4],name+":"+s[1]))
                else:
                    if previous_feature == "":
                        previous_feature= s[0]
                        vals.append((start-1,start+len(s[1]),float(s[2]),s[4],name+":"+s[1]))
                    elif previous_feature == s[0] and vals[index][2] < float(s[2]):
                            vals[index] = (start-1,start+len(s[1]),float(s[2]),s[4],name+":"+s[1])
                    elif previous_feature != s[0]:
                        previous_feature= s[0]
                        index           += 1
                        vals.append((start-1,start+len(s[1]),float(s[2]),s[4],name+":"+s[1]))
        if len(vals)>0:
            connection.executemany('insert into "'+cur_chr+'" (start,end,score,strand,name) values (?,?,?,?,?)',vals)
            connection.commit()
    connection.close()
    ex.add( sqlout, description="sql:"+description+"motif_scan.sql" )
    return sqlout

def false_discovery_rate_p_value(false_positive_list, true_positive_list, nb_false_positive_hypotesis=1.0):
    """
    Return false discovery rate
    """
    tp = 0
    fp = 0
    if len(true_positive_list) >0:
        tp = reduce(add, true_positive_list)* nb_false_positive_hypotesis
    if len(false_positive_list) >0:
        fp = reduce(add, false_positive_list)
    return fp / float(tp + fp)

def false_discovery_rate(false_positive, true_positive, alpha=1, nb_false_positive_hypotesis=1.0):
    """
    Return false discovery rate
    """
    result      = 0
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

    isRunning = True
    while isRunning:
        if index < 0:
            isRunning = False
        else:
            for x in  keys[index:]:
                if x in false_positive:
                    false_positive_list.append(false_positive[x])
                if x in true_positive:
                    true_positive_list.append(true_positive[x])
            p_value = false_discovery_rate_p_value(false_positive_list, true_positive_list, nb_false_positive_hypotesis=nb_false_positive_hypotesis)
            if p_value == alpha:
                isRunning = False
            elif p_value > alpha:
                isRunning = False
                index = (index + 1 < len(true_positive)) and index + 1 or index
            else:
                index -= 1
    return keys[index]

def sqlite_to_false_discovery_rate  (
                                        ex, motifs, background, genrep, chromosomes,
                                        description='',
                                        sqls=None,      beds=None,
                                        threshold=0,    via='lsf', keep_max_only=False,
                                        alpha=0.05, nb_false_positive_hypotesis=1.0
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
                                                        description     = description,
                                                        sql             = new_sql0, bed=None, threshold=threshold,
                                                        via             = via, keep_max_only = keep_max_only
                                                    )
        # random
        false_positive_result   = save_motif_profile(
                                                        ex, motifs, background, genrep, chromosomes,
                                                        description     = description,
                                                        sql             = new_sql1, bed=None, threshold=threshold,
                                                        via             = via, keep_max_only=keep_max_only
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
                                                        description     = description,
                                                        sql             = None, bed=new_bed0, threshold=threshold,
                                                        via             = via, keep_max_only=keep_max_only
                                                    )
        # random
        false_positive_result   = save_motif_profile(
                                                        ex, motifs, background, genrep, chromosomes,
                                                        description     = description,
                                                        sql             = None, bed=new_bed1, threshold=threshold,
                                                        via             = via, keep_max_only=keep_max_only
                                                    )

    fp_scores = get_score(ex.working_directory+"/"+false_positive_result)
    tp_scores = get_score(ex.working_directory+"/"+true_positive_result)

    return false_discovery_rate (
                                    fp_scores, tp_scores,
                                    alpha=alpha, nb_false_positive_hypotesis=nb_false_positive_hypotesis
                                )

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

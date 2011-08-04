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
from BeautifulSoup import BeautifulSoup

################################################################################
@program
def meme( fasta, maxsize = 10000000, args = None ):
    """Binding for the ``meme`` motif finder.
    """
    if args is None:
        args = []
    outname = unique_filename_in()
    call = ["meme", fasta, "-o", outname, "-dna", "-maxsize", str(maxsize), "-revcomp"]+args
    return {"arguments": call, "return_value": outname}

def parse_meme_html_output(ex, meme_file, fasta, chromosomes):
    """ Parse meme html file and convert to track """
    soup    = None
    with open(meme_file, "r") as meme_in:
        soup = BeautifulSoup("".join(meme_in.readlines()))
    pattern         = re.compile("^pspm(\d+)")
    html_matrices   = soup.html.body.form.findAll({"input" : True},  attrs = {"id" : pattern})

    for item in html_matrices:
        chrom_dict    = {}
        start               = 0
        end                 = 0
        result_id           = pattern.match(item["id"]).group(1)
        matrix_file         = unique_filename_in()
        sqlout              = unique_filename_in()
        # write martix
        raws                = item["value"].lstrip("\n").rstrip("\n\n").split("\n")
        raws[0]             = ">" + raws[0]
        raws[1:]            = [ "1 "+raws[i] for i in xrange(1, len(raws))]
        motif_length        = len(raws[1:])
        with open(matrix_file, "w") as file_out:
            file_out.write("\n".join(raws))
        # write motif position
        motif_position      = soup.html.body.form.find  (
                                                            {"table" : True},
                                                            attrs = {"id" : re.compile("tbl_sites_%s" %(result_id))}
                                                        )
        for tds in motif_position.tbody.findAll("tr"):
            is_searching_header   = True
            chrom_name          = ""
            strand_name         = tds.find({"td" : True}, { "class" : "strand_name"}).text
            strand_side         = tds.find({"td" : True}, { "class" : "strand_side"}).text
            strand_start        = tds.find({"td" : True}, { "class" : "strand_start"}).text
            strand_pvalue       = tds.find({"td" : True}, { "class" : "strand_pvalue"}).text
            # search fasta header
            with open(fasta, "r") as fasta_in:
                while is_searching_header:
                    line = fasta_in.readline()
                    if line is None:
                        is_searching_header = False
                    elif line.startswith(">"):
                        header = line[1:].split(" ")
                        if header[0] == strand_name:
                            is_searching_header               = False
                            chrom_name, feature_position    = header[1].split(":")
                            feature_start, feature_end      = feature_position.split("-")
            if strand_side == "+":
                start   = int(feature_start) + int(strand_start)
                end     = start + motif_length
            elif strand_side == "-":
                start   = int(feature_end) - (int(strand_start) + motif_length + 1)
                end     = int(feature_end) - int(strand_start)
            else:
                raise ValueError("Unknow strand side value: %s!" %(strand_side))
            if chrom_name in chrom_dict:
                chrom_dict[chrom_name] +=   [
                                                (int(start), int(end), strand_name, float(strand_pvalue), strand_side )
                                            ]
            else:
                chrom_dict[chrom_name] =    [
                                                (int(start), int(end), strand_name, float(strand_pvalue), strand_side)
                                            ]
        with track.new(sqlout,  format = "sql", datatype = "qualitative") as track_result:
            keys            = chromosomes.keys()
            chrom_used = {}
            track_result.attributes = {'datatype': 'qualitative', 'source': 'meme_file'}
            for chromosome in chrom_dict:
                is_searching_chromosome   = True
                index                   = 0
                while is_searching_chromosome:
                    if index >= len(keys):
                        raise ValueError("Chromosomes named: %s not found in select assembly!"%(chromosome))
                    elif chromosomes[keys[index]]["name"] == chromosome:
                        is_searching_chromosome   = False
                        chrom_used[chromosomes[keys[index]]["name"]] = {"length" : chromosomes[keys[index]]["length"] }
                    else:
                        index += 1
                track_result.write(chromosome, chrom_dict[chromosome], fields = track.Track.qualitative_fields)
            track_result.chrmeta = chrom_used
        #files.append((matrix_file, sqlout))
        ex.add( matrix_file, description = "matrix:"+item["name"] )
        ex.add( sqlout,  description = "sql:"+sqlout)
    #return files

def add_meme_files( ex, genrep, chromosomes, description = '',
                    bed = None, sql = None, meme_args = None, via = 'lsf' ):
    """Fetches sequences, then calls ``meme``
    on them and finally saves the results in the repository.
    """
    if meme_args is None:
        meme_args   = []
    if bed is not None:
        bed = os.path.expanduser(bed)
        if not os.path.isabs(bed):
            bed = os.path.normcase(bed)
    elif sql is not None:
        sql = os.path.expanduser(sql)
        if not os.path.isabs(sql):
            sql = os.path.normcase(sql)
    fasta, size = genrep.fasta_from_data( chromosomes, out = unique_filename_in(),
                                        bed = bed, sql = sql )
    meme_out    = meme.nonblocking( ex, fasta, maxsize = size*1.5, args = meme_args, via = via ).wait()
    html        = os.path.join(meme_out, "meme.html")
    parse_meme_html_output(ex, meme_out+"/meme.html", fasta, chromosomes)
    archive     = common.compress(ex, meme_out)
    ex.add( html, description = "html:"+description+"meme.html" )
    ex.add( archive, description = "archive:"+description )
    return meme_out

@program
def motif_scan( fasta, motif, background, threshold = 0 ):
    """Binding for the ``S1K`` motif scanner.
    """
    call = ["S1K", motif, background, str(threshold), fasta]
    return {"arguments": call, "return_value": None}

def save_motif_profile( ex, motifs, background, genrep, chromosomes, data_path,
                        description = '', threshold = 0, via = 'lsf', keep_max_only = False ):
    """Scan a set of motifs on a set of chromosomes and saves the results as an sql file.
    The 'motifs' argument is a dictionary with keys motif names and values PWM with 'n' rows like:
    "1 p(A) p(C) p(G) p(T)" where the sum of the 'p's is 1 and the first column allows to skip a position with a '0'.
    """
    if motifs is not None:
        for name in motifs:
            motifs[name] = os.path.normcase(os.path.expanduser(motifs[name]))
    if background is not None:
        background = os.path.normcase(os.path.expanduser(background))
    fasta, size      = genrep.fasta_from_data( chromosomes, data_path, out = unique_filename_in() )
    sqlout          = unique_filename_in()
    futures         = {}
    regions         = {}
    chromosomes_set = set()
    chrom_used = {}
    keys            = chromosomes.keys()
    if not(isinstance(motifs, dict)):
        raise ValueError("'Motifs' must be a dictionary with keys 'motif_names' and values the PWMs.")

    for name, pwm in motifs.iteritems():
        output = unique_filename_in()
        futures[name] = (
                            output,
                            motif_scan.nonblocking( ex, fasta, pwm, background, threshold, via = via, stdout = output )
                        )

    with track.load(data_path, chrmeta = chromosomes) as data:
        for value in chromosomes.values():
            for row in data.read(selection = value['name'], fields = ["start", "name"]):
                name = re.search(r'(\S+)\s*', row[1]).groups()[0]
                regions[name] = (value['name'], int(row[0]))

    with track.new(sqlout,  format = "sql", datatype = "qualitative", ) as track_result:
        track_result.attributes = {'source': 'S1K'}

    for name, future in futures.iteritems():
        vals            = []
        _               = future[1].wait()
        cur_chr         = ''
        index           = 0
        previous_feature = ""
        with open(future[0], 'r') as f_in:
            for line in f_in:
                string  = line.rstrip('\n').split('\t')
                reg     = regions[string[0]]
                start   = reg[1]+int(string[3])
                if cur_chr != '' and reg[0] != cur_chr:
                    chromosomes_set.add(cur_chr)
                    with track.load(sqlout,  format = "sql") as track_result:
                        track_result.write(cur_chr, vals)
                    index   = 0
                    vals    =   [
                                    (
                                        start-1,
                                        start+len(string[1]),
                                        float(string[2]),
                                        name,
                                        string[4],
                                        "sequence = "+string[1]
                                    )
                                ]
                    previous_feature = string[0]
                cur_chr = reg[0]
                if not keep_max_only:
                    vals.append (
                                    (
                                        start-1,
                                        start+len(string[1]),
                                        float(string[2]),
                                        name,
                                        string[4],
                                        "sequence = "+string[1]
                                    )
                                )
                else:
                    if previous_feature == "":
                        previous_feature = string[0]
                        vals.append(
                                        (
                                            start-1,
                                            start+len(string[1]),
                                            float(string[2]),
                                            name, string[4],
                                            "sequence = "+string[1]
                                        )
                                    )
                    elif previous_feature == string[0] and vals[index][2] < float(string[2]):
                        vals[index] = (
                                        start-1,
                                        start+len(string[1]),
                                        float(string[2]),
                                        name,
                                        string[4],
                                        "sequence = "+string[1]
                                      )
                    elif previous_feature != string[0]:
                        previous_feature = string[0]
                        index           += 1
                        vals.append(
                                    (
                                        start-1,
                                        start+len(string[1]),
                                        float(string[2]),
                                        name,
                                        string[4],
                                        "sequence = "+string[1]
                                    )
                                   )
        if len(vals) > 0:
            with track.load(sqlout,  format = "sql") as track_result:
                track_result.write(cur_chr, vals, fields = track.Track.qualitative_fields)

    for chromosome in chromosomes_set:
        is_searching_chromosome   = True
        i                       = 0
        while is_searching_chromosome:
            if i >= len(keys):
                raise ValueError("Chromosomes named: %s not found in selected assembly!"%(chromosome))
            elif chromosomes[keys[i]]["name"] == chromosome:
                is_searching_chromosome                           = False
                chrom_used[ chromosomes[keys[i]]["name"] ] = {"length" : chromosomes[keys[i]]["length"] }
            else:
                i += 1

    with track.load(sqlout,  format = "sql") as track_result:
        track_result.chrmeta = chrom_used

    ex.add( sqlout, description = "sql:"+description+"motif_scan.sql" )
    return sqlout

def false_discovery_rate_p_value(false_positive_list, true_positive_list, nb_false_positive_hypotesis = 1.0):
    """
    Return false discovery rate
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
    Return false discovery rate
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

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

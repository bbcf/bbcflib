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
def meme( fasta, maxsize=10000000, args=[] ):
    """Binding for the ``meme`` motif finder.
    """
    outname = unique_filename_in()
    call = ["meme",fasta,"-o",outname,"-dna","-maxsize",str(maxsize),"-revcomp"]+args
    return {"arguments": call, "return_value": outname}

def parse_meme_html_output(ex, meme, fasta, chromosomes):
    soup    = None
    files   = []
    with open(meme, "r") as meme_in:
        soup = BeautifulSoup("".join(meme_in.readlines()))
    pattern         = re.compile("^pspm(\d+)")
    html_matrices   = soup.html.body.form.findAll({"input" : True},  attrs={"id" : pattern})

    for item in html_matrices:
        dict_chromosomes    = {}
        start               = 0
        end                 = 0
        result_id           = pattern.match(item["id"]).group(1)
        matrix_file         = unique_filename_in()
        sqlout              = unique_filename_in()
        # write martix
        raws                = item["value"].lstrip("\n").rstrip("\n\n").split("\n")
        raws[0]             = ">" + raws[0]
        raws[1:]            = [ "1 "+raws[i] for i in xrange(1,len(raws))]
        motif_length        = len(raws[1:])
        with open(matrix_file, "w") as f:
            f.write("\n".join(raws))
        # write motif position
        motif_position  = soup.html.body.form.find( {"table" : True}, attrs={"id" : re.compile("tbl_sites_%s" %(result_id))} )
        #~ tds             = motif_position.findChildren({"td" : True}, limit=4) # strand_name, strand_side, strand_start, strand_pvalue
        for tds in motif_position.tbody.findAll("tr"):
            isSearchingHeader   = True
            chromosome_name     = ""
            strand_name         = tds.find({"td" : True}, { "class" : "strand_name"}).text
            strand_side         = tds.find({"td" : True}, { "class" : "strand_side"}).text
            strand_start        = tds.find({"td" : True}, { "class" : "strand_start"}).text
            strand_pvalue       = tds.find({"td" : True}, { "class" : "strand_pvalue"}).text
            # search fasta header
            with open(fasta, "r") as fasta_in:
                while isSearchingHeader:
                    line = fasta_in.readline()
                    if line is None:
                        isSearchingHeader = False
                    elif line.startswith(">"):
                        header = line[1:].split(" ")
                        if header[0] == strand_name:
                            isSearchingHeader                   = False
                            chromosome_name, feature_position   = header[1].split(":")
                            feature_start, feature_end          = feature_position.split("-")
            if strand_side == "+":
                start   = int(feature_start) + int(strand_start)
                end     = start + motif_length
            elif strand_side == "-":
                start   = int(feature_end) - (int(strand_start) + motif_length + 1)
                end     = int(feature_end) - int(strand_start)
            else:
                raise ValueError("Unknow strand side value: %s!" %(strand_side))
            if chromosome_name in dict_chromosomes:
                dict_chromosomes[chromosome_name]+= [
                                                        (int(start), int(end), strand_name, float(strand_pvalue), strand_side )
                                                    ]
            else:
                dict_chromosomes[chromosome_name]=  [
                                                        (int(start), int(end), strand_name, float(strand_pvalue), strand_side)
                                                    ]
        with track.new(sqlout,  format="sql", datatype="qualitative") as t:
            keys            = chromosomes.keys()
            chomosomes_used = []
            t.attributes = {'datatype': 'qualitative', 'source': 'meme'}
            for chromosome in dict_chromosomes:
                isSearchingChromosome   = True
                index                   = 0
                while isSearchingChromosome:
                    if index >= len(keys):
                        raise ValueError("Chromosomes named: %s not found in select assembly!"%(chromosome))
                    elif chromosomes[keys[index]]["name"] == chromosome:
                        isSearchingChromosome   = False
                        chomosomes_used.append( chromosomes[keys[index]] )
                    else:
                        index +=1
                t.write(chromosome, dict_chromosomes[chromosome], fields=track.Track.qualitative_fields)
            t.chrmeta = chomosomes_used
        #files.append((matrix_file, sqlout))
        ex.add( matrix_file, description="matrix:"+item["name"] )
        ex.add( sqlout,  description="sql:"+sqlout)
    #return files

def add_meme_files( ex, genrep, chromosomes, description='',
                    bed=None, sql=None, meme_args=[], via='lsf' ):
    """Fetches sequences, then calls ``meme``
    on them and finally saves the results in the repository.
    """
    if bed is not None:
        bed = os.path.expanduser(bed)
        if not os.path.isabs(bed):
            bed = os.path.normcase(bed)
    elif sql is not None:
        sql = os.path.expanduser(sql)
        if not os.path.isabs(sql):
            sql = os.path.normcase(sql)
    fasta,size = genrep.fasta_from_data( chromosomes, out=unique_filename_in(),
                                        bed=bed, sql=sql )
    meme_out    = meme.nonblocking( ex, fasta, maxsize=size*1.5, args=meme_args, via=via ).wait()
    html        = os.path.join(meme_out,"meme.html")
    files       = parse_meme_html_output(ex, meme_out+"/meme.html", fasta, chromosomes) # return or not?
    archive     = common.compress(ex, meme_out)
    ex.add( html, description="html:"+description+"meme.html" )
    ex.add( archive, description="archive:"+description )
    return meme_out

@program
def motif_scan( fasta, motif, background, description='', threshold=0 ):
    """Binding for the ``S1K`` motif scanner.
    """
    call = ["S1K",motif,background,str(threshold),fasta]
    return {"arguments": call, "return_value": None}

def save_motif_profile( ex, motifs, background, genrep, chromosomes, data_path,
                        description='', threshold=0, via='lsf', keep_max_only=False ):
    """Scan a set of motifs on a set of chromosomes and saves the results as an sql file.
    The 'motifs' argument is a dictionary with keys motif names and values PWM with 'n' rows like:
    "1 p(A) p(C) p(G) p(T)" where the sum of the 'p's is 1 and the first column allows to skip a position with a '0'.
    """
    if motifs is not None:
        for name in motifs:
            motifs[name] = os.path.normcase(os.path.expanduser(motifs[name]))
    if background is not None:
        background = os.path.normcase(os.path.expanduser(background))
    fasta,size      = genrep.fasta_from_data( chromosomes, data_path, out=unique_filename_in() )
    sqlout          = unique_filename_in()
    futures         = {}
    regions         = {}
    chromosomes_set = set()
    chomosomes_used = []
    keys            = chromosomes.keys()
    if not(isinstance(motifs,dict)):
        raise ValueError("'Motifs' must be a dictionary with keys 'motif_names' and values the PWMs.")

    for name,pwm in motifs.iteritems():
        output = unique_filename_in()
        futures[name] = (output,
                         motif_scan.nonblocking( ex, fasta, pwm, background, threshold,
                                                 via=via, stdout=output ))

    with track.load(data_path, chrmeta=chromosomes) as t:
        for v in chromosomes.values():
                for row in t.read(selection=v['name'], fields=["start","name"]):
                    name = re.search(r'(\S+)\s*',row[1]).groups()[0]
                    regions[name] = (v['name'],int(row[0]))

    with track.new(sqlout,  format="sql", datatype="qualitative",) as t:
        t.attributes= {'source': 'S1K'}

    for name,f in futures.iteritems():
        vals            = []
        _               = f[1].wait()
        cur_chr         = ''
        index           = 0
        previous_feature=""
        with open(f[0],'r') as f_in:
            for l in f_in:
                s       = l.rstrip('\n').split('\t')
                reg     = regions[s[0]]
                start   = reg[1]+int(s[3])
                if cur_chr != '' and reg[0] != cur_chr:
                    chromosomes_set.add(cur_chr)
                    with track.load(sqlout,  format="sql") as t:
                        t.write(cur_chr, vals)
                    index           = 0
                    vals            = [(start-1,start+len(s[1]),name+":"+s[1],float(s[2]),s[4])]
                    previous_feature= s[0]
                cur_chr = reg[0]
                if not keep_max_only:
                    vals.append((start-1,start+len(s[1]),name+":"+s[1],float(s[2]),s[4]))
                else:
                    if previous_feature == "":
                        previous_feature= s[0]
                        vals.append((start-1,start+len(s[1]),float(s[2]),s[4],name+":"+s[1]))
                    elif previous_feature == s[0] and vals[index][2] < float(s[2]):
                            vals[index] = (start-1,start+len(s[1]),name+":"+s[1],float(s[2]),s[4])
                    elif previous_feature != s[0]:
                        previous_feature= s[0]
                        index           += 1
                        vals.append((start-1,start+len(s[1]),name+":"+s[1],float(s[2]),s[4]))
        if len(vals)>0:
            with track.load(sqlout,  format="sql") as t:
                t.write(cur_chr, vals, fields=track.Track.qualitative_fields)

    for chromosome in chromosomes_set:
        isSearchingChromosome   = True
        i                       = 0
        while isSearchingChromosome:
            if i >= len(keys):
                raise ValueError("Chromosomes named: %s not found in selected assembly!"%(chromosome))
            elif chromosomes[keys[i]]["name"] == chromosome:
                isSearchingChromosome   = False
                chomosomes_used.append( chromosomes[keys[i]] )
            else:
                i +=1

    with track.load(sqlout,  format="sql") as t:
        t.chrmeta = chomosomes_used

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
                                        true_positive_data_path, false_positive_data_path,
                                        description='',
                                        threshold=0,    via='lsf', keep_max_only=False,
                                        alpha=0.05, nb_false_positive_hypotesis=1.0
                                    ):
    """
    sqls or beds take an array (list or tuple) like:
    sqls = (sql, sql_random)
    beds = (bed, bed_random)
    """
    if true_positive_data_path is None or false_positive_data_path is None:
        raise ValueError(u"Variables true_positive_data_path or false_positive_data_path is set to None!")
    true_positive_data_path = os.path.normcase(os.path.expanduser(true_positive_data_path))
    false_positive_data_path= os.path.normcase(os.path.expanduser(false_positive_data_path))
    true_positive_result    = None
    false_positive_result   = None

    # original
    true_positive_result    = save_motif_profile(
                                                    ex, motifs, background, genrep, chromosomes, true_positive_data_path,
                                                    description     = description, threshold=threshold,
                                                    via             = via, keep_max_only = keep_max_only
                                                )
    # random
    false_positive_result   = save_motif_profile(
                                                    ex, motifs, background, genrep, chromosomes, false_positive_data_path,
                                                    description     = description, threshold=threshold,
                                                    via             = via, keep_max_only=keep_max_only
                                                )

    fp_scores = None
    tp_scores = None
    with track.load(false_positive_result,  format="sql") as t:
        fp_scores = t.get_scores_frequencies()
    with track.load(true_positive_result,  format="sql") as t:
        tp_scores = t.get_scores_frequencies()
    fdr = false_discovery_rate  (
                                    fp_scores, tp_scores,
                                    alpha=alpha, nb_false_positive_hypotesis=nb_false_positive_hypotesis
                                )
    return true_positive_result,fdr

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

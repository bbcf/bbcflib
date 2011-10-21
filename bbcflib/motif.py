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
from . import track
from .common import set_file_descr, compress

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
    def _xmltree(_c,_t):
        seq_name = {}
        seq_chr = None
        for it in _t.getiterator():
            if it.tag == 'sequence':
                seq_name[it.attrib['id']] = it.attrib['name']
            if it.tag == 'scanned_sites':
                name = seq_name[it.attrib['sequence_id']]
                name,seq_chr,start,end = re.search(r'(.*)\|(.+):(\d+)-(\d+)',name).groups()
            if it.tag == 'scanned_site' and _c == seq_chr:
                start = int(start)+int(it.attrib['position'])-1
                end = str(start+ncol[it.attrib['motif_id']])
                start = str(start)
                strnd = it.attrib['strand'] == 'plus' and 1 or -1
                score = it.attrib['pvalue']
                yield (start,end,it.attrib['motif_id'],score,strnd)
    outsql = unique_filename_in()
    chrlist = dict((v['name'], {'length': v['length']}) for v in chromosomes.values())
    with track.new(outsql, format='sql', chrmeta=chrlist, datatype='qualitative') as t:
        track_result.attributes['source'] = 'Meme'
        for chrom in chrlist.keys():
            t.write(chrom,_xmltree(chrom,tree),fields=['start','end','name','score','strand'])
    return {'sql':outsql,'matrices':allmatrices}
                

def parallel_meme( ex, genrep, chromosomes, regions, 
                   name=None, meme_args=None, via='lsf' ):
    """Fetches sequences, then calls ``meme``
    on them and finally saves the results in the repository.
    """
    if meme_args is None:
        meme_args   = []
    if not(isinstance(regions,list)):
        regions = [regions]
    if not(isinstance(name,list)):
        if name == None:
            name = '_'
        name = [name]
    futures = {}
    for i,n in enumerate(name):
        (fasta, size) = genrep.fasta_from_regions( chromosomes, regions[i], out=unique_filename_in() )
        futures[n] = meme.nonblocking( ex, fasta, maxsize=size*1.5, args=meme_args, via=via )
    all_res = {}
    for n,f in futures.iteritems():
        meme_out = f.wait()
        archive = compress( ex, meme_out, 'gz' )
        meme_res = parse_meme_xml( ex, os.path.join(meme_out, "meme.xml"),
                                   chromosomes )
        ex.add( os.path.join(meme_out, "meme.html"),
                description=set_file_descr(n+"_meme.html",tag='meme',type='html') )
        ex.add( meme_res['sql'], description=set_file_descr(n+"_meme_sites.sql",tag='meme',type='sql') )
        ex.add( archive, description=set_file_descr(n+"_meme.tgz",tag='meme',type='tar') )
        for i,motif in enumerate(meme_res['matrices'].keys()):
            ex.add( meme_res['matrices'][motif], description=set_file_descr(n+"_meme_"+motif+".txt",tag='meme',type='txt') )
            ex.add( os.path.join(meme_out, "logo"+str(i+1)+".png"), description=set_file_descr(n+"_meme_"+motif+".png",tag='meme',type='png') )
        all_res[n] = meme_res
    return all_res

@program
def motif_scan( fasta, motif, background, threshold=0 ):
    """Binding for the ``S1K`` motif scanner.
    """
    call = ["S1K", motif, background, str(threshold), fasta]
    return {"arguments": call, "return_value": None}

def save_motif_profile( ex, motifs, background, genrep, chromosomes, regions,
                        keep_max_only=False, threshold=0, description='', via='lsf' ):
    """Scan a set of motifs on a set of regions and saves the results as an sql file.
    The 'motifs' argument is a single PWM file or a dictionary with keys motif names and values PWM files 
    with 'n' rows like:
    "1 p(A) p(C) p(G) p(T)" 
    where the sum of the 'p's is 1 and the first column allows to skip a position with a '0'.
    """
    fasta, size = genrep.fasta_from_regions( chromosomes, regions )
    sqlout = unique_filename_in()
    if not(isinstance(motifs, dict)):
        motifs = {"_": motifs}
#        raise ValueError("'Motifs' must be a dictionary with keys 'motif_names' and values the PWMs.")
    futures = {}
    for name, pwm in motifs.iteritems():
        output = unique_filename_in()
        futures[name] = ( output, motif_scan.nonblocking( ex, fasta, pwm, background,
                                                          threshold, stdout=output, via=via ) )
    chrlist = dict((v['name'], {'length': v['length']}) for v in chromosomes.values())
##############
    def _parse_s1k(_f,_c,_n):
        if keep_max_only:
            maxscore = {}
        with open(_f, 'r') as fin:
            for line in fin:
                row = line.split("\t")
                sname,seq_chr,start,end = re.search(r'(.*)\|(.+):(\d+)-(\d+)',row[0]).groups()
                if seq_chr != _c: 
                    continue
                start = int(row[3])+int(start)-1
                tag = row[1]
                end = start+len(tag)
                score = float(row[2])
                strnd = row[4]=='+' and 1 or -1
                if keep_max_only:
                    if not(sname in maxscore) or score>maxscore[sname][3]:
                        maxscore[sname] = (start,end,_n+":"+tag,score,strnd)
                else:
                    yield (start,end,_n+":"+tag,score,strnd)
        if keep_max_only:
            for x in maxscore.values():
                yield x
##############
    with track.new( sqlout, format="sql", datatype="qualitative", chrmeta=chrlist ) as track_result:
        track_result.attributes['source'] = 'S1K'
        for name, future in futures.iteritems():
            _ = future[1].wait()
            for chrom in chrlist.keys():
                track_result.write( chrom, _parse_s1k(future[0],chrom,name),
                                    fields=['start','end','name','score','strand'] )
    ex.add( sqlout, description=set_file_descr(description+"_motif_scan.sql",tag="sql",type="sql") )
    return sqlout

def FDR_threshold( ex, motif, background, genrep, chromosomes, regions, alpha=.1, nb_samples=1, via='lsf' ):
    """
    Computes a score threshold for 'motif' on 'regions' based on a false discovery rate < alpha and returns the 
    threshold or a dictionary with keys thresholds and values simulated FDRs when alpha < 0.
    """
    fasta, size = genrep.fasta_from_regions( chromosomes, regions )
    shuf_fasta, shuf_size = genrep.fasta_from_regions( chromosomes, regions, shuffled=True )
    output = unique_filename_in()
#### Threshold at -100 to get all scores!
    future = motif_scan.nonblocking( ex, fasta, motif, background, -100, stdout=output, via=via )
    shuf_futures = {}
    for i in range(nb_samples):
        out = unique_filename_in()
        shuf_futures[out] = motif_scan.nonblocking( ex, shuf_fasta, motif, background, -100, stdout=out, via=via )
    _ = future.wait()
    TP_scores = {}
    ntp = 0
    with open(output, 'r') as fin:
        for line in fin:
            row = line.split("\t")
            score = int(round(float(row[2])))
            if score in TP_scores: 
                TP_scores[score] += 1
            else: 
                TP_scores[score] = 1
            ntp += 1
    scores = sorted(TP_scores.keys(),reverse=True)
    scores = [scores[0]+1]+scores+[-101]
    FP_scores = dict((k,0) for k in scores)
    nfp = 0
    for file,fut in shuf_futures.iteritems():
        _ = fut.wait()
        with open(file, 'r') as fin:
            for line in fin:
                row = line.split("\t")
                fscore = int(round(float(row[2])))
                tscore = max([k for k in scores if k<=fscore])
                FP_scores[tscore] += 1
                nfp += 1
    TP_scores[scores[-1]] = ntp
    TP_scores[scores[0]] = 0
    FP_scores[scores[-1]] = nfp
    for i,sc in enumerate(scores[1:-1]):
        TP_scores[sc] = TP_scores[scores[i]]+TP_scores[sc]
        FP_scores[sc] = FP_scores[scores[i]]+FP_scores[sc]
    cur_fdr = 1.0
    threshold = scores[0]
    for k in sorted(FP_scores.keys()):
        if TP_scores[k] > 0 and FP_scores[k]/float(TP_scores[k]) < cur_fdr:
            cur_fdr = FP_scores[k]/float(TP_scores[k])
            if cur_fdr <= alpha:
                threshold = k
                break
        FP_scores[k] = cur_fdr
    if alpha < 0: return FP_scores
    return threshold

def sqlite_to_false_discovery_rate( ex, motif, background, genrep, chromosomes, regions,
                                    alpha=0.05, nb_samples=1,
                                    description='', via='lsf' ):
    """
    Computes a score threshold for 'motif' on 'regions' based on a false discovery rate < alpha and returns the 
    thresholded profile.
    """
    threshold = FDR_threshold( motif, background, genrep, chromosomes, regions, alpha=alpha, nb_samples=nb_samples, via=via )
    track = save_motif_profile( ex, motif, background, genrep, chromosomes, regions,
                                threshold=threshold, description=description, via=via )
    return track, threshold

def parse_meme_html_output(ex, meme_file, fasta, chromosomes):
    """Legacy signature"""
    meme_file = os.path.splitext(meme_file)[0]+".xml"
    return parse_meme_xml( ex, meme_file, chromosomes )
    
def add_meme_files( ex, genrep, chromosomes, description='',
                    bed=None, sql=None, meme_args=None, via='lsf' ):
    """Legacy signature"""
    regions = sql or bed
    return parallel_meme( ex, genrep, chromosomes, regions=regions, 
                          name=description, meme_args=meme_args, via=via )
    
#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

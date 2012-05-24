"""
=====================
Module: bbcflib.motif
=====================

No documentation
"""

# Built-in modules #
import re, os, tarfile
from operator import add

# Internal modules #
from bbcflib import btrack as track
from bbcflib.common import set_file_descr, unique_filename_in, gzipfile

# Other modules #
from bein import program
from bein.util import touch

################################################################################
@program
def meme( fasta, outdir, maxsize=10000000, args=None ):
    """Binding for the ``meme`` motif finder.
    """
    if args is None:
        args = []
    if not("minw" in args): args += ["-minw","6"]
    if not("maxw" in args): args += ["-maxw","16"]
    call = ["meme", fasta, "-o", outdir, "-dna", "-maxsize", str(maxsize)]+args
    return {"arguments": call, "return_value": None}

def parse_meme_xml( ex, meme_file, chrmeta ):
    """ Parse meme xml file and convert to track """
    from xml.etree import ElementTree as ET
    touch(ex,meme_file)
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
    def _xmltree(_t):#(_c,_t):
        seq_name = {}
        seq_chr = None
        for it in _t.getiterator():
            if it.tag == 'sequence':
                seq_name[it.attrib['id']] = it.attrib['name']
            if it.tag == 'scanned_sites':
                name = seq_name[it.attrib['sequence_id']]
                name,seq_chr,start,end = re.search(r'(.*)\|(.+):(\d+)-(\d+)',name).groups()
            if it.tag == 'scanned_site':# and _c == seq_chr:
                start = int(start)+int(it.attrib['position'])-1
                end = start+ncol[it.attrib['motif_id']]
                strnd = it.attrib['strand'] == 'plus' and 1 or -1
                score = it.attrib['pvalue']
                yield (seq_chr,str(start),str(end),it.attrib['motif_id'],score,strnd)
    outsql = unique_filename_in()+".sql"
    outtrack = track.track(outsql, chrmeta=chrmeta, info={'datatype':'qualitative'},
                           fields=['start','end','name','score','strand'])
    outtrack.write(track.FeatureStream(_xmltree(tree),fields=['chr']+outtrack.fields))
    outtrack.close()
    return {'sql':outsql,'matrices':allmatrices}


def parallel_meme( ex, assembly, regions, name=None, meme_args=None, via='lsf' ):
    """Fetches sequences, then calls ``meme`` on them and finally saves the results in the repository.
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
    fasta_files = []
    for i,n in enumerate(name):
        (fasta, size) = assembly.fasta_from_regions( regions[i], out=unique_filename_in() )
        tmpfile = unique_filename_in()
        outdir = unique_filename_in()
        futures[n] = (outdir, meme.nonblocking( ex, fasta, outdir, 
                                                maxsize=(size*3)/2, 
                                                args=meme_args, via=via, 
                                                stderr=tmpfile ))
        fasta_files.append(fasta)
    all_res = {}
    for n,f in futures.iteritems():
        f[1].wait()
        meme_out = f[0]
        archive = unique_filename_in()
        tgz = tarfile.open(archive, "w:gz")
        tgz.add( meme_out )
        tgz.close()
        meme_res = parse_meme_xml( ex, os.path.join(meme_out, "meme.xml"), assembly.chrmeta )
        if os.path.exists(os.path.join(meme_out, "meme.html")):
            ex.add( os.path.join(meme_out, "meme.html"),
                    description=set_file_descr(n+"_meme.html",step='meme',type='html',group=n) )
        ex.add( meme_res['sql'], description=set_file_descr(n+"_meme_sites.sql",step='meme',type='sql',group=n) )
        ex.add( archive, description=set_file_descr(n+"_meme.tgz",step='meme',type='tar',group=n) )
        gzipfile(ex,fasta_files[n])
        ex.add( fasta_files[n]+".gz", description=set_file_descr(n+"_sites.fa.gz",step='meme',type='fasta',group=n) )
        for i,motif in enumerate(meme_res['matrices'].keys()):
            ex.add( meme_res['matrices'][motif], description=set_file_descr(n+"_meme_"+motif+".txt",step='meme',type='txt',group=n) )
            ex.add( os.path.join(meme_out, "logo"+str(i+1)+".png"), description=set_file_descr(n+"_meme_"+motif+".png",step='meme',type='png',group=n) )
        all_res[n] = meme_res
    return all_res

@program
def motif_scan( fasta, motif, background, threshold=0 ):
    """Binding for the ``S1K`` motif scanner.
    """
    call = ["S1K", motif, background, str(threshold), fasta]
    return {"arguments": call, "return_value": None}

def save_motif_profile( ex, motifs, background, assembly, regions, keep_max_only=False, 
                        threshold=0, description='motif_scan.sql', via='lsf' ):
    """Scan a set of motifs on a set of regions and saves the results as an sql file.
    The 'motifs' argument is a single PWM file or a dictionary with keys motif names and values PWM files
    with 'n' rows like:
    "1 p(A) p(C) p(G) p(T)"
    where the sum of the 'p's is 1 and the first column allows to skip a position with a '0'.
    """
    fasta, size = assembly.fasta_from_regions( regions )
    sqlout = unique_filename_in()
    if not(isinstance(motifs, dict)):
        motifs = {"_": motifs}
#        raise ValueError("'Motifs' must be a dictionary with keys 'motif_names' and values the PWMs.")
    futures = {}
    for name, pwm in motifs.iteritems():
        output = unique_filename_in()
        futures[name] = ( output, motif_scan.nonblocking( ex, fasta, pwm, background,
                                                          threshold, stdout=output, via=via ) )
##############
    def _parse_s1k(_f,_n):
        if keep_max_only:
            maxscore = {}
        with open(_f, 'r') as fin:
            for line in fin:
                row = line.split("\t")
                sname,seq_chr,start,end = re.search(r'(.*)\|(.+):(\d+)-(\d+)',row[0]).groups()
                start = int(row[3])+int(start)-1
                tag = row[1]
                end = start+len(tag)
                score = float(row[2])
                strnd = row[4]=='+' and 1 or -1
                if keep_max_only:
                    if not((sname,seq_chr) in maxscore) or score>maxscore[(sname,seq_chr)][3]:
                        maxscore[(sname,seq_chr)] = (seq_chr,start,end,_n+":"+tag,score,strnd)
                else:
                    yield (seq_chr,start,end,_n+":"+tag,score,strnd)
        if keep_max_only:
            for x in maxscore.values():
                yield x
##############
    track_result = track.track( sqlout, chrmeta=assembly.chrmeta,
                                info={'datatype':'qualitative'},
                                fields=['start','end','name','score','strand'] )
    for name, future in futures.iteritems():
        future[1].wait()
        track_result.write( track.FeatureStream(_parse_s1k(future[0], name ),
                                                fields=['chr']+track_result.fields) )
    track_result.close()
    ex.add( sqlout, description=description )
    return sqlout

def FDR_threshold( ex, motif, background, assembly, regions, alpha=.1, nb_samples=1, via='lsf' ):
    """
    Computes a score threshold for 'motif' on 'regions' based on a false discovery rate < alpha and returns the
    threshold or a dictionary with keys thresholds and values simulated FDRs when alpha < 0.
    """
    fasta, size = assembly.fasta_from_regions( regions )
    shuf_fasta, shuf_size = assembly.fasta_from_regions( regions, shuffled=True )
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

def sqlite_to_false_discovery_rate( ex, motif, background, assembly, regions, alpha=0.05, nb_samples=1,
                                    description='', via='lsf' ):
    """
    Computes a score threshold for 'motif' on 'regions' based on a false discovery rate < alpha and returns the
    thresholded profile.
    """
    threshold = FDR_threshold( motif, background, assembly, regions, alpha=alpha, nb_samples=nb_samples, via=via )
    track = save_motif_profile( ex, motif, background, assembly, regions,
                                threshold=threshold, description=description, via=via )
    return track, threshold

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

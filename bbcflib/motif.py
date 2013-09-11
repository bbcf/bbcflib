"""
=====================
Module: bbcflib.motif
=====================

"""

# Built-in modules #
import re, os, tarfile
from operator import add

# Internal modules #
from bbcflib.track import track, FeatureStream
from bbcflib.common import set_file_descr, unique_filename_in, gzipfile, fasta_length

# Other modules #
from bein import program
from bein.util import touch

################################################################################
@program
def meme( fasta, outdir, background=None, maxsize=10000000, args=None ):
    """Binding of the ``meme`` motif finder."""
    if args is None: args = []
    if not("-minw" in args): args += ["-minw","6"]
    if not("-maxw" in args): args += ["-maxw","16"]
    if os.path.exists(background): args += ["-bfile",background]
    call = ["meme", fasta, "-o", outdir, "-dna", "-maxsize", str(maxsize)]+args
    return {"arguments": call, "return_value": None}

@program
def memechip( fasta, outdir, background=None, args=None ):
    """Binding of the ``meme-chip`` pipeline."""
    if args is None: args = []
    if not("-ccut" in args): args += ["-ccut","300"]
    if not("-minw" in args): args += ["-meme-minw","6"]
    if not("-maxw" in args): args += ["-meme-maxw","16"]
    if not("-nmeme" in args): args += ["-nmeme","300"]
    if os.path.exists(background): args += ["-bfile",background]
    call = ["meme-chip", fasta, "-o", outdir]+args
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
    outtrack = track(outsql, chrmeta=chrmeta, info={'datatype':'qualitative'},
                     fields=['start','end','name','score','strand'])
    outtrack.write(FeatureStream(_xmltree(tree),fields=['chr']+outtrack.fields))
    outtrack.close()
    return {'sql':outsql,'matrices':allmatrices}


def parallel_meme( ex, assembly, regions, name=None, chip=False, meme_args=None, via='lsf' ):
    """Fetches sequences, then calls ``meme`` on them and finally saves the results in the repository.
    
    """
    if meme_args is None: meme_args = []
    if not(isinstance(regions,list)): regions = [regions]
    if not(isinstance(name,list)): name = [name or '_']
    futures = {}
    fasta_files = {}
    background = assembly.statistics(unique_filename_in(),frequency=True)
    genomeRef = assembly.untar_genome_fasta()
    for i,n in enumerate(name):
        (fasta, size) = assembly.fasta_from_regions( regions[i], out=unique_filename_in(), path_to_ref=genomeRef )
        tmpfile = unique_filename_in()
        outdir = unique_filename_in()
        if chip:
            futures[n] = (outdir, memechip.nonblocking( ex, fasta, outdir, background,
                                                        args=meme_args, via=via, 
                                                        stderr=tmpfile, memory=6 ))
        else:
            futures[n] = (outdir, meme.nonblocking( ex, fasta, outdir, background,
                                                    maxsize=(size*3)/2, args=meme_args,
                                                    via=via, stderr=tmpfile, memory=6 ))
        fasta_files[n] = fasta
    all_res = {}
    for n,f in futures.iteritems():
        f[1].wait()
        meme_out = f[0]
        archive = unique_filename_in()
        tgz = tarfile.open(archive, "w:gz")
        tgz.add( meme_out, arcname=n[1]+"_meme",
                 exclude=lambda x: os.path.basename(x) in [fasta_files[n],background] )
        tgz.close()
        ex.add( archive, description=set_file_descr(n[1]+"_meme.tgz",
                                                    step='meme', type='tar',
                                                    groupId=n[0]) )
        gzipfile(ex,fasta_files[n],args=["-f"])
        ex.add( fasta_files[n]+".gz",
                description=set_file_descr(n[1]+"_sites.fa.gz",
                                           step='meme', type='fasta',
                                           groupId=n[0]) )
        if not(chip) and os.path.exists(os.path.join(meme_out, "meme.xml")):
            meme_res = parse_meme_xml( ex, os.path.join(meme_out, "meme.xml"),
                                       assembly.chrmeta )
            if os.path.exists(os.path.join(meme_out, "meme.html")):
                ex.add( os.path.join(meme_out, "meme.html"),
                        description=set_file_descr(n[1]+"_meme.html",
                                                   step='meme', type='html', 
                                                   groupId=n[0]) )
            ex.add( meme_res['sql'], description=set_file_descr(n[1]+"_meme_sites.sql",
                                                                step='meme', type='sql',
                                                                groupId=n[0]) )
            for i,motif in enumerate(meme_res['matrices'].keys()):
                ex.add( meme_res['matrices'][motif],
                        description=set_file_descr(n[1]+"_meme_"+motif+".txt",
                                                   step='meme', type='txt', 
                                                   groupId=n[0]) )
                ex.add( os.path.join(meme_out, "logo"+str(i+1)+".png"),
                        description=set_file_descr(n[1]+"_meme_"+motif+".png",
                                                   step='meme', type='png', 
                                                   groupId=n[0]) )
            all_res[n] = meme_res
    return all_res

@program
def motif_scan( fasta, motif, background, threshold=0 ):
    """Binding of the ``S1K`` motif scanner."""
    call = ["S1K", motif, background, str(threshold), fasta]
    return {"arguments": call, "return_value": None}

def save_motif_profile( ex, motifs, assembly, regions=None, fasta=None, background=None, keep_max_only=False,
                        threshold=0, output=None, description=None, via='lsf' ):
    """
    Scans a set of motifs on a set of regions and saves the results as a track file.
    The 'motifs' argument is a single PWM file or a dictionary with keys motif names and values PWM files
    with 'n' rows like:
    "1 p(A) p(C) p(G) p(T)"
    where the sum of the 'p's is 1 and the first column allows to skip a position with a '0'.
    """
    if regions is not None:
        fasta, size = assembly.fasta_from_regions( regions )
    if fasta is None:
        raise ValueError("Provide either a fasta file or a valid list of regions")
    if assembly:
        chrmeta = assembly.chrmeta
    else:
        chroms = fasta_length.nonblocking( ex, fasta, via=via ).wait()
        chrmeta = dict([(v['name'], {'length': v['length']}) for v in chroms.values()])
    if output is None:
        sqlout = unique_filename_in()
    else:
        sqlout = output
    if not(isinstance(motifs, dict)):
        motifs = {"_": motifs}
    futures = {}
    if background is None:
        background = assembly.statistics(unique_filename_in(),frequency=True,
                                         matrix_format=True)
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
                row = line.strip().split("\t")
                spatt = re.search(r'(.*)\|(.+):(\d+)-(\d+)',row[0])
                if spatt:
                    sname,seq_chr,start,end = spatt.groups()
                else:
                    sname = row[0]
                    seq_chr = row[0]
                    start = 0
                    end = 0
                start = int(row[3])+int(start)-1
                tag = row[1]
                end = start+len(tag)
                score = float(row[2])
                strnd = 1 if row[4][0] == '+' else -1
                if keep_max_only:
                    if not((sname,seq_chr) in maxscore) or score>maxscore[(sname,seq_chr)][3]:
                        maxscore[(sname,seq_chr)] = (seq_chr,start,end,_n+":"+tag,score,strnd)
                else:
                    yield (seq_chr,start,end,_n+":"+tag,score,strnd)
        if keep_max_only:
            for x in maxscore.values():
                yield x
##############
    track_result = track( sqlout, chrmeta=chrmeta, info={'datatype':'qualitative'},
                          fields=['chr','start','end','name','score','strand'] )
    for name, future in futures.iteritems():
        future[1].wait()
        track_result.write( FeatureStream(_parse_s1k(future[0], name), fields=track_result.fields) )
    track_result.close()
    if description is not None: ex.add( sqlout, description=description )
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
    outtrack = save_motif_profile( ex, motif, assembly, regions, background=background,
                                   threshold=threshold, description=description, via=via )
    return outtrack, threshold

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

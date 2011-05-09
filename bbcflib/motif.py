"""
===============
bbcflib.motif
===============
"""
import os
import sqlite3
from bein import *
from bein.util import *
from bbcflib import common

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

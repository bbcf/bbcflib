"""
===============
bbcflib.motif
===============
"""
import os
import sqlite3
from bein import *
from bein.util import *

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

def save_motif_profile( ex, motif, background, genrep, chromosomes,
                        description='', sql=None, bed=None, threshold=0, via='lsf' ):
    fasta,size = genrep.fasta_from_bed( chromosomes, out=unique_filename_in(), 
                                        bed=bed, sql=sql )
    output = unique_filename_in()
    future = motif_scan.nonblocking( ex, fasta, motif, background, threshold, 
                                     via=via, stdout=output )
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
        cur = connection.cursor()
        for v in chromosomes.values():
            cur.execute('select start,name from "'+v['name']+'"')
            connection.commit()
            for row in cur:
                name = re.search(r'(\S+)\s*',row[1]).groups()[0]
                regions[name] = (v['name'],int(row[0]))
            cur.close()
    else:
        raise TypeError("save_motif_profile requires either a 'sqlite' or a 'bed' file.")
    sqlout = unique_filename_in()
    connection = sqlite3.connect( sqlout )
    connection.execute('create table chrNames (name text, length integer)')
    connection.execute('create table attributes (key text, value text)')
    connection.execute('insert into attributes (key,value) values ("datatype","qualitative")')
    vals = [(v['name'],str(v['length'])) for v in chromosomes.values()]
    connection.executemany('insert into chrNames (name,length) values (?,?)',vals)
    [connection.execute('create table "'+v['name']
                        +'" (start integer, end integer, score real, strand char, name string)') 
     for v in chromosomes.values()]
    connection.commit()
    vals = []
    null = future.wait()
    cur_chr = ''
    with open(output,'r') as f:
        for l in f:
            s = l.rstrip('\n').split('\t')
            reg = regions[s[0]]
            start = reg[1]+int(s[3])
            vals.append((start-1,start+len(s[1]),float(s[2]),s[4],s[1]))
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

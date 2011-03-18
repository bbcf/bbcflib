"""
===============
bbcflib.motif
===============
"""
import os
from bein import *
from bein.util import *

@program
def meme( fasta, maxsize=10000000, args=[] ):
    """Binding for the ``meme`` motif finder.
    """
    outname = unique_filename_in()
    call = ["meme",fasta,"-o",outname,"-dna","-maxsize",str(maxsize),"-revcomp"]+args
    return {"arguments": call, "return_value": outname}

def add_meme_files( ex, genrep, chromosomes, description="",
                    bed=None, sql=None, meme_args=[], via='lsf' ):
    """Fetches sequences, then calls ``meme`` 
    on them and finally saves the results in the repository.
    """
    fasta,size = genrep.fasta_from_bed( chromosomes, out=unique_filename_in(), 
                                        bed=bed, sql=sql )
    meme_out = meme.nonblocking( ex, fasta, maxsize=size*1.5, via=via ).wait()
    html = os.path.join(meme_out,"meme.html")
    ex.add( html, description="html:"+description+"meme.html" )
    return meme_out

@program
def motif_scan( fasta, motif, background, threshold=0 ):
    """Binding for the ``S1K`` motif scanner.
    """
    #got to do something with output
    call = ["S1K",motif,background,threshold,fasta]
    return {"arguments": call, "return_value": ''}

"""
======================
Module: bbcflib.common
======================

Utility functions common to several pipelines.
"""

# Built-in modules #
import sys

###############################################################################
def normalize_url(url):
    """Produce a fixed form for an HTTP URL.

    Make sure the URL begins with http:// and does *not* end with a /.

    >>> normalize_url('http://bbcf.epfl.ch')
    'http://bbcf.epfl.ch'
    >>> normalize_url('http://bbcf.epfl.ch/')
    'http://bbcf.epfl.ch'
    >>> normalize_url('bbcf.epfl.ch/')
    'http://bbcf.epfl.ch'
    """
    url = url.lower()
    if not(url.startswith("http://")):
        url = "http://" + url
    if url.endswith("/"):
        url = url[:-1]
    return url

###############################################################################
def cat(files):
    """Concatenates files.
    """
    if len(files) > 1:
        out = unique_filename_in()
        with open(out,"w") as f1:
            for inf in files:
                with open(inf,"r") as f2:
                    [f1.write(l) for l in f2]
                    f2.close()
            f1.close()
    elif len(files) == 1:
        out = files[0]
    else:
        out = None
    return out

###############################################################################
def create_sql_track(path, chrmeta, datatype="quantitative", format='sql', name="Unamed"):
    ''' This function shouldn't be used'''
    from bbcflib import track
    with track.new(path, format, name, chrmeta, datatype) as t:
        for chrom in chrmeta: t.write(chrom, (), getattr(track.Track, datatype + '_fields'))

###############################################################################
###############################################################################
##################### BIG TRY STATEMENT STARTING... ###########################
###############################################################################
###############################################################################
try:
    from bein import *

    #-------------------------------------------------------------------------#
    def get_files( id_or_key, minilims ):
        """Retrieves a dictionary of files created by an htsstation job identified by its key or bein id in a MiniLIMS.

        The dictionary keys are the file types (e.g. 'pdf', 'bam', 'py' for python objects), the values are dictionaries with keys repository file names and values actual file descriptions (names to provide in the user interface).
        """
        import re
        if isinstance(id_or_key, str):
            try:
                exid = max(minilims.search_executions(with_text=id_or_key))
            except ValueError, v:
                raise ValueError("No execution with key "+id_or_key)
        else:
            exid = id_or_key
        file_dict = {}
        all = dict((y['repository_name'],y['description']) for y in
                   [minilims.fetch_file(x) for x in minilims.search_files(source=('execution',exid))])
        for f,d in all.iteritems():
            cat,name = re.search(r'([^:]+):(.*)$',d).groups()
            name = re.sub(r'\s+\(BAM INDEX\)','.bai',name)
            if cat in file_dict:
                file_dict[cat].update({f: name})
            else:
                file_dict[cat] = {f: name}
        return file_dict

    #-------------------------------------------------------------------------#
    @program
    def run_gMiner( job ):
        import pickle
        job_file = unique_filename_in()
        with open(job_file,'w') as f:
            pickle.dump(job,f)
        def get_output_files(p):
            with open(job_file,'r') as f:
                job = pickle.load(f)
            return job['job_output']
        return {"arguments": ["run_gminer.py",job_file],
                "return_value": get_output_files}

    #-------------------------------------------------------------------------#
    def merge_sql( ex, sqls, names, description="merged.sql", outdir=None, via='lsf' ):
        """Run ``gMiner``'s 'merge_score' function on a set of sql files
        """
        import os
        if outdir == None:
            outdir = unique_filename_in()
        if not(os.path.exists(outdir)):
            os.mkdir(outdir)
        if not(isinstance(names,list)):
            names = []
        if len(names) < len(sqls):
            n = sqls
            n[:len(names)] = names
            names = n
        gMiner_job = dict([('track'+str(i+1),f) for i,f in enumerate(sqls)]
                          +[('track'+str(i+1)+'_name',str(ni))
                            for i,ni in enumerate(names)])
        gMiner_job['operation_type'] = 'genomic_manip'
        gMiner_job['manipulation'] = 'merge_scores'
        gMiner_job['output_location'] = outdir
        files = run_gMiner.nonblocking(ex,gMiner_job,via=via).wait()
        ex.add( files[0], description=description )
        return files[0]

    #-------------------------------------------------------------------------#
    @program
    def merge_two_bed(file1,file2):
        """Binds ``intersectBed`` from the 'BedTools' suite.
        """
        return {"arguments": ['intersectBed','-a',file1,'-b',file2], "return_value": None}

    #-------------------------------------------------------------------------#
    def merge_many_bed(ex,files,via='lsf'):
        """Runs ``intersectBed`` iteratively over a list of bed files.
        """
        out = files[0]
        for f in files[1:]:
            next = unique_filename_in()
            _ = merge_two_bed.nonblocking( ex, out, f, via=via, stdout=next ).wait()
            out = next
        return out

    #-------------------------------------------------------------------------#
    @program
    def join_pdf(files):
        """Uses 'ghostscript' to join several pdf files into one.
        """
        out = unique_filename_in()
        gs_args = ['gs','-dBATCH','-dNOPAUSE','-q','-sDEVICE=pdfwrite',
                   '-sOutputFile=%s'%out]
        gs_args += files
        return {"arguments": gs_args, "return_value": out}

    #-------------------------------------------------------------------------#
    @program
    def compress(path, compression_type="lxzma"):
        """
        compression type allowed:
        - gunzip, gz
        - bzip2, bz
        - lxzma, xz
        """
        archive = unique_filename_in()
        call    = None
        if compression_type == "lxzma" or compression_type == "xz":
            call = ["tar", "cJvf", archive, path]
        elif compression_type == "bzip2" or compression_type == "bz2":
            call = ["tar", "cjvf", archive, path]
        elif compression_type == "gunzip" or compression_type == "gz":
            call = ["tar", "czvf", archive, path]
        else:
            raise ValueError("Compression type: %s not yet supported!" %(compression_type))
        return {"arguments": call, "return_value": archive}

    #-------------------------------------------------------------------------#
    @program
    def uncompress(path):
        """
        uncompress tar archive
        """
        output = unique_filename_in()
        call = ["tar", "xvf", path]
        return {"arguments": call, "return_value": output}

    #-------------------------------------------------------------------------#
    @program
    def scp(source, destination, args=[]):
        output = unique_filename_in()
        call = ["scp"] + args + [ source, destination ]
        return {"arguments": call, "return_value": output}

except:
    print >>sys.stderr, "Bein not found.  Skipping some common functions."

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

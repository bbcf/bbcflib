"""
======================
Module: bbcflib.common
======================

Utility functions common to several pipelines.
"""

# Built-in modules #
import os, sys, time, json, csv

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
    # url = url.lower()
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
        out = rstring()
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
def set_file_descr(filename,**kwargs):
    """Implements file naming compatible with the 'get_files' function.
    Examples: 
    >>>> set_file_descr("toto",**{'tag':'pdf','step':1,'type':'doc','comment':'ahaha'})
    'pdf:toto[step:1,type:doc] (ahaha)'
    if 'tag' and/or comment are ommitted:
    >>>> set_file_descr("toto",step=1,type='doc')
    'toto[step:1,type:doc]'
    """
    file_descr = filename
    argskeys = kwargs.keys()
    if 'tag' in kwargs:
        file_descr = kwargs['tag']+":"+filename
        argskeys.remove('tag')
    comment = ''
    if 'comment' in argskeys:
        comment = " ("+str(kwargs['comment'])+")"
        argskeys.remove('comment')
    plst = (",").join([str(k)+":"+str(kwargs[k]) for k in argskeys])
    file_descr += "["+plst+"]"
    file_descr += comment
    return file_descr
#    tag:filename[group:grpId,step:stepId,type:fileType,view:admin] (comment) 

#-------------------------------------------------------------------------#
def get_files_old( id_or_key, minilims ):
    """Retrieves a dictionary of files created by an htsstation job identified by its key
    or bein id in a MiniLIMS.

    The dictionary keys are the file types (e.g. 'pdf', 'bam', 'py' for python objects),
    the values are dictionaries with keys repository file names and values actual file
    descriptions (names to provide in the user interface).
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
def get_files( id_or_key, minilims, by_type=True ):
    """Retrieves a dictionary of files created by an htsstation job identified by its key
    or bein id in a MiniLIMS.

    The dictionary keys are the file types (e.g. 'pdf', 'bam', 'py' for python objects),
    the values are dictionaries with keys repository file names and values actual file
    descriptions (names to provide in the user interface).
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
        tag_d = d.split(':')
        tag = None
        if len(tag_d)>1 and not(re.search("\[",tag_d[0])):
            tag = tag_d[0]
            d = ":".join(tag_d[1:])
        pars_patt = re.search(r'\[(.*)\]',d)
        if not(pars_patt):
            pars = "group:0,step:0,type:none,view:admin"
            re.sub(r'([^\s\[]*)',r'\1['+pars+']',x,1)
        else:
            pars = pars_patt.groups()[0]
        par_dict = dict([x.split(":") for x in pars.split(",")])
	if re.search(r'\s+\(BAM INDEX\)',d):
            d = re.sub('[','.bai[',d)
        cat = (by_type and par_dict.get('type')) or tag or 'none'
        if cat in file_dict:
            file_dict[cat].update({f: d})
        else:
            file_dict[cat] = {f: d}
    return file_dict

#-------------------------------------------------------------------------#
def merge_sql( ex, sqls, names, description="merged.sql", outdir=None, via='lsf' ):
    """Run ``gMiner``'s 'merge_score' function on a set of sql files
    """
    import os
    if outdir == None:
        outdir = rstring()
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
def merge_many_bed(ex,files,via='lsf'):
    """Runs ``intersectBed`` iteratively over a list of bed files.
    """
    out = files[0]
    for f in files[1:]:
        next = rstring()
        _ = merge_two_bed.nonblocking( ex, out, f, via=via, stdout=next ).wait()
        out = next
    return out

#-------------------------------------------------------------------------#
def timer(function):
    """ A decorator that makes the decorated *function* return its execution time. """
    def wrapper(*args, **kwargs):
        t1 = time.time()
        result = function(*args, **kwargs)
        t2 = time.time()
        print "Execution time of function", function.__name__, ":", str(t2-t1), "s."
        return result
    return wrapper

#-------------------------------------------------------------------------#
def results_to_json(lims, exid):
    """Create a JSON string describing the results of execution *exid*.

    The execution is sought in *lims*, and all its output files and
    their descriptions are written to the string.
    """
    produced_file_ids = lims.search_files(source=('execution',exid))
    d = dict([(lims.fetch_file(i)['description'], lims.path_to_file(i))
              for i in produced_file_ids])
    j = json.dumps(d)
    return j

#-------------------------------------------------------------------------#
def rstring(len=20):
    """Generate a random string of length *len* (usually for filenames).
    Equivalent to bein's unique_filename_in(), without requiring the import.
    """
    import string, random
    return "".join([random.choice(string.letters+string.digits) for x in range(len)])

#-------------------------------------------------------------------------#

def writecols(file, cols, header=None, sep="\t"):
    """Write a list of iterables *cols* as columns in a *sep*-delimited text file.
    One can precise an array *header* to define column names.
    The parameter *sep* defines the delimiter character between columns.
    """
    ncols = len(cols)
    with open(file,"wb") as f:
        w = csv.writer(f, delimiter=sep)
        if header:
            if len(header) < len(cols):
                header = header + ["c"+str(i+1+len(header)) for i in range(ncols-len(header))]
                print >>sys.stderr, "Warning: header has less elements than there are columns."
            elif len(header) > len(cols):
                header = header[:ncols]
                print >>sys.stderr, "Warning: header has more elements than there are columns."
            w.writerow(header)
        maxlen = max([len(c) for c in cols])
        for c in cols:
            if len(c) < maxlen:
                c.extend( [""]*(maxlen-len(c)) )
        for row in zip(*cols):
            w.writerow(row)

def readcols(filename, header=False, sep="\t", skip=0):
    """Read a *sep*-delimited text file *filename* and stores it
    in a dictionary, which keys are column headers if present, and values
    are the columns of the file.
    If *header*=True, the first line is interpreted as the header.
    If *header* is an array, its elements become the column headers.
    One can skip the first *skip* lines of the file.
    """
    with open(filename,"rb") as f:
        r = csv.reader(f, delimiter=sep)
        for i in range(skip):
            f.next()
        if header:
            if not getattr(header, '__iter__', False): # if not iterable
                header = r.next() # then the first line is read as the header
        rows = [row for row in r]
        if not header:
            header = range(1,len(rows[0])+1)
        elif len(header) != len(rows[0]):
            header = header + [i+1+len(header) for i in range(len(rows[0])-len(header))]
            print >>sys.stderr, "Warning: header and columns don't have the same number of elements."
        parsedfile = dict(zip(header,zip(*rows)))
    return parsedfile

#-------------------------------------------------------------------------#

def isnum(s):
    """Return True if string *s* represents a number, False otherwise"""
    try:
        float(s)
        return True
    except ValueError:
        return False

#-------------------------------------------------------------------------#

def unique(seq, idfun=None):
    """
    Return all unique elements in *seq*, preserving order - unlike list(set(seq)),
    and almost as fast. Permits this sort of filter: unique(seq, lambda x: x.lower())
    """
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result



###############################################################################
###############################################################################
##################### BIG TRY STATEMENT STARTING... ###########################
###############################################################################
###############################################################################
try:
    from bein import program
    #-------------------------------------------------------------------------#
    @program
    def join_pdf(files):
        """Uses 'ghostscript' to join several pdf files into one.
        """
        out = rstring()
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
        archive = rstring()
        call    = None
        if compression_type == "lxzma" or compression_type == "xz":
            call = ["tar", "cJf", archive, path]
        elif compression_type == "bzip2" or compression_type == "bz2":
            call = ["tar", "cjf", archive, path]
        elif compression_type == "gunzip" or compression_type == "gz":
            call = ["tar", "czf", archive, path]
        else:
            raise ValueError("Compression type: %s not yet supported!" %(compression_type))
        return {"arguments": call, "return_value": archive}

    #-------------------------------------------------------------------------#
    @program
    def uncompress(path):
        """
        uncompress tar archive
        """
        output = rstring()
        call = ["tar", "xvf", path]
        return {"arguments": call, "return_value": output}

    #-------------------------------------------------------------------------#
    @program
    def scp(source, destination, args = None):
        if args is None:
            args = []
        output = rstring()
        call = ["scp"] + args + [ source, destination ]
        return {"arguments": call, "return_value": output}

    #-------------------------------------------------------------------------#
    @program
    def merge_two_bed(file1,file2):
        """Binds ``intersectBed`` from the 'BedTools' suite.
        """
        return {"arguments": ['intersectBed','-a',file1,'-b',file2], "return_value": None}

    #-------------------------------------------------------------------------#
    @program
    def run_gMiner( job ):
        import pickle
        job_file = rstring()
        with open(job_file,'w') as f:
            pickle.dump(job,f)
        def get_output_files(p):
            with open(job_file,'r') as f:
                job = pickle.load(f)
            return job['job_output']
        return {"arguments": ["run_gminer.py",job_file],
                "return_value": get_output_files}

except:
    print >>sys.stderr, "Bein not found.  Skipping some common functions."

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

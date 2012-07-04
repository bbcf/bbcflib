"""
======================
Module: bbcflib.common
======================

Utility functions common to several pipelines.
"""

# Built-in modules #
import os, sys, time, csv, string, random, pickle, re, json

#-------------------------------------------------------------------------#
def unique_filename_in(path=None):
    """Return a random filename unique in the given path.

    The filename returned is twenty alphanumeric characters which are
    not already serving as a filename in *path*.  If *path* is
    omitted, it defaults to the current working directory.
    """
    if path == None: path = os.getcwd()
    def random_string():
        return "".join([random.choice(string.letters + string.digits) for x in range(20)])
    while True:
        filename = random_string()
        files = [f for f in os.listdir(path) if f.startswith(filename)]
        if files == []: break
    return filename

#-------------------------------------------------------------------------#
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
    if not(url.startswith(("http://","https://","ftp://"))):
        url = "http://" + url
    return url.strip("/")

#-------------------------------------------------------------------------#
def cat(files,out=None,skip=0):
    """Concatenates files.
    """
    if len(files) > 1:
        out = out or unique_filename_in()
        if os.path.exists(out): mode = "a"
        else: mode = "w"
        with open(out,mode) as f1:
            for inf in files:
                with open(inf,"r") as f2:
                    for i in range(skip):
                        f2.readline()
                    [f1.write(l) for l in f2]
    elif len(files) == 1:
        out = files[0]
    return out

#-------------------------------------------------------------------------#
def track_header( descr, ftype, url, ffile ):
    header = "track type="+ftype+" name="+re.sub(r'\[.*','',descr)
    url += ffile
    header += " bigDataUrl="+url+" "
    style = " visibility=4 "
    if ftype=="bigWig":
        style = " visibility=2 windowingFunction=maximum"
        if re.search(r'_rev.bw ',descr): style+= " color=0,10,200"
        if re.search(r'_fwd.bw ',descr): style+= " color=200,10,0"
    return header+style+"\n"

#"track type="+type+" name='"+fname+"' bigDataUrl="+htsurl+style+"\n"

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
        for i in range(skip): f.next()
        if header: header = r.next()
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

def unique(seq, fun=None):
    """
    Return all unique elements in *seq*, preserving order - unlike list(set(seq)),
    and almost as fast. Permits this sort of filter: unique(seq, lambda x: x.lower())
    """
    if fun is None:
        def fun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = fun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result


##############################
#-------- PROGRAMS --------###
##############################


try:
    from bein import program

    @program
    def gzipfile(files,args=None):
        """Runs gzip on files."""
        if not(isinstance(files,list)):
            files=[files]
        gzcall = ['gzip']
        if args:
            gzcall += args
        gzcall += files
        return {"arguments": gzcall, "return_value": None}

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
    def coverageBed(file1,file2):
        """Binds ``coverageBed`` from the 'BedTools' suite.
        """
        return {"arguments": ['coverageBed','-a',file1,'-b',file2], "return_value": None}

        #-------------------------------------------------------------------------#
    @program
    def intersectBed(file1,file2):
        """Binds ``intersectBed`` from the 'BedTools' suite.
        """
        return {"arguments": ['intersectBed','-a',file1,'-b',file2], "return_value": None}

        #-------------------------------------------------------------------------#
    @program
    def track_convert( fileIn, fileOut=None, formatIn=None, formatOut=None,
                       assembly=None, chrmeta=None ):
        if fileOut is None: fileOut = unique_filename_in()
        args = [fileIn,fileOut]
        if formatIn:  args += ["--format1",str(formatIn)]
        if formatOut: args += ["--format2",str(formatOut)]
        if assembly:  args += ["--assembly",str(assembly)]
        if chrmeta:   args += ["--chrmeta",json.dumps(chrmeta)]
        return {"arguments": ["track_convert"]+args, "return_value": fileOut}
        #-------------------------------------------------------------------------#
    @program
    def gMiner_run( job ):
        gminer_args = ["--"+str(k)+"="+str(v) for k,v in job.iteritems()]
        def _split_output(p):
            files = ''.join(p.stdout).strip().split(',')
            return files
        return {"arguments": ["gminer_run"]+gminer_args,
                "return_value": _split_output}

        #-------------------------------------------------------------------------#
    def merge_sql( ex, sqls, outdir=None, datatype='quantitative', via='lsf' ):
        """Run ``gMiner``'s 'merge_scores' function on a set of sql files
        """
        if outdir == None:
            outdir = unique_filename_in()
        if not(os.path.exists(outdir)):
            os.mkdir(outdir)
        if len(sqls) == 1:
            return sqls[0]
        gMiner_job = {"operation": "merge_scores", "output": outdir,
                      "datatype": datatype,
                      "args": "'"+json.dumps({"trackList":",".join(sqls)})+"'"}
        files = gMiner_run.nonblocking(ex,gMiner_job,via=via).wait()
        return files[0]

        #-------------------------------------------------------------------------#
    def intersect_many_bed(ex,files,via='lsf'):
        """Runs ``intersectBed`` iteratively over a list of bed files.
        """
        out = files[0]
        for f in files[1:]:
            next = unique_filename_in()
            intersectBed.nonblocking( ex, out, f, via=via, stdout=next ).wait()
            out = next
        return out

    def set_file_descr(filename,**kwargs):
        """Implements file naming compatible with the 'get_files' function::

            >>> set_file_descr("toto",**{'tag':'pdf','step':1,'type':'doc','comment':'ahaha'})
            'pdf:toto[step:1,type:doc] (ahaha)'

        if 'tag' and/or comment are ommitted::

            >>> set_file_descr("toto",step=1,type='doc')
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
    def get_files( id_or_key, minilims, by_type=True, select_param=None ):
        """Retrieves a dictionary of files created by an htsstation job identified by its key
        or bein id in a MiniLIMS.

        The dictionary keys are the file types (e.g. 'pdf', 'bam', 'py' for python objects),
        the values are dictionaries with keys repository file names and values actual file
        descriptions (names to provide in the user interface).

        'select_param' can be used to select a subset of files: if it is a string or a list of strings,
        then only files containing these parameters will be returned,
        and if it is a dictionary, only files with parameters matching the key/value pairs will be returned.
        """
        if isinstance(select_param,str): select_param = [select_param]
        if isinstance(id_or_key, str):
            try:
                exid = max(minilims.search_executions(with_text=id_or_key))
            except ValueError, v:
                raise ValueError("No execution with key "+id_or_key)
        else:
            exid = id_or_key
        file_dict = {}
        allf = dict((y['repository_name'],y['description']) for y in
                    [minilims.fetch_file(x) for x in minilims.search_files(source=('execution',exid))])
        for f,d in allf.iteritems():
            tag_d = d.split(':')
            tag = None
            if len(tag_d)>1 and not(re.search("\[",tag_d[0])):
                tag = tag_d[0]
                d = ":".join(tag_d[1:])
            pars_patt = re.search(r'\[(.*)\]',d)
            if not(pars_patt):
                pars = "group:0,step:0,type:none,view:admin"
                re.sub(r'([^\s\[]*)',r'\1['+pars+']',d,1)
            else:
                pars = pars_patt.groups()[0]
            par_dict = dict(x.split(":") for x in pars.split(","))
            if select_param:
                if isinstance(select_param,dict):
                    if not(all([k in par_dict and par_dict[k]==v for k,v in select_param.iteritems()])):
                        continue
                else:
                    if not(all([k in par_dict for k in select_param])):
                        continue
            cat = (by_type and par_dict.get('type')) or tag or 'none'
            if cat in file_dict:
                file_dict[cat][f]=d
            else:
                file_dict[cat] = {f: d}
        return file_dict

        #-------------------------------------------------------------------------#
except ImportError: print "Warning: no module named 'bein'. @program imports skipped."

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

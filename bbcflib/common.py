"""
======================
Module: bbcflib.common
======================
"""

import pickle

###############################################################################
def normalize_url(url):
    """Produce a fixed form for an HTTP URL.

    Make sure the URL begins with http:// and does *not* end with a /.

    >>> normalize_url('http://www.google.com')
    'http://www.google.com'
    >>> normalize_url('http://www.google.com/')
    'http://www.google.com'
    >>> normalize_url('www.google.com/')
    'http://www.google.com'
    """
    url = url.lower()
    if not(url.startswith("http://")):
        url = "http://" + url
    if url.endswith("/"):
        url = url[:-1]
    return url

def natural_sort(item):
    """
    Will sort strings that contain numbers correctly

    >>> l = ['v1.3.12', 'v1.3.3', 'v1.2.5', 'v1.2.15', 'v1.2.3', 'v1.2.1']
    >>> l.sort(key=natural_sort)
    >>> l.__repr__()
    "['v1.2.1', 'v1.2.3', 'v1.2.5', 'v1.2.15', 'v1.3.3', 'v1.3.12']"
    """
    import re
    def try_int(s):
        try: return int(s)
        except: return s
    return map(try_int, re.findall(r'(\d+|\D+)', item))

###############################################################################
def named_temporary_path(suffix=''):
    ''' Often, one needs a new random and temporary file path
        instead of the random and temporary file object provided
        by the tempfile module'''
    import os, tempfile
    file = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
    path = file.name
    file.close()
    os.remove(path)
    return path

###############################################################################
class memoized_method(object):
    """
    Decorator that caches a function's return value the first time
    it is called. If called later, the cached value is returned, and
    not re-evaluated

    This class is meant to be used as a decorator of methods. The return value
    from a given method invocation will be cached on the instance whose method
    was invoked. All arguments passed to a method decorated with memoize must
    be hashable.

    If a memoized method is invoked directly on its class the result will not
    be cached. Instead the method will be invoked like a static method::

        class Obj(object):
            @memoized_method
            def add_to(self, arg):
                return self + arg
        Obj.add_to(1)    # not enough arguments
        Obj.add_to(1, 2) # returns 3, result is not cached

    """
    def __init__(self, func):
        self.func = func
    def __get__(self, obj, objtype=None):
        from functools import partial
        if obj is None:
            return self.func
        return partial(self, obj)
    def __call__(self, *args, **kw):
        obj = args[0]
        try:
            cache = obj.__cache
        except AttributeError:
            cache = obj.__cache = {}
        key = (self.func, args[1:], frozenset(kw.items()))
        try:
            res = cache[key]
        except KeyError:
            res = cache[key] = self.func(*args, **kw)
        return res

def sentinelize(iterable, sentinel):
    '''
    Add an item to the end of an iterable

    >>> list(sentinelize(range(4), 99))
    [0, 1, 2, 3, 99]
    '''
    for item in iterable: yield item
    yield sentinel

###############################################################################
terminal_colors = {
    'end':    '\033[0m',    # Text Reset
    'blink':  '\033[5m',    # Blink
    'txtblk': '\033[0;30m', # Black - Regular
    'txtred': '\033[0;31m', # Red
    'txtgrn': '\033[0;32m', # Green
    'txtylw': '\033[0;33m', # Yellow
    'txtblu': '\033[0;34m', # Blue
    'txtpur': '\033[0;35m', # Purple
    'txtcyn': '\033[0;36m', # Cyan
    'txtwht': '\033[0;37m', # White
    'bldblk': '\033[1;30m', # Black - Bold
    'bldred': '\033[1;31m', # Red
    'bldgrn': '\033[1;32m', # Green
    'bldylw': '\033[1;33m', # Yellow
    'bldblu': '\033[1;34m', # Blue
    'bldpur': '\033[1;35m', # Purple
    'bldcyn': '\033[1;36m', # Cyan
    'bldwht': '\033[1;37m', # White
    'unkblk': '\033[4;30m', # Black - Underline
    'undred': '\033[4;31m', # Red
    'undgrn': '\033[4;32m', # Green
    'undylw': '\033[4;33m', # Yellow
    'undblu': '\033[4;34m', # Blue
    'undpur': '\033[4;35m', # Purple
    'undcyn': '\033[4;36m', # Cyan
    'undwht': '\033[4;37m', # White
    'bnkblk': '\033[5;30m', # Black - Blinking
    'bnkred': '\033[5;31m', # Red
    'bnkgrn': '\033[5;32m', # Green
    'bnkylw': '\033[5;33m', # Yellow
    'bnkblu': '\033[5;34m', # Blue
    'bnkpur': '\033[5;35m', # Purple
    'bnkcyn': '\033[5;36m', # Cyan
    'bnkwht': '\033[5;37m', # White
    'bakblk': '\033[40m',   # Black - Background
    'bakred': '\033[41m',   # Red
    'bakgrn': '\033[42m',   # Green
    'bakylw': '\033[43m',   # Yellow
    'bakblu': '\033[44m',   # Blue
    'bakpur': '\033[45m',   # Purple
    'bakcyn': '\033[46m',   # Cyan
    'bakwht': '\033[47m',   # White
}

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

def create_sql_track( sql_name, chromosomes, datatype="quantitative" ):
    conn = sqlite3.connect( sql_name )
    conn.execute('CREATE TABLE chrNames (name text, length integer)')
    conn.execute('CREATE TABLE attributes (key text, value text)')
    conn.execute('INSERT INTO attributes (key,value) VALUES ("datatype","%s")' %datatype)
    vals = [(v['name'],str(v['length'])) for v in chromosomes.values()]
    conn.executemany('INSERT INTO chrNames (name,length) VALUES (?,?)',vals)
    if datatype=="quantitative":
        [conn.execute('CREATE TABLE "'+v['name']+'" (start integer, end integer, score real)')
         for v in chromosomes.values()]
        [conn.execute('CREATE INDEX "range_idx_'+v['name']+'" ON "'+v['name']+'"(start, end)')
         for v in chromosomes.values()]
    elif datatype=="qualitative":
        [conn.execute('CREATE TABLE "'+v['name']+'" (start integer,end integer,score real,name text,strand integer,attributes text)')
         for v in chromosomes.values()]
        [conn.execute('CREATE INDEX "name_idx_'+ v['name'] +'" ON "'+ v['name'] +'"(name)')
         for v in chromosomes.values()]
    else:
        raise ValueError("Supported datatypes are 'quantitative' and 'qualitative', not "+datatype)
    conn.commit()
    conn.close()
    return sql_name

try:
    from bein import MiniLIMS, unique_filename_in, ProgramOutput, program, execution

###############################################################################
    def get_files( id_or_key, minilims ):
        """Retrieves a dictionary of files created by htsstation job identified by its key or bein id in a MiniLIMS.

        The dictionarie's keys are the file types (e.g. 'pdf', 'bam', 'py' for python objects), the values are dictionaries with keys repository file names and values actual file descriptions (names to provide in the user interface).
        """
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

############ gMiner operations ############
    @program
    def run_gMiner( job ):
        job_file = unique_filename_in()
        with open(job_file,'w') as f:
            pickle.dump(job,f)
        def get_output_files(p):
            with open(job_file,'r') as f:
                job = pickle.load(f)
            return job['job_output']
        return {"arguments": ["run_gminer.py",job_file],
                "return_value": get_output_files}

    def merge_sql( ex, sqls, names, description="merged.sql", outdir=None, via='lsf' ):
        """Run ``gMiner``'s 'merge_score' function on a set of sql files
        """
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

############################################################

    @program
    def merge_two_bed(file1,file2):
        """Binds ``intersectBed`` from the 'BedTools' suite.
        """
        return {"arguments": ['intersectBed','-a',file1,'-b',file2], "return_value": None}

    def merge_many_bed(ex,files,via='lsf'):
        """Runs ``intersectBed`` iteratively over a list of bed files.
        """
        out = files[0]
        for f in files[1:]:
            next = unique_filename_in()
            _ = merge_two_bed.nonblocking( ex, out, f, via=via, stdout=next ).wait()
            out = next
        return out

    @program
    def join_pdf(files):
        """Uses 'ghostscript' to join several pdf files into one.
        """
        out = unique_filename_in()
        gs_args = ['gs','-dBATCH','-dNOPAUSE','-q','-sDEVICE=pdfwrite',
                   '-sOutputFile=%s'%out]
        gs_args += files
        return {"arguments": gs_args, "return_value": out}

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

    @program
    def uncompress(path):
        """
        uncompress tar archive
        """
        output = unique_filename_in()
        call = ["tar", "xvf", path]
        return {"arguments": call, "return_value": output}

    @program
    def scp(source, destination, user, host):
        output = unique_filename_in()
        call = ["scp", source, user+"@"+host, destination]
        return {"arguments": call, "return_value": output}

except:
    print >>sys.stderr, "Bein not found.  Skipping some common functions."


#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

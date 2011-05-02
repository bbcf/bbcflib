"""
===============
bbcflib.common
===============

"""

import pickle
from bein import *

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
    '''
    Will sort strings that contain numbers correctly

    >>> l = ['v1.3.12', 'v1.3.3', 'v1.2.5', 'v1.2.15', 'v1.2.3', 'v1.2.1']
    >>> l.sort(key=natural_sort)
    >>> l.__repr__()
    "['v1.2.1', 'v1.2.3', 'v1.2.5', 'v1.2.15', 'v1.3.3', 'v1.3.12']"
    '''
    import re
    def try_int(s):
        try: return int(s)
        except: return s
    return map(try_int, re.findall(r'(\d+|\D+)', item))

###############################################################################
def named_temporary_path(suffix=''):
    ''' Often, one needs a new random and temporary file path
        instead of the random and tempory file object provided
        by the tempfile module'''
    import os, tempfile
    file = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
    path = file.name
    file.close()
    os.remove(path)
    return path

###############################################################################
class memoize_once(object):
    '''Decorator that caches a function's return value the first time
    it is called. If called later, the cached value is returned, and
    not re-evaluated'''
    def __init__ (self, f):
        self.f = f
        self.mem = []
    def __call__ (self, *args, **kwargs):
        if len(self.mem) == 0:
            self.mem.append(self.f(*args, **kwargs))
        return self.mem[0]

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
        name = re.sub(r'\s+\(BAM index\)','.bai',name)
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
    conn.execute('create table chrNames (name text, length integer)')
    conn.execute('create table attributes (key text, value text)')
    conn.execute('insert into attributes (key,value) values ("datatype","%s")' %datatype)
    vals = [(v['name'],str(v['length'])) for v in chromosomes.values()]
    conn.executemany('insert into chrNames (name,length) values (?,?)',vals)
    if datatype=="quantitative":
        [conn.execute('create table "'+v['name']+'" (start integer, end integer, score real)') 
         for v in chromosomes.values()]
        [conn.execute('create index "range_idx_'+v['name']+'" on "'+v['name']+'" (start, end)') 
         for v in chromosomes.values()]
    elif datatype=="qualitative":
        [conn.execute('create table "'+v['name']+'" (start integer,end integer,score real,name text,strand integer,attributes text)') 
         for v in chromosomes.values()]
        [conn.execute('create index "name_idx_'+v['name']+'" on" "'+v['name']+'" (name)') 
         for v in chromosomes.values()] 
    else:
        raise ValueError("Supported datatypes are 'quantitative' and 'qualitative', not "+datatype)
    conn.commit()
    conn.close()
    return sql_name

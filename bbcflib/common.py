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

###############################################################################
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


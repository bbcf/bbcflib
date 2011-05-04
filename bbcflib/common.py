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
    '''Decorator that caches a function's return value the first time
    it is called. If called later, the cached value is returned, and
    not re-evaluated
    
    This class is meant to be used as a decorator of methods. The return value
    from a given method invocation will be cached on the instance whose method
    was invoked. All arguments passed to a method decorated with memoize must
    be hashable.
    
    If a memoized method is invoked directly on its class the result will not
    be cached. Instead the method will be invoked like a static method:
    class Obj(object):
        @memoized_method
        def add_to(self, arg):
            return self + arg
    Obj.add_to(1)    # not enough arguments
    Obj.add_to(1, 2) # returns 3, result is not cached
    '''
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

###############################################################################
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

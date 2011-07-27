"""
===============================
Submodule: bbcflib.track.common
===============================

Common stuff for the track package.
"""

###########################################################################
def check_path(path):
    """
    Raises an exception if the path *path* is already taken.
    """
    import os
    if os.path.exists(path): raise Exception("The location '" + path + "' is already taken")

###########################################################################
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
    """
    Often, one needs a new random and temporary file path
    instead of the random and temporary file object provided
    by the tempfile module
    """
    import os, tempfile
    file = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
    path = file.name
    file.close()
    os.remove(path)
    return path

###############################################################################
def sentinelize(iterable, sentinel):
    """
    Add an item to the end of an iterable

    >>> list(sentinelize(range(4), 99))
    [0, 1, 2, 3, 99]
    """
    for item in iterable: yield item
    yield sentinel

###############################################################################
def sqlcmp(file_a, file_b):
    """
    Compare two two sqlite3 databases via their dumps
    """
    import itertools, sqlite3
    A = sqlite3.connect(file_a)
    B = sqlite3.connect(file_b)
    for a,b in itertools.izip_longest(A.iterdump(), B.iterdump()):
        if a != b: return False
    return True

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
class ModifiedDict(object):
    ''' A dictionary like object that tracks modification
        and has a special boolean self.modified value'''

    def __init__(self, d=None):
        self.data = {}
        if d is not None: self.update(d)
        self.modified = False

    def __call__(self, *args, **kwargs):
        self.modified = True
        super(ModifiedDict, self)(*args, **kwargs)

    def __repr__(self):
        return repr(self.data)

    def __cmp__(self, d):
        if isinstance(d, ModifiedDict): return cmp(self.data, d.data)
        else: return cmp(self.data, d)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key):
        if key in self.data: return self.data[key]
        if hasattr(self.__class__, "__missing__"): return self.__class__.__missing__(self, key)
        raise KeyError(key)

    def __contains__(self, key):
        return key in self.data

    def __setitem__(self, key, item):
        self.modified = True
        self.data[key] = item

    def __delitem__(self, key):
        self.modified = True
        del self.data[key]

    def __iter__(self):
        return iter(self.data)

    def keys(self): return self.data.keys()
    def items(self): return self.data.items()
    def iteritems(self): return self.data.iteritems()
    def iterkeys(self): return self.data.iterkeys()
    def itervalues(self): return self.data.itervalues()
    def values(self): return self.data.values()
    def has_key(self, key): return key in self.data

    def get(self, key, failobj=None):
        if key not in self: return failobj
        return self[key]

    def clear(self):
        self.modified = True
        self.data.clear()

    def setdefault(self, key, failobj=None):
        if key not in self:
            self.modified = True
            self[key] = failobj
        return self[key]

    def pop(self, key, *args):
        self.modified = True
        return self.data.pop(key, *args)

    def popitem(self):
        self.modified = True
        return self.data.popitem()

    def copy(self):
        if self.__class__ is ModifiedDict: return ModifiedDict(self.data.copy())
        import copy
        data = self.data
        try:
            self.data = {}
            c = copy.copy(self)
        finally:
            self.data = data
        c.update(self)
        return c

    def update(self, d=None, **kwargs):
        if d is None: return
        self.modified = True
        if isinstance(d, ModifiedDict): self.data.update(d.data)
        elif isinstance(d, type({})) or not hasattr(d, 'items'): self.data.update(d)
        else:
            for k, v in d.items(): self[k] = v
        if len(kwargs): self.data.update(kwargs)

    @classmethod
    def fromkeys(cls, iterable, value=None):
        d = cls()
        for key in iterable: d[key] = value
        return d

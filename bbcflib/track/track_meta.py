"""
=======================================
Subpackage: bbcflib.track.track_chrmeta
=======================================

Takes care of managing chromosome meta data.
"""

# Built-in modules #
import os

# Internal modules #
from ..genrep import GenRep
from .common import named_temporary_path, ModifiedDict

###########################################################################
class ChromMetaData(ModifiedDict):
    def __call__(self, x):
        self.token = x
        self.modified = True
        self.data  = self.read_token(x)

    def read_token(self, x):
        # Null case #
        if not x: return {}
        # Dictionary case #
        elif isinstance(x, dict): return x
        # ModifiedDictionary case #
        elif isinstance(x, ModifiedDict): return x.data
        # Tuple case #
        elif isinstance(x, tuple): return self.read_tuple(x)
        # String case #
        elif isinstance(x, basestring):
        # File case #
            if os.path.exists(x) and not os.path.isdir(x): return self.read_file(x)
        # GenRep case #
            else: return self.read_genrep(x)

    #--------------------------------------------------------------------------#
    def read_tuple(self, t):
        column_names, rows = t
        rows = [dict([(k,r[i]) for i, k in enumerate(column_names)]) for r in rows]
        return dict([(r['name'], dict([(k, r[k]) for k in column_names if k != 'name'])) for r in rows])

    def read_genrep(self, name):
        g = GenRep()
        if g.is_down(): raise Exception("The assembly '" + name + "' could not be read because the Genrep server is down.")
        if not g.is_available(name): raise Exception("The Genrep server does not know about the assembly '" + name + "'.")
        return g.assembly(name).chrmeta

    def read_file(self, path):
        result = {}
        with open(path, 'r') as f:
            for line in f:
                line = line.strip('\n')
                if len(line) == 0:       continue
                if line.startswith("#"): continue
                if line.endswith(" \\"):
                    raise Exception("The file '" + path + "' includes linebreaks ('\\') which are not supported.")
                if '\t' in line: seperator = '\t'
                else:            seperator = ' '
                line = line.split(seperator)
                if len(line) != 2:
                    raise Exception("The file " + path + " does not seam to be a valid chromosome file.")
                name = line[0]
                try:
                    length = int(line[1])
                except ValueError:
                    raise Exception("The file '" + path + "' has invalid values.")
                result[name] = dict([('length', length)])
        if not result:
            raise Exception("The file '" + path + "' does not seam to contain any information.")
        return result

    def write_file(self, path=None):
        if not path: path = named_temporary_path()
        if os.path.exists(path): raise Exception("The location '" + path + "' is already taken")
        def lines():
            for k,v in self.items(): yield k + '\t' + str(v['length']) + '\n'
        with open(path, 'w') as f: f.writelines(lines())
        return path

    #--------------------------------------------------------------------------#
    @property
    def column_names(self):
        return ['name'] + self.keys()

    @property
    def rows(self):
        return [dict([['name', chrom]] + [(k,v) for k,v in self[chrom].items()]) for chrom in self]

    #--------------------------------------------------------------------------#
    def choose_max(self, chrmeta):
        for chrom in self.data:
            if chrom in chrmeta:
                if chrmeta[chrom] > self.data[chrom]:
                    self.data[chrom] = chrmeta[chrom]

###########################################################################
class TrackMetaData(ModifiedDict):
    def __call__(self, x):
        self.token = x
        self.modified = True
        self.data  = self.read_token(x)

    def read_token(self, x):
        if not x: return {}
        elif isinstance(x, dict): return x
        elif isinstance(x, ModifiedDict): return x.data

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

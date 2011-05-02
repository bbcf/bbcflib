"""
===================================
Submodule: bbcflib.track.format_wig
===================================

Implementation of the WIG format.
"""

from .track_proxy import ProxyTrack
from .track_text import TextTrack

###########################################################################
class GenomicFormat(ProxyTrack, TextTrack):
    @property
    def _type(self):
        return 'quantitative' 
   
    @property
    def _fields(self):
        return self.default_fields

    def _read(self):
        def all_features():
            params = {}
            for line in self.file:
                line = line.strip("\n").lstrip()
                if len(line) == 0:              continue
                if line.startswith("#"):        continue
                if line.startswith("track "):   continue
                if line.startswith("browser "): continue
                if line.endswith(" \\"):
                    raise Exception("The track " + self.location + " includes linebreaks ('\\') which are not supported.")
                if line.startswith("variableStep") or line.startswith("fixedStep"):
                    params = dict([p.split('=',1) for p in shlex.split('mode=' + line)])
                    if not params.get('chrom',False):
                        raise Exception("The track " + self.location + " does not specify a chromosome and is hence not valid.")
                    try:
                        params['span'] = int(params.get('span',1))
                    except ValueError:
                        raise Exception("The track " + self.location + " has a non integer as span value.")
                    if params['span'] < 1:
                        raise Exception("The track " + self.location + " has a negative or null span value.")
                    if line.startswith("fixedStep "):
                        if not 'start' in params:
                            raise Exception("The track " + self.location + " has a fixedStep directive without a start.")
                        try:
                            params['start'] = int(params['start'])
                        except ValueError:
                            raise Exception("The track " + self.location + " has a non integer as start value.")
                        try:
                            params['step'] = int(params.get('step',1))
                        except ValueError:
                            raise Exception("The track " + self.location + " has a non integer as step value.")
                        if params['step'] < 1:
                            raise Exception("The track " + self.location + " has a negative or null step value.")
                        params['count'] = 0
                    continue
                if not params:
                    raise Exception("The track " + self.location + " is missing a fixedStep or variableStep directive.")
                if params['mode'] == 'fixedStep':
                    try:
                        line = float(line)
                    except ValueError:
                        raise Exception("The track " + self.location + " has non floats as score values and is hence not valid.")
                    base = params['start'] + params['count'] * params['step'] 
                    yield [params['chrom'], base, base + params['span'], line]
                    params['count'] += 1
                if params['mode'] == 'variableStep':
                    line = line.split()
                    try:
                        line[0] = int(line[0])
                        line[1] = float(line[1])
                    except ValueError:
                        raise Exception("The track " + self.location + " has invalid values.")
                    except IndexError:
                        raise Exception("The track " + self.location + " has missing values.")
                    yield [params['chrom'], line[0], line[0] + params['span'], line[1]]
        def all_entries():
            sentinel = ('', sys.maxint, sys.maxint, 0.0)
            X = gm_com.sentinelize(all_features(), sentinel)
            x = X.next()
            if x == sentinel:
                yield x
                return 
            while True:
                x_next = X.next()
                if x_next == sentinel:
                    if x[3] != 0.0: yield tuple(x)
                    break
                if x[0] == x_next[0]:
                    if x[2] > x_next[1]:
                        raise Exception("The track " + self.location + " has a start larger than its end or a span larger than its step.") 
                    if x[2] == x_next[1] and x[3] == x_next[3]:
                        x[2] = x_next[2]
                        continue
                if x[3] != 0.0: yield tuple(x)
                x = x_next 

        global chr, entry, generator
        chr             = ''
        self._seen_chr  = []
        entry           = ['', '', '', '']
        generator       = all_entries()
        def get_next_entry():
            global entry, generator
            entry = generator.next()    
        def iter_until_different_chr():
            global chr, entry
            while True:
                if entry[0] != chr: break
                yield entry[1:]
                get_next_entry()
        get_next_entry()
        while True:
            if entry[0] == chr: break
            chr = entry[0]
            if chr in self._seen_chr:
                raise Exception("The track " + self.location + " is not sorted by chromosomes (" + chr + ").")
            if not chr in self.all_chrs:
                raise Exception("The track " + self.location + " has a value (" + chr + ") not specified in the chromosome file.")
            self._seen_chr.append(chr)
            yield chr, iter_until_different_chr()

    def dump(self):
        raise NotImplementedError

###########################################################################
def create(path):
    open(path, 'w').close()

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

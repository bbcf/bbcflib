"""
===================================
Submodule: bbcflib.track.format_wig
===================================

Implementation of the WIG format.
"""

# General modules #
import sys, shlex

# Specific module #
from .track_proxy import ProxyTrack
from .track_text import TextTrack
from ..common import sentinelize

###########################################################################
class GenomicFormat(TextTrack, ProxyTrack):
    def _all_features(self):
        self._file.seek(0)
        params = {}
        seen_track = False
        for line in self._file:
            line = line.strip("\n").lstrip()
            if len(line) == 0:              continue
            if line.startswith("#"):        continue
            if line.endswith(" \\"):
                raise Exception("The file '" + self._path + "' includes linebreaks ('\\') which are not supported.")
            if line.startswith("track "):
                if not seen_track:
                    seen_track = True
                    continue
                raise Exception("The file '" + self._path + "' contains a second 'track' directive. This is not supported.")
            if line.startswith("browser "): continue
            if line.startswith("variableStep") or line.startswith("fixedStep"):
                params = dict([p.split('=',1) for p in shlex.split('mode=' + line)])
                if not params.get('chrom',False):
                    raise Exception("The file '" + self._path + "' does not specify a chromosome and is hence not valid.")
                try:
                    params['span'] = int(params.get('span',1))
                except ValueError:
                    raise Exception("The file '" + self._path + "' has a non integer as span value.")
                if params['span'] < 1:
                    raise Exception("The file '" + self._path + "' has a negative or null span value.")
                if line.startswith("fixedStep "):
                    if not 'start' in params:
                        raise Exception("The file '" + self._path + "' has a fixedStep directive without a start.")
                    try:
                        params['start'] = int(params['start'])
                    except ValueError:
                        raise Exception("The file '" + self._path + "' has a non integer as start value.")
                    try:
                        params['step'] = int(params.get('step',1))
                    except ValueError:
                        raise Exception("The file '" + self._path + "' has a non integer as step value.")
                    if params['step'] < 1:
                        raise Exception("The file '" + self._path + "' has a negative or null step value.")
                    params['count'] = 0
                continue
            if not params:
                raise Exception("The file '" + self._path + "' is missing a fixedStep or variableStep directive.")
            if params['mode'] == 'fixedStep':
                try:
                    line = float(line)
                except ValueError:
                    raise Exception("The file '" + self._path + "' has non floats as score values and is hence not valid.")
                base = params['start'] + params['count'] * params['step'] 
                yield [params['chrom'], base, base + params['span'], line]
                params['count'] += 1
            if params['mode'] == 'variableStep':
                line = line.split()
                try:
                    line[0] = int(line[0])
                    line[1] = float(line[1])
                except ValueError:
                    raise Exception("The file '" + self._path + "' has invalid values.")
                except IndexError:
                    raise Exception("The file '" + self._path + "' has missing values.")
                yield [params['chrom'], line[0], line[0] + params['span'], line[1]]

    def _all_entries(self):
        # TODO: What if datatype is 'qualitative' and the wig has overlaps ?
        sentinel = ('', sys.maxint, sys.maxint, 0.0)
        X = sentinelize(self._all_features(), sentinel)
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
                    raise Exception("The file '" + self._path + "' has a start larger than its end or a span larger than its step.") 
                if x[2] == x_next[1] and x[3] == x_next[3]:
                    x[2] = x_next[2]
                    continue
            if x[3] != 0.0: yield tuple(x)
            x = x_next

    def _write(self):
        yield self._header_line
        for f in self.read(): yield "fixedStep chrom=%s start=%s span=%s\n%s\n" % (f[0],f[1],f[2]-f[1],f[3])

    #-----------------------------------------------------------------------------#
    @property
    def _datatype(self): 
        return self._datatype_

    @_datatype.setter
    def _datatype(self, datatype):
        if not datatype: datatype = 'quantitative'
        self._datatype_ = datatype

###########################################################################
def random_track(kind='fixed', number_of_features=32, size=100, range=1000, jump=10000, small_jump=200, orig_start=0, chrs=16):
    import random
    chr = 0
    if kind == 'fixed':
        yield 'track type=wiggle_0 name="Pol2 Signal" description="Chip-Seq" source="Random generator"\n'
        for i in xrange(number_of_features):
            if i % (number_of_features / chrs) == 0:
                chr += 1
                end  = orig_start
            start = end   + (random.randint(0,jump))
            end   = start + (random.randint(1,size))
            multiplier = random.randint(1,range) 
            sys.stdout.write('fixedStep chrom=chr' + str(chr) + ' start=' + str(start) + ' step=1' + '\n')
            for x in xrange(start, end):
                sys.stdout.write(str(multiplier + multiplier * random.random()) + '\n')
            for x in xrange(start, end):
                sys.stdout.write('0\n')
            start = end   + small_jump + (random.randint(0,jump))
            end   = start +              (random.randint(1,size))
            multiplier = random.randint(1,range) 
            for x in xrange(start,end):
                yield str(multiplier + multiplier * random.random()) + '\n'
    if kind == 'variable':
        yield 'track type=wiggle_0 name="Rap1 Peaks" description="Chip-Seq" source="Random generator"\n'
        for i in xrange(number_of_features*2):
            if i % ((number_of_features*2) / chrs) == 0:
                chr += 1
                end  = orig_start
                sys.stdout.write('variableStep chrom=chr' + str(chr) + '\n')
            start = end   + (random.randint(0,jump))
            end   = start + (random.randint(1,int(size/2)))
            multiplier = random.randint(1,range) 
            for x in xrange(start,end):
                yield str(x) + ' ' + str(multiplier + multiplier * random.random()) + '\n'

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

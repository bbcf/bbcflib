"""
=====================================
Submodule: bbcflib.track.formats.rand
=====================================


Implementation of a random track generator.
"""

# Built-in modules #
import sys, random, tempfile

# Internal modules #
from .. import Track

# Variables #
chrsuffix = 'Awfully super extra long chromosome denomination string '

###########################################################################
class TrackFormat(Track):
    def __init__(self, *args):
        # Essential parameters #
        self.all_chrs = [chrsuffix + str(x) for x in range(10)]
        # Other parameters #
        self.name_gen = tempfile._RandomNameSequence()
        self.name_gen.rng.seed(0)
        self.size = 500

    def read(self, selection=None, fields=None):
        if type(selection) != str: raise Exception(0, "You can't specify a region on a random track")
        if fields: raise Exception(0, "You can't specify fields on a random track")
        start = 0
        for feat in range(int(self.size + 4*self.size*random.random())):
            start = start + (random.randint(0,100))
            end = start + (random.randint(1,100) )
            score = random.gammavariate(1, 0.1) * 1000
            strand = map(lambda x: x==1 and 1 or -1, [random.randint(0,1)])[0]
            yield [start, end, self.name_gen.next(), score, strand]

    @property
    def datatype(self):
        return 'qualitative'

    @property
    def fields(self):
        return Track.qualitative_fields

    @property
    def chrmeta(self):
        return dict([(chr, dict([('length', sys.maxint)])) for chr in self.all_chrs])

    @property
    def attributes(self):
        return {'type':'random'}

    def unload(self, *args):
        pass

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

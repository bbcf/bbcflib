# General Modules #
import os

# Specific Modules #
from ..track import Track, new

# Specific Variables #
from .test_variables import track_collections, yeast_chr_file

# Unittesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

###################################################################################
class Test_Read(unittest.TestCase):
    def runTest(self):
        t = track_collections['Scores']['1']
        with Track(t['path'], chrfile=yeast_chr_file) as t['track']:
            # Just the first feature #
            data = t['track'].read()
            self.assertEqual(data.next(), ('chr1', 0, 10, 8.0))
            # Number of features #
            data = t['track'].read()
            self.assertEqual(len(list(data)), 4)

#-----------------------------------------------------------------------------#   
class Test_Write(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bed')
        with new(path, chrfile=yeast_chr_file) as t:
            features = {}
            features['chr1'] = [(10, 20, 'Lorem', 1.0, 1),
                                (30, 40, 'Ipsum', 2.0, 1)]
            features['chr2'] = [(10, 20, 'Dolor', 3.0, 1)]
            for chrom, data in sorted(features.items()):
                t.write(chrom, data)
        with open(path,'r') as f: A = f.read().split('\n')
        with open(track_collections['Validation']['4']['path'],'r') as f: B = f.read().split('\n')
        self.assertEqual(A[-4:], B[-4:])

Test_Write().runTest()

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

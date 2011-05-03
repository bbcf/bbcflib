# General Modules #
import os

# Specific Modules #
from ..track import Track, new
from ..common import named_temporary_path

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
        t = track_collections['Validation']['1']
        with Track(t['path'], chrfile=t['chrfile']) as t['track']:
            # Just the first feature #
            data = t['track'].read()
            self.assertEqual(data.next(), ('chr1', 0, 10, 'Validation feature 1', 10.0))
            # Number of features #
            data = t['track'].read()
            self.assertEqual(len(list(data)), 12)

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
        with open(path,                                        'r') as f: A = f.read().split('\n')
        with open(track_collections['Validation']['4']['path'],'r') as f: B = f.read().split('\n')
        self.assertEqual(A[1:], B)
        os.remove(path)

#-----------------------------------------------------------------------------#   
class Test_Roundtrips(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bed')
        for track_num, track in sorted(track_collections['Validation'].items()):
            with Track(track['path'], chrfile=track['chrfile']) as t:
                t.dump(path)
            with open(path,         'r') as f: A = f.read().split('\n')
            with open(track['path'],'r') as f: B = f.read().split('\n')
            self.assertEqual(A[1:], B)
            os.remove(path)

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

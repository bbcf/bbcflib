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
        t = track_collections['Signals'][1]
        with Track(t['path'], chromosomes_data=t['chromosomes_data']) as t['track']:
            # Just the first feature #
            data = t['track'].read()
            self.assertEqual(data.next(), ('chr1', 0, 10, -1.0))
            # Number of features #
            data = t['track'].read()
            self.assertEqual(len(list(data)), 3)

#-----------------------------------------------------------------------------#
class Test_Write(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bedGraph')
        with new(path, chromosomes_data=yeast_chr_file) as t:
            features = {}
            features['chr1'] = [(0,  10, -1.0),
                                (20, 30, -1.75),
                                (40, 50, -2.5)]
            for chrom, data in sorted(features.items()):
                t.write(chrom, data)
        with open(path,                                   'r') as f: A = f.read().split('\n')
        with open(track_collections['Signals'][1]['path'],'r') as f: B = f.read().split('\n')
        self.assertEqual(A[1:], B)
        os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Roundtrips(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bedGraph')
        for track_num, track in sorted(track_collections['Signals'].items()):
            with Track(track['path'], chromosomes_data=track['chromosomes_data']) as t:
                t.dump(path)
            with open(path,         'r') as f: A = f.read().split('\n')
            with open(track['path'],'r') as f: B = f.read().split('\n')
            self.assertEqual(A[1:], B)
            os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Format(unittest.TestCase):
    def runTest(self):
        t = track_collections['Signals'][1]
        with Track(t['path'], chrfile=t['chrfile']) as t:
            self.assertEqual(t.format, 'sql')
            self.assertEqual(t._format, 'bedGraph')

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

# Built-in modules #
import os, shutil

# Internal modules #
from ... import track
from ..common import named_temporary_path
from ..track_collection import track_collections

# Unittesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

###################################################################################
class Test_Read(unittest.TestCase):
    def runTest(self):
        t = track_collections['Signals'][1]
        with track.load(t['path']) as t['track']:
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
        with track.new(path) as t:
            self.assertEqual(t.datatype, 'quantitative')
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
        for track_num, d in sorted(track_collections['Signals'].items()):
            with track.load(d['path']) as t:
                t.dump(path)
            with open(path,              'r') as f: A = f.read().split('\n')
            with open(d['path'],'r') as f: B = f.read().split('\n')
            self.assertEqual(A[1:], B)
            os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Format(unittest.TestCase):
    def runTest(self):
        # Not specified #
        t = track_collections['Signals'][1]
        with track.load(t['path']) as t:
            self.assertEqual(t.format, 'bedGraph')
        # No extension #
        old = track_collections['Signals'][1]['path']
        new = named_temporary_path()
        shutil.copyfile(old, new)
        with track.load(new, 'bedGraph') as t:
            self.assertEqual(t.format, 'bedGraph')
        os.remove(new)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

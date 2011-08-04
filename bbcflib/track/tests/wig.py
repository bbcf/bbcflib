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
        t = track_collections['Scores'][1]
        with track.load(t['path']) as t['track']:
            # Just the first feature #
            data = t['track'].read()
            self.assertEqual(data.next(), ('chr1', 0, 10, 8.0))
            # Number of features #
            data = t['track'].read()
            self.assertEqual(len(list(data)), 4)

#-----------------------------------------------------------------------------#
class Test_Write(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.wig')
        with track.new(path) as t:
            features = {}
            features['chr1'] = [(20,   50,  20),
                                (50,   80, 300),
                                (120, 130,   0)]
            for chrom, data in sorted(features.items()):
                t.write(chrom, data)
        with open(path,                                  'r') as f: A = f.read().split('\n')
        with open(track_collections['Scores'][3]['path'],'r') as f: B = f.read().split('\n')
        self.assertEqual(A[1:], B)
        os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Roundtrips(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.wig')
        for track_num, d in sorted(track_collections['Scores'].items()):
            if track_num == 3: continue
            with track.load(d['path']) as t:
                t.dump(path)
            with open(path,              'r') as f: A = f.read().split('\n')
            with open(d['path'],'r') as f: B = f.read().split('\n')
            self.assertEqual(A[1:], B)
            os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Format(unittest.TestCase):
    def runTest(self):
        t = track_collections['Scores'][1]
        with track.load(t['path']) as t:
            self.assertEqual(t.format, 'wig')
            self.assertEqual(t.type_identifier, 'wiggle_0')

#------------------------------------------------------------------------------#
class Test_Format(unittest.TestCase):
    def runTest(self):
        # Not specified #
        t = track_collections['Scores'][1]
        with track.load(t['path']) as t:
            self.assertEqual(t.format, 'wig')
            self.assertEqual(t.type_identifier, 'wiggle_0')
        # No extension #
        old = track_collections['Scores'][1]['path']
        new = named_temporary_path()
        shutil.copyfile(old, new)
        with track.load(new, 'wig') as t:
            self.assertEqual(t.format, 'wig')
        os.remove(new)
        # Only track line #
        old = track_collections['Peaks']['Rap1']['path']
        new = named_temporary_path()
        shutil.copyfile(old, new)
        with track.load(new) as t:
            self.assertEqual(t.format, 'wig')
        os.remove(new)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

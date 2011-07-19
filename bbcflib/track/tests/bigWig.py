# Built-in modules #
import os

# Internal modules #
from ... import track
from ..common import named_temporary_path
from ..track_collection import track_collections, yeast_chr_file

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
        t = track_collections['Binary'][1]
        with track.load(t['path']) as t:
            # Just the first feature #
            data = t.read()
            self.assertEqual(data.next(), ('chr1', 0, 10, -1.0))
            # Number of features #
            data = t.read()
            self.assertEqual(len(list(data)), 3)

#-----------------------------------------------------------------------------#
class Test_Write(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bigWig')
        with track.new(path, chrmeta=yeast_chr_file) as t:
            self.assertEqual(t.datatype, 'quantitative')
            features = {}
            features['chr1'] = [(0,  10, -1.0),
                                (20, 30, -1.75),
                                (40, 50, -2.5)]
            for chrom, data in sorted(features.items()): t.write(chrom, data)
        import filecmp
        print path
        print  track_collections['Binary'][1]['path']
        self.assertTrue(filecmp.cmp(path, track_collections['Binary'][1]['path']))
        #os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Overwrite(unittest.TestCase):
    def runTest(self):
        pass

#-----------------------------------------------------------------------------#
class Test_Roundtrips(unittest.TestCase):
    def runTest(self):
        pass

#-----------------------------------------------------------------------------#
class Test_Format(unittest.TestCase):
    def runTest(self):
        pass

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

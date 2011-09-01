# Built-in modules #
import sys

# Internal modules #
from bbcflib import track
from bbcflib.track.track_collection import track_collections

# Unittesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = False

################################################################################
class Test_Read(unittest.TestCase):
    def runTest(self):
        if sys.version_info < (2, 7): self.skipTest("Gzip support only works in 2.7 for the moment")
        path = track_collections['Compressed'][1]['path']
        with track.load(path) as t:
            # Just the first feature #
            data = t['track'].read()
            self.assertEqual(data.next(), ('chr1', 0, 10, 'Validation feature 1', 10.0, 0))
            # Number of features #
            data = t['track'].read()
            self.assertEqual(len(list(data)), 12)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

# Genreral Modules #
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

Test_Read().runTest()

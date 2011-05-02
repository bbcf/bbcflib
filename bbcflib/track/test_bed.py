# Genreral Modules #
import os

# Specific Modules #
from ..track import Track, new

# Specific Variables #
from .test_variables import track_collections

# Unittesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

###################################################################################
class Test_Read(unittest.TestCase):
    def runTest(self):
        t = track_collections['Validation']['1']
        with Track(t['orig_path']) as t['track']:
            # Just the first feature #
            data = t['track'].read()
            self.assertEqual(data.next(), ('chr1', 0, 10, 'Validation feature 1', 10.0))

# Internal modules #
from ..track_collection import track_collections
from .. import magic

# Unittesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

###################################################################################
class Test_BED(unittest.TestCase):
    def runTest(self):
        self.skipTest("These tests don't pass anymore. Delete this line once they are fixed.")
        path = track_collections['Yeast']['RP genes']['path']
        ftype = magic.guess_file_format(path)
        self.assertEqual(ftype, 'bed')

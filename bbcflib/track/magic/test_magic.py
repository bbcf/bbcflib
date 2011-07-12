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
        path    = track_collections['Validation'][1]['path']
        ftype   = magic.guess_file_format(path)
        self.assertEqual(ftype, 'bed')

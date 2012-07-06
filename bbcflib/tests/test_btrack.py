# Built-in modules #
import cStringIO

# Internal modules #
from bbcflib.genrep import Assembly
from bbcflib.btrack import *

# Unitesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

# Path to testing files
path = "test_data/btrack/"


class Test_Btrack(unittest.TestCase):
    def setUp(self):
        pass
    def test_(self):
        pass

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

    def test_track(self):
        pass

    def test_convert(self):
        pass

    def test_concat_fields(self):
        pass

    def test_map_chromosomes(self):
        pass

    def test_split_field(self):
        pass

    def test_score_threshold(self):
        pass

    def test_map_chromosomes(self):
        pass

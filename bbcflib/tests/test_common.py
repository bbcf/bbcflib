# Internal modules #
from bbcflib.common import *

# Built-in modules #
import os

# Unitesting module #
try:
    import unittest2 as unittest
    assert unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

#Path to testing files
path = "test_data/common/"


class Test_Common(unittest.TestCase):
    def setUp(self):
        pass

    def test_unique_filename_in(self):
        fn = unique_filename_in(path)
        f = open(os.path.join(path,fn),'wb') # create a new file
        self.assertIn(fn,os.listdir(path))
        f.close(); os.remove(os.path.join(path,fn))

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#



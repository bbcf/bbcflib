# Unitesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Tested module #
from .. import rnaseq

# Other modules #
import os

#--------------------------------------------------------------------------------#
class TestWorkflow(unittest.TestCase):
    def runTest(self):
        self.skipTest("No tests yet")

# Unitesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Tested module #
from ..rnaseq.rnaseq import rnaseq_workflow

# Other modules #
import re, sys, os, random

#--------------------------------------------------------------------------------#
class TestWorkflow(unittest.TestCase):
    def runTest(self):
        pass

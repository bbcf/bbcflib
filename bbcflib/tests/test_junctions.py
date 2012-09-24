
# Built-in modules #

# Internal modules #
from bbcflib.rnaseq import *
from bein import execution

# Unitesting modules #
try: import unittest2 as unittest
except ImportError: import unittest
from numpy.testing import assert_almost_equal

# Other modules #
from numpy import array

# Nosetest flag #
__test__ = True

# Local test:
# run_unctions.py -v local -c /archive/epfl/bbcf/jdelafon/test_rnaseq/config/junc.txt -d junctions

class fakejob(object):
    def __init__(self,assembly):
        self.assembly = genrep.Assembly(assembly)
        self.groups = {1:{1:{}}}


class Test_Junctions(unittest.TestCase):
    def setUp(self):
        self.assembly = genrep.Assembly('hg19')
        self.index = "/archive/epfl/bbcf/jdelafon/soapsplice_index/hg19.index"
        self.unmapped = "/archive/epfl/bbcf/jdelafon/test_junctions/unmapped.fastq"
        self.bam_files = {1:{1:{'unmapped_fastq':'unmapped.fastq'}}}

    def test_junctions_workflow(self):
        with execution(None) as ex:
            junctions_workflow(ex,job=1,bam_files=self.bam_files,index=self.index)


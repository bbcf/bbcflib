
# Built-in modules #

# Internal modules #
from bbcflib.junctions import *
from bbcflib.frontend import Job
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

def fakejob(assembly):
        job = Job(1,'today','key',assembly.id,'test','email',{})
        job.add_group(id=1,name='group1')
        job.add_run(group_id=1, id=1)
        return job


class Test_Junctions(unittest.TestCase):
    def setUp(self):
        self.assembly = genrep.Assembly('hg19')
        self.index = "/archive/epfl/bbcf/jdelafon/soapsplice_index/hg19.index"
        self.path = "/archive/epfl/bbcf/jdelafon/bin/SOAPsplice-v1.9/bin/soapsplice"
        if not os.path.exists(self.path):
            self.path = "/Users/delafont/bin/SOAPsplice-v1.9/bin/soapsplice"
        unmapped = "test_data/junctions/junc_reads_chr1-100k"
        self.bam_files = {1:{1:{'unmapped_fastq':(unmapped+'_R1',unmapped+'_R2')}}}
        self.job = fakejob(self.assembly)

    def test_junctions_workflow(self):
        with execution(None) as ex:
            junctions_workflow(ex,job=self.job,bam_files=self.bam_files,index=self.index,
                               path_to_soapsplice=self.path)


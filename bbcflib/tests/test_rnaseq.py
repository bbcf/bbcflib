
# Built-in modules #

# Internal modules #
from bbcflib.rnaseq import *
from bbcflib import genrep
from bbcflib.frontend import Job
from bein import execution

# Unitesting modules #
try:
    import unittest2 as unittest
    assert unittest
except ImportError:
    import unittest


# Nosetest flag #
__test__ = True

# Local test:
# run_rnaseq.py -v local -c /archive/epfl/bbcf/jdelafon/test_rnaseq/config/gapkowt.txt -d rnaseq -p genes,transcripts

testpath = os.path.abspath("test_data/rnaseq/")

def fakejob(assembly):
        job = Job(1,'today','key',assembly.id,'test','email',{})
        job.add_group(id=1,name='group1')
        job.add_run(group_id=1, id=1)
        return job


class Test_Counting(unittest.TestCase):
    def setUp(self):
        self.assembly = genrep.Assembly('mm9')
        self.job = fakejob(self.assembly)

    def test_count_reads(self):
        stranded = False
        with execution(None) as ex:
            C = Counter(ex,"local",self.job,self.assembly,["KO.1","KO.2"],sys.stderr,sys.stdout,"genes",False,stranded)
            bamfiles = [os.path.join(testpath,"gapdhKO.bam")]*2
            gtf = os.path.join(testpath,"mm9_3genes_renamed.gtf")
            count_files = C.count_reads(bamfiles, gtf)
            with open(count_files["genes"]) as genes:
                self.assertEqual(len(genes.readlines()), 4)
            with open(count_files["transcripts"]) as trans:
                self.assertEqual(len(trans.readlines()), 18)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

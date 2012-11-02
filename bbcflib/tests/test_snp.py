# Internal modules #
from bbcflib.snp import *
from bbcflib import genrep

# Unitesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

#Path to testing files
path = "test_data/snp/"

# Local test:
# cd /scratch/cluster/monthly/jdelafon/snp
# run_snp.py -v local -c /archive/epfl/bbcf/jdelafon/test_snp/config/snp.config -d /scratch/cluster/monthly/jdelafon/snp/snp_minilims -w wdir -f /archive/epfl/bbcf/jdelafon/test_snp/reference/chrV.tar.gz
# run_snp.py -v local -c /archive/epfl/bbcf/jdelafon/test_snp/config/hart.config -d /scratch/cluster/monthly/jdelafon/snp/snp_minilims -w wdir


class Test_SNP(unittest.TestCase):
    def setUp(self):
        pass

    def test_annotate_snps(self):
        assembly = genrep.Assembly('sacCer2')
        filedict = {'chrV':path+'chrV'}
        self.outall, self.outexons = annotate_snps(filedict, ["s1","s2"], assembly)
        with open(self.outall,'r') as f: print "\nAll SNPs ('outall'):\n",f.read()
        with open(self.outexons,'r') as g: print "\nExonic SNPs ('outexons'):\n",g.read()
        raise IOError("Error raised voluntarily to print test outputs.")

    def tearDown(self):
        os.remove(self.outall)
        os.remove(self.outexons)


#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#



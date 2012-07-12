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
# run_snp.py -v local -c /scratch/cluster/monthly/jdelafon/snp/config/snp.config -d /scratch/cluster/monthly/jdelafon/snp/snp_minilims -w wdir -f /scratch/cluster/monthly/jdelafon/snp/reference/sacCer2/chrV.tar.gz


class Test_SNP(unittest.TestCase):
    def setUp(self):
        pass

    def test_annotate_snps(self):
        assembly = genrep.Assembly('sacCer2')
        filedict = {'chrV':path+'chrV'}
        outall, outexons = annotate_snps(filedict, ["s1","s2"], assembly)
        with open(outall,'r') as f: print '\noutall\n',f.read()
        with open(outexons,'r') as g: print '\noutexons\n',g.read()
        os.remove(outall)
        os.remove(outexons)
        raise

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#



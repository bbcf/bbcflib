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

#Path tpo testing files
path = "test_data/snp/"

class Test_SNP(unittest.TestCase):
    def setUp(self):
        self.sample_names = ["sample1"]
        self.assembly = genrep.Assembly('sacCer2')

    def test_annotate_snps(self):
        filedict = {'chrV':path+"chrV.txt", 'chrIV':path+"chrIV.txt",
                    'chrVI':path+"chrVI.txt", 'chrII':path+"chrII.txt"}
        outall, outexons = annotate_snps(filedict, self.sample_names, self.assembly)
        #with open(outall,'r') as f: print 'outall\n',f.read()
        #with open(outexons,'r') as g: print 'outexons\n',g.read()
        os.remove(outall)
        os.remove(outexons)
        #raise


#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#



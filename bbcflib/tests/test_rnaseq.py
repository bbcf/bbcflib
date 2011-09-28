# Built-in modules #
import datetime, ConfigParser, cStringIO

# Internal modules #
from bbcflib.rnaseq import *

# Unitesting modules #
try: import unittest2 as unittest
except ImportError: import unittest
from numpy.testing import assert_almost_equal, assert_equal
from numpy import array

# Nosetest flag #
__test__ = True


class Test_Rnaseq(unittest.TestCase):
    def setUp(self):
	# add some that are not in the mappings from Ensembl
        e1="e1"; e2="e2";
        t1="t1"; t2="t2";
        g1="g1";
        self.M = numpy.matrix('1 1; 0 1')
        self.counts = numpy.array([[27,12],[3,3]]) # [[cond1],[cond2]]
        self.dexons = dict(zip([e1,e2], zip(*self.counts)))
        self.gene_ids = [g1]
        self.exon_mapping = {e1:([t1,t2],g1), e2:([t2],g1)}
        self.transcript_mapping = {t1:g1, t2:g1}
        self.trans_in_gene = {g1:[t1,t2]}
        self.exons_in_trans = {t1:[e1], t2:[e1,e2]}
        """
            Cond1         Cond2
        |=====g1====| |=====g1====|
        |-t1-|        |-t1-| 
        |-----t2----| |-----t2----|
          v      v      v      v  
          27     12     3      3 
        |.e1.| |.e2.| |.e1.| |.e2.|
	"""
    def test_transcripts_expression(self):
        texp = transcripts_expression(self.gene_ids, self.transcript_mapping,
		self.trans_in_gene, self.exons_in_trans, self.dexons)
        assert_almost_equal(texp["t1"], array([15., 0.]))
        assert_almost_equal(texp["t2"], array([12., 3.]))
        
    def test_genes_expression(self):
        gexp = genes_expression(self.gene_ids, self.exon_mapping, self.dexons)
        assert_almost_equal(gexp["g1"], array([39., 6.]))
        
    def test_estimate_size_factors(self):
        res, size_factors = estimate_size_factors(self.counts)
        self.assertIsInstance(self.counts, numpy.ndarray)
        self.assertEqual(type(self.counts), type(res))
        self.assertEqual(self.counts.shape, res.shape)
        assert_almost_equal(size_factors, array([2.5, 0.41666]), decimal=3)
        assert_almost_equal(res, array([[10.8, 4.8],[7.2, 7.2]]))
        
    def test_lsqnonneg(self):
        C = numpy.array([[0.0372, 0.2869],
                         [0.6861, 0.7071],
                         [0.6233, 0.6245],
                         [0.6344, 0.6170]])
        C1 = numpy.array([[0.0372, 0.2869, 0.4],
                          [0.6861, 0.7071, 0.3],
                          [0.6233, 0.6245, 0.1],
                          [0.6344, 0.6170, 0.5]])
        C2 = numpy.array([[0.0372, 0.2869, 0.4],
                          [0.6861, 0.7071,-0.3],
                          [0.6233, 0.6245,-0.1],
                          [0.6344, 0.6170, 0.5]])
        d = numpy.array([0.8587, 0.1781, 0.0747, 0.8405])

        [x, resnorm, residual] = lsqnonneg(C, d)
        dres = abs(resnorm - 0.8315)          # compare with matlab result
        self.assertLess(dres, 0.001)

        [x, resnorm, residual] = lsqnonneg(C1, d)
        dres = abs(resnorm - 0.1477)          # compare with matlab result
        self.assertLess(dres, 0.01)

        [x, resnorm, residual] = lsqnonneg(C2, d)
        dres = abs(resnorm - 0.1027)          # compare with matlab result
        self.assertLess(dres, 0.01)

        k = numpy.array([[0.1210, 0.2319, 0.4398, 0.9342, 0.1370],
                         [0.4508, 0.2393, 0.3400, 0.2644, 0.8188],
                         [0.7159, 0.0498, 0.3142, 0.1603, 0.4302],
                         [0.8928, 0.0784, 0.3651, 0.8729, 0.8903],
                         [0.2731, 0.6408, 0.3932, 0.2379, 0.7349],
                         [0.2548, 0.1909, 0.5915, 0.6458, 0.6873],
                         [0.8656, 0.8439, 0.1197, 0.9669, 0.3461],
                         [0.2324, 0.1739, 0.0381, 0.6649, 0.1660],
                         [0.8049, 0.1708, 0.4586, 0.8704, 0.1556],
                         [0.9084, 0.9943, 0.8699, 0.0099, 0.1911]])
        k1 = k-0.5
        l = numpy.array([0.4225, 0.8560, 0.4902, 0.8159, 0.4608, 0.4574, 0.4507, 0.4122, 0.9016, 0.0056])

        [x, resnorm, residual] = lsqnonneg(k, l)
        dres = abs(resnorm - 0.3695)          # compare with matlab result
        self.assertLess(dres, 0.01)

        [x, resnorm, residual] = lsqnonneg(k1, l)
        dres = abs(resnorm - 2.8639)          # compare with matlab result
        self.assertLess(dres, 0.01)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

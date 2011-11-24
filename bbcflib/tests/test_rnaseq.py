# coding: utf-8

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


class Test_Expressions1(unittest.TestCase):
    """Two conditions, different lengths, invertible"""
    def setUp(self):
        e1="e1"; e2="e2"; t1="t1"; t2="t2"; g1="g1";
        self.ncond = 2
        self.counts = numpy.array([[27,12],[3,3]]) # [[cond1],[cond2]]
        self.rpkms = numpy.array([[27/3.,12/6.],[3/3.,3/6.]])
        self.exons_data = [[e1,e2]]+list(self.counts)+list(self.rpkms)+[[0.,3.],[3.,9.],["g1","g1"],["gg1","gg1"]]
        self.transcript_lengths = {t1:3., t2:9.}
        self.exon_lengths = {e1:3., e2:6.}
        self.exon_to_gene = {e1:g1, e2:g1}
        self.trans_in_gene = {g1:[t1,t2]}
        self.exons_in_trans = {t1:[e1], t2:[e1,e2]}
        """
           *Counts Cond1*             *Rpkm Cond1*
        g1 |===.======|               |===.===|
        t1 |---|         21           |---|      7
        t2 |---.------|  6  +  12     |---.---| (2) 2

             27    12                   9   2
           |.e1.||.e2.|               |e1| |e2|
        """
        """
           *Counts Cond2*             *Rpkm Cond2*
        g1 |===.======|               |===.===|
        t1 |---|         1.5          |---|      0.5
        t2 |---.------|  1.5 + 3      |---.---| (0.5) 0.5

             3     3                    1   0.5
           |.e1.||.e2.|               |e1| |e2|
        """
    def test_transcripts_expression(self):
        trpk, tcounts, err = transcripts_expression(self.exons_data, self.exon_lengths, self.transcript_lengths,
                 self.trans_in_gene, self.exons_in_trans, self.ncond)
        assert_almost_equal(trpk["t1"], array([7., 0.5]))
        assert_almost_equal(trpk["t2"], array([2., 0.5]))
        assert_almost_equal(tcounts["t1"], array([21., 1.5]))
        assert_almost_equal(tcounts["t2"], array([18., 4.5]))
        self.assertEqual(sum(sum(tcounts.values())), sum(sum(self.counts)))

    def test_genes_expression(self):
        gcounts, grpkms = genes_expression(self.exons_data, self.exon_to_gene, self.ncond)
        assert_almost_equal(gcounts["g1"], array([39., 6.]))


class Test_Expressions2(unittest.TestCase):
    """One condition, equal lengths, invertible"""
    def setUp(self):
        e1="e1"; e2="e2"; e3="e3"; t1="t1"; t2="t2"; t3="t3"; g1="g1";
        self.ncond = 1
        self.counts = numpy.array([[10,15,10]]) # [[cond1]]
        self.rpkms = numpy.array([[10/5.,15/5.,10/5.]])
        self.exons_data = [[e1,e2,e3]]+list(self.counts)+list(self.rpkms)+\
               [[0.,5.,10],[5.,10.,15.],["g1","g1","g1"],["gg1","gg1","gg1"]]
        self.transcript_lengths = {t1:10., t2:10., t3:15.}
        self.exon_lengths = {e1:5., e2:5., e3:5.}
        self.exon_to_gene = {e1:g1, e2:g1, e3:g1}
        self.trans_in_gene = {g1:[t1,t2,t3]}
        self.exons_in_trans = {t1:[e1,e2], t2:[e2,e3], t3:[e1,e2,e3]}
        """
           *Counts Cond1*
        g1 |===.===.===|
        t1 |-------|      5 + 5
        t2     |-------|  5 +     5
        t3 |-----------|  5 + 5 + 5

            10  15   10
           |e1||e2| |e3|
        """
    def test_transcripts_expression(self):
        trpk, tcounts, err = transcripts_expression(self.exons_data, self.exon_lengths, self.transcript_lengths,
                 self.trans_in_gene, self.exons_in_trans, self.ncond)
        assert_almost_equal(tcounts["t1"], array([10.]))
        assert_almost_equal(tcounts["t2"], array([10.]))
        assert_almost_equal(tcounts["t3"], array([15.]))
        self.assertEqual(sum(sum(tcounts.values())), sum(sum(self.counts)))


class Test_Expressions3(unittest.TestCase):
    """Underdetermined system"""
    def setUp(self):
        e1="e1"; e2="e2"; e3="e3"; t1="t1"; t2="t2"; t3="t3"; t4="t4"; g1="g1";
        self.ncond = 1
        self.counts = numpy.array([[10,15,10]]) # [[cond1]]
        self.rpkms = numpy.array([[10/5.,15/5.,10/5.]])
        self.exons_data = [[e1,e2,e3]]+list(self.counts)+list(self.rpkms)+\
               [[0.,5.,10],[5.,10.,15.],["g1"]*3,["gg1"]*3]
        self.transcript_lengths = {t1:10., t2:10., t3:15., t4:10.}
        self.exon_lengths = {e1:5., e2:5., e3:5.}
        self.exon_to_gene = {e1:g1, e2:g1, e3:g1}
        self.trans_in_gene = {g1:[t1,t2,t3,t4]}
        self.exons_in_trans = {t1:[e1,e2], t2:[e2,e3], t3:[e1,e2,e3], t4:[e1,e3]}
        """
           *Counts Cond1*
        g1 |===.===.===|
        t1 |-------|      7.5 + 7.5
        t2     |-------|        7.5 + 7.5
        t3 |-----------|  0     0     0
        t4 |---|   |---|  2.5 +       2.5

            10  15   10
           |e1||e2| |e3|
        """
    def test_transcripts_expression(self):
        trpk, tcounts, err = transcripts_expression(self.exons_data, self.exon_lengths, self.transcript_lengths,
                 self.trans_in_gene, self.exons_in_trans, self.ncond)
        self.assertEqual(sum(sum(tcounts.values())), sum(sum(self.counts)))


#@unittest.skip("fix")
class Test_Expressions4(unittest.TestCase):
    """Overdetermined system"""
    def setUp(self):
        e1="e1"; e2="e2"; e3="e3"; e4="e4"; t1="t1"; t2="t2"; t3="t3"; g1="g1";
        self.ncond = 1
        self.counts = numpy.array([[10.,10.,10.,10.]]) # [[cond1]]
        self.rpkms = numpy.array([[10/5.,10/5.,10/5.,10/5.]])
        self.exons_data = [[e1,e2,e3,e4]]+list(self.counts)+list(self.rpkms)+\
               [[0,5,10,15],[5,10,15,20],["g1"]*4,["gg1"]*4]
        self.transcript_lengths = {t1:10., t2:15., t3:10.}
        self.exon_lengths = {e1:5., e2:5., e3:5., e4:5.}
        self.exon_to_gene = {e1:g1, e2:g1, e3:g1, e4:g1}
        self.trans_in_gene = {g1:[t1,t2,t3]}
        self.exons_in_trans = {t1:[e1,e2], t2:[e2,e3,e4], t3:[e2,e4]}
        """
           *Counts Cond1*
        g1 |===.===.===.===|
        t1 |-------|          12 + 12
        t2     |-----------|       4  + 4  + 4
        t3     |---|   |---|  0    0    0    0

            10  10  10   10
           |e1||e2||e3| |e4|
        """
    def test_transcripts_expression(self):
        trpk, tcounts, err = transcripts_expression(self.exons_data, self.exon_lengths, self.transcript_lengths,
                 self.trans_in_gene, self.exons_in_trans, self.ncond)
        print tcounts
        self.assertGreater(sum(sum(tcounts.values()))/sum(sum(self.counts)), 0.9)


#@unittest.skip("fix")
class Test_Expressions5(unittest.TestCase):
    """Even more underdetermined system"""
    def setUp(self):
        e1="e1"; e2="e2"; e3="e3"; e4="e4"; e5="e5"; t1="t1"; t2="t2"; t3="t3"; g1="g1";
        self.ncond = 1
        self.counts = numpy.array([[10.,10.,10.,10.,10.]]) # [[cond1]]
        self.rpkms = numpy.array([[10/5.,10/5.,10/5.,10/5.,10/5.]])
        self.exons_data = [[e1,e2,e3,e4,e5]]+list(self.counts)+list(self.rpkms)+\
               [[0,5,10,15,20],[5,10,15,20,25],["g1"]*5,["gg1"]*5]
        self.transcript_lengths = {t1:10., t2:15., t3:10.}
        self.exon_lengths = {e1:5., e2:5., e3:5., e4:5., e5:5.}
        self.exon_to_gene = {e1:g1, e2:g1, e3:g1, e4:g1, e5:g1}
        self.trans_in_gene = {g1:[t1,t2,t3]}
        self.exons_in_trans = {t1:[e1,e2,e3], t2:[e3,e4,e5], t3:[e2,e4,e5]}
        """
           *Counts Cond1*
        g1 |===.===.===.===.===|
        t1 |-----------|          4.6 + 4.6 + 4.6
        t2         |-----------|              4.6 + 4.6 + 4.6
        t3     |---|   |-------|        3.0 +       3.0 + 3.0

            10  10  10   10  10
           |e1||e2||e3| |e4||e5|
        """
    def test_transcripts_expression(self):
        trpk, tcounts, err = transcripts_expression(self.exons_data, self.exon_lengths, self.transcript_lengths,
                 self.trans_in_gene, self.exons_in_trans, self.ncond)
        print tcounts
        self.assertGreater(sum(sum(tcounts.values()))/sum(sum(self.counts)), 0.9)


class Test_others(unittest.TestCase):
    def setUp(self):
        self.counts = numpy.array([[27,12],[3,3]]) # [[cond1],[cond2]]

    def test_estimate_size_factors(self):
        res, size_factors = estimate_size_factors(self.counts)
        self.assertIsInstance(self.counts, numpy.ndarray)
        self.assertEqual(type(self.counts), type(res))
        self.assertEqual(self.counts.shape, res.shape)
        assert_almost_equal(size_factors, array([2.5, 0.41666]), decimal=3)
        assert_almost_equal(res, array([[10.8, 4.8],[7.2, 7.2]]))

    def test_save_results(self):
        with execution(None) as ex:
            save_results(ex,cols=[[1,2,3],('a','b','c','d')],conditions=["num","char"],
                         header=["num","char"], desc="test")


class Test_NNLS(unittest.TestCase):
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


class Test_Pileup(unittest.TestCase):
    def setUp(self):
        self.gene_name = "Gapdh" # Control gene - always highly expressed
        self.gene_id = "ENSG00000111640"
        self.chr = 12
        self.start = 6643093
        self.end = 6647537
        self.exons = []
        self.counts = []
        self.bam = "test_data/Gapdh.bam"

    def test_exon_labels(self):
        self.exons = exons_labels(self.bam)
        exons = [("ENSMUSE00000569415|ENSMUSG00000057666|125115289|125115329|-1", 41),
                 ("ENSMUSE00000709315|ENSMUSG00000057666|125114615|125115329|-1", 715)]
        self.assertIn(exons[0],self.exons); self.assertIn(exons[1],self.exons)
        self.assertEqual(len(self.exons),19)

    def test_pileup_file(self):
        self.exons = exons_labels(self.bam)
        self.counts = pileup_file(self.bam, self.exons)
        counts = [0, 35, 0, 0, 0, 0, 0, 3679, 3707, 0, 0, 0, 149, 3, 0, 0, 55, 0, 161]
        self.assertItemsEqual(self.counts,counts)


#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

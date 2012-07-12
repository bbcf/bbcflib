# Built-in modules #
import cStringIO
import os, sys, math

# Internal modules #
from bbcflib import btrack, genrep
from bbcflib.btrack import FeatureStream as fstream
from bbcflib.bFlatMajor.common import sentinelize, reorder, unroll, sorted_stream, shuffled, fusion, cobble
from bbcflib.bFlatMajor.stream.annotate import getNearestFeature
from bbcflib.bFlatMajor.stream.intervals import concatenate, neighborhood, combine, segment_features
from bbcflib.bFlatMajor.stream.intervals import exclude, require, disjunction, intersection, union
from bbcflib.bFlatMajor.stream.scores import merge_scores, mean_score_by_feature, window_smoothing
from bbcflib.bFlatMajor.numeric.regions import feature_matrix, average_feature_matrix
from bbcflib.bFlatMajor.numeric.signal import normalize, correlation
#from bbcflib.bFlatMajor.figure.rplots import scatterplot, lineplot, boxplot, heatmap

# Other modules #
import numpy

# Unitesting modules #
try:
    import unittest2 as unittest
except ImportError:
    import unittest
from numpy.testing import assert_almost_equal

# Nosetest flag #
__test__ = True

# Path to testing files
path = "test_data/bFlatMajor/"

# Numpy print options #
numpy.set_printoptions(precision=3,suppress=True)


class Test_Common(unittest.TestCase):
    def setUp(self):
        self.a = genrep.Assembly('sacCer2')

    def test_reorder(self):
        stream = fstream([(10,12,0.5), (14,15,1.2)], fields=['start','end','score'])
        expected = [(12,0.5,10), (15,1.2,14)]
        rstream = list(reorder(stream,['end','score','start']))
        self.assertListEqual(rstream,expected)

    def test_unroll(self):
        stream = fstream([(10,12,0.5), (14,15,1.2)], fields=['start','end','score'])
        expected = [(0,),(0.5,),(0.5,),(0,),(0,),(1.2,),(0,)]
        ustream = list(unroll(stream,(9,16)))
        self.assertListEqual(ustream, expected)

        stream = fstream([(0,1,5,'n'),(1,2,9,'n'),(2,3,11,'n')], fields=['start','end','score'])
        expected = [(5,'n'),(9,'n'),(11,'n')]
        ustream = list(unroll(stream,(0,3)))
        self.assertListEqual(ustream, expected)

    def test_sorted_stream(self):
        stream = fstream([(10,0.8),(15,2.8),(12,19.5),(12,1.4),(13,0.1)], fields=['start','score'])
        sstream = list(sorted_stream(stream,fields=['start']))
        expected = [(10,0.8),(12,19.5),(12,1.4),(13,0.1),(15,2.8)]
        self.assertListEqual(sstream,expected)

        stream = fstream([(10,0.8),(15,2.8),(12,1.4),(13,0.1),(12,19.5)], fields=['start','score'])
        sstream = list(sorted_stream(stream,fields=['start','score']))
        expected = [(10,0.8),(12,1.4),(12,19.5),(13,0.1),(15,2.8)]
        self.assertListEqual(sstream,expected)

        stream = fstream([('chrX',0,1,0.8),('chrIX',3,5,2.8),('chrIX',3,9,1.4),('chrIX',2,10,0.1),('chrIX',7,10,0.8)],
                         fields=['chr','start','end','score'])
        sstream = list(sorted_stream(stream, fields=['start','chr']))
        expected = [('chrX',0,1,0.8),('chrIX',2,10,0.1),('chrIX',3,5,2.8),('chrIX',3,9,1.4),('chrIX',7,10,0.8)]
        self.assertListEqual(sstream,expected)

        stream = fstream([('chrX',0,1,0.8),('chrIX',3,5,2.8),('chrIX',3,9,1.4),('chrIX',2,10,0.1),('chrIX',7,10,0.8)],
                         fields=['chr','start','end','score'])
        sstream = list(sorted_stream(stream, fields=['chr','start','score'], chrnames=self.a.chrnames))
        expected = [('chrIX',2,10,0.1),('chrIX',3,9,1.4),('chrIX',3,5,2.8),('chrIX',7,10,0.8),('chrX',0,1,0.8)]
        self.assertListEqual(sstream,expected)

    def test_shuffled(self):
        pass

    def test_fusion(self):
        stream = fstream([('chr1',10,15,'A',1),('chr1',13,18,'B',-1),('chr1',18,25,'C',-1)],
                         fields = ['chr','start','end','name','score'])
        expected = [(10, 18, 'chr1', 'A|B', 0),(18, 25, 'chr1', 'C', -1)]
        fused = list(fusion(stream))
        self.assertEqual(fused,expected)

    def test_cobble(self):
        stream = fstream([('chr1',10,15,'A',1),('chr1',13,18,'B',-1),('chr1',18,25,'C',-1)],
                         fields = ['chr','start','end','name','score'])
        expected = [(10, 13, 'chr1', 'A',   1),
                    (13, 15, 'chr1', 'A|B', 0),
                    (15, 18, 'chr1', 'B',  -1),
                    (18, 25, 'chr1', 'C',  -1)]
        cobbled = list(cobble(stream))
        self.assertEqual(cobbled,expected)


################### STREAM ######################


class Test_Annotate(unittest.TestCase):
    def setUp(self):
        self.assembly = genrep.Assembly('ce6')
        """
        ----- 14,795,328 ---- 14,798,158 - 14,798,396 ---- 14,800,829 -----
              |                            |
               ->     Y54E2A.11             ->     Y54E2A.12
        """
    def test_getNearestFeature(self):
        features = fstream([('chrII',14795327,14798367)], fields=['chr','start','end'])
        expected = [(14795327, 14798367, 'chrII', 'Y54E2A.12|tbc-20_Y54E2A.11|eif-3.B', 'Promot_Included', '28_0')]
        annotations = self.assembly.gene_track(chromlist=['chrII'])
        result = list(getNearestFeature(features,annotations))
        self.assertItemsEqual(result,expected)


class Test_Intervals(unittest.TestCase):
    def setUp(self):
        pass

    def test_concatenate(self):
        pass

    def test_neighborhood(self):
        pass

    def test_combine(self):
        pass

    def test_segment_features(self):
        pass

    def test_exclude(self):
        pass

    def test_disjunction(self):
        pass

    @unittest.skip("fix the deal with fusion and cobble first")
    def test_intersection(self):
        expected = (91143, 91144,'chr', ('C','*A','0','|EBMYCG00000002479|Rv0083',1,0))
        a = genrep.Assembly('mycoTube_H37RV')
        c = btrack.concat_fields(a.annot_track('CDS','chr'),infields=['name','strand','frame'], as_tuple=True)
        feat = btrack.FeatureStream(iter([('chr',91143,91144,('C','*A','0'))]), fields=['chr','start','end','rest'])
        g = combine([feat,c], intersection, win_size=10000)
        self.assertEqual(g.next(),expected)

    def test_union(self):
        pass


class Test_Scores(unittest.TestCase):
    def setUp(self):
        pass

    def test_merge_scores(self):
        pass

    def test_mean_score_by_feature(self):
        pass

    def test_window_smoothing(self):
        pass


################### NUMERIC ######################


class Test_Regions(unittest.TestCase):
    def setUp(self):
        pass

    def test_feature_matrix(self):
        pass

    def test_average_feature_matrix(self):
        pass


class Test_Signal(unittest.TestCase):
    def setUp(self):
        pass

    def test_normalize(self):
        x = [1,2,3,4,5] # mean=15/5=3, var=(1/5)*(4+1+0+1+4)=2
        assert_almost_equal(normalize(x), numpy.array([-2,-1,0,1,2])*(1/math.sqrt(2)))

    #@unittest.skip("")
    def test_correlation(self):
        numpy.set_printoptions(precision=3,suppress=True)
        # Create 2 vectors of scores, zero everywhere except a random position
        N = 10
        x = numpy.zeros(N)
        y = numpy.zeros(N)
        xpeak = numpy.random.randint(0,N)
        ypeak = numpy.random.randint(0,N)
        x[xpeak] = 10
        y[ypeak] = 10
        x = (x-numpy.mean(x))/numpy.std(x)
        y = (y-numpy.mean(y))/numpy.std(y)

        # Make tracks out of them and compute cross-correlation with our own function
        X = [('chr',k,k+1,s) for k,s in enumerate(x)]
        Y = [('chr',k,k+1,s) for k,s in enumerate(y)]
        X = btrack.FeatureStream(iter(X),fields=['chr','start','end','score'])
        Y = btrack.FeatureStream(iter(Y),fields=['chr','start','end','score'])
        corr = correlation([X,Y], regions=(0,N))#, limits=[-N+1,N-1])

        # Compute cross-correlation "by hand" and using numpy.correlate(mode='valid')
        raw = []
        np_corr_valid = []
        for k in range(N):
            """
            X         |- - - - -|          k=0
            Y              <- |- - - - -|
            up to
            X         |- - - - -|          k=4
            Y         |- - - - -|
            """
            raw.append(numpy.dot(x[-k-1:],y[:k+1]) / N)
            np_corr_valid.extend(numpy.correlate(x[-k-1:],y[:k+1],mode='valid'))
        for k in range(N-1,0,-1):
            """
            X         |- - - - -|          k=4
            Y    <- |- - - - -|
            up to
            X         |- - - - -|          k=1
            Y |- - - - -|
            """
            raw.append(numpy.dot(x[:k],y[-k:]) / N)
            np_corr_valid.extend(numpy.correlate(x[:k],y[-k:],mode='valid'))

        # Compute cross-correlation using numpy.correlate(mode='full')
        np_corr_full = numpy.correlate(x,y,mode="full")[::-1] / N
        np_corr_valid = numpy.asarray(np_corr_valid) / N

        # Test if all methods yield the same result
        assert_almost_equal(corr, numpy.asarray(raw))
        assert_almost_equal(corr, np_corr_full)
        assert_almost_equal(corr, np_corr_valid)

        # Test if the lag between the two tracks is correcty detected
        self.assertEqual(numpy.argmax(corr)-(N-1), ypeak-xpeak)


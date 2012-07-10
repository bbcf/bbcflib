# Built-in modules #
import cStringIO
import os, sys, math

# Internal modules #
from bbcflib import btrack, genrep
from bbcflib.btrack import FeatureStream as fstream
from bbcflib.bFlatMajor.common import sentinelize, reorder, unroll, sorted_stream, shuffled, fusion
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
        pass

    def test_reorder(self):
        expected = [(12,0.5,10), (15,1.2,14)]
        stream = fstream([(10,12,0.5), (14,15,1.2)], fields=['start','end','score'])
        rstream = list(reorder(stream,['end','score','start']))
        self.assertListEqual(rstream,expected)

    def test_unroll(self):
        expected = [(0,),(0.5,),(0.5,),(0,),(0,),(1.2,),(0,)]
        stream = fstream([(10,12,0.5), (14,15,1.2)], fields=['start','end','score'])
        ustream = list(unroll(stream,9,16))
        self.assertListEqual(ustream, expected)

        expected = [(5,),(9,),(11,)]
        stream = fstream([(0,1,5),(1,2,9),(2,3,11)], fields=['start','end','score'])
        ustream = list(unroll(stream,0,3))
        self.assertListEqual(ustream, expected)

    def test_sorted_stream(self):
        pass

    def test_shuffled(self):
        pass

    def test_fusion(self):
        pass



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

    def test_intersection(self):
        pass

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
        corr = correlation([X,Y], start=0, end=N)#, limits=[-N+1,N-1])

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
        np_corr_full = numpy.correlate(x,y,mode="full") / N
        np_corr_valid = numpy.asarray(np_corr_valid[::-1]) / N
        raw = numpy.asarray(raw[::-1])

        # Test if all methods yield the same result
        assert_almost_equal(corr, raw)
        assert_almost_equal(corr, np_corr_full)
        assert_almost_equal(corr, np_corr_valid)

        # Test if the lag between the two tracks is correcty detected
        self.assertEqual(numpy.argmax(corr)-(N-1), xpeak-ypeak)

        #print 'x = ',x
        #print 'y = ',y
        #print 'corr', corr
        #print 'raw ', raw
        #print 'lag      :',xpeak-ypeak
        #print 'Py idx   :',numpy.argmax(corr)
        #print 'human idx:',numpy.argmax(corr)+1
        #raise


################### FIGURE ######################






# Built-in modules #
import cStringIO
import os, sys

# Internal modules #
from bbcflib import btrack
from bbcflib.genrep import Assembly
from bbcflib.bFlatMajor.stream.annotate import getNearestFeature
from bbcflib.bFlatMajor.stream.intervals import concatenate, neighborhood, combine, segment_features
from bbcflib.bFlatMajor.stream.intervals import exclude, require, disjunction, intersection, union
from bbcflib.bFlatMajor.stream.scores import merge_scores, mean_score_by_feature, window_smoothing
from bbcflib.bFlatMajor.numeric.regions import feature_matrix, average_feature_matrix
from bbcflib.bFlatMajor.numeric.signal import normalize, correlation
#from bbcflib.bFlatMajor.figure.rplots import scatterplot, lineplot, boxplot, heatmap

# Other modules #
import numpy
from numpy.random import randint

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

################### STREAM ######################


class Test_Annotate(unittest.TestCase):
    def setUp(self):
        pass

    def test_getNearestFeature(self):
        pass

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
        pass

    def test_correlation(self):
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
        Q1 = 'track1.bed'
        Q2 = 'track2.bed'
        if os.path.exists(Q1): os.remove(Q1)
        if os.path.exists(Q2): os.remove(Q2)
        X = btrack.FeatureStream(iter(x),fields=['chr','start','end','score'])
        Y = btrack.FeatureStream(iter(y),fields=['chr','start','end','score'])
        with btrack.track(Q1, assembly='sacCer2', fields=['chr','start','end','score']) as q1:
            for k in range(N):
                q1.write(X,fields=['chr','start','end','score'])
                #q1.write('chrI', [(k,k+1,'',x[k],0)])
        with btrack.track(Q2, assembly='sacCer2', fields=['chr','start','end','score']) as q2:
            for k in range(N):
                q2.write(Y,fields=['chr','start','end','score'])
                #q2.write('chrI', [(k,k+1,'',y[k],0)])
        fig, corr = correlation()(Q1, Q2, 'chrI')

        # Compute cross-correlation "by hand" and using numpy.correlate(mode='valid')
        raw = []; np_corr_valid = []
        for k in range(N):
            """
            X         |- - - - -|
            Y              <- |- - - - -|
            up to
            X         |- - - - -|
            Y         |- - - - -|
            """
            raw.append(numpy.dot(x[-k-1:],y[:k+1]))
            np_corr_valid.extend(numpy.correlate(x[-k-1:],y[:k+1],mode='valid'))
        for k in range(N-1,0,-1):
            """
            X         |- - - - -|
            Y      <- |- - - - -|
            up to
            X         |- - - - -|
            Y |- - - - -|
            """
            raw.append(numpy.dot(x[:k],y[-k:]))
            np_corr_valid.extend(numpy.correlate(x[:k],y[-k:],mode='valid'))

        # Compute cross-correlation using numpy.correlate(mode='full')
        np_corr_full = numpy.correlate(x,y,mode="full")[::-1]

        # Test if all methods yield the same result
        assert_almost_equal(corr, numpy.array(raw))
        assert_almost_equal(corr, np_corr_full)
        # Test if the lag between the two tracks is correcty detected
        self.assertEqual(numpy.argmax(corr)-(N-1), ypeak-xpeak)


################### FIGURE ######################






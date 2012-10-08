from bbcflib.bFlatMajor.stream import score_by_feature, segment_features
from bbcflib.bFlatMajor.common import sorted_stream
import numpy

def feature_matrix(trackScores,trackFeatures,segment=False,fn='mean',**kw):
    """
    Return an array with as many lines as there are features in *trackFeatures*, and as many columns
    as there are score tracks in *trackScores*. Each element in the matrix thus corresponds to the
    (average) score of some genomic feature.

    If *segment* is True, each feature will be segmented into bins using
    bbcflib.bFlatMajor.stream.intervals.segment_features (additional parameters in *\*\*kw* will be passed to this function).
    Then each element of the array is itself an array with *nbin* lines and one column for each track in *trackScores*.

    If *segment* is False, then each element of the array is an array with one element for each track in *trackScores*.

    Example::

                      gene1                 gene2
        X: -----#####|#####|#####--------###|###|###-----  (features)
        Y: _____________666|66666________666|666|666_____  (scores1)
        Z: _____22222|22222|22222________________________  (scores2)

        With segment=True, nbins=3:

              Y   Z
        R: [[[0.  2.],    # bin0 - gene1
             [2.  2.],    # bin1
             [6.  2.]],   # bin2
            [[6.  0.],    # bin0 - gene2
             [6.  0.],    # bin1
             [6.  0.]]]   # bin2

        With segment=False:

              Y   Z
        R:  [[3.  2.]
             [6.  0.]]

    Note: the whole segmented features track will be loaded in memory.

    :param trackScores: (FeatureStream, or list of FeatureStream objects) score track(s).
    :param trackFeatures: (FeatureStream) feature track.
    :param segment: (bool) segment each feature into bins.[False]
    :param fn: (str) Operation applied to the list of scores for one feature.
        It is the `fn` argument to `stream.score_by_feature` - one of 'sum','mean','median','min','max'.
    :param **kw: arguments to pass to segment_features (`nbins`,`upstream`,`downstream`).
    :rtype: tuple (numpy.ndarray of strings, numpy.ndarray of floats)
    """
    nbins = 1
    nscores = 1
    if segment:
        trackFeatures = sorted_stream(segment_features(trackFeatures,**kw))
        nbins = kw.get('nbins',segment_features.func_defaults[0]) \
                + kw.get('upstream',(0,0))[1] \
                + kw.get('downstream',(0,0))[1]
    all_means = score_by_feature(trackScores,trackFeatures,fn=fn)
    nfields = len(trackFeatures.fields)
    if isinstance(trackScores,(list,tuple)):
        nscores = len(trackScores)
    scores_dict = {}
    if segment:
        empty_mat = numpy.zeros(shape=(nbins,nscores))
    else:
        empty_mat = numpy.zeros(nscores)
    name_idx = all_means.fields.index('name')
    for t in all_means:
        _n = t[name_idx]
        scores_dict.setdefault(_n, empty_mat.copy())
        if segment:
            scores_dict[_n][t[nfields-1]] = t[nfields:]
        else:
            scores_dict[_n] = t[nfields:]
    feat_names = numpy.array(scores_dict.keys())
    scores_mat = numpy.array(scores_dict.values())
    return (feat_names,scores_mat)

def summed_feature_matrix(trackScores,trackFeatures,**kw):
    """
    Each feature in *trackFeatures* is segmented into bins using bbcflib.bFlatMajor.stream.segment_features
    (with parameters passed from *\*\*kw*).
    This creates a matrix with a column for each track in *trackScores* and a row for each bin in the segmented features.
    The values of a matrix entry is the score from one track in *trackScores* in one bin summed over all features.


    Example::

                      gene1                 gene2
        X: -----#####|#####|#####--------###|###|###-----  (features, nbins=3)
        Y: _____________666|66666________666|666|666_____
        Z: _____22222|22222|22222________________________

             Y   Z
        R: [[3.  1.],   # bin 0
            [4.  1.],   # bin 1
            [6.  1.]]   # bin 2

    Note: the whole segmented features track will be loaded in memory.

    :param trackScores: (FeatureStream, or list of FeatureStream objects) score track(s).
    :param trackFeatures: (FeatureStream) feature track.
    :param **kw: arguments to pass to segment_features (`nbins`,`upstream`,`downstream`).
    :rtype: numpy.ndarray, int (number of features)
    """
    nfields = len(trackFeatures.fields)
    trackFeatures = sorted_stream(segment_features(trackFeatures,**kw))
    all_means = score_by_feature(trackScores,trackFeatures)
    if isinstance(trackScores,(list,tuple)):
        nscores = len(trackScores)
    else:
        nscores = 1
    nbins = kw.get('nbins',segment_features.func_defaults[0]) \
            + kw.get('upstream',(0,0))[1] \
            + kw.get('downstream',(0,0))[1]
    averages = numpy.zeros(shape=(nbins,nscores))
    ntot = -1
    for ntot,x in enumerate(all_means):
        averages[x[nfields]] += x[(nfields+1):]
    return averages, (ntot+1)/nbins



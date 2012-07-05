from bbcflib.bFlatMajor.stream import mean_score_by_feature, segment_features
from bbcflib.bFlatMajor.common import sorted_stream
import numpy

def feature_matrix(trackScores,trackFeatures,segment=False,**kw):
    """
    ???

    :param trackScores: (FeatureStream) score track.
    :param trackFeatures: (FeatureStream) feature track.
    :param segment: (bool) ??? .[False]
    :rtype: FeatureStream
    """
    nbins = 1
    nscores = 1
    if segment:
        trackFeatures = sorted_stream(segment_features(trackFeatures,**kw))
        nbins = kw.get('nbins',segment_features.__defaults__[0])\
                +kw.get('upstream',(0,0))[1]\
                +kw.get('downstream',(0,0))[1]
    all_means = mean_score_by_feature(trackScores,trackFeatures)
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
        if not(_n in scores_dict): scores_dict[_n] = empty_mat.copy()
        if segment:
            scores_dict[_n][t[nfields-1]] = t[nfields:]
        else:
            scores_dict[_n] = t[nfields:]
    feat_names = numpy.array(scores_dict.keys())
    scores_mat = numpy.array(scores_dict.values())
    return (feat_names,scores_mat)

def average_feature_matrix(trackScores,trackFeatures,**kw):
    """
    ???

    :param trackScores: (FeatureStream) score track.
    :param trackFeatures: (FeatureStream) feature track.
    :rtype: FeatureStream
    """
    trackFeatures = sorted_stream(segment_features(trackFeatures,**kw))
    all_means = mean_score_by_feature(trackScores,trackFeatures)
    if isinstance(trackScores,(list,tuple)):
        nscores = len(trackScores)
    else:
        nscores = 1
    nfields = len(trackFeatures.fields)
    nbins = kw.get('nbins',segment_features.__defaults__[0])\
            +kw.get('upstream',(0,0))[1]\
            +kw.get('downstream',(0,0))[1]
    averages = numpy.zeros(shape=(nbins,nscores))
    for ntot,x in enumerate(all_means):
        averages[x[nfields]] += x[(nfields+1):]
    averages *= 1.0*(nbins/(ntot+1))
    return averages



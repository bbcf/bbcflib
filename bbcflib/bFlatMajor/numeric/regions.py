from bbcflib.bFlatMajor.stream import mean_score_by_feature, segment_features
import numpy

def feature_matrix(trackScores,trackFeatures):
    all_means = mean_score_by_feature(trackScores,trackFeatures)
    nfields = len(trackFeatures.fields)
    if 'name' in all_means.fields:
        name_field = all_means.fields.index('name')
    else:
        name_field = all_means.fields.index('start')
    scores_dict = dict((str(x[name_field]),x[nfields:]) for x in all_means)
    feat_names = numpy.array(scores_dict.keys())
    scores_mat = numpy.array(scores_dict.values())
    return (feat_names,scores_mat)

def average_feature_matrix(trackScores,trackFeatures,nbins=20,):
    nfields = len(trackFeatures.fields)
    trackFeatures = segment_features(trackFeatures,nbins=nbins)
    all_means = mean_score_by_feature(trackScores,trackFeatures)
    nscores = len(trackScores)
    averages = numpy.zeros(shape=(nbins,nscores))
    for ntot,x in enumerate(all_means):
        averages[x[-1]] += x[nfields:-1]
    averages *= 1.0/(ntot+1)
    return averages
    
    

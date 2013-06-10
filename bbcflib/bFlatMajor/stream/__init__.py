"""
The *stream* module contains algorithms which produce :func:`FeatureStream <bbcflib.track.FeatureStream>` out of one or more :func:`FeatureStream <bbcflib.track.FeatureStream>` objects.
"""
from bbcflib.bFlatMajor import *
###############################################################################
_members = {'merge_scores': ['trackList'],
            'filter_scores': ['trackFeatures','trackScores'],
            'score_by_feature': ['trackFeatures', 'trackScores'],
            'window_smoothing': ['trackList'],
            'normalize': ['trackList'],
            'concatenate': ['trackList'],
            'selection': ['trackList'],
            'neighborhood': ['trackList'],
            'combine': ['trackList'],
            'segment_features': ['trackList'],
            'getNearestFeature': ['features','annotations'],
            }

class stream(bFlatMajorGroup):
    def __init__(self):
        bFlatMajorGroup.__init__(self,_members)
###############################################################################

from .intervals import *
from .scores import *
from .annotate import *



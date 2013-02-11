"""
The *stream* module contains algorithms which produce :func:`FeatureStream <bbcflib.btrack.FeatureStream>` out of one or more :func:`FeatureStream <bbcflib.btrack.FeatureStream>` objects.
"""
from bbcflib.bFlatMajor import *
###############################################################################
_members = {'filter_scores': ['trackFeatures','trackScores'],
            'score_by_feature': ['trackFeatures', 'trackScores'],
            'getNearestFeature': ['features','annotations'],
            }

class stream(bFlatMajorGroup):
    def __init__(self):
        bFlatMajorGroup.__init__(self,_members)
###############################################################################

from .intervals import *
from .scores import *
from .annotate import *



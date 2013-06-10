"""
The *numeric* module contains algorithms which return numeric objects (typically `numpy.arrays`) from one or more :func:`FeatureStream <bbcflib.track.FeatureStream>` objects.
"""
from bbcflib.gfminer import *

###############################################################################
_members = {'score_array': ['trackList'],
            'correlation': ['trackList'],
            'feature_matrix': ['trackScores','trackFeatures'],
            'summed_feature_matrix': ['trackScores','trackFeatures']
            }
class numeric(gfminerGroup):
    def __init__(self):
        gfminerGroup.__init__(self,_members)
###############################################################################

from .signal import *
from .regions import *

"""
The *numeric* module contains algorithms which return numeric objects (typically `numpy.arrays`) from one or more :func:`FeatureStream <bbcflib.btrack.FeatureStream>` objects.
"""
from bbcflib.bFlatMajor import *

###############################################################################
_members = {'score_array': ['trackList'],
            'correlation': ['trackList'],
            'feature_matrix': ['trackScores','trackFeatures'],
            'summed_feature_matrix': ['trackScores','trackFeatures']
            }
class numeric(bFlatMajorGroup):
    def __init__(self):
        bFlatMajorGroup.__init__(self,_members)
###############################################################################

from .signal import *
from .regions import *

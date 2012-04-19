from bbcflib.bFlatMajor import *

###############################################################################
_members = {'correlation': ['trackList'],
            'feature_matrix': ['trackScores','trackFeatures'],
            'scaled_feature_matrix': ['trackScores','trackFeatures']
            }
class numeric(bFlatMajorGroup):
    def __init__(self):
        bFlatMajorGroup.__init__(self,_members)
###############################################################################

from .signal import *
from .regions import *

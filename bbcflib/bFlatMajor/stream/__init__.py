from bbcflib.bFlatMajor import *
###############################################################################
_members = {'merge_scores': ['trackList'],
            'mean_score_by_feature': ['trackScores', 'trackFeatures'],
            'window_smoothing': ['trackList'],
            'neighborhood': ['trackList'],
            'combine': ['trackList'],
            'getNearestFeature': ['features','annotations']
            }
class stream(bFlatMajorGroup):
    def __init__(self):
        bFlatMajorGroup.__init__(self,_members)
###############################################################################

from .intervals import *
from .scores import *
from .annotate import *

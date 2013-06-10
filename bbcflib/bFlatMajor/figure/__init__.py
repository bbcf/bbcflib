"""
The *figure* module contains algorithms which produce figure files (such as PDF or PNG) from one or more :func:`FeatureStream <bbcflib.track.FeatureStream>` objects.
"""
from bbcflib.bFlatMajor import *
###############################################################################
_members = {}
class figure(bFlatMajorGroup):
    def __init__(self):
        bFlatMajorGroup.__init__(self,_members)
###############################################################################

from .rplots import *

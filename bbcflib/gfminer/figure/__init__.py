"""
The *figure* module contains algorithms which produce figure files (such as PDF or PNG) from one or more :func:`FeatureStream <bbcflib.track.FeatureStream>` objects.
"""
from bbcflib.gfminer import *
###############################################################################
_members = {}
class figure(gfminerGroup):
    def __init__(self):
        gfminerGroup.__init__(self,_members)
###############################################################################

from .rplots import *

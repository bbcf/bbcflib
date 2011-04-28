"""
============================
Submodule: bbcflib.track.bed
============================

Implementation of the BED format.
"""

from ..track import *

###########################################################################
class GenomicFormat(Track):
    all_fields_possible = ['start', 'end', 'name', 'score', 'strand', 'thick_start', 'thick_end',
                           'item_rgb', 'block_count', 'block_sizes', 'block_starts']

    def load(self):
        self.file = open(self.path)

    def unload(self, type, value, trackback):
        self.file.close()

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

"""
=======================================
Submodule: bbcflib.track.formats.bigWig
=======================================

Implementation of the bigWig format.
"""

# Internal modules #
from ..track_binary import TrackBinary
from .bedGraph import TrackFormat as TrackBedgraph

###########################################################################
class TrackFormat(TrackBinary, TrackBedgraph):
    type_identifier = 'bigWig'

    def binary_to_text(self, source, dest):
        self.run_tool('bigWigToBedGraph', [source, dest])

    def text_to_binary(self, source, dest):
        self.run_tool('bedGraphToBigWig', [source, self.generate_chr_file(), dest])

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

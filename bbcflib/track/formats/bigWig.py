"""
=======================================
Submodule: bbcflib.track.formats.bigWig
=======================================

Implementation of the bigWig format.
"""

# Internal modules #
from bbcflib.track.track_binary import TrackBinary
from bbcflib.track.formats.wig  import TrackFormat as TrackWIG

###########################################################################
class TrackFormat(TrackBinary, TrackWIG):
    type_identifier = 'bigWig'

    def binary_to_text(self, source, dest):
        self.run_tool('bigWigToWig', [source, dest])

    def text_to_binary(self, source, dest):
        self.run_tool('wigToBigWig', [source, self.generate_chr_file(), dest])

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

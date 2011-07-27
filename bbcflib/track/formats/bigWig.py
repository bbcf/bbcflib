"""
=======================================
Submodule: bbcflib.track.formats.bigWig
=======================================

Implementation of the bigWig format. Requires the command line utilities "bigWigToBedGraph" and "bedGraphToBigWig" to be installed and present in the $PATH environment variable.
"""

# Built-in modules #
import os

# Internal modules #
from ..common import check_executable
from ..track_binary import TrackBinary
from .bedGraph import TrackFormat as TrackBedgraph

###########################################################################
class TrackFormat(TrackBinary, TrackBedgraph):
    backend_format = 'bedGraph'

    def __init__(self, *args, **kwargs):
        check_executable('bigWigToBedGraph')
        check_executable('bedGraphToBigWig')
        super(TrackFormat, self).__init__(*args, **kwargs)

    def binary_to_text(self, source, dest):
        self.run_tool('bigWigToBedGraph', [source, dest])

    def text_to_binary(self, source, dest):
        chrfile = self.chrmeta.write_file()
        self.run_tool('bedGraphToBigWig', [source, chrfile, dest])
        os.remove(chrfile)

    #-----------------------------------------------------------------------------#
    @property
    def _datatype(self):
        return 'quantitative'

    @_datatype.setter
    def _datatype(self, datatype):
        if datatype and datatype != 'quantitative':
            raise Exception("The track '" + self._path + "' cannot be loaded as a '" + datatype + "' datatype.")

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

"""
=======================================
Submodule: bbcflib.btrack.formats.bigWig
=======================================

Implementation of the bigWig format. Requires the command line utilities "bigWigToBedGraph" and "bedGraphToBigWig" to be installed and present in the $PATH environment variable.
"""

###########################################################################
###########################################################################
## WARNING: The bbcflib.track package is deprecated.                     ##
##          A new project simply called 'track' replaces it.             ##
###########################################################################
###########################################################################

# Built-in modules #
import os

# Internal modules #
from bbcflib.btrack.common import check_executable
from bbcflib.btrack.track_binary import TrackBinary
from bbcflib.btrack.formats.bedGraph import TrackFormat as TrackBedgraph

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
        chrfile = self.chrmeta.write_file(sep=' ')
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

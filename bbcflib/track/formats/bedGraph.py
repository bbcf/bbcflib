"""
=========================================
Submodule: bbcflib.track.formats.bedgraph
=========================================

Implementation of the bedGraph format.
"""

# Internal modules #
from ..track_proxy import TrackProxy
from ..track_text import TrackText






###########################################################################
class TrackFormat(TrackText, TrackProxy):
    type_identifier   = 'bedGraph'

    def _read_entries(self):
        self._file.seek(0)
        seen_track = False
        for line in self._file:
            line = line.strip("\n").lstrip()
            if len(line) == 0:              continue
            if line.startswith("#"):        continue
            if line.endswith(" \\"):
                raise Exception("The file '" + self._path + "' includes linebreaks ('\\') which is not supported.")
            if line.startswith("track "):
                if not seen_track:
                    seen_track = True
                    continue
                raise Exception("The file '" + self._path + "' contains a second 'track' directive. This is not supported.")
            if line.startswith("browser "): continue
            line = line.split()
            if len(line) != 4:
                raise Exception("The file '" + self._path + "' must have 4 columns on every line.")
            try:
                line[1] = int(line[1])
                line[2] = int(line[2])
            except ValueError:
                raise Exception("The file '" + self._path + "' has non integers as interval bounds and is hence not valid.")
            if line[2] <= line[1]:
                raise Exception("The file '" + self._path + "' has negative or null intervals and is hence not valid.")
            try:
                line[3] = float(line[3])
            except ValueError:
                raise Exception("The file '" + self._path + "' has non floats as score values and is hence not valid.")
            yield line

    def _write_entries(self):
        for f in self.read(): yield ' '.join([str(x) for x in f]) + '\n'

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

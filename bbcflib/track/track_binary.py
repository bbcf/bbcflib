"""
=====================================
Submodule: bbcflib.track.track_binary
=====================================

Methods that create a text file upon opening a binary file in the temporary directory and reconverts everything to a binary file file upon closing, all transparently.
"""

# Built-in modules #
import os, subprocess

# Internal modules #
from .common import named_temporary_path

###########################################################################
class TrackBinary(object):
    def __init__(self, path, format=None, name=None, chrmeta=None, datatype=None, readonly=False):
        # Parameters with a double underscore refer to the underlying binary track #
        self.__path    = path
        self.__format  = self.type_identifier
        # Create the text track #
        tmp_path = named_temporary_path()
        self.binary_to_text(self.__path, tmp_path)
        # Load the new text track as self #
        super(TrackBinary, self).__init__(tmp_path, format, name, chrmeta, datatype, readonly)

    def dump(self, path=None):
        super(TrackBinary, self).dump()
        if not path: path = self.__path
        elif os.path.exists(path): raise Exception("The location '" + path + "' is already taken.")
        self.text_to_binary(self._path, self.__path)

    #-----------------------------------------------------------------------------#
    def run_tool(self, tool_name, args):
        return_code = subprocess.call([tool_name] + args)
        if return_code != 1:
            raise Exception("The tool '" + tool_name + "' did not exit with success code 1.")

    def generate_chr_file(self):
        tmp_path = named_temporary_path()
        with open(tmp_path, 'w') as f: f.writelines(self.chromosome_file)
        return tmp_path

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#


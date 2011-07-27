"""
=====================================
Submodule: bbcflib.track.track_binary
=====================================

Methods that create a text file upon opening a binary file in the temporary directory and reconverts everything to a binary file file upon closing, all transparently.
"""

# Built-in modules #
import os, shutil, subprocess

# Internal modules #
from . import Track
from .common import named_temporary_path
from .track_proxy import TrackProxy, TrackBackend, backend_format

###########################################################################
class TrackBinary(object):
    def __init__(self, path, format=None, name=None, chrmeta=None, datatype=None, readonly=False, empty=False):
        # Parameters with a double underscore refer to the underlying binary track #
        self.__path  = path
        self._format = self.backend_format
        # Create the text track #
        tmp_path = named_temporary_path('.' + self.backend_format)
        if empty: open(tmp_path, 'w').close()
        else: self.binary_to_text(self.__path, tmp_path)
        # Load the new text track as self #
        super(TrackBinary, self).__init__(tmp_path, format, name, chrmeta, datatype, readonly, empty)

    def dump(self, path=None):
        super(TrackBinary, self).dump()
        if not path: path = self.__path
        elif os.path.exists(path): raise Exception("The location '" + path + "' is already taken.")
        self.text_to_binary(self._path, path)

    def unload(self, datatype=None, value=None, traceback=None):
        if self.modified and not self.readonly: self.dump()
        if os.path.exists(self.path): os.remove(self.path)
        if os.path.exists(self._path): os.remove(self._path)

    #-----------------------------------------------------------------------------#
    def run_tool(self, tool_name, args):
        proc = subprocess.Popen([tool_name] + args, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if stderr: raise Exception("The tool '" + tool_name + "' exited with message: " + '"' + stderr.strip('\n') + '"')

    #--------------------------------------------------------------------------#
    @classmethod
    def mutate_destination(cls, self, path, format):
        # Maybe use the same temporary SQL (Case 1)
        if issubclass(self.__class__, TrackProxy):
            self._path = named_temporary_path('.' + cls.backend_format)
            self.__path = path
            self.__class__ = cls
            self.modified = True
        # Or create a new one by copying it (Case 3)
        elif issubclass(self.__class__, TrackBackend):
            tmp_path = named_temporary_path('.' + backend_format)
            shutil.copy(self.path, tmp_path)
            self.__init__(tmp_path)
            self._path = named_temporary_path('.' + cls.backend_format)
            self.__path = path
            self.modified = True
            self.__class__ = cls
        else: Track.mutate_destination(path, format)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

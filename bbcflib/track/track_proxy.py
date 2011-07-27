"""
====================================
Submodule: bbcflib.track.track_proxy
====================================

Methods that create an SQL file upon opening in the temporary directory and reconverts everything to a text file file upon closing, all transparently.
"""

# Built-in modules #
import os, shutil

# Internal modules #
from . import Track, new
from .common import named_temporary_path, check_path
from .formats.sql import TrackFormat as TrackBackend

# Globals #
backend_format = 'sql'

###########################################################################
class TrackProxy(TrackBackend):
    def __init__(self, path, format=None, name=None, chrmeta=None, datatype=None, readonly=False, empty=False):
        # Parameters with underscore refer to the overlying track #
        # Parameters without the underscore refer to the underlying track #
        self._path     = path
        self._datatype = datatype
        # Create the SQL track #
        tmp_path = named_temporary_path('.' + backend_format)
        with new(tmp_path, backend_format, name=name, datatype=self._datatype) as t:
            if not empty:
                with open(self._path, 'r') as self._file:
                    t.attributes.update(self._read_header())
                    fields = self._fields
                    for chrom, data in self._read(): t.write(chrom, data, fields)
        # Load the new SQL track as self #
        super(TrackProxy, self).__init__(tmp_path, backend_format, name=None, chrmeta=chrmeta, datatype=None, readonly=readonly)

    def unload(self, datatype=None, value=None, traceback=None):
        if self.modified and not self.readonly: self.dump()
        if os.path.exists(self.path): os.remove(self.path)

    def commit(self):
        super(TrackProxy, self).commit()
        self.dump()

    def dump(self, path=None):
        if not path: path = self._path
        else: check_path(path)
        with open(path, 'w') as file: file.writelines(self._write())

    #--------------------------------------------------------------------------#
    def mutate_source(self, path, format):
        # Move the SQL temporary file (Case 2)
        if format == backend_format:
            self.format = format
            check_path(path)
            super(TrackProxy, self).unload()
            shutil.move(self.path, path)
            self.format = backend_format
            self.__class__ = TrackBackend
            self.__init__(path, format)
        else: super(TrackProxy, self).mutate_source(path, format)

    @classmethod
    def mutate_destination(cls, self, path, format):
        # Maybe use the same temporary SQL (Case 1)
        if issubclass(self.__class__, TrackProxy):
            self._path = path
            self.__class__ = cls
            self.modified = True
        # Or create a new one by copying it (Case 3)
        elif issubclass(self.__class__, TrackBackend):
            self.unload()
            tmp_path = named_temporary_path('.' + backend_format)
            shutil.copy(self.path, tmp_path)
            self.__init__(tmp_path)
            self._path = path
            self.modified = True
            self.__class__ = cls
        else: Track.mutate_destination(path, format)

    #--------------------------------------------------------------------------#
    @staticmethod
    def create(path):
        open(path, 'w').close()

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

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
from .common import named_temporary_path
from .formats.sql import TrackFormat as TrackBackend

# Globals #
backend_format = 'sql'

###########################################################################
class TrackProxy(TrackBackend):
    def __init__(self, path, format=None, name=None, chrmeta=None, datatype=None, readonly=False, empty=False):
        # Parameters with underscore refer to the overlying track #
        self._path     = path
        self._datatype = datatype
        # Parameters without the underscore refer to the underlying track #
        self.modified = False
        # Create the SQL track #
        tmp_path = named_temporary_path('.' + backend_format)
        with new(tmp_path, backend_format, name) as t:
            if not empty:
                with open(self._path, 'r') as self._file:
                    t.attributes.update(self._read_header())
                    fields = self._fields
                    for chrom, data in self._read(): t.write(chrom, data, fields)
        # Load the new SQL track as self #
        super(TrackProxy, self).__init__(tmp_path, backend_format, name, chrmeta, self._datatype, readonly)

    def unload(self, datatype=None, value=None, traceback=None):
        if self.modified and not self.readonly: self.dump()
        if os.path.exists(self.path): os.remove(self.path)

    def commit(self):
        super(TrackProxy, self).commit()
        self.dump()

    def dump(self, path=None):
        if not path: path = self._path
        elif os.path.exists(path): raise Exception("The location '" + path + "' is already taken")
        with open(path, 'w') as file: file.writelines(self._write())

    def convert(self, path, format=backend_format):
        if os.path.exists(path): raise Exception("The location '" + path + "' is already taken")
        if format == backend_format: shutil.move(self.path, path)
        else: super(TrackProxy, self).convert(path, format)

    @classmethod
    def create(cls, path, format, name, chrmeta, datatype):
        open(path, 'w').close()
        return Track(path, format, name, chrmeta, datatype, empty=True)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

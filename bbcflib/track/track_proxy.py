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
    def __init__(self, path, format=None, name=None, chrmeta=None, datatype=None, readonly=False):
        # Parameters with underscore refer to the underlying track #
        self._path     = path
        self._format   = self.type_identifier
        self._datatype = datatype
        # Parameters without the underscore refer to the exposed track #
        self.chrmeta  = chrmeta
        self.modified = False
        # Create the SQL track #
        tmp_path = named_temporary_path()
        with open(self._path, 'r') as self._file:
            with new(tmp_path, backend_format, self._datatype, name) as t:
                # Prepare meta data #
                self._fields
                self._meta_chr
                self._meta_track
                # Copy data #
                for chrom, data in self._read():
                    t.write(chrom, data, self._fields)
                # Copy meta track #
                self._meta_track_dict = self._meta_track
                self._meta_track_dict['converted_by']   = __package__
                self._meta_track_dict['converted_from'] = self._path
                t.meta_track = self._meta_track_dict
                # Copy meta chr #
                t.meta_chr   = [chr for chr in self._meta_chr if chr['name'] in self._seen_chr]
        # Load the new SQL track as self #
        super(TrackProxy, self).__init__(tmp_path, backend_format, name, chrmeta=None, datatype=None, readonly=readonly)

    #-----------------------------------------------------------------------------#
    def unload(self, datatype=None, value=None, trackback=None):
        if self.modified and not self.readonly: self.dump()
        super(TrackProxy, self).unload(datatype, value, trackback)
        if os.path.exists(self.path): os.remove(self.path)

    def commit(self):
        super(TrackProxy, self).commit()
        self.dump()

    def dump(self, path=None):
        if not path: path = self._path
        elif os.path.exists(path): raise Exception("The location '" + path + "' is already taken")
        with open(path, 'w') as file: file.writelines(self._write())

    def convert(self, path, format=backend_format):
        if os.path.exists(path):
            raise Exception("The location '" + path + "' is already taken")
        if format == backend_format:
            shutil.move(self.path, path)
        else:
            super(TrackProxy, self).convert(path, format)

    #-----------------------------------------------------------------------------#
    def _read(self):
        raise NotImplementedError

    def _write(self):
        raise NotImplementedError

    @property
    def _fields(self):
        return getattr(Track, self._datatype + '_fields')

    @property
    def _datatype(self):
        raise NotImplementedError

    @_datatype.setter
    def _datatype(self, datatype):
        raise NotImplementedError

    #-----------------------------------------------------------------------------#
    @property
    def meta_chr(self):
        return super(TrackProxy, self).meta_chr

    @meta_chr.setter
    def meta_chr(self, data):
        self.modified = True
        super(TrackProxy, self).set_meta_chr(data)

    @property
    def meta_track(self):
        return super(TrackProxy, self).meta_track

    @meta_track.setter
    def meta_track(self, data):
        self.modified = True
        super(TrackProxy, self).set_meta_track(data)

    def write(self, chrom, data, fields=None):
        self.modified = True
        super(TrackProxy, self).write(chrom, data, fields)

    def remove(self, chrom=None):
        self.modified = True
        super(TrackProxy, self).remove(chrom)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

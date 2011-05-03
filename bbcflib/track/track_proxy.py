"""
====================================
Submodule: bbcflib.track.track_proxy
====================================

Methods that create an SQL file upon opening in the temporary directory and reconverts everything to a text file file upon closing, all transparently.
"""

# General modules #
import os, shutil

# Specific module #
from ..common import named_temporary_path
from ..track import new
from .format_sql import GenomicFormat as SQLTrack

###########################################################################
class ProxyTrack(SQLTrack):
    def __init__(self, path, format=None, name=None, chrfile=None, type=None):
        # Parameters with underscore refer to the underlying track #
        self._path    = path
        self._format  = format
        self._type    = type 
        # Parameters without the underscore refer to the exposed track #
        self.chrfile  = chrfile
        self.modified = False
        # Create the SQL track #
        tmp_path = named_temporary_path()
        with open(self._path, 'r') as self._file:
            with new(tmp_path, 'sql', self._type, name) as t:
                # Prepare meta data #
                self._meta_chr
                self._meta_track
                self._fields
                # Copy data #
                for chrom, data in self._read():
                    t.write(chrom, data, self._fields)
                # Copy meta track #
                self._meta_track_dict = self._meta_track
                self._meta_track_dict['datatype']       = self._type
                self._meta_track_dict['converted_by']   = __package__
                self._meta_track_dict['converted_from'] = self._path
                t.meta_track = self._meta_track_dict
                # Copy meta chr #
                t.meta_chr   = [chr for chr in self._meta_chr if chr['name'] in self._seen_chr]
        # Load the new SQL track as self #
        super(ProxyTrack, self).__init__(tmp_path, 'sql', name)

    #-----------------------------------------------------------------------------#
    def unload(self, type, value, trackback):
        if self.modified: self.dump()
        super(ProxyTrack, self).unload(type, value, trackback)
        if os.path.exists(self.path): os.remove(self.path)

    def commit(self):
        super(ProxyTrack, self).commit()
        self.dump()
    
    def dump(self, path=None):
        if not path: path = self._path
        elif os.path.exists(path): raise Exception("The location '" + path + "' is already taken")
        with open(path, 'w') as file: file.writelines(self._write())

    def convert(self, path, format='sql'):
        if os.path.exists(path):
            raise Exception("The location '" + path + "' is already taken")
        if format == 'sql':
            shutil.move(self.path, path)
        else:
            super(ProxyTrack, self).convert(path, format, type, mean_scores)

    #-----------------------------------------------------------------------------#
    def _read(self):
        raise NotImplementedError

    def _write(self):
        raise NotImplementedError

    @property
    def _fields(self):
        return self.default_fields
    
    @property
    def _type(self): 
        raise NotImplementedError

    @_type.setter
    def _type(self, datatype):
        raise NotImplementedError

    #-----------------------------------------------------------------------------#
    @property
    def meta_chr(self): 
        return super(ProxyTrack, self).meta_chr

    @meta_chr.setter
    def meta_chr(self, data):
        self.modified = True
        super(ProxyTrack, self).set_meta_chr(data)

    @property
    def meta_track(self): 
        return super(ProxyTrack, self).meta_track
    
    @meta_track.setter
    def meta_track(self, data): 
        self.modified = True
        super(ProxyTrack, self).set_meta_track(data)

    def write(self, chrom, data, fields=None):
        self.modified = True
        super(ProxyTrack, self).write(chrom, data, fields)

    def remove(self, chrom=None):
        self.modified = True
        super(ProxyTrack, self).remove(chrom)
   
#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

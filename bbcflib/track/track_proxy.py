"""
====================================
Submodule: bbcflib.track.track_proxy
====================================

Methods that create an SQL file upon opening in the temporary directory and reconverts everything to a text file file upon closing, all transparently.
"""

# General modules #
import os

# Specific module #
from ..common import named_temporary_path
from ..track import new
from .format_sql import GenomicFormat as SQLTrack

###########################################################################
class ProxyTrack(SQLTrack):
    def __init__(self, path, format=None, name=None, chrfile=None):
        # Parameters with underscore refer to the underlying track #
        self._path    = path
        self._format  = format
        # Parameters without the underscore refer to the exposed track #
        self.chrfile = chrfile
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
                self._file.seek(0)
                for chrom, data in self._read():
                    t.write(chrom, data, self._fields)
                # Copy meta track #
                self._meta_track['datatype']       = self._type
                self._meta_track['converted_by']   = __package__
                self._meta_track['converted_from'] = self._path
                t.meta_track = self._meta_track
                # Copy meta chr #
                t.meta_chr   = [chr for chr in self._meta_chr if chr['name'] in self._seen_chr]
                t.meta_chr   = [chr for chr in self._meta_chr if chr['name'] in self._seen_chr]
                t.meta_chr   = [chr for chr in self._meta_chr if chr['name'] in self._seen_chr]
        # Load the new SQL track as self #
        super(ProxyTrack, self).__init__(tmp_path, 'sql', name)

    def unload(self, type, value, trackback):
        super(ProxyTrack, self).unload(type, value, trackback)
        if self.modified: self.commit()

    def commit(self, path=None):
        self.dump()
    
    def dump(self, path=None):
        if not path:
            path = self._path
        elif os.path.exists(path):
            raise Exception("The location '" + path + "' is already taken")
        with open(path, 'w') as file:
            file.writelines(self._ouput())

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
    
    #-----------------------------------------------------------------------------#
    def convert(self, path, format='sql'):
        if os.path.exists(path):
            raise Exception("The location '" + path + "' is already taken")
        if format == 'sql':
            os.rename(self.path, path)
        else:
            super(ProxyTrack, self).convert(path, format, type, mean_scores)

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

"""
====================================
Submodule: bbcflib.track.proxy_track
====================================

Methods that create an SQL file upon opening in the temporary directory and reconverts everything to a text file file upon closing, all transparently.
"""

from ..common import named_temporary_path
from ..track import new

###########################################################################
class ProxyTrack(object):
    def __init__(self, path, format=None, name=None, chrfile=None):
        # Parameters with underscore refer to the underlying track #
        self._path    = path
        self._format  = format
        self._chrfile = chrfile
        # Parameters without the underscore refer to the exposed track #
        self.modified = False
        # Create the SQL track #
        tmp_path = named_temporary_path()
        with open(self._path, 'r') as self._file:
            with new(tmp_path, 'sql', self._type, name) as t:
                # Prepare meta data #
                self._meta_chr   = self._get_meta_chr()
                self._meta_track = self._get_meta_track()
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
                t.meta_chr   = [chr for chr in self._chr_meta if chr['name'] in self._seen_chr]
        # Load the new SQL track as self #
        super(ProxyTrack, self).__init__(tmp_path, 'sql', name)

    def unload(self, type, value, trackback):
        super(ProxyTrack, self).unload(type, value, trackback)
        if self.modified: self.commit()

    def commit(self):
        self.dump()

    #-----------------------------------------------------------------------------#
    @meta_chr.setter
    def set_meta_chr(self, data):
        self.modified = True
        super(ProxyTrack, self).set_meta_chr(data)

    @meta_track.setter
    def set_meta_track(self, data): 
        self.modified = True
        super(ProxyTrack, self).set_meta_track(data)

    def write(self, chrom, data, fields=None):
        self.modified = True
        super(ProxyTrack, self).write(chrom, data, fields)

    def remove(self, chrom=None):
        self.modified = True
        super(ProxyTrack, self).remove(chrom)
    
    #-----------------------------------------------------------------------------#
    def convert(self, path, format='sql', type=None, mean_scores=False):
        if format == 'sql' and type == self.type:
            pass
        else:
            super(ProxyTrack, self).convert(path, format, type, mean_scores)

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

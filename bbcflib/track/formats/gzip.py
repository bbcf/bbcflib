"""
====================================
Submodule: bbcflib.track.formats.sql
====================================

Implementation of a compression layer GNU zip format.
"""

# Built-in modules #
import gzip

# Internal modules #
from bbcflib.track.track_util import determine_format, import_implementation

################################################################################
def TrackFormat(path, format=None, name=None, chrmeta=None, datatype=None, readonly=False, empty=False):
    if not format: format = determine_format(path)
    implementation = import_implementation(format)
    base = implementation.TrackFormat
    #----------------------------------#
    class TrackCompressed(base):
        @property
        def file_obj(self): return gzip.open(self.path, 'r')
    #----------------------------------#
    return TrackCompressed(path, format, name, chrmeta, datatype, readonly, empty)

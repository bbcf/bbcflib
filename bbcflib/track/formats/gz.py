"""
====================================
Submodule: bbcflib.track.formats.sql
====================================

Implementation of a compression layer GNU zip format.
"""

# Built-in modules #
import os, gzip

# Internal modules #
from bbcflib.track.track_util import import_implementation

################################################################################
def TrackFormat(path, format=None, name=None, chrmeta=None, datatype=None, readonly=False, empty=False):
    # TODO: recognize the format of a compressed file if no extension in the path
    format = os.path.splitext(os.path.splitext(path)[0])[1][1:]
    implementation = import_implementation(format)
    base = implementation.TrackFormat
    #----------------------------------#
    class TrackCompressed(base):
        def open(self, path, mode):
            return gzip.open(path, mode)
    #----------------------------------#
    return TrackCompressed(path, format, name, chrmeta, datatype, readonly, empty)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

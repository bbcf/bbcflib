"""
===================================
Submodule: bbcflib.track.track_text
===================================

Methods common to the text formats.
"""

# Genreral Modules #
import os, shlex

# Specific Modules #
from ..common import memoize_once

#-----------------------------------------------------------------------------#
def strand_to_int(strand):
    if strand == '+': return 1
    if strand == '-': return -1
    return 0

def int_to_strand(num):
    if num == 1: return  '+'
    if num == -1: return '-'
    return '.'

###########################################################################
class TextTrack(object):
    @property
    @memoize_once
    def _meta_chr(self):
        if not self.chrfile:
            raise Exception("The track '" + self._path + "' does not have a chromosome file associated.")
        if not os.path.exists(self.chrfile):
            raise Exception("The file '" + self.chrfile + "' cannot be found")
        if os.path.isdir(self.chrfile):
            raise Exception("The location '" + self.chrfile + "' is a directory (a file was expected).")
        result = []
        with open(self.chrfile, 'r') as f:
            for line in f:
                line = line.strip('\n')            
                if len(line) == 0:       continue
                if line.startswith("#"): continue
                if line.endswith(" \\"):
                    raise Exception("The file '" + self.chrfile + "' includes linebreaks ('\\') which are not supported.")
                if '\t' in line: seperator = '\t'
                else:            seperator = ' '
                line = line.split(seperator)
                if len(line) != 2:
                    raise Exception("The file " + self.chrfile + " does not seam to be a valid chromosome file.")
                name = line[0]
                try:
                    length = int(line[1])
                except ValueError as err:
                    raise Exception("The file '" + self.chrfile + "' has invalid values.", err)
                result.append(dict([('name', name),('length', length)]))
        if not result:
            raise Exception("The file '" + self.chrfile + "' does not seam to contain any information.", err)
        return result

    @property
    @memoize_once
    def _meta_track(self):
        self._file.seek(0)
        result = {}    
        for line in self._file:
            line = line.strip("\n").lstrip()
            if len(line) == 0:              continue
            if line.startswith("#"):        continue
            if line.endswith(" \\"):
                raise Exception("The track '" + self._path + "' includes linebreaks ('\\') which are not supported.")
            if line.startswith("browser "): continue
            if line.startswith("track "):
                try:
                    result = dict([p.split('=',1) for p in shlex.split(line[6:])])
                except ValueError as err:
                    raise Exception("The <track> header line for the file '" + self._path + "' seams to be invalid", err)
            return result

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

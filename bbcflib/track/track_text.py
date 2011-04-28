"""
====================================
Submodule: bbcflib.track.text_format
====================================

Methods common to the text formats.
"""

import os, shlex

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
    def _get_meta_chr(self):
        if not self._chr_file:
            raise Exception("The track '" + self._path + "' does not have a chromosome file associated.")
        if not os.path.exists(self._chr_file):
            raise Exception("The file '" + self._chr_file + "' cannot be found")
        if os.path.isdir(self._chr_file):
            raise Exception("The location '" + self._chr_file + "' is a directory (a file was expected).")
        result = []
        with open(self._chr_file, 'r') as f:
            for line in f:
                line = line.strip('\n')            
                if len(line) == 0:       continue
                if line.startswith("#"): continue
                if line.endswith(" \\"):
                    raise Exception("The file '" + self._chr_file + "' includes linebreaks ('\\') which are not supported.")
                if '\t' in line: seperator = '\t'
                else:            seperator = ' '
                line = line.split(seperator)
                if len(line) != 2:
                    raise Exception("The file " + self._chr_file + " does not seam to be a valid chromosome file.")
                name = line[0]
                try:
                    length = int(line[1])
                except ValueError as err:
                    raise Exception("The file '" + self._chr_file + "' has invalid values.", err)
                result.append(dict([('name', name),('length', length)]))
        if not result:
            raise Exception("The file '" + self._chr_file + "' does not seam to contain any information.", err)
        return result

    @property
    def _get_meta_track(self):
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

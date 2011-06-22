"""
===================================
Submodule: bbcflib.track.track_text
===================================

Methods common to the text formats.
"""

# General Modules #
import os, shlex

# Internal Modules #
from ..common import memoized_method
from ..genrep import GenRep

# Functions #
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
    def _read(self):
        global chrom, entry, generator
        chrom           = ''
        entry           = ['', '', '', '']
        generator       = self._all_entries()
        self._seen_chr  = set()
        def get_next_entry():
            global entry, generator
            entry = generator.next()
        def iter_until_different_chr():
            global chrom, entry
            while True:
                if entry[0] != chrom: break
                yield entry[1:]
                get_next_entry()
        get_next_entry()
        while True:
            if entry[0] == chrom: break
            chrom = entry[0]
            if not chrom in self._all_chrs:
                raise Exception("The file '" + self._path + "' has a value (" + chrom + ") not specified in the chromosome file.")
            self._seen_chr.add(chrom)
            yield chrom, iter_until_different_chr()

    #-----------------------------------------------------------------------------#
    @property
    @memoized_method
    def _meta_chr(self):
        def parse_dict(info):
            return [dict([("name", self.chrmeta[chr]["name"]),("length", self.chrmeta[chr]["length"])]) for chr in self.chrmeta]
        # Is a dictionary #
        if isinstance(self.chrmeta, dict):
            return parse_dict(self.chrmeta)
        # Is a path #
        elif os.path.exists(self.chrmeta):
            if not self.chrmeta:
                raise Exception("The file '" + self._path + "' does not have a chromosome file associated.")
            if not os.path.exists(self.chrmeta):
                raise Exception("The file '" + self.chrmeta + "' cannot be found")
            if os.path.isdir(self.chrmeta):
                raise Exception("The location '" + self.chrmeta + "' is a directory (a file was expected).")
            result = []
            with open(self.chrmeta, 'r') as f:
                for line in f:
                    line = line.strip('\n')
                    if len(line) == 0:       continue
                    if line.startswith("#"): continue
                    if line.endswith(" \\"):
                        raise Exception("The file '" + self.chrmeta + "' includes linebreaks ('\\') which are not supported.")
                    if '\t' in line: seperator = '\t'
                    else:            seperator = ' '
                    line = line.split(seperator)
                    if len(line) != 2:
                        raise Exception("The file " + self.chrmeta + " does not seam to be a valid chromosome file.")
                    name = line[0]
                    try:
                        length = int(line[1])
                    except ValueError:
                        raise Exception("The file '" + self.chrmeta + "' has invalid values.")
                    result.append(dict([('name', name),('length', length)]))
            if not result:
                raise Exception("The file '" + self.chrmeta + "' does not seam to contain any information.")
            return result
        # Is a string assembly #
        else:
            g = GenRep()
            if not g.is_available(self.chrmeta):
                raise Exception("The genrep server does not know about the assembly '" + self.chrmeta + "'.")
            return parse_dict(g.assembly(self.chrmeta).chromosomes())

    @property
    @memoized_method
    def _meta_track(self):
        self._file.seek(0)
        result = {}
        for line in self._file:
            line = line.strip("\n").lstrip()
            if len(line) == 0:              continue
            if line.startswith("#"):        continue
            if line.endswith(" \\"):
                raise Exception("The file '" + self._path + "' includes linebreaks ('\\') which are not supported.")
            if line.startswith("browser "): continue
            if line.startswith("track "):
                try:
                    result = dict([p.split('=',1) for p in shlex.split(line[6:])])
                except ValueError:
                    raise Exception("The <track> header line for the file '" + self._path + "' seams to be invalid")
        return result

    @property
    def _all_chrs(self):
       return [x['name'] for x in self._meta_chr]

    @property
    def _header_line(self):
        self.meta_track_dict = self.meta_track
        self.meta_track_dict['type']           = self.type_identifier
        self.meta_track_dict['converted_by']   = __package__
        self.meta_track_dict['converted_from'] = self.path
        return "track " + ' '.join([key + '="' + value + '"' for key, value in self.meta_track_dict.items()]) + '\n'

    #-----------------------------------------------------------------------------#
    @staticmethod
    def create(path, datatype, name):
        open(path, 'w').close()

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

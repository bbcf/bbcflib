"""
===================================
Submodule: bbcflib.track.track_text
===================================

Methods common to the text formats.
"""

# Built-in modules #
import shlex

# Internal modules #
from . import Track

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
class TrackText(object):
    def _read(self):
        global chrom, entry, generator
        chrom           = ''
        entry           = ['', '', '', '']
        generator       = self._read_entries()
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
            self._seen_chr.add(chrom)
            yield chrom, iter_until_different_chr()

    def _write(self):
        yield self._write_header()
        for l in self._write_entries(): yield l

    #--------------------------------------------------------------------------#
    def _read_header(self):
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

    def _write_header(self):
        d = self.attributes
        d['type'] = self.type_identifier
        d['converted_by'] = __package__
        return "track " + ' '.join([k + '="' + v + '"' for k, v in d.items()]) + '\n'

    #--------------------------------------------------------------------------#
    @property
    def _fields(self):
        return getattr(Track, self._datatype + '_fields')

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

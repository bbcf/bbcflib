"""
===================================
Submodule: bbcflib.track.format_bed
===================================

Implementation of the BED format.
"""

# Specific Modules #
from .track_proxy import ProxyTrack
from .track_text import TextTrack, strand_to_int, int_to_strand
from ..common import memoized_method

# Global variables #
all_fields_possible = ['start', 'end', 'name', 'score', 'strand', 'thick_start', 'thick_end',
                       'item_rgb', 'block_count', 'block_sizes', 'block_starts']

###########################################################################
class GenomicFormat(TextTrack, ProxyTrack):
    @property
    @memoized_method
    def _fields(self):
        self._file.seek(0)
        for line in self._file:
            line = line.strip("\n").lstrip()
            if len(line) == 0:              continue
            if line.startswith("#"):        continue
            if line.endswith(" \\"):
                raise Exception("The track '" + self._path + "' includes linebreaks ('\\') which are not supported.")
            if line.startswith("track "):   continue
            if line.startswith("browser "): continue
            else:
                if '\t' in line:  self._seperator = '\t'
                else:             self._seperator = ' '
                line = line.split(self._seperator)
            self.num_fields = len(line) - 1
            if self.num_fields < 2:
                raise Exception("The track '" + self._path + "' has less than three columns and is hence not a valid BED file.")
            if self.num_fields > len(all_fields_possible):
                raise Exception("The track '" + self._path + "' has too many columns and is hence not a valid BED file.")
            return all_fields_possible[0:self.num_fields]

    def _read(self):
        global line, chrom
        line = ''
        chrom  = ''
        self._seen_chr = []
        def get_next_line():
            global line, chrom
            for line in self._file:
                line = line.strip("\n").lstrip()
                if len(line) == 0:              continue
                if line.startswith("#"):        continue
                if line.endswith(" \\"):
                    raise Exception("The track '" + self._path + "' includes linebreaks ('\\') which is not supported.")
                if line.startswith("track "):
                    if not chrom: continue
                    raise Exception("The file '" + self._path + "' contains a second 'track' directive. This is not supported.")
                if line.startswith("browser "):
                    if not chrom: continue
                    raise Exception("The file '" + self._path + "' contains a 'browser' directive. This is not supported.")
                line = line.split(self._seperator)
                if len(line) != self.num_fields + 1:
                    raise Exception("The track '" + self._path + "' has a varying number of columns. This is not supported.")
                try:
                    line[1] = int(line[1])
                    line[2] = int(line[2])
                except ValueError:
                    raise Exception("The track '" + self._path + "' has non integers as interval bounds and is hence not valid.")
                if line[2] <= line[1]:
                    raise Exception("The track '" + self._path + "' has negative or null intervals and is hence not valid.")
                if len(line) > 4:
                    if line[4] == '.': line[4] = 0.0
                    try:
                        line[4] = float(line[4])
                    except ValueError:
                        raise Exception("The track '" + self._path + "' has non floats as score values and is hence not valid.")
                if len(line) > 5:
                    line[5] = strand_to_int(line[5])
                if len(line) > 6:
                    try:
                        line[6] = float(line[6])
                    except ValueError:
                        raise Exception("The track '" + self._path + "' has non integers as thick starts and is hence not valid.")
                if len(line) > 7:
                    try:
                        line[7] = float(line[7])
                    except ValueError:
                        raise Exception("The track '" + self._path + "' has non integers as thick ends and is hence not valid.")
                break
        def iter_until_different_chr():
            global line
            while True:
                if line[0] != chrom: break
                yield line[1:]
                get_next_line()
        self._file.seek(0)
        get_next_line()
        if line:
            while True:
                if line[0] == chrom: break
                chrom = line[0]
                if chrom in self._seen_chr:
                    raise Exception("The track '" + self._path + "' is not sorted by chromosomes (" + chrom + ").")
                if not chrom in self._all_chrs:
                    raise Exception("The track '" + self._path + "' has a value (" + chrom + ") not specified in the chromosome file.")
                self._seen_chr.append(chrom)
                yield chrom, iter_until_different_chr()

    def _write(self):
        # Get fields #
        fields = ['start', 'end'] 
        for f in all_fields_possible[2:]:
            if f in self.fields: fields.append(f)
            else: break
        # Write everything #
        yield self._header_line
        for feature in self.read(fields=fields):
            f = list(feature)
            try:
                f[5] = int_to_strand(feature[5])
            except IndexError:
                pass
            yield '\t'.join([str(x) for x in f]) + '\n'

    #-----------------------------------------------------------------------------#
    @property
    def _type(self): 
        return 'qualitative' 

    @_type.setter
    def _type(self, datatype):
        if not datatype: datatype = 'qualitative'
        if datatype != 'qualitative':
            raise Exception("The track '" + self._path + "' cannot be loaded as a '" + datatype + "' data type.")

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

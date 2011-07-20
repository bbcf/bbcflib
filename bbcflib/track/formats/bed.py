"""
====================================
Submodule: bbcflib.track.formats.bed
====================================

Implementation of the BED format.
"""

# Internal modules #
from ..track_proxy import TrackProxy
from ..track_text import TrackText, strand_to_int, int_to_strand

# Global variables #
all_fields_possible = ['start', 'end', 'name', 'score', 'strand', 'thick_start', 'thick_end',
                       'item_rgb', 'block_count', 'block_sizes', 'block_starts']

###########################################################################
class TrackFormat(TrackText, TrackProxy):
    type_identifier   = 'bed'

    def _read_entries(self):
        self._file.seek(0)
        seen_track = False
        for line in self._file:
            line = line.strip("\n").lstrip()
            if len(line) == 0:              continue
            if line.startswith("#"):        continue
            if line.endswith(" \\"):
                raise Exception("The file '" + self._path + "' includes linebreaks ('\\') which is not supported.")
            if line.startswith("track "):
                if not seen_track:
                    seen_track = True
                    continue
                raise Exception("The file '" + self._path + "' contains a second 'track' directive. This is not supported.")
            if line.startswith("browser "): continue
            line = line.split(self._seperator)
            if len(line) != self.num_fields + 1:
                raise Exception("The file '" + self._path + "' has a varying number of columns. This is not supported.")
            try:
                line[1] = int(line[1])
                line[2] = int(line[2])
            except ValueError:
                raise Exception("The file '" + self._path + "' has non integers as interval bounds and is hence not valid.")
            if line[2] <= line[1]:
                raise Exception("The file '" + self._path + "' has negative or null intervals and is hence not valid.")
            try:
                if line[4] == '.' or line[4] == '': line[4] = 0.0
            except IndexError:
                line.append(0.0)
            try:
                line[4] = float(line[4])
            except ValueError:
                raise Exception("The file '" + self._path + "' has non floats as score values and is hence not valid.")
            try:
                line[5] = strand_to_int(line[5])
            except IndexError:
                line.append(0)
            if len(line) > 6:
                try:
                    line[6] = float(line[6])
                except ValueError:
                    raise Exception("The file '" + self._path + "' has non integers as thick starts and is hence not valid.")
                if len(line) > 7:
                    try:
                        line[7] = float(line[7])
                    except ValueError:
                        raise Exception("The file '" + self._path + "' has non integers as thick ends and is hence not valid.")
            yield line

    def _write_entries(self):
        # Get fields #
        fields = ['start', 'end']
        for f in all_fields_possible[2:]:
            if f in self.fields: fields.append(f)
            else: break
        # Write everything #
        for feature in self.read(fields=fields):
            f = list(feature)
            try:
                f[5] = int_to_strand(feature[5])
            except IndexError:
                pass
            yield '\t'.join([str(x) for x in f]) + '\n'

    #-----------------------------------------------------------------------------#
    @property
    def _fields(self):
        self._file.seek(0)
        for line in self._file:
            line = line.strip("\n").lstrip()
            if len(line) == 0:              continue
            if line.startswith("#"):        continue
            if line.endswith(" \\"):
                raise Exception("The file '" + self._path + "' includes linebreaks ('\\') which are not supported.")
            if line.startswith("track "):   continue
            if line.startswith("browser "): continue
            else:
                if '\t' in line:  self._seperator = '\t'
                else:             self._seperator = ' '
                line = line.split(self._seperator)
            self.num_fields = len(line) - 1
            if self.num_fields < 2:
                raise Exception("The file '" + self._path + "' has less than three columns and is hence not a valid BED file.")
            if self.num_fields > len(all_fields_possible):
                raise Exception("The file '" + self._path + "' has too many columns and is hence not a valid BED file.")
            return all_fields_possible[0:max(5,self.num_fields)]

    #-----------------------------------------------------------------------------#
    @property
    def _datatype(self):
        return 'qualitative'

    @_datatype.setter
    def _datatype(self, datatype):
        if datatype and datatype != 'qualitative':
            raise Exception("The file '" + self._path + "' cannot be loaded as a '" + datatype + "' datatype.")

###########################################################################
def random_track(number_of_features=15000000, size=1000, jump=1000, orig_start=0, chrs=20):
    import random, tempfile
    yield 'track type=bed name="Features" description="Intervals" source="Random generator"\n'
    name_gen = tempfile._RandomNameSequence()
    chr = 0
    for i in range(number_of_features):
        if i % (number_of_features / chrs) == 0:
            chr += 1
            start = orig_start
        start       =   start + (random.randint(0,jump))
        end         =   start + (random.randint(1,size))
        thick_start =   start + (random.randint(-size*0.25,size*0.25))
        thick_end   =   end   + (random.randint(-size*0.25,size*0.25))
        name        = name_gen.next() + name_gen.next()
        strand      = random.random() < 0.5 and '+' or '-'
        score       = random.random()
        line        = ['chr' + str(chr), str(start), str(end), name, score, strand, str(thick_start), str(thick_end)]
        yield ('\t'.join(line) + '\n')

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

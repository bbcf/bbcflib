"""
===================================
Submodule: bbcflib.track.extras.sql
===================================

Implementation of special methods of the Track object specific
to the SQL format and not part of the standard functionality
that all formats should have.
"""

# Built-in modules #
import random

# Internal modules #
from ... import track

###################################################################################
class TrackExtras(object):
    def get_scores_frequencies(self):
        scores = {}
        for chrom in self:
            self.cursor.execute("SELECT count (*), score FROM '" + chrom + "' GROUP BY score")
            for result in self.cursor:
                if result[1] not in scores: scores[result[1]] =  result[0]
                else:                       scores[result[1]] += result[0]
        return scores

    def shuffle_track(self, new_path, repeat_number=1):
        ''' Makes a new track with the number of features as the
        original track.'''
        with track.new(new_path, "sql", name="Random track") as t:
            # Copy meta data #
            t.meta_chr   = self.meta_chr
            t.meta_track = self.meta_track
            # Define generator #
            def shuffle_features(length, data):
                for feature in data:
                    distance     = feature[0] - feature[1]
                    random_start = random.randint(0, length - distance)
                    random_end   = random_start + distance
                    yield (random_start, random_end, "Random feature", 0.0, 0)
            # Iterate #
            for chrom in self:
                for i in range(repeat_number):
                    t.write(chrom, shuffle_features(self.chr_length(chrom), self.read(chrom)))

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
        '''Example output:
        scores = {100.0: 1, 10.0: 3, 1.0: 2}
        '''
        scores = {}
        for chrom in self:
            self.cursor.execute("SELECT count (*), score FROM '" + chrom + "' GROUP BY score")
            for result in self.cursor:
                if result[1] not in scores: scores[result[1]] =  result[0]
                else:                       scores[result[1]] += result[0]
        return scores

    def shuffle_track(self, new_path, repeat_number=1):
        ''' Makes a new track with the same number of features as the
        original track or multiplied with a factor.'''
        with track.new(new_path, "sql", name="Random track") as t:
            # Copy meta data #
            t.chrmeta    = self.chrmeta
            t.attributes = self.attributes
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
                    t.write(chrom, shuffle_features(self.chrmeta[chrom]['length'], self.read(chrom)))

"""
===================================
Submodule: bbcflib.track.extras.sql
===================================

Implementation of special methods of the Track object specific
to the SQL format and not part of the standard functionality
that all formats should have.
"""

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

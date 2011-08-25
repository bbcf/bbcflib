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
from bbcflib import track

###################################################################################
class TrackExtras(object):
    def get_cumul_binned_scores(self, roundoff=1):
        '''Computes cumulative counts for binned scores. 'roundoff' gives the accuracy of the bins, namely
        score*roundoff is an integer
        '''
        counts = {}
        subqry = "(SELECT ROUND(score*%i)/%i AS s, count(*) as c from '%s' group by s)"
        for chrom in self:
            self.cursor.execute("SELECT y.s,sum(x.c) FROM "+subqry %(roundoff,roundoff,chrom)+" AS x, "+subqry %(roundoff,roundoff,chrom)+" AS y WHERE x.s>y.s GROUP BY y.s")
            temp = dict(tuple(x) for x in self.cursor)
            if not(len(temp)): continue
            kmax = max(temp.keys())
            if not(kmax in counts): counts[kmax] = 0
            ckeys = counts.keys()
            k2 = ckeys.pop(0)
            for k in temp.keys():
                while k2<k:
                    k2 = ckeys.pop(0)
                counts[k] = counts[k2]+temp[k]
        return counts

    def get_scores_frequencies(self):
        '''Example output:
        scores = {100.0: 1, 10.0: 3, 1.0: 2}
        '''
        scores = {}
        for chrom in self:
            self.cursor.execute("SELECT count(*), ROUND(score) AS s FROM '" + chrom + "' GROUP BY s")
            for result in self.cursor:
                if result[1] not in scores: scores[result[1]] =  result[0]
                else:                       scores[result[1]] += result[0]
        return scores

    def read_shuffled(self, chrom, repeat_number=1, **kwargs):
        ''' Yields randomly located features of the same length as the original track.'''
        chr_len = self.chrmeta[chrom]['length']
        for i in range(repeat_number):
            for feat in self.read(chrom, **kwargs):
                feat_len = feat[1]-feat[0]
                f0 = random.randint(0, chr_len-feat_len)
                f1 = f0+feat_len
                yield (f0,f1)+feat[2:]

    def shuffle_track(self, new_path, repeat_number=1, **kwargs):
        ''' Makes a new track with 'repeat_number' features for every feature of the
        original track.'''
        with track.new(new_path, "sql", name="Random track") as t:
            t.chrmeta    = self.chrmeta
            t.attributes = self.attributes
            for chrom in self:
                t.write( chrom, self.read_shuffled( chrom, repeat_number, **kwargs ) )


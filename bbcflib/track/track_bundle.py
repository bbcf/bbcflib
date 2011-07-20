"""
======================================
Subpackage: bbcflib.track.track_bundle
======================================

Provides the TrackBundle class that can be useful for grouping many Track objects together and represent them as one.

"""

#############################################################################################
class TrackBundle(object):
    def __init__(self, tracks):
        self.tracks = tracks
        self.name = 'Collection of ' + str(len(tracks)) + ' track' + ((len(tracks) > 1) and 's' or '')
        self.chrs = list(reduce(set.intersection, map(set,[t.chrs for t in tracks])))
        self.chrmeta = tracks[0].chrmeta

    def read(self, selection, fields, cursor): return [t.read(selection, fields, cursor=cursor) for t in self.tracks]

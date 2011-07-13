# Built-in modules #
import os

# Internal modules #
from ..track_collection import track_collections
from ..common import named_temporary_path
from ... import track

# Unittesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

###################################################################################
class Test_Shuffle(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.sql')
        d = track_collections['Validation'][1]
        with track.load(d['path'], chrmeta=d['chrmeta']) as t:
            t.shuffle_track(path, 3)
        with track.load(path, chrmeta=d['chrmeta']) as t:
            self.assertEqual(t.count(), 36)
        os.remove(path)

#-----------------------------------------------------------------------------#
class Test_ScoreFrequencies(unittest.TestCase):
    def runTest(self):
        t = track_collections['Scores'][4]
        with track.load(t['path'], chrmeta=t['chrmeta']) as t:
            freq = t.get_scores_frequencies()
        expected = {100.0: 1, 8.0: 1, 10.0: 3, 9000.0: 1, 40.0: 2, 50.0: 1, 20.0: 3}
        self.assertEqual(freq, expected)

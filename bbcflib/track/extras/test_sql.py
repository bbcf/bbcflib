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
        new = named_temporary_path('.sql')
        old = track_collections['Validation'][1]['path_sql']
        with track.load(old) as t:
            t.shuffle_track(new, 3)
        with track.load(new) as t:
            self.assertEqual(t.count(), 36)
        os.remove(new)

#-----------------------------------------------------------------------------#
class Test_ScoreFrequencies(unittest.TestCase):
    def runTest(self):
        path = track_collections['Scores'][4]['path_sql']
        with track.load(path) as t:
            freq = t.get_scores_frequencies()
        expected = {100.0: 1, 8.0: 1, 10.0: 3, 9000.0: 1, 40.0: 2, 50.0: 1, 20.0: 3}
        self.assertEqual(freq, expected)

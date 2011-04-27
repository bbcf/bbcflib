"""
===============================
Submodule: bbcflib.track._test_
===============================

Tests for the bbcflib.track subpackage.
"""

# Genreral Modules #
import sys, os

# Specific Modules #
from ..track import Track, new

# Unittesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Tracks path #
tracks_path = 'extras/tracks/' 
try:
    tracks_path = os.path.abspath('/'.join(__file__.split('/')[:-1]) + '../../../' + tracks_path) + '/'
except NameError:
    pass

# Tracks collection #
yeast_chr_file = tracks_path + 'chr/yeast.chr' 
track_collections = {
'Validation': {
  '1': {'path':tracks_path+'qual/bed/validation1.bed', 'type':'qualitative', 'fields':Track.qualitative_fields[:4], 'chrfile':yeast_chr_file},
  '2': {'path':tracks_path+'qual/bed/validation2.bed', 'type':'qualitative', 'fields':Track.qualitative_fields,     'chrfile':yeast_chr_file},
  '3': {'path':tracks_path+'qual/bed/validation3.bed', 'type':'qualitative', 'fields':Track.qualitative_fields,     'chrfile':yeast_chr_file},
  },
}

# Modify tracks collection #
for col_name, col in sorted(track_collections.items()):
    for track_name, track in sorted(col.items()):
        # Make the SQL path #
        old_path = track['path']
        old_name = old_path.split('/')[-1]
        new_dir  = '/'.join(old_path.split('/')[:-2]) + '/' + 'sql' + '/'
        new_name = '.'.join(old_name.split('.')[:-1]  + ['sql'])
        track['path_sql'] = new_dir + new_name        

################################################################################### 
class Test_SQL_Read(unittest.TestCase):
    def runTest(self):
        t = track_collections['Validation']['1']
        with Track(t['path_sql']) as t['track']:
            # Just the first feature #
            data = t['track'].read()
            self.assertEqual(data.next(), ('chr1', 0, 10, 'Validation feature 1', 10.0))
            # Number of features #
            data = t['track'].read()
            self.assertEqual(len(list(data)), 12)
            # Different fields #
            data = t['track'].read('chr1', fields=['score'])
            expected = [(10.0,), (0.0,), (10.0,), (0.0,), (0.0,), (10.0,), (10.0,), (10.0,), (10.0,), (10.0,), (10.0,), (5.0,)]
            self.assertEqual(list(data), expected)
            # Region #
            expected =[(20, 30, u'Validation feature 3', 10.0),
                       (25, 30, u'Validation feature 4',  0.0),
                       (40, 45, u'Validation feature 6',  0.0),
                       (40, 50, u'Validation feature 5', 10.0)]
            data = t['track'].read({'chr':'chr1','start':22,'end':46})
            self.assertEqual(list(data), expected)
            # Strict region #
            data = t['track'].read({'chr':'chr1','start':22,'end':46,'inclusion':'strict'})
            self.assertEqual(list(data), expected[1:3])
            # Empty result #
            data = t['track'].read({'chr':'chr2','start':0,'end':10})
            self.assertEqual(list(data), [])

################################################################################### 
class Test_SQL_Write(unittest.TestCase):
    def runTest(self):
        self.assertEqual(1, 1)

################################################################################### 
class Test_SQL_Conversion(unittest.TestCase):
    def runTest(self):
        self.assertEqual(1, 1)

################################################################################### 
class Test_SQL_Creation(unittest.TestCase):
    def runTest(self):
        self.assertEqual(1, 1)

################################################################################### 
Test_SQL_Read().runTest()

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#


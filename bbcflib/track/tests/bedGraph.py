# Built-in modules #
import os, shutil, filecmp

# Internal modules #
from ... import track
from ..common import named_temporary_path, sqlcmp
from ..track_collection import track_collections, yeast_chr_file

# Unittesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

###################################################################################
class Test_Read(unittest.TestCase):
    def runTest(self):
        t = track_collections['Signals']['A']
        with track.load(t['path']) as t['track']:
            # Just the first feature #
            data = t['track'].read()
            self.assertEqual(data.next(), ('chr1', 0, 10, -1.0))
            # Number of features #
            data = t['track'].read()
            self.assertEqual(len(list(data)), 3)

#-----------------------------------------------------------------------------#
class Test_Write(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bedGraph')
        with track.new(path) as t:
            self.assertEqual(t.datatype, 'quantitative')
            features = {}
            features['chr1'] = [(0,  10, -1.0),
                                (20, 30, -1.75),
                                (40, 50, -2.5)]
            for chrom, data in sorted(features.items()):
                t.write(chrom, data)
        with open(path,                                     'r') as f: A = f.read().split('\n')
        with open(track_collections['Signals']['A']['path'],'r') as f: B = f.read().split('\n')
        self.assertEqual(A[1:], B)
        os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Roundtrips(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bedGraph')
        for track_num, d in sorted(track_collections['Signals'].items()):
            with track.load(d['path']) as t:
                t.dump(path)
            with open(path,              'r') as f: A = f.read().split('\n')
            with open(d['path'],         'r') as f: B = f.read().split('\n')
            self.assertEqual(A[1:], B)
            os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Format(unittest.TestCase):
    def runTest(self):
        # Not specified #
        t = track_collections['Signals']['A']
        with track.load(t['path']) as t:
            self.assertEqual(t.format, 'bedGraph')
        # No extension #
        old = track_collections['Signals']['A']['path']
        new = named_temporary_path()
        shutil.copyfile(old, new)
        with track.load(new, 'bedGraph') as t:
            self.assertEqual(t.format, 'bedGraph')
        os.remove(new)

#-------------------------------------------------------------------------------#
class Test_Convert(unittest.TestCase):
    def runTest(self):
        # Paths #
        path_ref_bg  = track_collections['Signals'][1]['path']
        path_new_bg  = named_temporary_path('.bedGraph')
        path_ref_wig = track_collections['Scores'][1]['path']
        path_new_wig = named_temporary_path('.wig')
        path_ref_sql = track_collections['Scores'][1]['path_sql']
        path_new_sql = named_temporary_path('.sql')
        # Case 1: BEDGRAPH to WIG #
        with track.load(path_ref_bg) as t:
            t.convert(path_new_wig)
            self.assertEqual(t.format, 'wig')
        with open(path_new_wig, 'r') as f: A = f.read().split('\n')
        with open(path_ref_wig, 'r') as f: B = f.read().split('\n')
        # Case 2: BEDGRAPH to SQL #
        with track.load(path_ref_bg, chrmeta=yeast_chr_file) as t:
            t.convert(path_new_sql)
            self.assertEqual(t.format, 'sql')
        self.assertTrue(sqlcmp(path_new_sql, path_ref_sql))
        # Case 3: SQL to BEDGRAPH #
        with track.load(path_ref_sql) as t:
            t.convert(path_new_bg)
            self.assertEqual(t.format, 'bedGraph')
        with open(path_new_bg, 'r') as f: A = f.read().split('\n')
        with open(path_ref_bg, 'r') as f: B = f.read().split('\n')
        self.assertEqual(A[1:], B)
        # Case 4: WIG to BEDGRAPH #
        os.remove(path_new_bg)
        with track.load(path_ref_wig) as t:
            t.convert(path_new_bg)
            self.assertEqual(t.format, 'bedGraph')
        with open(path_new_bg, 'r') as f: A = f.read().split('\n')
        with open(path_ref_bg, 'r') as f: B = f.read().split('\n')
        self.assertEqual(A[1:], B)
        # Cleanup #
        os.remove(path_new_bg)
        os.remove(path_new_wig)
        os.remove(path_new_sql)

#-------------------------------------------------------------------------------#
class Test_Export(unittest.TestCase):
    def runTest(self):
        # Paths #
        path_ref_bg  = track_collections['Signals'][1]['path']
        path_new_bg  = named_temporary_path('.bedGraph')
        path_ref_wig = track_collections['Scores'][1]['path']
        path_new_wig = named_temporary_path('.wig')
        path_ref_sql = track_collections['Scores'][1]['path_sql']
        path_new_sql = named_temporary_path('.sql')
        # Case 1: BEDGRAPH to WIG #
        with track.load(path_ref_bg) as t:
            t.export(path_new_wig)
        with open(path_new_wig, 'r') as f: A = f.read().split('\n')
        with open(path_ref_wig, 'r') as f: B = f.read().split('\n')
        # Case 2: BEDGRAPH to SQL #
        with track.load(path_ref_bg, chrmeta=yeast_chr_file) as t:
            t.export(path_new_sql)
        self.assertTrue(sqlcmp(path_new_sql, path_ref_sql))
        # Case 3: SQL to BEDGRAPH #
        with track.load(path_ref_sql) as t:
            t.export(path_new_bg)
        with open(path_new_bg, 'r') as f: A = f.read().split('\n')
        with open(path_ref_bg, 'r') as f: B = f.read().split('\n')
        self.assertEqual(A[1:], B)
        # Case 4: WIG to BEDGRAPH #
        os.remove(path_new_bg)
        with track.load(path_ref_wig) as t:
            t.export(path_new_bg)
        with open(path_new_bg, 'r') as f: A = f.read().split('\n')
        with open(path_ref_bg, 'r') as f: B = f.read().split('\n')
        self.assertEqual(A[1:], B)
        # Cleanup #
        os.remove(path_new_bg)
        os.remove(path_new_wig)
        os.remove(path_new_sql)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

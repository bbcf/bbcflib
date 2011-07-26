# Built-in modules #
import os, filecmp, shutil

# Internal modules #
from ... import track
from ..common import named_temporary_path
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
        t = track_collections['Binary']['A']
        with track.load(t['path']) as t:
            # Just the first feature #
            data = t.read()
            self.assertEqual(data.next(), ('chr1', 0, 10, -1.0))
            # Number of features #
            data = t.read()
            self.assertEqual(len(list(data)), 3)

#-----------------------------------------------------------------------------#
class Test_Write(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bigWig')
        with track.new(path, chrmeta=yeast_chr_file) as t:
            self.assertEqual(t.datatype, 'quantitative')
            features = {}
            features['chr1'] = [(0,  10, -1.0),
                                (20, 30, -1.75),
                                (40, 50, -2.5)]
            for chrom, data in sorted(features.items()): t.write(chrom, data)
        self.assertTrue(filecmp.cmp(path, track_collections['Binary']['A']['path']))
        os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Overwrite(unittest.TestCase):
    def runTest(self):
        old_path = track_collections['Binary']['B']['path']
        new_path = named_temporary_path('.bigWig')
        shutil.copyfile(old_path, new_path)
        feature = (10, 20, 9999.0)
        chrom = 'chr2'
        with track.load(new_path, chrmeta=yeast_chr_file) as t:
            t.write(chrom, (feature,))
        with track.load(new_path, chrmeta=yeast_chr_file) as t:
            self.assertEqual(feature, t.read(chrom).next())
        os.remove(new_path)

#-----------------------------------------------------------------------------#
class Test_Roundtrips(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bigWig')
        for i in ('A','B'):
            d = track_collections['Binary'][i]
            with track.load(d['path'], chrmeta=yeast_chr_file) as t: t.dump(path)
            self.assertTrue(filecmp.cmp(d['path'], path))
            os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Format(unittest.TestCase):
    def runTest(self):
        # Not specified #
        t = track_collections['Binary']['A']
        with track.load(t['path']) as t:
            self.assertEqual(t.format, 'bigWig')
        # No extension #
        old = track_collections['Binary']['B']['path']
        new = named_temporary_path()
        shutil.copyfile(old, new)
        with track.load(new, 'bigWig') as t:
            self.assertEqual(t.format, 'bigWig')
        os.remove(new)

#-------------------------------------------------------------------------------#
class Test_Convert(unittest.TestCase):
    def runTest(self):
        # Paths #
        path_ref_bw  = track_collections['Binary'][1]['path']
        path_new_bw  = named_temporary_path('.bigWig')
        path_ref_wig = track_collections['Scores'][1]['path']
        path_new_wig = named_temporary_path('.wig')
        path_ref_sql = track_collections['Scores'][1]['path_sql']
        path_new_sql = named_temporary_path('.sql')
        # Case 1b: BIGWIG to WIG #
        with track.load(path_ref_bw) as t:
            t.convert(path_new_wig)
            self.assertEqual(t.format, 'wig')
        with open(path_new_wig, 'r') as f: A = f.read().split('\n')
        with open(path_ref_wig, 'r') as f: B = f.read().split('\n')
        # Case 2b: BIGWIG to SQL #
        with track.load(path_ref_bw, chrmeta=yeast_chr_file) as t:
            t.convert(path_new_sql)
            self.assertEqual(t.format, 'sql')
        self.assertTrue(filecmp.cmp(path_new_sql, path_ref_sql))
        # Case 3b: SQL to BIGWIG #
        with track.load(path_ref_sql) as t:
            t.convert(path_new_bw)
            self.assertEqual(t.format, 'bigWig')
        self.assertTrue(filecmp.cmp(path_new_bw, path_ref_bw))
        # Case 4b: WIG to BIGWIG #
        os.remove(path_new_bw)
        with track.load(path_ref_wig, chrmeta=yeast_chr_file) as t:
            t.convert(path_new_bw)
            self.assertEqual(t.format, 'bigWig')
        self.assertTrue(filecmp.cmp(path_new_bw, path_ref_bw))
        # Cleanup #
        os.remove(path_new_bw)
        os.remove(path_new_wig)
        os.remove(path_new_sql)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

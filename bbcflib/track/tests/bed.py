# Built-in modules #
import os, shutil

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

################################################################################
class Test_Read(unittest.TestCase):
    def runTest(self):
        t = track_collections['Validation'][1]
        with track.load(t['path']) as t['track']:
            # Just the first feature #
            data = t['track'].read()
            self.assertEqual(data.next(), ('chr1', 0, 10, 'Validation feature 1', 10.0, 0))
            # Number of features #
            data = t['track'].read()
            self.assertEqual(len(list(data)), 12)

#------------------------------------------------------------------------------#
class Test_Write(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bed')
        with track.new(path) as t:
            self.assertEqual(t.datatype, 'qualitative')
            features = {}
            features['chr1'] = [(10, 20, 'Lorem', 1.0, 1),
                                (30, 40, 'Ipsum', 2.0, 1)]
            features['chr2'] = [(10, 20, 'Dolor', 3.0, 1)]
            for chrom, data in sorted(features.items()):
                t.write(chrom, data)
        with open(path,                                      'r') as f: A = f.read().split('\n')
        with open(track_collections['Validation'][4]['path'],'r') as f: B = f.read().split('\n')
        self.assertEqual(A[1:], B)
        os.remove(path)

#------------------------------------------------------------------------------#
class Test_Overwrite(unittest.TestCase):
    def runTest(self):
        old_path = track_collections['Validation'][1]['path']
        new_path = named_temporary_path('.bed')
        shutil.copyfile(old_path, new_path)
        feature = (10, 20, 'Dolor', 3.0, 1)
        chrom = 'chr2'
        with track.load(new_path) as t:
            t.write(chrom, (feature,))
        with track.load(new_path) as t:
            self.assertEqual(feature, t.read(chrom).next())
        os.remove(new_path)

#------------------------------------------------------------------------------#
class Test_Roundtrips(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bed')
        for i in (2,3,4):
            track_dict = track_collections['Validation'][i]
            with track.load(track_dict['path']) as t: t.dump(path)
            with open(path,              'r') as f: A = f.read().split('\n')
            with open(track_dict['path'],'r') as f: B = f.read().split('\n')
            self.assertEqual(A[1:], B)
            os.remove(path)

#------------------------------------------------------------------------------#
class Test_Chrmeta(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.bed')
        info = {'chr1': {'length': 197195432}, 'chr2': {'length': 129993255}}
        # Setting #
        with track.new(path) as t:
            t.chrmeta = info
            self.assertEqual(t.chrmeta, info)
        os.remove(path)
        # Dictionary #
        d = track_collections['Validation'][1]
        with track.load(d['path'], chrmeta=info, readonly=True) as t:
            self.assertEqual(t.chrmeta['chr1']['length'], 197195432)
        # Genrep #
        with track.load(d['path'], chrmeta='hg19', readonly=True) as t:
            self.assertEqual(t.chrmeta['chr1']['length'], 249250621)
        # File #
        with track.load(d['path'], chrmeta=yeast_chr_file, readonly=True) as t:
            self.assertEqual(t.chrmeta['chr1']['length'], 230208)

#------------------------------------------------------------------------------#
class Test_Format(unittest.TestCase):
    def runTest(self):
        # Not specified #
        t = track_collections['Validation'][1]
        with track.load(t['path']) as t:
            self.assertEqual(t.format, 'bed')
        # No extension #
        old = track_collections['Validation'][1]['path']
        new = named_temporary_path()
        shutil.copyfile(old, new)
        with track.load(new, 'bed') as t:
            self.assertEqual(t.format, 'bed')
        os.remove(new)
        # Only track line #
        old = track_collections['Yeast']['RP genes']['path']
        new = named_temporary_path()
        shutil.copyfile(old, new)
        with track.load(new) as t:
            self.assertEqual(t.format, 'bed')
        os.remove(new)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

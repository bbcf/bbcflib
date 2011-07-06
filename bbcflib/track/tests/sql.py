# Built-in modules #
import os

# Internal modules #
from ... import track
from ..common import named_temporary_path
from ..track_collection import track_collections

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
        t = track_collections['Validation'][1]
        with track.load(t['path_sql']) as t['track']:
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
                       (40, 45, u'Validation feature 5',  0.0),
                       (40, 50, u'Validation feature 6', 10.0)]
            data = t['track'].read({'chr':'chr1','start':22,'end':46})
            self.assertEqual(list(data), expected)
            # Strict region #
            data = t['track'].read({'chr':'chr1','start':22,'end':46,'inclusion':'strict'})
            self.assertEqual(list(data), expected[1:3])
            # Empty result #
            data = t['track'].read({'chr':'chr2','start':0,'end':10})
            self.assertEqual(list(data), [])

#-----------------------------------------------------------------------------#
class Test_Creation(unittest.TestCase):
    def runTest(self):
        format = 'sql'
        path = named_temporary_path('.' + format)
        # With format #
        with track.new(path=path, format=format) as t:
            self.assertEqual(t.meta_track, {'datatype':'qualitative'})
        os.remove(path)
        # Without format #
        with track.new(path=path) as t:
            self.assertEqual(t.format, format)
        os.remove(path)
        # Different datatype #
        with track.new(path=path, format=format, datatype='quantitative') as t:
            self.assertEqual(t.meta_track, {'datatype':'quantitative'})
        os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Write(unittest.TestCase):
    def runTest(self):
        format = 'sql'
        path = named_temporary_path('.' + format)
        with track.new(path) as t:
            # Single feature #
            chrom = '9toto'
            features = [(10, 20, 'A', 0.0, 1)]
            t.write(chrom, features)
            self.assertEqual(list(t.read(chrom)), features)
            # Addition #
            features.append((10, 20, 'B', 0.0, 1))
            t.write(chrom, [features[1]])
            self.assertEqual(list(t.read(chrom)), features)
            # Multi feature #
            chrom = 'chr2'
            features = [(i*10, i*10+5, 'X', 0.0, 0) for i in xrange(5)]
            t.write(chrom, features)
            self.assertEqual(list(t.read(chrom)), features)
            # Memory Copy #
            t.write('chrX', list(t.read(chrom)))
            self.assertEqual(list(t.read('chrX')), list(t.read(chrom)))
            # Live copy #
            t.write('chrY', t.read(chrom, cursor=True))
            self.assertEqual(list(t.read('chrY')), list(t.read(chrom)))
            # More fields #
            chrom = 'chr3'
            features = [(10, 20, 'A', 0.0, 1, 8, 22)]
            t.write(chrom, features, fields=track.Track.qualitative_fields + ['thick_start', 'thick_end'])
            self.assertEqual(list(t.read(chrom)), features)
        os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Count(unittest.TestCase):
    def runTest(self):
        t = track_collections['Yeast']['All genes']
        with track.load(t['path_sql']) as t['track']:
            num = t['track'].count()
            self.assertEqual(num, 6717)

#-----------------------------------------------------------------------------#
class Test_Meta(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.sql')
        with track.new(path) as t:
            # Chromosome #
            info = [{'name': 'chr1', 'length': 1500}, {'name': 'chr2', 'length': 2000}]
            t.meta_chr = info
            self.assertEqual(t.meta_chr, info)
            # Track attributes #
            info = {'datatype': 'quantitative', 'source': 'SGD'}
            t.meta_track = info
            self.assertEqual(t.meta_track, info)
            t.meta_track.update({'k':'v'})
            self.assertEqual(t.meta_track, info)
        os.remove(path)

#-----------------------------------------------------------------------------#
class Test_Remove(unittest.TestCase):
    def runTest(self):
        path = named_temporary_path('.sql')
        with track.new(path) as t:
            chrom = 'chr1'
            features = [(i*10, i*10+5, 'X', 0.0, 0) for i in xrange(5)]
            t.write(chrom, features)
            t.remove()
            self.assertEqual(list(t.read()), [])
        os.remove(path)

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

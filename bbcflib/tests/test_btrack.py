# Built-in modules #
import cStringIO, tempfile, os, shutil

# Internal modules #
from bbcflib.genrep import Assembly
from bbcflib.btrack import track, convert, map_chromosomes, score_threshold
from bbcflib.btrack import FeatureStream, text

# Unitesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

# Path to testing files
path = "test_data/btrack/"


class Test_Track(unittest.TestCase):
    def setUp(self):
        self.assembly_id = 'sacCer2'
        self.bed = os.path.join(path,"yeast_genes.bed")

    def test_track(self):
        t = track(self.bed)
        self.assertIsInstance(t, text.BedTrack)
        # Check all attributes
        for attr in ['path','filehandle','format','fields','assembly','chrmeta','info']:
            self.assertTrue(hasattr(t,attr))
        for attr in ['read','write','open','close','save']:
            self.assertTrue(hasattr(t,attr))
        self.assertEqual(t.path, self.bed)
        self.assertIsInstance(t.filehandle, file)
        self.assertEqual(t.format, 'bed')
        self.assertIsInstance(t.fields, list)
        self.assertEqual(t.assembly, None)
        self.assertIsInstance(t.chrmeta, dict)
        self.assertIsInstance(t.info, dict)

    def test_read(self):
        t = track(self.bed)
        content = t.read()
        self.assertIsInstance(content, FeatureStream)


class Test_Operations(unittest.TestCase):
    def setUp(self):
        pass

    def test_map_chromosomes(self):
        pass

    def test_score_threshold(self):
        pass


class Test_Formats(unittest.TestCase):
    def setUp(self):
        self.assembly_id = 'sacCer2'
        self.bed = os.path.join(path,"yeast_genes.bed")

    def test_bed(self):
        no_extension = self.bed.split('.bed')[0]
        shutil.copy(self.bed, no_extension)
        t = track(no_extension, format='bed')
        self.assertIsInstance(t, text.BedTrack)
        self.assertEqual(t.format,'bed')
        self.assertListEqual(t.fields, ['chr','start','end','name','score','strand'])
        os.remove(no_extension)

    def test_bedgraph(self):
        bg = os.path.join(path,'test.bedgraph')
        convert(self.bed,bg)
        print bg
        os.remove(bg)
        raise

    def test_wig(self):
        pass

    def test_bigwig(self):
        pass

    def test_gff(self):
        pass

    def test_sga(self):
        pass


# Built-in modules #
import cStringIO, tempfile, os, shutil

# Internal modules #
from bbcflib.genrep import Assembly
from bbcflib.btrack import track, convert, map_chromosomes, score_threshold
from bbcflib.btrack import FeatureStream
from bbcflib.btrack.text import BedTrack, BedGraphTrack, WigTrack, SgaTrack, GffTrack
from bbcflib.btrack.bin import BigWigTrack, BamTrack
from bbcflib.btrack.sql import SqlTrack

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
        self.assertIsInstance(t, BedTrack)
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
        self.assertIsInstance(t, BedTrack)
        self.assertEqual(t.format,'bed')
        self.assertListEqual(t.fields, ['chr','start','end','name','score','strand'])
        os.remove(no_extension)

    def test_bedgraph(self):
        bg = os.path.join(path,'test.bedgraph')
        bg_track = convert(self.bed, bg)
        self.assertIsInstance(bg_track, BedGraphTrack)
        os.remove(bg)

    def test_wig(self):
        wig = os.path.join(path,'test.wig')
        wig_track = convert(self.bed, wig)
        self.assertIsInstance(wig_track, WigTrack)
        os.remove(wig)

    #@unittest.skip('Format not implemented yet')
    def test_bigwig(self):
        bw = os.path.join(path,'test.bw')
        bw_track = convert(self.bed, bw)
        self.assertIsInstance(bw_track, BigWigTrack)
        os.remove(bw)

    def test_gff(self):
        gff = os.path.join(path,'test.gff')
        gff_track = convert(self.bed, gff)
        self.assertIsInstance(gff_track, GffTrack)
        os.remove(gff)

    @unittest.skip('Fails')
    def test_sga(self):
        sga = os.path.join(path,'test.sga')
        sga_track = convert(self.bed, sga)
        self.assertIsInstance(sga_track, SgaTrack)
        os.remove(sga)


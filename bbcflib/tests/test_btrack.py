# Built-in modules #
import os, shutil, time

# Internal modules #
from bbcflib.btrack import track, convert, FeatureStream
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
        self.assembly = 'sacCer2'
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

    def test_index(self):
        # Temp file containing only chrI
        tempfile = "temp.txt"
        g = open(tempfile,'wb')
        g.writelines(["chrI\t1\n"]*200)
        g.close()

        # Never read chrI from this track: the whole file is read at each iteration.
        t = track(tempfile,fields=['chr','end'],chrmeta=self.assembly)
        tnorm = tskip = 0
        for chr in ['2micron','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX']:
            t1 = time.time()
            s = t.read(chr, skip=False)
            for x in s: pass
            t2 = time.time()
            tnorm += t2-t1
        t.close()

        # Read chrI, then read others: chrI (the whole file) is skipped at each iteration.
        t = track(tempfile,fields=['chr','end'],chrmeta=self.assembly)
        for chr in ['chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX']:
            t1 = time.time()
            s = t.read(chr, skip=True)
            for x in s: pass
            t2 = time.time()
            tskip += t2-t1
        t.close()

        self.assertLess(tskip/tnorm,0.5) # Second case it at least twice faster
        os.remove(tempfile)

class Test_Formats(unittest.TestCase):
    """Converting from bed to every other available format, using btrack.convert."""
    def setUp(self):
        self.assembly = 'sacCer2'
        self.bed = os.path.join(path,"yeast_genes.bed")

    def test_bed(self): # as general TextTrack
        no_extension = self.bed.split('.bed')[0] # guess extension from header
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
        s = bg_track.read(); s.next()
        os.remove(bg)

    def test_wig(self):
        wig = os.path.join(path,'test.wig')
        wig_track = convert(self.bed, wig)
        self.assertIsInstance(wig_track, WigTrack)
        s = wig_track.read(); s.next()
        os.remove(wig)

    @unittest.skip('Works but creates temp files when testing (tempfile)')
    def test_bigwig(self):
        try:
            bw = os.path.join(path,'test.bw')
            bw_track = convert(self.bed, bw)
            self.assertIsInstance(bw_track, BigWigTrack)
            s = bw_track.read(); s.next()
            os.remove(bw)
        except OSError: pass

    @unittest.skip('Converting to bam is not implemented yet.')
    def test_bam(self):
        bam = os.path.join(path,'test.bam')
        bam_track = convert(self.bed, bam)
        self.assertIsInstance(bam_track, BamTrack)
        s = bam_track.read(); s.next()
        os.remove(bam)

    def test_gff(self):
        gff = os.path.join(path,'test.gff')
        gff_track = convert(self.bed, gff)
        self.assertIsInstance(gff_track, GffTrack)
        #s = gff_track.read(); s.next() # problems with empty fields after conversion
        os.remove(gff)

    def test_sql(self):
        sql = os.path.join(path,'test.sql')
        sql_track = convert(self.bed, sql, chrmeta=self.assembly)
        self.assertIsInstance(sql_track, SqlTrack)
        s = sql_track.read(cursor=True); s.next()
        os.remove(sql)

    def test_sga(self):
        sga = os.path.join(path,'test.sga')
        sga_track = convert(self.bed, sga)
        self.assertIsInstance(sga_track, SgaTrack)
        s = sga_track.read(); s.next()
        os.remove(sga)


class Test_Bam(unittest.TestCase):
    def setUp(self):
        self.assembly = 'sacCer2'
        self.bam = os.path.join(path,'yeast3_chrV_150k-175k.bam')

    def test_coverage(self):
        t = track(self.bam)
        res = t.coverage(region=('chrV',160000,160002))
        expected = [('chrV',160000,160001,3),('chrV',160001,160002,3)]
        self.assertListEqual(list(res),expected)

    def test_count(self):
        t = track(self.bam)
        res = t.count(regions=[('chrV',150000,175000)])
        expected = [('chrV',150000,175000,2514)]
        self.assertListEqual(list(res),expected)


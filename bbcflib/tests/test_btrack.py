# Built-in modules #
import os, shutil, time

# Internal modules #
from bbcflib.btrack import track, convert, FeatureStream, check_ordered
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
        s = t.read(); s.next()
        self.assertIsInstance(s, FeatureStream)

        # zipped file
        t = track(self.bed+'.gz')
        s = t.read(); s.next()
        self.assertIsInstance(s, FeatureStream)

    def test_index(self):
        # Temp file containing only chrX
        tempfile = os.path.join(path,"temp1.txt")
        with open(tempfile,'wb') as g:
            g.writelines(["chrX\t1\n"]*200)

        # Read chrX at last: the whole file is read at each iteration.
        t = track(tempfile,fields=['chr','end'],chrmeta=self.assembly)
        tnorm = tskip = nnorm = nskip = 0
        for chr in ['chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX']:
            t1 = time.time()
            s = t.read(chr, skip=False)
            for x in s: nnorm += 1
            t2 = time.time()
            tnorm += t2-t1
        t.close()

        # Read chrX, then read others: chrX (the whole file) is skipped at each iteration.
        t = track(tempfile,fields=['chr','end'],chrmeta=self.assembly)
        for chr in ['chrX','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrI']:
            t1 = time.time()
            s = t.read(chr, skip=True)
            for x in s: nskip += 1
            t2 = time.time()
            tskip += t2-t1
        t.close()

        self.assertEqual(nnorm,nskip)
        self.assertLess(tskip/tnorm,0.5) # Second case it at least twice faster

    def test_check_ordered(self):
        # All sorted
        tempfile = os.path.join(path,'temp2.txt')
        with open(tempfile,'wb') as g:
            g.writelines(["chrX\t1\t2\n", "chrX\t3\t4\n", "chrV\t1\t2\n"])
        self.assertEqual(check_ordered(tempfile),True)

        # Chromosomes unsorted
        tempfile = os.path.join(path,'temp3.txt')
        with open(tempfile,'wb') as g:
            g.writelines(["chrX\t1\t2\n", "chrV\t1\t2\n", "chrX\t3\t4\n"])
        self.assertEqual(check_ordered(tempfile),False)

        # Starts unsorted
        tempfile = os.path.join(path,'temp4.txt')
        with open(tempfile,'wb') as g:
            g.writelines(["chrX\t3\t4\n", "chrX\t1\t2\n", "chrV\t1\t2\n"])
        self.assertEqual(check_ordered(tempfile),False)

        # Ends unsorted
        tempfile = os.path.join(path,'temp5.txt')
        with open(tempfile,'wb') as g:
            g.writelines(["chrX\t1\t4\n", "chrX\t1\t3\n", "chrV\t1\t2\n"])
        self.assertEqual(check_ordered(tempfile),False)

    def test_chr_loop(self):
        tempfile = os.path.join(path,'temp6.txt')
        t = track(self.bed)
        out = track(tempfile, fields=t.fields)
        print out.fields
        for chr in ['chrII','chrIII','chrIV']:
            s = t.read(chr)
            out.write(s)

    def tearDown(self):
        for test_file in ['temp1','temp2','temp3','temp4','temp5','temp6']:
            test_file = os.path.join(path,test_file)+'.txt'
            if os.path.exists(test_file): os.remove(test_file)

class Test_Formats(unittest.TestCase):
    """Converting from bed to every other available format, using btrack.convert."""
    def setUp(self):
        self.assembly = 'sacCer2'
        self.bed = os.path.join(path,"yeast_genes.bed")
        self.fields = ['chr','start','end','name','score','strand']

    def test_bed(self): # as general TextTrack
        shutil.copy(self.bed, os.path.join(path,'test')) # guess extension from header
        t = track('test', format='bed', fields=self.fields)
        self.assertIsInstance(t, BedTrack)
        self.assertEqual(t.format,'bed')
        self.assertListEqual(t.fields, self.fields)

    def test_bedgraph(self):
        bg = os.path.join(path,'test.bedgraph')
        t = convert(self.bed, bg)
        self.assertIsInstance(t, BedGraphTrack)
        s = t.read(); s.next()
        self.assertListEqual(t.fields, ['chr','start','end','score'])

    def test_wig(self):
        wig = os.path.join(path,'test.wig')
        t = convert(self.bed, wig)
        self.assertIsInstance(t, WigTrack)
        s = t.read(); s.next()
        self.assertListEqual(t.fields, ['chr','start','end','score'])

    @unittest.skip('Works but creates temp files when testing (tempfile)')
    def test_bigwig(self):
        try:
            bw = os.path.join(path,'test.bw')
            t = convert(self.bed, bw)
            self.assertIsInstance(t, BigWigTrack)
            s = t.read(); s.next()
            self.assertListEqual(t.fields, self.fields)
        except OSError: pass

    @unittest.skip('Converting to bam is not implemented yet.')
    def test_bam(self):
        bam = os.path.join(path,'test.bam')
        t = convert(self.bed, bam)
        self.assertIsInstance(t, BamTrack)
        s = t.read(); s.next()
        self.assertListEqual(t.fields, self.fields)

    def test_gff(self):
        gff = os.path.join(path,'test.gff')
        t = convert(self.bed, gff)
        self.assertIsInstance(t, GffTrack)
        self.assertListEqual(t.fields, ['chr','source','name','start','end','score','strand','frame','attributes'])

    def test_sql(self):
        sql = os.path.join(path,'test.sql')
        t = convert(self.bed, sql, chrmeta=self.assembly)
        self.assertIsInstance(t, SqlTrack)
        s = t.read(); s.next()
        self.assertListEqual(t.fields, self.fields)

    def test_sga(self):
        sga = os.path.join(path,'test.sga')
        t = convert(self.bed, sga)
        self.assertIsInstance(t, SgaTrack)
        s = t.read(); s.next()
        self.assertListEqual(t.fields, ['chr','start','end','name','strand','counts'])

    def tearDown(self):
        for ext in ['','.bed','.bw','.wig','.bedGraph','.bam','.sql','.sga','.gff']:
            test_file = os.path.join(path,'test'+ext)
            if os.path.exists(test_file): os.remove(test_file)

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


class Test_Conversions(unittest.TestCase):
    def setUp(self):
        self.assembly = 'sacCer2'
        self.bam = os.path.join(path,'yeast3_chrV_150k-175k.bam')

    def test_to_sql(self):
        pass

# Built-in modules #
import os, shutil, time

# Internal modules #
from bbcflib.track import track, convert, FeatureStream, check_ordered, check_format
from bbcflib.track.text import BedTrack, BedGraphTrack, WigTrack, SgaTrack, GffTrack
from bbcflib.track.bin import BigWigTrack, BamTrack
from bbcflib.track.sql import SqlTrack

# Unitesting module #
try:
    import unittest2 as unittest
    assert unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

# Path to testing files
path = "test_data/track/"


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
        s = t.read()
        self.assertIsInstance(s, FeatureStream)
        x = s.next()
        y = t.readline()
        self.assertEqual(x,y)

        # zipped file
        t = track(self.bed+'.gz')
        s = t.read(); s.next()
        self.assertIsInstance(s, FeatureStream)

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

    def test_check_format(self):
        self.assertTrue(check_format(self.bed))

    def test_chr_loop(self):
        tempfile = os.path.join(path,'temp6.txt')
        t = track(self.bed)
        out = track(tempfile, fields=t.fields)
        for chr in ['chrII','chrIII','chrIV']:
            s = t.read(chr)
            out.write(s)

    def test_get_chrmeta(self):
        t = track(self.bed,chrmeta=self.assembly)
        self.assertEqual(t.chrmeta['chrV'],{'length':576869, 'ac':'2508_NC_001137.2'})
        # "guess"
        t = track(self.bed,chrmeta="guess")
        self.assertEqual(t.chrmeta, {'chrII':{'length':607135}, 'chrIII':{'length':178216},
                                     'chrIV':{'length':1402556}} )
    def tearDown(self):
        for test_file in ['temp1','temp2','temp3','temp4','temp5','temp6']:
            test_file = os.path.join(path,test_file)+'.txt'
            if os.path.exists(test_file): os.remove(test_file)

class Test_Formats(unittest.TestCase):
    """Converting from bed to every other available format, using track.convert."""
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
        bg = os.path.join(path,'test.bedGraph')
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
        self.assertListEqual(t.fields, ['chr','start','end','name','strand','score'])

    def tearDown(self):
        for ext in ['','.bed','.bw','.wig','.bedGraph','.bam','.sql','.sga','.gff']:
            test_file = os.path.join(path,'test'+ext)
            if os.path.exists(test_file): os.remove(test_file)


class Test_Index(unittest.TestCase):
    def setUp(self):
        self.assembly = 'sacCer2'
        self.bed = os.path.join(path,"yeast_genes.bed")

    def _test_index(self, filename):
        t = track(filename,fields=['chr','end'],chrmeta=self.assembly)
        tnorm = tskip = 0
        for chr in t.chrmeta:
            t1 = time.time()
            s = t.read(chr, skip=False)
            norm = list(s)
            t2 = time.time()
            tnorm += t2-t1
        for chr in t.chrmeta:
            t1 = time.time()
            s = t.read(chr, skip=True)
            skip = list(s)
            t2 = time.time()
            tskip += t2-t1
        self.assertListEqual(norm,skip)
        self.assertLess(tskip,tnorm) # Second case it faster
        return t

    def test_index_bed(self):
        t = self._test_index(self.bed)
        self.assertEqual(t.index, {'chrII':[41, 360],'chrIII':[360, 393],'chrIV':[393, 994]})

    def test_index_wig(self):
        wig = os.path.join(path,'test.wig')
        convert(self.bed, wig)
        t = self._test_index(wig)

    def test_index_sga(self):
        sga = os.path.join(path,'test.sga')
        convert(self.bed, sga)
        t = self._test_index(sga)

    def tearDown(self):
        for ext in ['.wig','.sga']:
            test_file = os.path.join(path,'test'+ext)
            if os.path.exists(test_file): os.remove(test_file)


class Test_Bam(unittest.TestCase):
    def setUp(self):
        self.assembly = 'sacCer2'
        self.bam = os.path.join(path,'yeast3_chrV_150k-175k.bam')

    def test_coverage(self):
        t = track(self.bam)
        res = t.coverage(region=('chrV',160000,160002))
        expected = [('chrV',160000,160002,3)]
        self.assertListEqual(list(res),expected)

    def test_count(self):
        t = track(self.bam)
        res = t.count(regions=[('chrV',150000,175000)])
        expected = [('chrV',150000,175000,2512)]
        self.assertListEqual(list(res),expected)


class Test_Conversions(unittest.TestCase):
    def setUp(self):
        self.assembly = 'sacCer2'
        self.bam = os.path.join(path,'yeast3_chrV_150k-175k.bam')

    def test_to_sql(self):
        pass


class Test_Header(unittest.TestCase):
    def setUp(self):
        self.assembly = 'sacCer2'
        self.bed = os.path.join(path,"yeast_genes.bed")

    def test_skip_header(self):
        t = track(self.bed) # skips the first line by default (header=None)
        L1 = len([line for line in t.read()])
        t = track(self.bed, header=None) # same
        L11 = len([line for line in t.read()])
        t = track(self.bed, header=True) # skips the first line
        L111 = len([line for line in t.read()])
        t = track(self.bed, header=5) # skips 5 lines
        L2 = len([line for line in t.read()])
        t = track(self.bed, header='track') # skips lines starting with 'track'
        L3 = len([line for line in t.read()])
        t = track(self.bed, header=['track','chr']) # skips lines starting with 'track' or 'chr'
        L4 = len([line for line in t.read()])
        t = track(self.bed, header=['track','chrII']) # skips lines starting with 'track' or 'chrII'
        L5 = len([line for line in t.read()])
        self.assertEqual(L1,L11)
        self.assertEqual(L1,L111)
        self.assertEqual(L1-4,L2)
        self.assertEqual(L1,L3)
        self.assertEqual(L4,0)
        self.assertEqual(L1-11,L5)
        t.close()

    def test_make_header(self):
        t = track(self.bed)
        o = track("temp.bed",fields=t.fields)
        t.open()
        o.make_header(t.filehandle.readline()) # copy the header
        t.filehandle.seek(0)
        o.write(t.read())
        o.close()
        self.assertEqual(list(t.read()),list(o.read()))

    def test_intermediate_header(self):
        # header lines in the middle
        t = track(os.path.join(path,"yeast_genes_tracklines.bed"))
        s = t.read()
        for x in s: pass

    def tearDown(self):
        for test_file in ['temp.bed']:
            if os.path.exists(test_file): os.remove(test_file)

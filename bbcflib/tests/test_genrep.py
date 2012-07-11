# Built-in modules #
import cStringIO, os

# Internal modules #
from bbcflib.genrep import Assembly, GenRep

# Unitesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

# Path to testing files
path = "test_data/genrep/"

###################################################################################

class Test_Assembly(unittest.TestCase):
    def setUp(self):
        self.assembly = Assembly('ce6')
        self.root = self.assembly.genrep.root
        """
        Gene Y54E2A.11 (-1 to each Ensembl start by convention)
        4 transcripts, Y54E2A.11a.1, Y54E2A.11a.2, Y54E2A.11b.1, Y54E2A.11b.2

               5327    5434      5503  5697     5742   6075      6128      6836      6906   8367
        a.1.1   |--------|   b.2.2 |----|   a.1.3 |------|   a.1.4 |---------|   a.1.5 |------|
        b.1.1     |------|                  b.1.3 |----|     b.1.4    |------|   a.2.5 |----|
        b.2.4       |----|                                   b.2.4    |----|

                |========|         |====|         |======|         |=========|         |======|
        2863   =   107        +     194       +     333        +       708       +       1461
        """
    def test_fasta_from_regions(self):
        expected = ({'chrI':['G','C','AAGCCTAAGCCTAAGCCTAA','TTTTTGAAAT']}, 32)
        # GenRep request, list form
        regions = [('chrI',0,1),('chrI',1,2),('chrI',10,30),('chrI',1010,1020)]
        url_seq = self.assembly.fasta_from_regions(regions=regions,out={})
        self.assertEqual(url_seq, expected)
        # Custom fasta file, dict form
        regions = {'chrI':[(0,1),(1,2),(10,30),(1010,1020)]}
        custom_seq = self.assembly.fasta_from_regions(regions=regions,out={},
                            path_to_ref=os.path.join(path,"chrI_ce6_30lines.fa"))
        self.assertEqual(custom_seq, expected)

    def test_get_features_from_gtf(self):
        expected = {'eif-3.B': [[14795327, 14795434, 1, 'chrII'], [14795331, 14795434, 1, 'chrII'],
                                [14795333, 14795434, 1, 'chrII'], [14795503, 14795697, 1, 'chrII'],
                                [14795742, 14795907, 1, 'chrII'], [14795742, 14796075, 1, 'chrII'],
                                [14796128, 14796836, 1, 'chrII'], [14796213, 14796354, 1, 'chrII'],
                                [14796213, 14796836, 1, 'chrII'], [14796906, 14797767, 1, 'chrII'],
                                [14796906, 14798367, 1, 'chrII']]}
        h = {'keys':'gene_name', 'values':'start,end,strand',
             'conditions':'gene_id:Y54E2A.11,type:exon', 'uniq':'1'}
        # Test with local database request
        zc = self.assembly.get_features_from_gtf(h, chr='chrII')
        self.assertItemsEqual(zc,expected)
        # Test with url request via GenRep
        self.assembly.genrep.root = ''
        zc = self.assembly.get_features_from_gtf(h, chr='chrII')
        self.assertItemsEqual(zc['eif-3.B'],expected['eif-3.B'])
        self.assembly.genrep.root = self.root

    @unittest.skip('long')
    def test_get_gene_mapping(self):
        expected = ('eif-3.B',14795327,14798367,2803,1,'chrII')
        # Test with local database request
        map = self.assembly.get_gene_mapping()
        zc = map['Y54E2A.11']
        self.assertEqual(zc,expected)
        # Test with url request via GenRep
        self.assembly.genrep.root = ''
        map = self.assembly.get_gene_mapping()
        zc = map['Y54E2A.11']
        self.assertEqual(zc,expected)
        self.assembly.genrep.root = self.root

    @unittest.skip('long')
    def test_get_transcript_mapping(self):
        expected = ('Y54E2A.11',14795327,14798367,2803,1,'chrII')
        # Test with local database request
        map = self.assembly.get_transcript_mapping()
        zc = map['Y54E2A.11a.1']
        self.assertEqual(zc,expected)
        # Test with url request via GenRep
        self.assembly.genrep.root = ''
        map = self.assembly.get_transcript_mapping()
        zc = map['Y54E2A.11a.1']
        self.assertEqual(zc,expected)
        self.assembly.genrep.root = self.root

    @unittest.skip('long')
    def test_get_exon_mapping(self):
        expected = (['Y54E2A.11a.1'],'Y54E2A.11',14795327,14795434,1,'chrII')
        # Test with local database request
        map = self.assembly.get_exon_mapping()
        zc = map['Y54E2A.11a.1.1']
        self.assertEqual(zc,expected)
        # Test with url request via GenRep
        self.assembly.genrep.root = ''
        map = self.assembly.get_exon_mapping()
        zc = map['Y54E2A.11a.1.1']
        print "zc",zc
        print "exp",expected
        self.assertEqual(zc,expected)
        self.assembly.genrep.root = self.root

    @unittest.skip('long')
    def test_get_exons_in_trans(self):
        expected = ['Y54E2A.11a.1.1','Y54E2A.11b.2.2','Y54E2A.11a.1.3',
                    'Y54E2A.11a.1.4','Y54E2A.11b.1.5'] # Y54E2A.11a.1.5 = Y54E2A.11b.1.5
        # Test with local database request
        map = self.assembly.get_exons_in_trans()
        zc = map['Y54E2A.11a.1']
        self.assertItemsEqual(zc,expected)
        # Test with url request via GenRep
        self.assembly.genrep.root = ''
        map = self.assembly.get_exons_in_trans()
        zc = map['Y54E2A.11a.1']
        self.assertItemsEqual(zc,expected)
        self.assembly.genrep.root = self.root

    @unittest.skip('long')
    def test_trans_in_gene(self):
        expected = ['Y54E2A.11a.1','Y54E2A.11a.2','Y54E2A.11b.1','Y54E2A.11b.2']
        # Test with local database request
        map = self.assembly.get_trans_in_gene()
        zc = map['Y54E2A.11']
        self.assertItemsEqual(zc,expected)
        # Test with url request via GenRep
        self.assembly.genrep.root = ''
        map = self.assembly.get_trans_in_gene()
        zc = map['Y54E2A.11']
        self.assertItemsEqual(zc,expected)
        self.assembly.genrep.root = self.root


class Test_GenRep(unittest.TestCase):
    def setUp(self):
        self.assembly = Assembly('ce6')
        self.assembly.genrep = GenRep('http://bbcftools.vital-it.ch/genrep/','/db/genrep')
        self.assembly.intype = '0'
        self.chromosomes = {(3066, u'NC_003279', 6): {'length': 15072421, 'name': u'chrI'},
                            (3067, u'NC_003280', 7): {'length': 15279323, 'name': u'chrII'},
                            (3068, u'NC_003281', 8): {'length': 13783681, 'name': u'chrIII'},
                            (3069, u'NC_003282', 5): {'length': 17493785, 'name': u'chrIV'},
                            (3070, u'NC_003283', 8): {'length': 20919568, 'name': u'chrV'},
                            (3071, u'NC_003284', 7): {'length': 17718854, 'name': u'chrX'},
                            (2948, u'NC_001328', 1): {'length': 13794,    'name': u'chrM'}}

    def test_attributes(self):
        self.assertItemsEqual(self.assembly.chromosomes,self.chromosomes)
        self.assertEqual(self.assembly.id,14)

    def test_config_correctly_loaded(self):
        self.assertEqual(self.assembly.genrep.url, 'http://bbcftools.vital-it.ch/genrep')
        self.assertEqual(self.assembly.genrep.root, '/db/genrep')

    def test_get_sequence(self):
        expected = ['G','C','AAGCCTAAGCCTAAGCCTAA','TTTTTGAAAT']
        coord_list = [(0,1),(1,2),(10,30),(1010,1020)]
        url_seq = self.assembly.genrep.get_sequence(chr_id=3066, coord_list=coord_list)
        custom_seq = self.assembly.genrep.get_sequence(chr_id=0, coord_list=coord_list,
                                              path_to_ref=os.path.join(path,"chrI_ce6_30lines.fa"))
        self.assertEqual(url_seq,expected)
        self.assertEqual(custom_seq,expected)


#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#



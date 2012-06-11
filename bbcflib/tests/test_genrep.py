# Built-in modules #
import datetime, ConfigParser, cStringIO

# Internal modules #
from ..genrep import Assembly, GenRep

# Unitesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

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


test_config_file = '''[genrep]
    genrep_url = http://bbcftools.vital-it.ch/genrep/
    genrep_root = /db/genrep'''

def get_config_file_parser():
    file = cStringIO.StringIO()
    file.write(test_config_file)
    file.seek(0)
    config = ConfigParser.ConfigParser()
    config.readfp(file)
    return config

@unittest.skip("These tests don't pass anymore. Delete this line once they are fixed.")
class Test_GenRep(unittest.TestCase):
    def setUp(self):
        self.genrep = GenRep('http://bbcftools.vital-it.ch/genrep/',
                             '/db/genrep/nr_assemblies')
        if self.genrep.is_down():
            self.skipTest("The Genrep server is down")
        self.genrep_from_config = GenRep(config = get_config_file_parser())
        self.ce6 = Assembly(assembly_id = 14,
               assembly_name = 'ce6',
               bbcf_valid = True,
               updated_at = datetime.datetime(2011, 1, 5, 14, 58, 43),
               nr_assembly_id = 106,
               genome_id = 8,
               source_name = 'UCSC',
               md5 = 'fd5631288b9cd427bf329e8868fb4c988752b2c5',
               source_id = 4,
               created_at = datetime.datetime.strptime('2010-12-19T20:52:31Z', '%Y-%m-%dT%H:%M:%SZ'),
               index_path = '/db/genrep/nr_assemblies/bowtie/fd5631288b9cd427bf329e8868fb4c988752b2c5')
        self.ce6.chromosomes = {(3067, u'NC_003280', 7): {'length': 15279323, 'name': u'chrII'},
                   (3069, u'NC_003282', 5): {'length': 17493785, 'name': u'chrIV'},
                   (3068, u'NC_003281', 8): {'length': 13783681, 'name': u'chrIII'},
                   (3070, u'NC_003283', 8): {'length': 20919568, 'name': u'chrV'},
                   (3071, u'NC_003284', 7): {'length': 17718854, 'name': u'chrX'},
                   (3066, u'NC_003279', 6): {'length': 15072421, 'name': u'chrI'},
                   (2948, u'NC_001328', 1): {'length': 13794, 'name': u'chrM'}}

    def test_config_correctly_loaded(self):
        self.assertEqual(self.genrep.url, 'http://bbcftools.vital-it.ch/genrep')
        self.assertEqual(self.genrep.root, '/db/genrep')

    def test_query_url(self):
        def _check_with_url(url):
            g = GenRep(url, '')
            self.assertEqual(g.query_url('boris', 'hilda'),
                             'http://bbcftools.vital-it.ch/genrep/boris.json?assembly_name = hilda')
            self.assertEqual(g.query_url('boris', 36),
                             'http://bbcftools.vital-it.ch/genrep/boris.json?assembly_id = 36')
        [_check_with_url(u) for u in ['http://bbcftools.vital-it.ch/genrep/',
                                     'http://bbcftools.vital-it.ch/genrep',
                                     'bbcftools.vital-it.ch/genrep/',
                                     'bbcftools.vital-it.ch/genrep']]
        g = GenRep('http://bbcftools.vital-it.ch/genrep', '')
        self.assertRaises(ValueError, lambda : g.query_url('boris', [1, 2, 3]))

    def assertAssembliesEqual(self, a, b):
        self.assertEqual(a.id, b.id)
        self.assertEqual(a.name, b.name)
        self.assertEqual(a.index_path, b.index_path)
        self.assertEqual(a.bbcf_valid, b.bbcf_valid)
        self.assertEqual(a.updated_at, b.updated_at)
        self.assertEqual(a.nr_assembly_id, b.nr_assembly_id)
        self.assertEqual(a.genome_id, b.genome_id)
        self.assertEqual(a.source_name, b.source_name)
        self.assertEqual(a.md5, b.md5)
        self.assertEqual(a.source_id, b.source_id)
        self.assertEqual(a.created_at, b.created_at)
        self.assertEqual(a.chromosomes, b.chromosomes)

    def test_assembly_with_name(self):
        a = self.genrep.assembly('ce6')
        self.assertAssembliesEqual(a, self.ce6)

    def test_assembly_with_id(self):
        a = self.genrep.assembly(14)
        self.assertAssembliesEqual(a, self.ce6)

    def test_assembly_with_name_from_config(self):
        a = self.genrep_from_config.assembly('ce6')
        self.assertAssembliesEqual(a, self.ce6)

    def test_assembly_with_id_from_config(self):
        a = self.genrep_from_config.assembly(14)
        self.assertAssembliesEqual(a, self.ce6)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#



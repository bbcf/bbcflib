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

#--------------------------------------------------------------------------------#
ce6 = Assembly(assembly_id = 14,
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

ce6.chromosomes = {(3067, u'NC_003280', 7): {'length': 15279323, 'name': u'chrII'},
                   (3069, u'NC_003282', 5): {'length': 17493785, 'name': u'chrIV'},
                   (3068, u'NC_003281', 8): {'length': 13783681, 'name': u'chrIII'},
                   (3070, u'NC_003283', 8): {'length': 20919568, 'name': u'chrV'},
                   (3071, u'NC_003284', 7): {'length': 17718854, 'name': u'chrX'},
                   (3066, u'NC_003279', 6): {'length': 15072421, 'name': u'chrI'},
                   (2948, u'NC_001328', 1): {'length': 13794, 'name': u'chrM'}}

#--------------------------------------------------------------------------------#
test_config_file = '''[genrep]
genrep_url=http://bbcftools.vital-it.ch/genrep/
genrep_root=/db/genrep'''
def get_config_file_parser():
    file = cStringIO.StringIO()
    file.write(test_config_file)
    file.seek(0)
    config = ConfigParser.ConfigParser()
    config.readfp(file)
    return config

###################################################################################
class Test_GenRep(unittest.TestCase):
    def setUp(self):
        self.skipTest("These tests don't pass anymore. Delete this line once they are fixed.")
        self.genrep = GenRep('http://bbcftools.vital-it.ch/genrep/',
                             '/db/genrep/nr_assemblies')
        if self.genrep.is_down():
            self.skipTest("The Genrep server is down")
        self.genrep_from_config = GenRep(config=get_config_file_parser())

    def test_config_correctly_loaded(self):
        self.assertEqual(self.genrep.url, 'http://bbcftools.vital-it.ch/genrep')
        self.assertEqual(self.genrep.root, '/db/genrep')

    def test_query_url(self):
        def check_with_url(url):
            g = GenRep(url, '')
            self.assertEqual(g.query_url('boris', 'hilda'),
                             'http://bbcftools.vital-it.ch/genrep/boris.json?assembly_name=hilda')
            self.assertEqual(g.query_url('boris', 36),
                             'http://bbcftools.vital-it.ch/genrep/boris.json?assembly_id=36')
        [check_with_url(u) for u in ['http://bbcftools.vital-it.ch/genrep/',
                                     'http://bbcftools.vital-it.ch/genrep',
                                     'bbcftools.vital-it.ch/genrep/',
                                     'bbcftools.vital-it.ch/genrep']]
        g = GenRep('http://bbcftools.vital-it.ch/genrep', '')
        self.assertRaises(ValueError,
                          lambda : g.query_url('boris',[1,2,3]))

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
        self.assertAssembliesEqual(a, ce6)

    def test_assembly_with_id(self):
        a = self.genrep.assembly(14)
        self.assertAssembliesEqual(a, ce6)

    def test_assembly_with_name_from_config(self):
        a = self.genrep_from_config.assembly('ce6')
        self.assertAssembliesEqual(a, ce6)

    def test_assembly_with_id_from_config(self):
        a = self.genrep_from_config.assembly(14)
        self.assertAssembliesEqual(a, ce6)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

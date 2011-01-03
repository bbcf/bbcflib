from bbcflib import *
from unittest import TestCase, TestSuite, main
from datetime import datetime

ce6 = Assembly(assembly_id = 14,
               assembly_name = 'ce6',
               bbcf_valid = True,
               updated_at = datetime.strptime('2010-12-20T07:47:55Z', '%Y-%m-%dT%H:%M:%SZ'),
               nr_assembly_id = 103,
               genome_id = 8,
               source_name = 'UCSC',
               md5 = '75d3de3127a40c1aa5fd835ba35984d40f3405a2',
               source_id = 4,
               created_at = datetime.strptime('2010-12-19T20:52:31Z', '%Y-%m-%dT%H:%M:%SZ'),
               index_path = '/scratch/frt/yearly/genrep/nr_assemblies/bowtie/75d3de3127a40c1aa5fd835ba35984d40f3405a2')
ce6.chromosomes = {(3067, u'NC_003280', 7): {'length': 15279323, 'name': u'chrII'},
                   (3069, u'NC_003282', 5): {'length': 17493785, 'name': u'chrIV'}, 
                   (3068, u'NC_003281', 8): {'length': 13783681, 'name': u'chrIII'}, 
                   (3070, u'NC_003283', 8): {'length': 20919568, 'name': u'chrV'},
                   (3071, u'NC_003284', 7): {'length': 17718854, 'name': u'chrX'}, 
                   (3066, u'NC_003279', 6): {'length': 15072421, 'name': u'chrI'}, 
                   (2948, u'NC_001328', 1): {'length': 13794, 'name': u'chrM'}}


class TestGenRep(TestCase):
    def setUp(self):
        self.genrep = GenRep('http://bbcftools.vital-it.ch/genrep/',
                             '/scratch/frt/yearly/genrep/nr_assemblies/bowtie')

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

    def test_assembly_with_name(self):
        a = self.genrep.assembly('ce6')
        self.assertEqual(a.id, ce6.id)
        self.assertEqual(a.name, ce6.name)
        self.assertEqual(a.index_path, ce6.index_path)
        self.assertEqual(a.bbcf_valid, ce6.bbcf_valid)
        self.assertEqual(a.updated_at, ce6.updated_at)
        self.assertEqual(a.nr_assembly_id, ce6.nr_assembly_id)
        self.assertEqual(a.genome_id, ce6.genome_id)
        self.assertEqual(a.source_name, ce6.source_name)
        self.assertEqual(a.md5, ce6.md5)
        self.assertEqual(a.source_id, ce6.source_id)
        self.assertEqual(a.created_at, ce6.created_at)
        self.assertEqual(a.chromosomes, ce6.chromosomes)

    def test_assembly_with_id(self):
        a = self.genrep.assembly(14)
        self.assertEqual(a.id, ce6.id)
        self.assertEqual(a.name, ce6.name)
        self.assertEqual(a.index_path, ce6.index_path)
        self.assertEqual(a.bbcf_valid, ce6.bbcf_valid)
        self.assertEqual(a.updated_at, ce6.updated_at)
        self.assertEqual(a.nr_assembly_id, ce6.nr_assembly_id)
        self.assertEqual(a.genome_id, ce6.genome_id)
        self.assertEqual(a.source_name, ce6.source_name)
        self.assertEqual(a.md5, ce6.md5)
        self.assertEqual(a.source_id, ce6.source_id)
        self.assertEqual(a.created_at, ce6.created_at)
        self.assertEqual(a.chromosomes, ce6.chromosomes)


class TestEmailReport(TestCase):
    pass

class TestDAFLIMS(TestCase):
    pass

class TestFrontend(TestCase):
    pass

def test_all():
    main()

if __name__ == '__main__':
    test_all()

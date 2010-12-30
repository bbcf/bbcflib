from bbcflib import *
from unittest import TestCase, TestSuite, main
from datetime import datetime

ce6 = Assembly(assembly_id = 14,
               assembly_name = 'ce6',
               bbcf_valid = True
               updated_at = datetime.strptime('2010-12-20T07:47:55Z', '%Y-%m-%dT%H:%M:%SZ'),
               nr_assembly_id = 103,
               genome_id = 8,
               source_name = 'UCSC',
               md5 = '75d3de3127a40c1aa5fd835ba35984d40f3405a2',
               source_id = 4,
               created_at = datetime.strptime('2010-12-19T20:52:31Z', '%Y-%m-%dT%H:%M:%SZ'),
               index_path = '/scratch/frt/yearly/genrep/75d3de3127a40c1aa5fd835ba35984d40f3405a2')
# ce6.add_chromosome(id = 2948, 
#                    refseq_locus = 'NC_001328', 
#                    refseq_version = 1,
#                    name = 
#                    length = 13784,


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

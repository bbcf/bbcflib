# Built-in modules #
import os
from functools import wraps

# Internal modules #
from bbcflib.genrep import Assembly, GenRep

# Unitesting module #
try:
    import unittest2 as unittest
    assert unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

# Path to testing files
path = "test_data/genrep/"


class Test_Assembly(unittest.TestCase):
    def setUp(self):
        self.assembly = Assembly('ce6')
        self.root = self.assembly.genrep.root
        self.intype = 0
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

    def with_without_genrep(test):
        """Decorator. Runs *test* with genrep.root successively activated (via /db/)
        and disabled (via URL to GenRep)."""
        @wraps(test) # gives to the wrapper the original function name
        def wrapper(self):
            root = self.assembly.genrep.root
            test(self)
            self.assembly.genrep.root = ''
            test(self)
            self.assembly.genrep.root = root
        return wrapper

    def test_fasta_from_regions(self):
        expected = ({'chrI':['G','C','AAGCCTAAGCCTAAGCCTAA','CTAAGCCTAAGCCTAAGCCT','TTTTTGAAAT']}, 52)
        # GenRep request, list form
        regions = [('chrI',0,1),('chrI',1,2),('chrI',10,30),('chrI',20,40),('chrI',1010,1020)]
        url_seq = self.assembly.fasta_from_regions(regions=regions,out={})
        self.assertEqual(url_seq, expected)
        # Custom fasta file, dict form
        regions = {'chrI':[(0,1),(1,2),(10,30),(20,40),(1010,1020)]}
        custom_seq = self.assembly.fasta_from_regions(regions=regions,out={},
                            path_to_ref=os.path.join(path,"chrI_ce6_30lines.fa"))
        self.assertEqual(custom_seq, expected)

    def test_transcripts_fasta(self):
        import sqlite3,itertools
        # cDNA fasta headers are of the form ">trans_id (...)"
        def _intersect_transcripts(regions):
            transcript_regions = {}
            dbpath = self.assembly.sqlite_path()
            db = sqlite3.connect(dbpath)
            cursor = db.cursor()
            for chrom,regs in regions.iteritems():
                regs = sorted(regs)
                for reg in regs:
                    request = """SELECT DISTINCT transcript_id,exon_id,start,end from '%s'
                                 WHERE end>=%d AND start<%d AND type='exon';""" \
                              % (chrom,reg[0],reg[1])
                    cursor.execute(request)
                    transcripts = dict((key,sorted((g[2],g[3]) for g in group)) \
                                       for key,group in itertools.groupby(cursor, lambda x: x[0]))
                    for tid,exon_coords in transcripts.iteritems():
                        tstart = exon_coords[0][0]
                        for exon in exon_coords:
                            st = max(exon[0],reg[0])
                            en = min(exon[1],reg[1])
                            if en > st:
                                transcript_regions.setdefault(tid, []).append((st-tstart,en-tstart))
            return transcript_regions

        #import json
        #exon_mapping = self.assembly.get_exon_mapping()
        #exons_in_trans = self.assembly.get_exons_in_trans()
        #with open('exons_in_trans_ce6.json','wb') as f:
        #    json.dump(exons_in_trans,f)
        #with open('exon_mapping_ce6.json','wb') as f:
        #    json.dump(exon_mapping,f)
        #with open('exons_in_trans_ce6.json') as f:
        #    exons_in_trans = json.load(f)
        #with open('exon_mapping_ce6.json') as f:
        #    exon_mapping = json.load(f)

        intype = self.assembly.intype
        chromosomes = self.assembly.chromosomes
        self.assembly.intype = 2
        self.assembly.chromosomes = dict((t[3].split('|')[0],{'length':t[2]-t[1], 'name':t[3].split('|')[0]}) \
                                          for t in self.assembly.transcript_track())
        regions = {'chrI':[(4100,10250)]}
        transcript_regions = _intersect_transcripts(regions)
        print "Regions:",transcript_regions
        print "Fasta_path:",self.assembly.fasta_path()
        seqs = self.assembly.fasta_from_regions(regions=transcript_regions, path_to_ref=self.assembly.fasta_path(), out={})
        self.assembly.intype = intype
        self.assembly.chromosomes = chromosomes
        print "Answer:",seqs
        raise

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

    ################
    # Annot tracks #
    ################

    @unittest.skip('slow')
    @with_without_genrep
    def test_gene_track(self):
        expected = ('chrI',4118,10232,'Y74C9A.3|Y74C9A.3',-1)
        track = self.assembly.gene_track()
        self.assertEqual(track.next(),expected)

    @unittest.skip('slow')
    @with_without_genrep
    def test_exon_track(self):
        expected = ('chrI',4118,4358,'Y74C9A.3.1.5|Y74C9A.3|Y74C9A.3',-1,-1)
        track = self.assembly.exon_track()
        self.assertEqual(track.next(),expected)

    @unittest.skip('slow')
    @with_without_genrep
    def test_transcript_track(self):
        expected = ('chrI',4118,10232,'Y74C9A.3.1|Y74C9A.3',-1)
        track = self.assembly.transcript_track()
        self.assertEqual(track.next(),expected)

    ############
    # Mappings #
    ############

    @unittest.skip('slow')
    @with_without_genrep
    def test_get_gene_mapping(self):
        expected = ('eif-3.B',14795327,14798367,2803,1,'chrII')
        map = self.assembly.get_gene_mapping()
        zc = map['Y54E2A.11']
        self.assertEqual(zc,expected)

    @with_without_genrep
    def test_get_gene_mapping_myco(self):
        expected = ('gdh',2777387,2782262,4875,-1,'chr')
        assembly = Assembly('mycoTube_H37RV')
        # Test with local database request
        map = assembly.get_gene_mapping()
        zc = map['EBMYCG00000000868']
        self.assertEqual(zc,expected)
        # Test with url request via GenRep
        self.assembly.genrep.root = ''
        map = assembly.get_gene_mapping()
            # First item
        zc = map['EBMYCG00000000868']
        self.assertEqual(zc,expected)
            # Second item
        expected = ('xseA',1234212,1235460,1248,-1,'chr')
        zc = map['EBMYCG00000001689']
        self.assertEqual(zc,expected)
            # Last item
        expected = ('Rv3085', 3450919, 3451750, 831, 1, 'chr')
        zc = map['EBMYCG00000002832']
        self.assertEqual(zc,expected)
        self.assembly.genrep.root = self.root
            # Check all lengths are set
        for v in map.itervalues():
            self.assertGreaterEqual(v[3],1)

    @unittest.skip('slow')
    @with_without_genrep
    def test_get_transcript_mapping(self):
        expected = ('Y54E2A.11',14795327,14798367,2803,1,'chrII')
        map = self.assembly.get_transcript_mapping()
        zc = map['Y54E2A.11a.1']
        self.assertEqual(zc,expected)

    @unittest.skip('slow')
    @with_without_genrep
    def test_get_exon_mapping(self):
        expected = (['Y54E2A.11a.1'],'Y54E2A.11','eif-3.B',14795327,14795434,1,'chrII')
        map = self.assembly.get_exon_mapping()
        zc = map['Y54E2A.11a.1.1']
        self.assertEqual(zc,expected)

    @unittest.skip('slow')
    @with_without_genrep
    def test_get_exons_in_trans(self):
        expected = ['Y54E2A.11a.1.1','Y54E2A.11b.2.2','Y54E2A.11a.1.3',
                    'Y54E2A.11a.1.4','Y54E2A.11b.1.5'] # Y54E2A.11a.1.5 = Y54E2A.11b.1.5
        map = self.assembly.get_exons_in_trans()
        zc = map['Y54E2A.11a.1']
        self.assertItemsEqual(zc,expected)

    @unittest.skip('slow')
    @with_without_genrep
    def test_trans_in_gene(self):
        expected = ['Y54E2A.11a.1','Y54E2A.11a.2','Y54E2A.11b.1','Y54E2A.11b.2']
        map = self.assembly.get_trans_in_gene()
        zc = map['Y54E2A.11']
        self.assertItemsEqual(zc,expected)


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
        # Sequence existing in genrep
        url_seq = self.assembly.genrep.get_sequence(chr_id=(3066,u'NC_003279',6), coord_list=coord_list)
        # Custom fasta
        custom_seq = self.assembly.genrep.get_sequence(chr_id=None, chr_name='chrI', coord_list=coord_list,
                                              path_to_ref=os.path.join(path,"chrI_ce6_30lines.fa"))
        self.assertEqual(url_seq,expected)
        self.assertEqual(custom_seq,expected)


#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#




# Built-in modules #

# Internal modules #
from bbcflib.rnaseq import *
from bbcflib.frontend import Job
from bein import execution, MiniLIMS

# Unitesting modules #
try: import unittest2 as unittest
except ImportError: import unittest
from numpy.testing import assert_almost_equal

# Other modules #
from numpy import array

# Nosetest flag #
__test__ = True

# Local test:
# run_rnaseq.py -v local -c /archive/epfl/bbcf/jdelafon/test_rnaseq/config/gapkowt.txt -d rnaseq -p genes,transcripts

class Assem(object):
    def __init__(self):
        self.chrmeta = None

def fakejob(assembly):
        job = Job(1,'today','key',assembly.id,'test','email',{})
        job.add_group(id=1,name='group1')
        job.add_run(group_id=1, id=1)
        return job

e1="e1"; e2="e2"; e3="e3"; e4="e4"; e5="e5"; t1="t1"; t2="t2"; t3="t3"; t4="t4"; g1="g1"; gname="gg1"; c="c"

class Test_Expressions1(unittest.TestCase):
    """Two conditions, different lengths, invertible"""
    def setUp(self):
        self.ncond = 2
        self.nreads = array([1e9]*self.ncond)
        self.counts = array([[27,12],[3,3]]) # [[cond1],[cond2]]
        self.rpkms = array([[27/3.,12/6.],[3/3.,3/6.]])
        self.exon_counts = ([e1,e2],self.counts.T)
        self.transcript_mapping = {t1:(g1,gname,0,3,3.,1,c), t2:(g1,gname,0,9,9.,1,c)}
        self.gene_mapping = {g1:(gname,0,9,9.,1,c)}
        self.exon_mapping = {e1:([t1,t2],g1,gname,0,3,1,c), e2:([t1],g1,gname,3,9,1,c)}
        self.trans_in_gene = {g1:[t1,t2]}
        self.exons_in_trans = {t1:[e1], t2:[e1,e2]}
        """
           *Counts Cond1*             *Rpkm Cond1*
        g1 |===.======|               |===.===|
        t1 |---|         21           |---|      7
        t2 |---.------|  6  +  12     |---.---| (2) 2

             27    12                   9   2
           |.e1.||.e2.|               |e1| |e2|
        """
        """
           *Counts Cond2*             *Rpkm Cond2*
        g1 |===.======|               |===.===|
        t1 |---|         1.5          |---|      0.5
        t2 |---.------|  1.5 + 3      |---.---| (0.5) 0.5

             3     3                    1   0.5
           |.e1.||.e2.|               |e1| |e2|
        """
    def test_transcripts_expression(self):
        tcounts, trpkm = transcripts_expression(self.exon_counts, self.exon_mapping, self.transcript_mapping,
                 self.trans_in_gene, self.exons_in_trans, self.ncond, self.nreads)
        assert_almost_equal(trpkm["t1"], array([7., 0.5]))
        assert_almost_equal(trpkm["t2"], array([2., 0.5]))
        assert_almost_equal(tcounts["t1"], array([21., 1.5]))
        assert_almost_equal(tcounts["t2"], array([18., 4.5]))
        self.assertEqual(sum(sum(tcounts.values())), sum(sum(self.counts)))

    def test_genes_expression(self):
        gcounts, grpkm = genes_expression(self.exon_counts, self.gene_mapping,
                self.exon_mapping, self.ncond, self.nreads)
        assert_almost_equal(gcounts["g1"], array([39., 6.]))
        assert_almost_equal(grpkm["g1"],array([39./9,6./9]))

    def test_coherence(self):
        # Test if sum of transcript counts equals gene counts
        tcounts, trpkmm = transcripts_expression(self.exon_counts, self.exon_mapping, self.transcript_mapping,
                     self.trans_in_gene, self.exons_in_trans, self.ncond, self.nreads)
        gcounts, grpkms = genes_expression(self.exon_counts, self.gene_mapping,
                     self.exon_mapping, self.ncond, self.nreads)
        self.assertEqual(sum(gcounts["g1"]), sum([sum(tcounts[t]) for t in self.trans_in_gene["g1"]]))


class Test_Expressions2(unittest.TestCase):
    """One condition, equal lengths, invertible"""
    def setUp(self):
        self.ncond = 1
        self.nreads = array([1e9])
        self.counts = array([[10,15,10]]) # [[cond1]]
        self.rpkms = array([[10/5.,15/5.,10/5.]])
        self.exon_counts = ([e1,e2,e3],self.counts.T)
        self.transcript_mapping = {t1:(g1,gname,0,10,10.,1,c), t2:(g1,gname,5,15,10.,1,c), t3:(g1,gname,0,15,15.,1,c)}
        self.gene_mapping = {g1:(gname,0,9,9.,1,c)}
        self.exon_mapping = {e1:([t1,t2],g1,gname,0,5,1,c), e2:([t1],g1,gname,5,10,1,c),
                             e3:([t2,t3],g1,gname,10,15,1,c)}
        self.trans_in_gene = {g1:[t1,t2,t3]}
        self.exons_in_trans = {t1:[e1,e2], t2:[e2,e3], t3:[e1,e2,e3]}
        """
           *Counts Cond1*
        g1 |===.===.===|
        t1 |-------|      5 + 5
        t2     |-------|  5 +     5
        t3 |-----------|  5 + 5 + 5

            10  15   10
           |e1||e2| |e3|
        """
    def test_transcripts_expression(self):
        tcounts, trpkm = transcripts_expression(self.exon_counts, self.exon_mapping, self.transcript_mapping,
                 self.trans_in_gene, self.exons_in_trans, self.ncond, self.nreads)
        assert_almost_equal(tcounts["t1"], array([10.]))
        assert_almost_equal(tcounts["t2"], array([10.]))
        assert_almost_equal(tcounts["t3"], array([15.]))
        assert_almost_equal(trpkm["t1"], array([1.]))
        assert_almost_equal(trpkm["t2"], array([1.]))
        assert_almost_equal(trpkm["t3"], array([1.]))
        self.assertEqual(sum(sum(tcounts.values())), sum(sum(self.counts)))

    def test_genes_expression(self):
        gcounts, grpkm = genes_expression(self.exon_counts, self.gene_mapping,
                self.exon_mapping, self.ncond, self.nreads)
        assert_almost_equal(gcounts["g1"], array([35.]))
        assert_almost_equal(grpkm["g1"],array([35./9]))


class Test_Expressions3(unittest.TestCase):
    """Underdetermined system"""
    def setUp(self):
        self.ncond = 1
        self.nreads = array([1e9])
        self.counts = array([[10,15,10]]) # [[cond1]]
        self.rpkms = array([[10/5.,15/5.,10/5.]])
        self.exon_counts = ([e1,e2,e3],self.counts.T)
        self.transcript_mapping = {t1:(g1,gname,0,10,10.,1,c), t2:(g1,gname,5,15,10.,1,c), \
                                   t3:(g1,gname,0,15,15.,1,c), t4:(g1,gname,0,15,10.,1,c)}
        self.gene_mapping = {g1:(gname,0,9,9.,1,c)}
        self.exon_mapping = {e1:([t1,t3,t4],g1,gname,0,5,1,c), e2:([t1,t2,t3],g1,gname,5,10,1,c),
                             e3:([t2,t3,t4],g1,gname,10,15,1,c)}
        self.trans_in_gene = {g1:[t1,t2,t3,t4]}
        self.exons_in_trans = {t1:[e1,e2], t2:[e2,e3], t3:[e1,e2,e3], t4:[e1,e3]}
        """
           *Counts Cond1*
        g1 |===.===.===|
        t1 |-------|      7.5 + 7.5
        t2     |-------|        7.5 + 7.5
        t3 |-----------|  0     0     0
        t4 |---|   |---|  2.5 +       2.5

            10  15   10
           |e1||e2| |e3|
        """
    def test_transcripts_expression(self):
        tcounts, trpkm = transcripts_expression(self.exon_counts, self.exon_mapping, self.transcript_mapping,
                 self.trans_in_gene, self.exons_in_trans, self.ncond, self.nreads)
        self.assertEqual(sum(sum(tcounts.values())), sum(sum(self.counts)))


class Test_Expressions4(unittest.TestCase):
    """Overdetermined system"""
    def setUp(self):
        self.ncond = 1
        self.nreads = array([1e9])
        self.counts = array([[10.,10.,10.,10.]]) # [[cond1]]
        self.rpkms = array([[10/5.,10/5.,10/5.,10/5.]])
        self.exon_counts = ([e1,e2,e3,e4],self.counts.T)
        self.transcript_mapping = {t1:(g1,gname,0,10,10.,1,c), t2:(g1,gname,5,20,15.,1,c), t3:(g1,gname,5,20,10.,1,c)}
        self.gene_mapping = {g1:(gname,0,12,12.,1,c)}
        self.exon_mapping = {e1:([t1],g1,gname,0,5,1,c), e2:([t1,t2,t3],g1,gname,5,10,1,c),
                             e3:([t2],g1,gname,10,15,1,c), e4:([t2,t3],g1,gname,15,20,1,c)}
        self.trans_in_gene = {g1:[t1,t2,t3]}
        self.exons_in_trans = {t1:[e1,e2], t2:[e2,e3,e4], t3:[e2,e4]}
        """
           *Counts Cond1*
        g1 |===.===.===.===|
        t1 |-------|          12 + 12
        t2     |-----------|       4  + 4  + 4
        t3     |---|   |---|  0    0    0    0

            10  10  10   10
           |e1||e2||e3| |e4|
        """
    def test_transcripts_expression(self):
        tcounts, trpkm = transcripts_expression(self.exon_counts, self.exon_mapping, self.transcript_mapping,
                 self.trans_in_gene, self.exons_in_trans, self.ncond, self.nreads)
        self.assertGreater(sum(sum(tcounts.values()))/sum(sum(self.counts)), 0.9)


class Test_Expressions5(unittest.TestCase):
    """Even more underdetermined system"""
    def setUp(self):
        self.ncond = 1
        self.nreads = array([1e9])
        self.counts = array([[10.,10.,10.,10.,10.]]) # [[cond1]]
        self.rpkms = array([[10/5.,10/5.,10/5.,10/5.,10/5.]])
        self.exon_counts = ([e1,e2,e3,e4,e5],self.counts.T)
        self.transcript_mapping = {t1:(g1,gname,0,15,10.,1,c), t2:(g1,gname,10,25,15.,1,c), t3:(g1,gname,5,25,10.,1,c)}
        self.gene_mapping = {g1:(gname,0,15,15.,1,c)}
        self.exon_mapping = {e1:([t1],g1,gname,0,5,1,c), e2:([t1,t3],g1,gname,5,10,1,c),
                             e3:([t1,t2],g1,gname,10,15,1,c), e4:([t2,t3],g1,gname,15,20,1,c),
                             e5:([t2,t3],g1,gname,20,25,1,c)}
        self.trans_in_gene = {g1:[t1,t2,t3]}
        self.exons_in_trans = {t1:[e1,e2,e3], t2:[e3,e4,e5], t3:[e2,e4,e5]}
        """
           *Counts Cond1*
        g1 |===.===.===.===.===|
        t1 |-----------|          4.6 + 4.6 + 4.6
        t2         |-----------|              4.6 + 4.6 + 4.6
        t3     |---|   |-------|        3.0 +       3.0 + 3.0

            10  10  10   10  10
           |e1||e2||e3| |e4||e5|
        """
    def test_transcripts_expression(self):
        tcounts, trpkm = transcripts_expression(self.exon_counts, self.exon_mapping, self.transcript_mapping,
                 self.trans_in_gene, self.exons_in_trans, self.ncond, self.nreads)
        self.assertGreater(sum(sum(tcounts.values()))/sum(sum(self.counts)), 0.9)

@unittest.skip("buld exon_mapping first")
class Test_Expressions_Solenne(unittest.TestCase):
    def setUp(self):
        e1="ENSMUSE00000139423"; e2="ENSMUSE00000139425"; t1="ENSMUST00000015723";
        g1="ENSMUSG00000015579"; c="chr17"; gname="Nkx2-5"
        self.ncond = 6
        self.nreads = array([1e9]*self.ncond)
        self.counts = array([[0,6.8],[0,10.6],[0,28.9],
                                   [4.4,160.9],[1.9,87.2],[2.7,58.5]])
        self.rpkms = array([[0,0.00691],[0,0.01077],[0,0.02937],
                                  [0.008133,0.16352],[0.00351,0.08862],[0.00499,0.05945]])
        self.exon_counts = ([e1,e2],self.counts.T)
        self.transcript_mapping = {t1:(g1,gname,26975609,26978510,1525.,1,c)}
        self.gene_mapping = {g1:(gname,26975609,26978510,2901.,1,c)}
        #self.exon_mapping = {e1:([t1],g1,gname,0,5,1,c), e2:([t1,t3],g1,gname,5,10,1,c),
        #self.exon_lengths = {e1:541., e2:984.}
        self.trans_in_gene = {g1:[t1]}
        self.exons_in_trans = {t1:[e1,e2]}
        """
        g1 |===.===.===|  6.8
        t1 |-----------|  6.8
             0      6.8
           |e1-|   |-e2|
        """
    def test_transcripts_expression(self):
        tcounts, trpkm = transcripts_expression(self.exon_counts, self.exon_mapping, self.transcript_mapping,
                 self.trans_in_gene, self.exons_in_trans, self.ncond, self.nreads)
        # Exact solution, special case len(tg)==1
        assert_almost_equal(tcounts["ENSMUST00000015723"],
                            array([6.8,10.6,28.9,165.3,89.1,61.2]))
        assert_almost_equal(trpkm["ENSMUST00000015723"],
                            array([0.00446,0.00695,0.01895,0.10839,0.05843,0.04013]), 5)
        self.assertAlmostEqual(sum(sum(tcounts.values())), sum(sum(self.counts)))
        # Pseudo-inverse solution
        #assert_almost_equal(tcounts["ENSMUST00000015723"], array([7.25,11.31,30.83,175.52,94.7,64.78]))
        #self.assertAlmostEqual(sum(sum(tcounts.values())), 384.39) #pinv solution


class Test_others(unittest.TestCase):
    def setUp(self):
        self.counts = array([[27,12],[3,3]]) # [[cond1],[cond2]]
        self.assembly = Assem()

    @unittest.skip("check - unused yet")
    def test_estimate_size_factors(self):
        res, size_factors = estimate_size_factors(self.counts)
        self.assertIsInstance(self.counts, numpy.ndarray)
        self.assertEqual(type(self.counts), type(res))
        self.assertEqual(self.counts.shape, res.shape)
        assert_almost_equal(size_factors, array([2.5, 0.41666]), decimal=3)
        assert_almost_equal(res, array([[10.8, 4.8],[7.2, 7.2]]))

    @unittest.skip("")
    def test_save_results(self):
        with execution(None) as ex:
            save_results(ex,cols=[[1,2,3],('a','b','c','d'),[],[],[],[],[]],
                     conditions=["num.1","char.1"],
                     group_ids={'num':1,'char':2}, assembly=self.assembly,
                     header=["num","char","e","e","e","e","e"], feature_type="test")

    def test_to_rpkm(self):
        lengths = array([3.,6.])
        nreads = array([10e6])
        # array form
        rpkm = to_rpkm(self.counts, lengths, nreads)
        assert_almost_equal(rpkm, array([[900.,200.],[100.,50.]]))
        # numeric form
        rpkm = to_rpkm(35., 7, 10e6)
        self.assertEqual(rpkm, 500)
        # dict form
        counts = dict(zip(["a","b"],self.counts))
        rpkm = to_rpkm(counts, lengths, nreads)
        assert_almost_equal(rpkm["a"], array([900.,200.]))
        assert_almost_equal(rpkm["b"], array([100.,50.]))

class Test_NNLS(unittest.TestCase):
    def test_lsqnonneg(self):
        C = array([[0.0372, 0.2869],
                         [0.6861, 0.7071],
                         [0.6233, 0.6245],
                         [0.6344, 0.6170]])
        C1 = array([[0.0372, 0.2869, 0.4],
                          [0.6861, 0.7071, 0.3],
                          [0.6233, 0.6245, 0.1],
                          [0.6344, 0.6170, 0.5]])
        C2 = array([[0.0372, 0.2869, 0.4],
                          [0.6861, 0.7071,-0.3],
                          [0.6233, 0.6245,-0.1],
                          [0.6344, 0.6170, 0.5]])
        d = array([0.8587, 0.1781, 0.0747, 0.8405])

        [x, resnorm, residual] = lsqnonneg(C, d)
        dres = abs(resnorm - 0.8315)          # compare with matlab result
        self.assertLess(dres, 0.001)

        [x, resnorm, residual] = lsqnonneg(C1, d)
        dres = abs(resnorm - 0.1477)          # compare with matlab result
        self.assertLess(dres, 0.01)

        [x, resnorm, residual] = lsqnonneg(C2, d)
        dres = abs(resnorm - 0.1027)          # compare with matlab result
        self.assertLess(dres, 0.01)

        k = array([[0.1210, 0.2319, 0.4398, 0.9342, 0.1370],
                         [0.4508, 0.2393, 0.3400, 0.2644, 0.8188],
                         [0.7159, 0.0498, 0.3142, 0.1603, 0.4302],
                         [0.8928, 0.0784, 0.3651, 0.8729, 0.8903],
                         [0.2731, 0.6408, 0.3932, 0.2379, 0.7349],
                         [0.2548, 0.1909, 0.5915, 0.6458, 0.6873],
                         [0.8656, 0.8439, 0.1197, 0.9669, 0.3461],
                         [0.2324, 0.1739, 0.0381, 0.6649, 0.1660],
                         [0.8049, 0.1708, 0.4586, 0.8704, 0.1556],
                         [0.9084, 0.9943, 0.8699, 0.0099, 0.1911]])
        k1 = k-0.5
        l = array([0.4225, 0.8560, 0.4902, 0.8159, 0.4608, 0.4574, 0.4507, 0.4122, 0.9016, 0.0056])

        [x, resnorm, residual] = lsqnonneg(k, l)
        dres = abs(resnorm - 0.3695)          # compare with matlab result
        self.assertLess(dres, 0.01)

        [x, resnorm, residual] = lsqnonneg(k1, l)
        dres = abs(resnorm - 2.8639)          # compare with matlab result
        self.assertLess(dres, 0.01)


class Test_Pileup(unittest.TestCase):
    def setUp(self):
        self.gene_name = "Gapdh"
        self.gene_id = "ENSG00000111640"
        self.bam = "test_data/Gapdh.bam"

    @unittest.skip("deprecated")
    def test_build_pileup(self):
        self.exons = fetch_labels(self.bam)
        self.counts = build_pileup(self.bam, self.exon, self.assembly)
        self.assertEqual([self.counts[e[0].split('|')[0]] for e in self.exons],\
                         [0, 35, 0, 0, 0, 0, 0, 3679, 3707, 0, 0, 0, 149, 3, 0, 0, 55, 0, 161])
        counts = {'ENSMUSE00000745781': 3, 'ENSMUSE00000886744': 0, 'ENSMUSE00000719706': 0, 'ENSMUSE00000511617': 0, 'ENSMUSE00000873350': 0, 'ENSMUSE00000775454': 0, 'ENSMUSE00000709315': 35, 'ENSMUSE00000822668': 149, 'ENSMUSE00000472146': 0, 'ENSMUSE00000487077': 3679, 'ENSMUSE00000751942': 3707, 'ENSMUSE00000740765': 0, 'ENSMUSE00000512722': 0, 'ENSMUSE00000569415': 0, 'ENSMUSE00000879959': 0, 'ENSMUSE00000815727': 55, 'ENSMUSE00000491786': 0, 'ENSMUSE00000721118': 161, 'ENSMUSE00000881196': 0}
        self.assertItemsEqual(self.counts,counts)


class Test_Junctions(unittest.TestCase):
    def setUp(self):
        # Example from Lingner TERRA study on hg19, reads unmapped to the genome
        self.assembly = genrep.Assembly('hg19')
        self.path = "/archive/epfl/bbcf/jdelafon/bin/SOAPsplice-v1.9/bin/soapsplice"
        self.index = os.path.abspath("test_data/rnaseq/index_chr1_150k/chr1_150k.index")
        unmapped = os.path.abspath("test_data/rnaseq/junc_reads_chr1-100k")
        self.bam_files = {1:{1:{'unmapped_fastq':(unmapped+'_R1',unmapped+'_R2')}}}
        self.job = fakejob(self.assembly)
        self.exon_mapping={}
        self.transcript_mapping={}
        self.exons_in_trans={}

    def test_unmapped(self):
        with execution(None) as ex:
            unmapped(ex,self.job,self.bam_files,self.assembly,group_names={1:'group1'}, \
                     exon_mapping=self.exon_mapping,transcript_mapping=self.transcript_mapping, \
                     exons_in_trans=self.exons_in_trans, via='local')

    def test_find_junctions(self):
        #minilims = MiniLIMS('test_junc_lims')
        minilims = None
        with execution(minilims) as ex:
            options = {'-q':1} # Sanger read quality format
            find_junctions(ex,self.job,self.bam_files,self.assembly,self.index,
                           path_to_soapsplice=self.path,soapsplice_options=options)


#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

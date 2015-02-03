
# Built-in modules #

# Internal modules #
from bbcflib.rnaseq import *
from bbcflib import genrep
from bbcflib.frontend import Job
from bein import execution
import shutil

# Unitesting modules #
try:
    import unittest2 as unittest
    assert unittest
except ImportError:
    import unittest


# Nosetest flag #
__test__ = True

# Local test:
# run_rnaseq.py -v local -c /archive/epfl/bbcf/jdelafon/test_rnaseq/config/gapkowt.txt -d rnaseq -p genes,transcripts

testpath = os.path.abspath("test_data/rnaseq/")
gene_counts = os.path.join(testpath,"gene_counts.txt")
gene_counts_clean = os.path.join(testpath,"gene_counts_clean.txt")
genes_diff = os.path.join(testpath,"genes_diff.txt")
genes_diff_clean = os.path.join(testpath,"genes_diff_clean.txt")

def fakejob(assembly):
        job = Job(1,'today','key',assembly.id,'test','email',{})
        job.add_group(id=1,name='group1')
        job.add_run(group_id=1, id=1)
        return job


class Test_Counting(unittest.TestCase):
    def setUp(self):
        self.assembly = genrep.Assembly('mm9')
        self.job = fakejob(self.assembly)
        stranded = False
        self.args = ("local",self.job,self.assembly,["KO.1","KO.2"],sys.stderr,sys.stdout,"genes",False,stranded)

    def test_count_reads(self):
        with execution(None) as ex:
            C = Counter(ex, *self.args)
            bamfiles = [os.path.join(testpath,"gapdhKO.bam")]*2
            gtf = os.path.join(testpath,"mm9_3genes_renamed.gtf")
            count_files = C.count_reads(bamfiles, gtf)
            with open(count_files["transcripts"]) as trans:
                self.assertEqual(len(trans.readlines()), 18)
            with open(count_files["genes"]) as genes:
                self.assertEqual(len(genes.readlines()), 4)
            # Edit the numbers so that DE analysis makes sense
            genes = [x.split('\t') for x in open(count_files["genes"]).readlines()]
            genes[0][1] = 'Count.KO.1'; genes[0][2] = 'Count.WT.1'
            genes[1][1] = '2.0'; genes[1][2] = '3.0'
            genes[3][1] = '800.0'; genes[3][2] = '900.0'
            genes = ['\t'.join(x) for x in genes]
            with open(gene_counts,'wb') as g:
                g.writelines(genes)


class Test_DE(unittest.TestCase):
    def setUp(self):
        self.assembly = genrep.Assembly('mm9')
        self.job = fakejob(self.assembly)
        stranded = False
        self.args = ("local",self.job,self.assembly,["KO.1","KO.2"],sys.stderr,sys.stdout,"genes",False,stranded)

    def test_clean_before_deseq(self):
        with execution(None) as ex:
            DE = DE_Analysis(ex, *self.args)
            clean = DE.clean_before_deseq(gene_counts)
            shutil.copy(clean, gene_counts_clean)
            with open(clean) as c:
                self.assertEqual(len(c.readlines()), 4)

    @unittest.skip('')
    def test_differential_analysis(self):
        with execution(None) as ex:
            DE = DE_Analysis(ex, *self.args)
            diffs = DE.differential_analysis(gene_counts_clean)
            for diff in diffs:
                shutil.copy(diff, genes_diff)
                with open(diff) as d:
                    self.assertEqual(len(d.readlines()), 5)

    def test_clean_deseq_output(self):
        with execution(None) as ex:
            DE = DE_Analysis(ex, *self.args)
            clean = DE.clean_deseq_output(genes_diff)
            shutil.copy(clean, genes_diff_clean)
            with open(clean) as c:
                self.assertEqual(len(c.readlines()), 5)


#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

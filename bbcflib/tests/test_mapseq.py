import socket
import re
from unittest2 import TestCase, main, skipIf, skip

from bbcflib.mapseq import *
from bein import execution
from bein.util import md5sum, remove_lines_matching

try:
    no_pysam = False
except:
    no_pysam = True

def hostname_contains(pattern):
    hostname = socket.gethostbyaddr(socket.gethostname())[0]
    if re.search(pattern, hostname) == None:
        return False
    else:
        return True

if hostname_contains('vital-it.ch'):
    not_vital_it = False
else:
    not_vital_it = True



class TestBowtie(TestCase):
    @skip("no @SQ line in the header (?)")
    def test_parallel_bowtie_local(self):
        with execution(None) as ex:
            bam = parallel_bowtie(ex, '../test_data/selected_transcripts',
                                  '../test_data/reads.raw', n_lines=250,
                                  via='local')
            sam = bam_to_sam(ex, bam)
            new_sam = remove_lines_matching(ex, '@PG', sam)
            new_bam = sam_to_bam(ex, new_sam)
            self.assertEqual(md5sum(ex, new_bam), '2e6bd8ce814949075715b8ffddd1dcd5')

    @skip("install add_nh_flags")
    @skipIf(no_pysam, "Test requires pysam to run.")
    def test_parallel_bowtie_local_with_nh_flags(self):
        with execution(None) as ex:
            bam = parallel_bowtie(ex, '../test_data/selected_transcripts',
                                  '../test_data/reads.raw', n_lines=250,
                                  add_nh_flags=True, via='local')
            sam = bam_to_sam(ex, bam)
            new_sam = remove_lines_matching(ex, '@PG', sam)
            new_bam = sam_to_bam(ex, new_sam)
            self.assertEqual(md5sum(ex, new_bam), '529cd218ec0a35d5d0a23fd7b842ee20')

    @skipIf(not_vital_it, "Not running on VITAL-IT.")
    def test_parallel_bowtie_lsf(self):
        with execution(None) as ex:
            bam = parallel_bowtie(ex, '../test_data/selected_transcripts',
                                  '../test_data/reads.raw', n_lines=250,
                                  via='lsf')
            sam = bam_to_sam(ex, bam)
            new_sam = remove_lines_matching(ex, '@PG', sam)
            new_bam = sam_to_bam(ex, new_sam)
            m = md5sum(ex, new_bam)
        self.assertEqual(m, 'find right md5sum')

    @skip("install add_nh_flags")
    @skipIf(not_vital_it, "Not running on VITAL-IT.")
    def test_parallel_bowtie_lsf_with_nh_flags(self):
        with execution(None) as ex:
            bam = parallel_bowtie(ex, '../test_data/selected_transcripts',
                                  '../test_data/reads.raw', n_lines=250,
                                  add_nh_flags=True, via='lsf')
            sam = bam_to_sam(ex, bam)
            new_sam = remove_lines_matching(ex, '@PG', sam)
            new_bam = sam_to_bam(ex, new_sam)
            m = md5sum(ex, new_bam)
        self.assertEqual(m, '7b7c270a3980492e82591a785d87538f')


class TestAddNhFlag(TestCase):
    @skip("install add_nh_flags")
    @skipIf(no_pysam, "No PySam")
    def test_internal_add_nh_flag(self):
        with execution(None) as ex:
            f = add_nh_flag('../test_data/mapped.sam')
            m = md5sum(ex, f)
        self.assertEqual(m, '50798b19517575533b8ccae5b1369a3e')

    @skip("install add_nh_flags")
    @skipIf(no_pysam, "No PySam")
    def test_external_add_nh_flag(self):
        with execution(None) as ex:
            f = external_add_nh_flag(ex, '../test_data/mapped.sam')
            g = add_nh_flag('../test_data/mapped.sam')
            m = md5sum(ex, f)
            m2 = md5sum(ex, g)
        self.assertEqual(m, m2)

if __name__ == '__main__':
    main()

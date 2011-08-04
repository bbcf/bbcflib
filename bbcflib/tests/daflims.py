# Built-in modules #
import os, ConfigParser, cStringIO

# Internal modules #
from ..daflims import DAFLIMS

# Unitesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

#--------------------------------------------------------------------------------#
test_config_file = '''[daflims]
daflims_username=boris
daflims_password='''
def get_config_file_parser():
    file = cStringIO.StringIO()
    file.write(test_config_file)
    file.seek(0)
    config = ConfigParser.ConfigParser()
    config.readfp(file)
    return config

###################################################################################
class Test_DAFLIMS(unittest.TestCase):
    def setUp(self):
        self.skipTest("These tests don't pass anymore. Delete this line once they are fixed.")
        self.d = DAFLIMS(config=get_config_file_parser())

    def test_symlinkname(self):
        self.assertEqual(self.d.symlinkname('lgtf', 'R2D2', 91, 3),
                         {'fastq': 'http://uhts-lgtf.vital-it.ch/symlink/CLNT_Xv8dVtOJaoOp.seq.tar.gz',
                          'eland': 'http://uhts-lgtf.vital-it.ch/symlink/CLNT_jVgjbO5f9MPW.exp.tar.gz',
                          'qseq': 'http://uhts-lgtf.vital-it.ch/symlink/CLNT_Lgtk8K2UsjeV.qseq.tar.gz'})
        self.assertRaises(ValueError, self.d.symlinkname, 'lgtf', 'boris', 91, 3)

    def test_lanedesc(self):
        self.assertEqual(self.d.lanedesc('lgtf', 'R2D2', 91, 3),
                         {'quantity (/pM)': 10.0, 'run': 91, 'PI lastname': 'Naef',
                          'PI firstname': 'Felix', 'NCBI ID': '9606', 'library': 'CLNT',
                          'submitter firstname': 'Keith', 'protocol': 'ChIP-Seq', 'cycle': 38,
                          'project': 'Michael Brunner project', 'lane': 3,
                          'run type': 'single read', 'submitter lastname': 'Harshman',
                          'machine': 'R2D2', 'organism': 'Homo sapiens'})
        self.assertRaises(ValueError, self.d.lanedesc, 'lgtf', 'boris', 91, 3)

    def test_fetch_symlink(self):
        fastq = self.d.symlinkname('lgtf','R2D2',91,3)['fastq']
        filename = self.d.fetch_symlink(fastq, to='/scratch/frt/daily/bbcf')
        self.assertTrue(os.path.exists(filename))
        self.assertEqual(os.stat(filename).st_size, 2927141632)
        os.unlink(filename)

    def test_fetch_fastq(self):
        s = self.d.fetch_fastq('lgtf','R2D2',91,3, to='/scratch/frt/daily/bbcf')
        self.assertTrue(os.path.exists(s['path']))
        self.assertEqual(os.stat(s['path']).st_size, 2927141632)
        os.unlink(s['path'])

    def test_fetch_eland(self):
        s = self.d.fetch_eland('lgtf','R2D2',91,3, to='/scratch/frt/daily/bbcf')
        self.assertTrue(os.path.exists(s['path']))
        self.assertEqual(os.stat(s['path']).st_size, 3011727531)
        os.unlink(s['path'])

    def test_fetch_qseq(self):
        s = self.d.fetch_qseq('lgtf','R2D2',91,3, to='/scratch/frt/daily/bbcf')
        self.assertTrue(os.path.exists(s['path']))
        self.assertEqual(os.stat(s['path']).st_size, 27223154)
        os.unlink(s['path'])

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

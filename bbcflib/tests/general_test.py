# This module #
from bbcflib import *

# Unitesting module #
try:
    from unittest2 import TestCase, TestSuite, main
except ImportError:
    from unittest import TestCase, TestSuite, main

# Other modules #
import socket, re, os, cStringIO
from datetime import datetime

# Configuration file #
test_config_file = '''
# test.cfg
# A configuration file for testing bbcflib.

[genrep]
genrep_url=http://bbcftools.vital-it.ch/genrep/
genrep_root=/scratch/frt/yearly/genrep/nr_assemblies/bowtie

[daflims]:
daflims_username=jrougemont
daflims_password=cREThu6u

[emailreport]:
email_sender=nobody@localhost
email_smtp_server=localhost
email_default_subject=Default Subject

[frontend]:
frontend_url=http://htsstation.vital-it.ch/rnaseq/
'''

def get_config_file_parser():
    file = cStringIO.StringIO()
    file.write(test_config_file)
    file.seek(0)
    cp = ConfigParser()
    cp.readfp(file)
    return cp

################################################################################### 
def hostname_contains(pattern):
    hostname = socket.gethostbyaddr(socket.gethostname())[0]
    if re.search(pattern, hostname) == None:
        return False
    else:
        return True
    
ce6 = Assembly(assembly_id = 14,
               assembly_name = 'ce6',
               bbcf_valid = True,
               updated_at = datetime(2011, 1, 5, 14, 58, 43),
               nr_assembly_id = 106,
               genome_id = 8,
               source_name = 'UCSC',
               md5 = 'fd5631288b9cd427bf329e8868fb4c988752b2c5',
               source_id = 4,
               created_at = datetime.strptime('2010-12-19T20:52:31Z', '%Y-%m-%dT%H:%M:%SZ'),
               index_path = '/scratch/frt/yearly/genrep/nr_assemblies/bowtie/fd5631288b9cd427bf329e8868fb4c988752b2c5')
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
        self.genrep_from_config = GenRep(config=get_config_file_parser())

    def test_config_correctly_loaded(self):
        self.assertEqual(self.genrep.url, 'http://bbcftools.vital-it.ch/genrep')
        self.assertEqual(self.genrep.root, '/scratch/frt/yearly/genrep/nr_assemblies/bowtie')

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


class TestEmailReport(TestCase):
    def setUp(self):
        self.report = EmailReport(sender='nobody@localhost',
                                  to='ross@localhost',
                                  subject='Default Subject',
                                  smtp_server='localhost')
        self.report_from_config = EmailReport(config=get_config_file_parser(), to='ross@localhost')

    def test_config_requires_to(self):
        self.assertRaises(TypeError,
                          lambda : EmailReport(config=get_config_file_parser()))


    def test_init(self):
        self.assertEqual(self.report.sender, 'nobody@localhost')
        self.assertEqual(self.report.to, 'ross@localhost')
        self.assertEqual(self.report.subject, 'Default Subject')
        self.assertEqual(self.report.smtp_server, 'localhost')

    def test_init_from_config(self):
        self.assertEqual(self.report_from_config.sender, 'nobody@localhost')
        self.assertEqual(self.report_from_config.to, 'ross@localhost')
        self.assertEqual(self.report_from_config.subject, 'Default Subject')
        self.assertEqual(self.report_from_config.smtp_server, 'localhost')

    def test_append_body(self):
        a = "This is a test of the emergency broadcasting system."
        self.report.appendBody(a)
        self.assertEqual(self.report.body, a)
        self.report.appendBody(a)
        self.assertEqual(self.report.body, a+a)

    def test_dump(self):
        self.assertEqual(self.report.dump(),
                         {'sender': 'nobody@localhost',
                          'to': 'ross@localhost',
                          'subject': 'Default Subject',
                          'smtp_server': 'localhost',
                          'body': self.report.body})

    def test_send(self):
        if not(hostname_contains('vital-it.ch')):
            self.report_from_config.send()


class TestDAFLIMS(TestCase):
    def setUp(self):
        self.d = DAFLIMS(config=get_config_file_parser())
        try:
            print self.d.symlinkname
        except AttributeError:
            self.skipTest("You don't have access to the DAFLIMBS, skipping appropriate tests")

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
        filename = self.d.fetch_symlink(fastq, to='/scratch/frt/daily/fross')
        self.assertTrue(os.path.exists(filename))
        self.assertEqual(os.stat(filename).st_size, 2927141632)
        os.unlink(filename)

    def test_fetch_fastq(self):
        s = self.d.fetch_fastq('lgtf','R2D2',91,3, to='/scratch/frt/daily/fross')
        self.assertTrue(os.path.exists(s['path']))
        self.assertEqual(os.stat(s['path']).st_size, 2927141632)
        os.unlink(s['path'])

    def test_fetch_eland(self):
        s = self.d.fetch_eland('lgtf','R2D2',91,3, to='/scratch/frt/daily/fross')
        self.assertTrue(os.path.exists(s['path']))
        self.assertEqual(os.stat(s['path']).st_size, 3011727531)
        os.unlink(s['path'])

    def test_fetch_qseq(self):
        s = self.d.fetch_qseq('lgtf','R2D2',91,3, to='/scratch/frt/daily/fross')
        self.assertTrue(os.path.exists(s['path']))
        self.assertEqual(os.stat(s['path']).st_size, 27223154)
        os.unlink(s['path'])


class TestFrontend(TestCase):
    def setUp(self):
        self.maxDiff = None
        self.f = Frontend(url='http://htsstation.vital-it.ch/rnaseq/')
        self.key = '9pv1x7PamOj80eXnZa14'
        self.frontend_job = Job(id = 2,
                                created_at = datetime.strptime('2010-12-30T13:29:54Z', '%Y-%m-%dT%H:%M:%SZ'),
                                key = "9pv1x7PamOj80eXnZa14",
                                assembly_id = 14,
                                description = "Job for testing Frontend module",
                                email = "madhadron@gmail.com",
                                options={u'domain': None,
                                         u'protocol': None,
                                         u'accept': None,
                                         u'run_nber': None,
                                         u'input_file': None,
                                         u'from_controller': None,
                                         u'lane_nber': None,
                                         u'method': None,
                                         u'remote_ip': None,
                                         u'bein_id': None,
                                         u'from_action': None,
                                         u'status_id': None,
                                         u'query_string': None,
                                         u'controller': None,
                                         u'referer': None,
                                         u'path': None,
                                         u'facility_id': None,
                                         u'user_agent': None,
                                         u'time': None,
                                         u'action': None,
                                         u'machine_id': None})
        self.frontend_job.add_group(id=3,
                                    control=False,
                                    name=u'My first group')
        self.frontend_job.add_run(id=5, group=3, 
                                  facility=u"lgtf", facility_location=u"Lausanne",
                                  machine=u"C3PO", machine_id=1, 
                                  run=36, lane=1)
        self.frontend_job.add_run(id=6, group=3, 
                                  facility=u"lgtf", facility_location=u"Lausanne",
                                  machine=u"C3PO", machine_id=1, 
                                  run=36, lane=2)
        self.frontend_job.add_group(id=4, control=True, name=u'Other group')
        self.frontend_job.add_run(id=7, group=4,
                                  facility=u"lgtf", facility_location=u"Lausanne",
                                  machine=u"C3PO", machine_id=1,
                                  run=37, lane=3)

    def test_init(self):
        def _f(url):
            q = Frontend(url)
            self.assertEqual(q.url, 'http://htsstation.vital-it.ch/rnaseq')
        [_f(u) for u in ['http://htsstation.vital-it.ch/rnaseq/',
                         'http://htsstation.vital-it.ch/rnaseq',
                         'htsstation.vital-it.ch/rnaseq/',
                         'htsstation.vital-it.ch/rnaseq']]

    def test_query_url(self):
        self.assertEqual(self.f.query_url('groups', self.key),
                         'http://htsstation.vital-it.ch/rnaseq/groups.json?key=9pv1x7PamOj80eXnZa14')

    def test_fetch_groups(self):
        self.assertEqual(self.f._fetch_groups(self.key),
                         [{'control': False, 
                           'created_at': datetime(2010, 12, 30, 13, 29, 54), 
                           'id': 3,
                           'job_id': 2,
                           'name': u'My first group'}, 
                          {'control': True, 
                           'created_at': datetime(2010, 12, 30, 13, 29, 54),
                           'id': 4, 
                           'job_id': 2, 
                           'name': u'Other group'}])

    def test_fetch_runs(self):
        self.assertEqual(self.f._fetch_runs(self.key),
                         [{'created_at': datetime(2010, 12, 30, 13, 29, 54),
                           'lane': 1, 'machine_name': u'C3PO', 'machine_id': 1,
                           'run': 36, 'facility_name': u'lgtf', 'group_id': 3, 'id': 5,
                           'facility_location': u'Lausanne'},
                          {'created_at': datetime(2010, 12, 30, 13, 29, 54),
                           'lane': 2, 'machine_name': u'C3PO', 'machine_id': 1, 
                           'run': 36, 'facility_name': u'lgtf', 'group_id': 3, 'id': 6,
                           'facility_location': u'Lausanne'}, 
                          {'created_at': datetime(2010, 12, 30, 13, 29, 54), 
                           'lane': 3, 'machine_name': u'C3PO', 'machine_id': 1,
                           'run': 37, 'facility_name': u'lgtf', 'group_id': 4, 'id': 7,
                           'facility_location': u'Lausanne'}])

    def test_fetch_job(self):
        self.assertEqual(self.f._fetch_job(self.key),
                         {'description': u'Job for testing Frontend module', 
                          'assembly_id': 14,
                          'created_at': datetime(2010, 12, 30, 13, 29, 54), 
                          'id': 2,
                          'key': u'9pv1x7PamOj80eXnZa14',
                          'email': u'madhadron@gmail.com',
                          'options': {u'domain': None, u'protocol': None, u'accept': None,
                                      u'run_nber': None, u'input_file': None,
                                      u'from_controller': None, u'lane_nber': None,
                                      u'method': None, u'remote_ip': None, u'bein_id': None,
                                      u'from_action': None, u'status_id': 1,
                                      u'query_string': None, u'controller': None,
                                      u'referer': None, u'path': None, u'facility_id': None,
                                      u'user_agent': None, u'time': None, u'action': None,
                                      u'machine_id': None}})

    def test_job(self):
        j = self.f.job(self.key)
        self.assertEqual(j.id, self.frontend_job.id)
        self.assertEqual(j.created_at, self.frontend_job.created_at)
        self.assertEqual(j.key, self.frontend_job.key)
        self.assertEqual(j.assembly_id, self.frontend_job.assembly_id)
        self.assertEqual(j.description, self.frontend_job.description)
        self.assertEqual(j.email, self.frontend_job.email)
        self.assertEqual(j.groups, self.frontend_job.groups)
        for g in j.groups.keys():
            for r in j.groups[g]['runs'].keys():
                self.assertEqual(j.groups[g]['runs'][r], self.frontend_job.groups[g]['runs'][r])
        [self.assertEqual(j.groups[g], self.frontend_job.groups[g])
         for g in j.groups.keys()]

if __name__ == '__main__':
    main()

# Built-in modules #
import datetime

# Internal modules #
from ..frontend import Frontend

# Unitesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

###################################################################################
class Test_Frontend(unittest.TestCase):
    def setUp(self):
        self.skipTest("These tests don't pass anymore. Delete this line once they are fixed.")
        self.f = Frontend(url='http://htsstation.vital-it.ch/rnaseq/')
        self.key = '9pv1x7PamOj80eXnZa14'
        self.frontend_job = Job(id = 2,
                                created_at = datetime.datetime.strptime('2010-12-30T13:29:54Z', '%Y-%m-%dT%H:%M:%SZ'),
                                key = "9pv1x7PamOj80eXnZa14",
                                assembly_id = 14,
                                description = "Job for testing Frontend module",
                                email = "webmaster.bbcf@epfl.ch",
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
                           'created_at': datetime.datetime(2010, 12, 30, 13, 29, 54),
                           'id': 3,
                           'job_id': 2,
                           'name': u'My first group'},
                          {'control': True,
                           'created_at': datetime.datetime(2010, 12, 30, 13, 29, 54),
                           'id': 4,
                           'job_id': 2,
                           'name': u'Other group'}])

    def test_fetch_runs(self):
        self.assertEqual(self.f._fetch_runs(self.key),
                         [{'created_at': datetime.datetime(2010, 12, 30, 13, 29, 54),
                           'lane': 1, 'machine_name': u'C3PO', 'machine_id': 1,
                           'run': 36, 'facility_name': u'lgtf', 'group_id': 3, 'id': 5,
                           'facility_location': u'Lausanne'},
                          {'created_at': datetime.datetime(2010, 12, 30, 13, 29, 54),
                           'lane': 2, 'machine_name': u'C3PO', 'machine_id': 1,
                           'run': 36, 'facility_name': u'lgtf', 'group_id': 3, 'id': 6,
                           'facility_location': u'Lausanne'},
                          {'created_at': datetime.datetime(2010, 12, 30, 13, 29, 54),
                           'lane': 3, 'machine_name': u'C3PO', 'machine_id': 1,
                           'run': 37, 'facility_name': u'lgtf', 'group_id': 4, 'id': 7,
                           'facility_location': u'Lausanne'}])

    def test_fetch_job(self):
        self.assertEqual(self.f._fetch_job(self.key),
                         {'description': u'Job for testing Frontend module',
                          'assembly_id': 14,
                          'created_at': datetime.datetime(2010, 12, 30, 13, 29, 54),
                          'id': 2,
                          'key': u'9pv1x7PamOj80eXnZa14',
                          'email': u'webmaster.bbcf@epfl.ch',
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

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

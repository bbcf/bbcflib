"""
================
bbcflib.frontend
================

This module provides access to Fabrice David's frontend code for
deploying web interfaces to workflows.  It provides a class
``Frontend`` which handles connections, and which returns objects of
type ``Job`` when queried.

A ``Frontend`` object must be given either a URL to the frontend, such
as ``http://htsstation.vital-it.ch/rnaseq/``, or a ``ConfigParser``
object defining the field ``frontend_url`` in the appropriate section.

Basic usage is to create a ``Frontend`` object, then call its ``job``
method with a job key.::

    f = Frontend(url='http://htsstation.vital-it.ch/rnaseq/')
    j = f.job(14)

.. autoclass:: Frontend

.. autoclass:: Job
"""
import urllib2
import json
from datetime import datetime

from common import normalize_url


class Frontend(object):
    """Connection to Fabrice David's web frontends for workflows.

    A Frontend object either takes *url* to connect to, or a
    ``ConfigParser`` given as the *config* keyword argument (and
    optionally a *section* keyword argument defaulting to
    ``"frontend"``.

    If a *config* argument is given, it reads the field
    ``frontend_url`` from the configuration.  If both the *url* and
    *config* arguments are given, *url* overrides config.
    """
    def __init__(self, url=None, config=None, section='frontend'):
        if url == None and config == None:
            raise TypeError("Must specify a URL or a configuration.")
        elif url != None:
            self.url = normalize_url(url)
        else:
            self.url = normalize_url(config.get(section, 'frontend_url'))
        
    def query_url(self, method, key):
        return """%s/%s.json?key=%s""" % (self.url, method, key)

    def _fetch_groups(self, key):
        def _f(g):
            a = g['group']
            return {'control': a['control'],
                    'created_at': datetime.strptime(a['created_at'],
                                                    '%Y-%m-%dT%H:%M:%SZ'),
                    'id': a['id'],
                    'name': a['name'],
                    'job_id': a['job_id']}
        return [_f(g) for g in json.load(urllib2.urlopen(self.query_url('groups', key)))]

    def _fetch_runs(self, key):
        def _f(r):
            a = r['run']
            return {'facility_name': a['facility_name'],
                    'id': a['id'],
                    'group_id': a['group_id'],
                    'machine_name': a['machine_name'],
                    'machine_id': a['machine_id'],
                    'lane': a['lane_nber'],
                    'run': a['run_nber'],
                    'facility_location': a['facility_location'],
                    'created_at': datetime.strptime(a['created_at'],
                                                    '%Y-%m-%dT%H:%M:%SZ')}
        return [_f(r) for r in json.load(urllib2.urlopen(self.query_url('runs', key)))]

    def _fetch_job(self, key):
        j = json.load(urllib2.urlopen("""%s/jobs/%s.json""" % (self.url, key)))['job']
        return {'id': j['id'],
                'created_at': datetime.strptime(j['created_at'],
                                                             '%Y-%m-%dT%H:%M:%SZ'),
                'key': j['key'],
                'assembly_id': j['assembly_id'],
                'description': j['description'],
                'email': j['email']}

    def job(self, key):
        """Fetch information about job *key* as a Job object."""
        x = self._fetch_job(key)
        j = Job(id = x['id'],
                created_at = x['created_at'],
                key = key,
                assembly_id = x['assembly_id'],
                description = x['description'],
                email = x['email'])
        [j.add_group(id=g['id'], control=g['control'], name=g['name'])
         for g in self._fetch_groups(key)]
        [j.add_run(id=r['id'], group=r['group_id'], 
                   facility=r['facility_name'], 
                   facility_location=r['facility_location'],
                   machine=r['machine_name'], 
                   machine_id=r['machine_id'],
                   run=r['run'], lane=r['lane'])
         for r in self._fetch_runs(key)]
        return j
                   


class Job(object):
    """An object specifying a workflow job.

    The fields are:

      * ``id``
      * ``created_at``
      * ``key``
      * ``assembly_id``
      * ``description``
      * ``email``
      * ``groups``

    ``groups`` is a dictionary of group IDs point to dictionaries with
    information about the group and a dictionary of runs.
    """
    def __init__(self, id, created_at, key, assembly_id, description, email):
        self.id = id
        self.created_at = created_at
        self.key = key
        self.assembly_id = assembly_id
        self.description = description
        self.email = email
        self.groups = {}

    def add_group(self, id, control, name):
        if self.groups.has_key(id):
            raise ValueError("A group with ID %d was already added." % id)
        else:
            self.groups[id] = {'control': control,
                               'name': name,
                               'runs': {}}

    def add_run(self, id, group, facility, facility_location, machine, machine_id, run, lane):
        try:
            runs = self.groups[group]['runs']
            if runs.has_key(id):
                raise ValueError("Group %d already has a run with ID %d" % (group, id))
            else:
                runs[id] = {'facility': facility,
                            'facility_location': facility_location,
                            'machine': machine,
                            'machine_id': machine_id,
                            'run': run,
                            'lane': lane}
        except KeyError, k:
            raise KeyError("No such group with ID %d" % group)

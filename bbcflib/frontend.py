"""
========================
Module: bbcflib.frontend
========================

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

# Built-in modules #
import json, urllib2
from datetime import datetime

# Internal modules #
from .common import normalize_url

################################################################################
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
                    'name': str(a['name']),
                    'job_id': a['job_id']}
        return [_f(g) for g in json.load(urllib2.urlopen(self.query_url('groups', key)))]

    def _fetch_runs(self, key):
        def _f(r):
            a = r['run']
            if not(a.get('url') == None):
                a['url'] = str(a['url'])
            else:
                a['url'] = None
            if not(a.get('key') == None):
                a['key'] = str(a['key'])
            else:
                a['key'] = None
            return {'facility_name': str(a['facility_name']),
                    'id': a['id'],
                    'group_id': a['group_id'],
                    'machine_name': str(a['machine_name']),
                    'machine_id': a['machine_id'],
                    'lane': a['lane_nber'],
                    'run': a['run_nber'],
                    'facility_location': str(a['facility_location']),
                    'url': a['url'],
                    'key': a['key'],
                    'created_at': datetime.strptime(a['created_at'],
                                                    '%Y-%m-%dT%H:%M:%SZ')}
        return [_f(r) for r in json.load(urllib2.urlopen(self.query_url('runs', key)))]

    def _fetch_job(self, key):
        j = json.load(urllib2.urlopen("""%s/jobs/%s.json""" % (self.url, key)))['job']
        ret_val = {'id': j.pop('id'),
                   'created_at': datetime.strptime(j.pop('created_at'),
                                                   '%Y-%m-%dT%H:%M:%SZ'),
                   'key': j.pop('key'),
                   'assembly_id': j.pop('assembly_id'),
                   'description': j.pop('description'),
                   'email': j.pop('email')}
        ret_val.update({'options': j})
        return ret_val

    def job(self, key):
        """Fetch information about job *key* as a Job object."""
        x = self._fetch_job(key)
        j = Job(id = x['id'],
                created_at = x['created_at'],
                key = key,
                assembly_id = x['assembly_id'],
                description = str(x['description']),
                email = str(x['email']),
                options = x['options'])
        [j.add_group(id=g['id'], control=g['control'], name=g['name'])
         for g in self._fetch_groups(key)]
        [j.add_run(id=r['id'], group=r['group_id'],
                   facility=r['facility_name'],
                   facility_location=r['facility_location'],
                   machine=r['machine_name'],
                   machine_id=r['machine_id'],
                   run=r['run'], lane=r['lane'],
                   url=r['url'],
                   key=r['key'])
         for r in self._fetch_runs(key)]
        return j

################################################################################
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
      * ``options``

    ``groups`` is a dictionary of group IDs point to dictionaries with
    information about the group and a dictionary of runs.
    """
    def __init__(self, id, created_at, key, assembly_id, description, email, options):
        self.id = id
        self.created_at = created_at
        self.key = key
        self.assembly_id = assembly_id
        self.description = description
        self.email = email
        self.groups = {}
        self.options = options

    def add_group(self, id, control, name):
        if self.groups.has_key(id):
            raise ValueError("A group with ID %d was already added." % id)
        else:
            self.groups[id] = {'control': control,
                               'name': name,
                               'runs': {}}

    def add_run(self, id, group, facility, facility_location, machine, machine_id, run, lane, url, key):
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
                            'lane': lane,
                            'url': url,
                            'key': key}
        except KeyError, k:
            raise KeyError("No such group with ID %d" % group)


def parseConfig( file ):
    """Constructs a Job object from parsing a text config file with ConfigObj.
    """
    from configobj import ConfigObj
    import time

    config = ConfigObj( file )
    if not('Job' in config and 'Groups' in config and 'Runs' in config):
        raise ValueError("Need 'Job', 'Groups' and 'Runs' sections in the configuration, only had: "+", ".join(config.keys()))
    job = Job(
        id = int(config['Job'].get('id') or 0),
        created_at = int(time.time()),
        key = config['Job'].get('key') or 'localconfig',
        assembly_id = config['Job']['assembly_id'],
        description = str(config['Job'].get('description')),
        email = str(config['Job'].get('email')),
        options = config.get('Options') or {})
    for gid, group in config['Groups'].iteritems():
        if not('name' in group):
            raise ValueError("Each entry in 'Groups' must have a 'name'")
        job.add_group(id=int(gid),
                      control=group.get('control').lower() in ['1','true','t'],
                      name=str(group['name']))

    for rid, run in config['Runs'].iteritems():
        if not('group_id' in run):
            raise ValueError("Each entry in 'Runs' must have a 'group_id'")
        job.add_run(id=int(rid),
                    group=int(run['group_id']),
                    facility=str(run.get('facility_name')),
                    facility_location=str(run.get('facility_location')),
                    machine=str(run.get('machine_name')),
                    machine_id=int((run.get('machine_id') is None) and "0" or run.get('machine_id')),
                    run =int((run.get('run') is None) and "0" or run.get('run')),
                    lane=int((run.get('lane') is None) and "0" or run.get('lane')),
                    url=run.get('url'),
                    key=run.get('key'))
    globals = config.get('Global variables') or {}
    return (job,globals)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

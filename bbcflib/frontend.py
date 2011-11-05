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

    def _fix_dict(self,a):
        for k,v in a.iteritems():
            if isinstance(v,unicode): a[k] = str(v)
            if isinstance(k,unicode): a[str(k)] = a.pop(k)
        if isinstance(a.get('created_at'),str):
            a['created_at'] = datetime.strptime(a['created_at'], '%Y-%m-%dT%H:%M:%SZ')
        return a

    def _fetch_groups(self, key):
        return [self._fix_dict(g['group']) for g in json.load(urllib2.urlopen(self.query_url('groups', key)))]

    def _fetch_runs(self, key):
        def _f(a):
            a['lane'] = a.pop('lane_nber')
            a['run'] = a.pop('run_nber')
            a['facility'] = a.pop('facility_name')
            a['machine']=a.pop('machine_name')
            return self._fix_dict(a)
        return [_f(r['run']) for r in json.load(urllib2.urlopen(self.query_url('runs', key)))]

    def _fetch_job(self, key):
        j = json.load(urllib2.urlopen("""%s/jobs/%s.json""" % (self.url, key)))['job']
        j = self._fix_dict(j)
        ret_val = {'id': j.pop('id'),
                   'created_at': j.pop('created_at'),
                   'key': j.pop('key'),
                   'assembly_id': j.pop('assembly_id'),
                   'description': j.pop('description'),
                   'email': j.pop('email')}
        ret_val.update({'options': j})
        return ret_val

    def job(self, key):
        """Fetch information about job *key* as a Job object."""
        x = self._fetch_job(key)
        x = self._fix_dict(x)
        j = Job(**x)
        [j.add_group(id=g.pop('id'), name=g.pop('name'), group=g)
         for g in self._fetch_groups(key)]
        [j.add_run(**r) for r in self._fetch_runs(key)]
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

    def add_group(self, id, name, group = None):
        if self.groups.has_key(id):
            if group is None:
                raise ValueError("A group with ID %d was already added." % id)
        else:
            self.groups[id] = {'name': name, 'runs': {}}
        if "library_file_type_id" in group:
            group.update({"library_file_url": a.get("library_file_url" or ""),
                          "library_id": a.get("library_id") or 0,
                          "library_param_file": a.get("library_param_file") or ""})
        self.groups[id].update(group)

    def add_run(self, **kwargs):
        group = kwargs.pop('group_id')
        try:
            id = kwargs.pop('id')
            runs = self.groups[group]['runs']
            if runs.has_key(id):
                raise ValueError("Group %d already has a run with ID %d" % (group, id))
            else:
                runs[id] = kwargs
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
                      name=str(group['name']),
                      group={'control': (group.get('control').lower() in ['1','true','t'])})

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
                    sequencing_library=(run.get('sequencing_library') is None) and None or str(run.get('sequencing_library')),
                    url=run.get('url'),
                    key=run.get('key'))
    globals = config.get('Global variables') or {}
    return (job,globals)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

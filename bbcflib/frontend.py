"""
========================
Module: bbcflib.frontend
========================

This module provides access to Fabrice David's frontend code for
deploying web interfaces to workflows.  It provides a class
``Frontend`` which handles connections, and which returns objects of
type ``Job`` when queried.

A ``Frontend`` object must be given either a URL to the frontend, such
as ``http://htsstation.epfl.ch/rnaseq/``, or a ``ConfigParser``
object defining the field ``frontend_url`` in the appropriate section.

Basic usage is to create a ``Frontend`` object, then call its ``job``
method with a job key.::

    f = Frontend(url='http://htsstation.epfl.ch/rnaseq/')
    j = f.job(14)

.. autoclass:: Frontend

.. autoclass:: Job
"""

# Built-in modules #
import json, urllib2
from datetime import datetime

# Internal modules #
from bbcflib.common import normalize_url

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
#        if isinstance(a.get('created_at'),str):
#            a['created_at'] = datetime.strptime(a['created_at'], '%Y-%m-%dT%H:%M:%SZ')
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

      * ``id`` (int)
      * ``created_at`` (str)
      * ``key`` (str)
      * ``assembly_id`` (str)
      * ``description`` (str)
      * ``email`` (str)
      * ``groups`` (dict)
      * ``options`` (dict)

    ``groups`` is a dictionary of the form ``{'id':group_id, 'name':group_name, 'runs':{runs_info}}``.
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
        self.files = {}

    def add_group(self, id, name, group=None):
        """Add info on the group to `self.groups`. If a dictionary *group* is given,
        `self.groups` is updated with the info it contains."""
        if self.groups.has_key(id):
            if group is None:
                raise ValueError("A group with ID "+str(id)+" was already added.")
        else:
            self.groups[id] = {'name': name or str(id), 'runs': {}}
        if group:
            if "library_file_type_id" in group:
                group.setdefault("library_file_url","")
                group.setdefault("library_id",0)
                group.setdefault("library_param_file","")
            self.groups[id].update(group)

    def add_run(self, **kwargs):
        """Add info on the run to a specified group info in `self.groups`.
        Mandatory keyword args are: *group_id*, *id* (run ID). All others will be added
        as complementary info to the dictionary describing the run."""
        group = kwargs.pop('group_id')
        try:
            id = kwargs.pop('id')
            runs = self.groups[group]['runs']
            if runs.has_key(id):
                raise ValueError("Group "+str(group)+" already has a run with ID "+str(id))
            else:
                runs[id] = kwargs
        except KeyError, k:
            raise KeyError("No such group with ID "+str(group))


def parseConfig( file, job=None, gl=None ):
    """Constructs or updates a Job object from parsing a text config file with ConfigObj.
    """
    from configobj import ConfigObj

    config = ConfigObj( file, unrepr=True )

    if not(job or ('Job' in config and 'Groups' in config and 'Runs' in config)):
        raise ValueError("Need 'Job', 'Groups' and 'Runs' sections in the configuration, only had: "+", ".join(config.keys()))
    id = 0
    created_at = datetime.now()
    key = 'userconfig'
    assembly_id = None
    description = None
    email = None
    options = {}
    if isinstance(job,Job):
        id = job.id
        created_at = job.created_at
        key = job.key
        assembly_id = job.assembly_id
        description = job.description
        email = job.email
        options = job.options
    if 'Job' in config:
        id = int(config['Job'].get('id',id))
        key = config['Job'].get('key',key)
        assembly_id = config['Job'].get('assembly_id',assembly_id)
        description = str(config['Job'].get('description',description))
        email = str(config['Job'].get('email',email))
    if isinstance(config.get('Options'),dict):
        for k,v in config['Options'].iteritems():
            if k in options:
                if not v: continue
                if isinstance(v,dict):
                    for k2, v2 in v.iteritems():
                        if v2: options[k][k2] = v2
                    continue
            options[k] = v
    newjob = Job( id, created_at, key, assembly_id, description, email, options )
    if isinstance(job,Job):
        newjob.groups = job.groups
        newjob.files = job.files
    for gid, group in config.get('Groups',{}).iteritems():
        if not('name' in group):
            raise ValueError("Each entry in 'Groups' must have a 'name'")
        if isinstance(group.get('control'),str):
            group['control'] = group['control'].lower() in ['1','true','t']
        else:
            group['control'] = False
        newjob.add_group(id=int(gid),name=group.pop('name'),group=group)
        newjob.files.setdefault(int(gid),{})
    for rid, run in config.get('Runs',{}).iteritems():
        if not('group_id' in run):
            raise ValueError("Each entry in 'Runs' must have a 'group_id'")
        run['group_id'] = int(run['group_id'])
        newjob.add_run(id=int(rid),**run)
    for rid, run in config.get('Files',{}).iteritems():
        gid = int(run['group_id'])
        run.pop('group_id')
        newjob.files[gid].setdefault(int(rid),{}).update(run)
    newgl = config.get('Global variables',{})
    if gl is None: gl = {}
    for k,v in gl.iteritems():
        if not(k in newgl): newgl[k] = v
    return (newjob,newgl)


#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

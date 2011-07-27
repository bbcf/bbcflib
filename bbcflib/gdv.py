"""
===================
Module: bbcflib.gdv
===================

Python API for the GDV genome viewer.
"""

# Built-in modules #
import json, urllib, urllib2

# Internal modules #
from .common import normalize_url

################################################################################
# GDV requests #

def create_gdv_project( gdv_key, gdv_email,
                        name, nr_assembly_id,
                        gdv_url="http://svitsrv25.epfl.ch/gdv", public=False ):
    '''
    Create a new project on GDV interface
    :param gdv_email: your login in TEQUILA
    :param gdv_key: your user key (get it from your GDV account)
    :param nr_assembly_id: the nrAssembly identifier of the species in Genrep
    :param name: name of the project
    :param public: 'true' to make the project public -optionnal-
    :rtype: a json : {'project_id':<the id>,'public_url':<the public url>} or {'project_id':<the id>} if you didn't make the
    project public
    '''
    request = { "id": "gdv_post",
                "mail": gdv_email,
                "key": gdv_key,
                "command": "new_project",
                "name": str(name),
                "seq_id": str(nr_assembly_id),
                "public": str(public).lower() }
    gdv_url = normalize_url(gdv_url)+"/post"
    return json.load(urllib2.urlopen( gdv_url, urllib.urlencode(request)))

def get_project_id(json):
    return json['project_id']

def get_public_url(json):
    return json['public_url']

def add_gdv_track( gdv_key, gdv_email,
                   project_id,
                   url,
                   name=None,
                   gdv_url="http://svitsrv25.epfl.ch/gdv"):
    '''
    Add a new track on a project on GDV
    :param gdv_email: your login in TEQUILA
    :param gdv_key: your user key (get it from your GDV account)
    :param name: name of the track -optionnal- (will take the file name by default)
    :param project_id: the project id to add the track
    :param url : the URL where to fetch the file
    :rtype: a json {job_id:<the job id>}
    '''
    request = { "id": "gdv_post",
                "mail": gdv_email,
                "key": gdv_key,
                "command": "new_track",
                "project_id": str(project_id),
                "url": str(url) }
    if name != None:
        request['name']=name
    gdv_url = normalize_url(gdv_url)+"/post"
    return urllib2.urlopen( gdv_url, urllib.urlencode(request) ).read()

def add_gdv_sqlite( gdv_key, gdv_email,
                    project_id,
                    url,
                    name=None,
                    gdv_url="http://svitsrv25.epfl.ch/gdv"):
    '''
    Deprecated :  use add_gdv_track instead
    '''
    return add_gdv_track(gdv_key,gdv_email,project_id,url,name,gdv_url)


def add_sql_files( gdv_key, gdv_email,
                   project_id,
                   files, names,
                   serv_url="http://htsstation.vital-it.ch/lims/chipseq/chipseq_minilims.files",
                   gdv_url="http://svitsrv25.epfl.ch/gdv"):
    '''
    Run `add_gdv_sqlite` on a list of files
    '''
    serv_url = normalize_url(serv_url)
    return [add_gdv_track( gdv_key, gdv_email, project_id,
                           serv_url+"/"+f, names[i],
                           gdv_url)
            for i,f in enumerate(files)]

def get_job_status(gdv_key,gdv_email,job_id):
    '''
    Get the status of a job in GDV
    :rtype: a json {job_id:<the job id>, status:<running,error or success>}
    '''
    request = { "id": "gdv_post",
                "mail": gdv_email,
                "key": gdv_key,
                "command": "status",
                "job_id": job_id,
                "gdv_url":"http://svitsrv25.epfl.ch/gdv"}
    gdv_url = normalize_url(gdv_url)+"/post"
    return urllib2.urlopen( gdv_url, urllib.urlencode(request) ).read()


def get_assemblies(gdv_key,gdv_email):
    '''
    Get all assemblies that are used in GDV
    :rtype: a JSON list [{id:<assembly id>,name:<assembly name>}
    '''
    request = { "id": "gdv_post",
                "mail": gdv_email,
                "key": gdv_key,
                "command":"assemblies" }
    gdv_url = normalize_url(gdv_url)+"/post"
    return urllib2.urlopen( gdv_url, urllib.urlencode(request) ).read()

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

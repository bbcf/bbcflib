"""
===================
Module: bbcflib.gdv
===================

Python API for the GDV genome viewer.
"""
import urllib, urllib2
import json
from bbcflib.common import normalize_url
############ GDV requests ############
                                                                                                                                  
def create_gdv_project( gdv_key, gdv_email,
                        name, run_key, nr_assembly_id,
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
                   gdv_url="http://svitsrv25.epfl.ch/gdv",
                   datatype="quantitative",
                   sqlite = False ):
    '''
    Add a new track on a project on GDV
    :param gdv_email: your login in TEQUILA
    :param gdv_key: your user key (get it from your GDV account)
    :param name: name of the track -optionnal- (will take the file name by default)
    :param project_id: the project id to add the track
    :param url : the URL where to fetch the file
    '''
    request = { "id": "gdv_post",
                "mail": gdv_email,
                "key": gdv_key,
                "command": "add_track",
                "project_id": str(project_id),
                "url": normalize_url(str(url)) }
    if sqlite:
        request["command"] = "add_sqlite"
        request["datatype"] = datatype
    if name != None: 
        request['name']=name
    gdv_url = normalize_url(gdv_url)+"/post"
    return urllib2.urlopen( gdv_url, urllib.urlencode(request) ).read()
    
def add_gdv_sqlite( gdv_key, gdv_email,
                    project_id,
                    url,
                    name=None,
                    gdv_url="http://svitsrv25.epfl.ch/gdv",
                    datatype="quantitative" ):
    '''
    Backward compatible interface to `add_gdv_track` for sqlite files.
    '''
    return add_gdv_track(gdv_key,gdv_email,project_id,url,name,gdv_url,datatype,True)


def add_sql_files( gdv_key, gdv_email,
                   project_id,
                   files, names,
                   serv_url="http://htsstation.vital-it.ch/lims/chipseq/chipseq_minilims.files",
                   gdv_url="http://svitsrv25.epfl.ch/gdv",
                   datatype="quantitative" ):
    '''
    Run `add_gdv_sqlite` on a list of files
    '''
    serv_url = normalize_url(serv_url)
    return [add_gdv_track( gdv_key, gdv_email, project_id,
                           serv_url+"/"+f, names[i],
                           gdv_url, dataype, True ) 
            for i,f in enumerate(files)]

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

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


import warnings

################################################################################
# GDV requests #
def new_project(mail, key, name, assembly_id, serv_url='http://gdv.epfl.ch/pygdv'):
    '''
    Create a new project on GDV.
    :param mail : login in TEQUILA
    :param key : an user-specific key (ask it to GDV admin)
    :param name : name of the project
    :param assembly_id : the assembly identifier in GenRep (must be BBCF_VALID)
    :return a JSON
    '''
    query_url = '%s/%s' % (normalize_url(serv_url), 'projects/create')
    request = {'mail':mail, 
               'key':key,
               'name':name,
               'assembly':assembly_id
               }
    return send_it(query_url, request)

def new_track(mail, key, assembly_id=None, project_id=None, urls=None, url=None, fsys=None, fsys_list=None, serv_url='http://gdv.epfl.ch/pygdv', file_names=None):
    '''
    Create a new track on GDV.
    :param mail : login in TEQUILA
    :param key : an user-specific key (ask it to GDV admin)
    :param name : name of the project
    :param assembly_id : the assembly identifier in GenRep (must be BBCF_VALID). Could be optional if a project_id is specified.
    :param project_id : the project identifier to add the track to.
    :param url : an url pointing to a file.
    :param urls : a list of urls separated by whitespaces. 
    :param fsys : if the file is on the same file system, a filesystem path
    :param fsys_list :  a list of fsys separated by whitespaces.
    :param file_names : a list of file name, in the same order than the files uploaded.
    If there is differents parameters given, the first file uploaded will be file_upload, 
    then urls, url, fsys and finally fsys_list. The list is separated by whitespaces.
    :return a JSON
    '''
    query_url = '%s/%s' % (normalize_url(serv_url), 'tracks/create')
    request = {'mail' : mail, 
               'key' : key}
    if assembly_id :
        request['assembly'] = assembly_id
    if project_id :
        request['project_id'] = project_id
    if urls :
        request['urls'] = urls
    if url :
        request['url'] = url
    if fsys :
        request['fsys'] = fsys
    if fsys_list :
        request['fsys_list'] = fsys_list
    if file_names :
        request['file_names'] = file_names
    return send_it(query_url, request)

def send_it(url, request, return_type='json'):
    '''
    Send the request to GDV and return the result. As JSON or as a request.read().
    '''
    req = urllib2.urlopen(url, urllib.urlencode(request))
    if return_type == 'json':
        return json.load(req)
    else :
        return req.read()



def transfert(prefix, file_path, delimiter, mail, key, serv_url, datatype, nr_assemblies=False):
    '''
    Tranfert files from old GDV to new one. A file contains all information on the files to transfert.
    The file is a simple tsv file with columns : (sql file name, track name, assembly identifier)
    :param prefix : where old files are
    :param file_path : file containing the list of files to transfert
    :param delimiter : the delimiter in the file
    :param mail : login in TEQUILA
    :param key : an user-specific key (ask it to GDV admin)
    :param serv_url : the server url of the new GDV
    :param nr_assemblies : True if the list of ids is nr_assembly and not an assembly id.
    '''
    import csv, os, sqlite3
    # read the .psv file
    with open(file_path) as f:
        sha1s = []
        correspond = None
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            sql_file = row[0].strip()
            file_name = row[1].strip()
            assembly_id = int(row[2])
            
            # change nr_assembly_id to an assembly id
            if nr_assemblies:
                if correspond is None:
                    from . import genrep
                    assemblies = genrep.GenRep().get_genrep_objects('assemblies', 'assembly')
                    correspond = dict([ (ass.nr_assembly_id, ass.id) for ass in assemblies if ass.bbcf_valid])
                assembly_id = correspond[assembly_id]
            # look if it's not already done
            if sql_file not in sha1s:
                sha1s.append(sql_file)
                # ensure that the datatype is set
                conn = sqlite3.connect(str(os.path.join(prefix, sql_file)))
                c = conn.cursor()
                dt = c.execute('select value from attributes where key = "datatype" limit 1;').fetchone()[0]
                if dt is None or dt != datatype:
                    c.execute('insert into attributes values (?, ?)', ('datatype', datatype))
                conn.commit()
                c.close()
                conn.close()
                # send the request
                new_track(mail, key, assembly_id=assembly_id, fsys=os.path.join(prefix, sql_file), serv_url=serv_url, file_names=file_name)














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
    warnings.simplefilter('always')
    request = { "id": "gdv_post",
                "mail": gdv_email,
                "key": gdv_key,
                "command": "new_project",
                "name": str(name),
                "seq_id": str(nr_assembly_id),
                "public": str(public).lower() }
    gdv_url = normalize_url(gdv_url)+"/post"
    warnings.warn('This method is used by the old version of GDV', DeprecationWarning)
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
    warnings.warn('This method is used by the old version of GDV', DeprecationWarning)
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
    warnings.warn('This method is used by the old version of GDV', DeprecationWarning)
    return add_gdv_track(gdv_key,gdv_email,project_id,url,name,gdv_url)


def add_sql_files( gdv_key, gdv_email,
                   project_id,
                   files, names,
                   serv_url="http://htsstation.vital-it.ch/lims/chipseq/chipseq_minilims.files",
                   gdv_url="http://svitsrv25.epfl.ch/gdv"):
    '''
    Run `add_gdv_sqlite` on a list of files
    '''
    warnings.warn('This method is used by the old version of GDV', DeprecationWarning)
    serv_url = normalize_url(serv_url)
    return [add_gdv_track( gdv_key, gdv_email, project_id,
                           serv_url+"/"+f, names[i],
                           gdv_url)
            for i,f in enumerate(files)]

def get_job_status(gdv_key, gdv_email, job_id, gdv_url):
    '''
    Get the status of a job in GDV
    :rtype: a json {job_id:<the job id>, status:<running,error or success>}
    '''
    request = { "id":       "gdv_post",
                "mail":     gdv_email,
                "key":      gdv_key,
                "command":  "status",
                "job_id":   job_id,
                "gdv_url":  gdv_url}
    gdv_url = normalize_url(gdv_url)+"/post"
    warnings.warn('This method is used by the old version of GDV', DeprecationWarning)
    return urllib2.urlopen( gdv_url, urllib.urlencode(request) ).read()


def get_assemblies(gdv_key, gdv_email, gdv_url):
    '''
    Get all assemblies that are used in GDV
    :rtype: a JSON list [{id:<assembly id>, name:<assembly name>}
    '''
    request = { "id": "gdv_post",
                "mail": gdv_email,
                "key": gdv_key,
                "command":"assemblies" }
    gdv_url = normalize_url(gdv_url)+"/post"
    warnings.warn('This method is used by the old version of GDV', DeprecationWarning)
    return urllib2.urlopen( gdv_url, urllib.urlencode(request) ).read()

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

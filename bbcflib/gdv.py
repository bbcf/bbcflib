"""
===================
Module: bbcflib.gdv
===================

Python API for the GDV genome viewer.
"""

# Built-in modules #
import json, urllib, urllib2

# Internal modules #
from bbcflib.common import normalize_url

default_url = 'http://gdv.epfl.ch/pygdv'

################################################################################
def get_project(mail, key, project_key, serv_url=default_url):
    '''
    Get a project by it's key.
    :param mail : login in TEQUILA
    :param key : an user-specific key (ask it to GDV admin)
    :param project_key : the project key
    '''
    query_url = '%s/%s' % (normalize_url(serv_url), 'projects/get')
    request = {'mail':mail, 
               'key':key,
               'project_key':project_key
               }
    return send_it(query_url, request)


def delete_project(mail, key, project_id, serv_url=default_url):
    query_url = '%s/%s/%s' % (normalize_url(serv_url), 'projects/delete', project_id)
    request = {'mail':mail,
               'key':key}
    return send_it(query_url, request)

def delete_track(mail, key, track_id, serv_url=default_url):
    query_url = '%s/%s/%s' % (normalize_url(serv_url), 'tracks/delete', track_id)
    request = {'mail':mail,
               'key':key}
    return send_it(query_url, request)

def new_project(mail, key, name, assembly_id, serv_url=default_url):
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


def single_track(mail, key, serv_url=default_url,
                 assembly_id=None, project_id=None, 
                 url=None, fsys=None, trackname=None, extension=None, 
                 force=False, delfile=False):
    '''
    Create a new track on GDV.
    :param mail : login in TEQUILA
    :param key : an user-specific key (ask it to GDV admin)
    :param name : name of the project
    :param assembly_id : the assembly identifier in GenRep (must be BBCF_VALID). Could be optional if a project_id is specified.
    :param project_id : the project identifier to add the track to.
    :param url : an url pointing to a file.
    :param extension : extension of the file provided.
    :param fsys : if the file is on the same file system, a filesystem path
    :param trackname : the name to give to the track 
    :param force : A boolean. Force the file to be recomputed.
    :param delfile : If true and file comming from fsys, the original file will be removed after job success
    :return a JSON
    '''
    query_url = '%s/%s' % (normalize_url(serv_url), 'tracks/create')
    request = {'mail' : mail, 'key' : key}
    if assembly_id: request['assembly'] = assembly_id
    if project_id:  request['project_id'] = project_id
    if url:         request['url'] = url
    if fsys:        request['fsys'] = fsys
    if trackname:   request['trackname'] = trackname
    if force:       request['force'] = force
    if extension:   request['extension'] = extension
    if delfile:     request['delfile'] = delfile
    return send_it(query_url, request)


def multiple_tracks(mail, key, serv_url=default_url,
                    assembly_id=None, project_id=None, 
                    urls=[], fsys_list=[], tracknames=[], extensions=[], 
                    force=False):
    '''
    Create tracks on GDV
    :param extensions : a list of extensions, separated by whitespace.
    :param urls : a list of urls separated by whitespaces. 
    :param fsys : if the file is on the same file system, a filesystem path
    :param fsys_list :  a list of fsys separated by whitespaces.
    :param tracknames : a list of file name, in the same order than the files uploaded.
    If there is differents parameters given, the first file uploaded will be file_upload, 
    then urls, url, fsys and finally fsys_list. The list is separated by whitespaces.
    For other params :see single_track
    '''
    tracks = []
    if isinstance(tracknames,basestring): tracknames = tracknames.split()
    if isinstance(extensions,basestring): extensions = extensions.split()
    if isinstance(urls,basestring): urls = urls.split()
    if isinstance(fsys_list,basestring): fsys_list = fsys_list.split()
    ntracks = len(urls)+len(fsys_list)
    tr = tracknames+[None]*(ntracks-len(tracknames))
    ex = extensions+[None]*(ntracks-len(extensions))
    tracks = [single_track(mail, key, serv_url=serv_url, 
                           assembly_id=assembly_id, project_id=project_id, 
                           url=u, trackname=tr[n], extension=ex[n], force=force) 
              for n,u in enumerate(urls)] \
              + [single_track(mail, key, serv_url=serv_url, 
                              assembly_id=assembly_id, project_id=project_id, 
                              fsys=f, trackname=tr[len(urls)+n], extension=ex[len(urls)+n], force=force) 
                 for n,f in enumerate(fsys_list)]
    return tracks

def send_it(url, request, return_type='json'):
    '''
    Send the request to GDV and return the result. As JSON or as a request.read().
    '''
    req = urllib2.urlopen(url, urllib.urlencode(request))
    if return_type == 'json': return json.load(req)
    return req.read()



#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

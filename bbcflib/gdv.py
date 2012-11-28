"""
===================
Module: bbcflib.gdv
===================

Python API for the GDV genome data viewer.
"""

import json, urllib, urllib2
from bbcflib.common import normalize_url

GDV_SERVER_URL = 'http://gdv.epfl.ch/pygdv'


def new_project(mail, key, name, assembly_id, serv_url=GDV_SERVER_URL):
    '''
    Create a new project on GDV.
    mail: your contact email (same as in GDV)
    key: your secret key
    name: name of the project
    assembly_id: the assembly identifier in GenRep (must be BBCF_VALID)
    serv_url: GDV url.
    return a JSON
    '''
    try:
        int(assembly_id)
    except:
        raise
    query_url = '%s/%s' % (serv_url, 'projects/create')
    request = {'mail': mail,
               'key': key,
               'name': name,
               'assembly': assembly_id
               }
    return send_request(query_url, request)


def get_project(mail, key, project_key='', project_id='', serv_url=GDV_SERVER_URL):
    '''
    Get a project by it's key, or id.
    If you don't specify an id or a key, it will fetch all projects.
    mail: your contact email (same as in GDV)
    key: your secret key
    project_key: the project key
    project_id: the project_id
    return a JSON
    '''
    query_url = '%s/%s' % (serv_url, 'projects/get')
    request = {'mail': mail,
               'key': key,
               'project_key': project_key,
               'project_id': project_id
               }
    return send_request(query_url, request)


def new_track(mail, key, serv_url=GDV_SERVER_URL,
                 assembly_id='', project_id='',
                 url='', fsys='', trackname='', extension='',
                 force='', delfile=False):
    '''
    Create a new track on GDV.
    mail: your contact email (same as in GDV)
    key: your secret key
    assembly_id: the assembly identifier in GenRep. Could be optional if a project_id is specified.
    project_id: the project identifier to add the track to.
    url: an url pointing to a file.
    extension: extension of the file provided.
    fsys: if the file is on the same file system, a filesystem path
    trackname: the name to give to the track.
    force: A boolean. Force the file to be recomputed.
    delfile: If true and file comming from fsys, the original file will be removed after job success.
    return a JSON
    '''
    query_url = '%s/%s' % (serv_url, 'tracks/create')
    request = {'mail': mail, 'key': key}
    if assembly_id:
        request['assembly'] = assembly_id
    if project_id:
        request['project_id'] = project_id
    if url:
        request['url'] = url
    if fsys:
        request['fsys'] = fsys
    if trackname:
        request['trackname'] = trackname
    if force:
        request['force'] = force
    if extension:
        request['extension'] = extension
    if delfile:
        request['delfile'] = delfile
    return send_request(query_url, request)


import json
import urllib
import urllib2


def send_request(url, request, return_type='json'):
    '''
    Send the request to GDV and return the result.
    As JSON or as a request.read().
    :param url: the url of GDV.
    :param request: the request.
    :return: GDV response as json or a stream.
    '''
    print '[pygdv] Sending request on %s with body : %s.' % (url, request)
    req = urllib2.urlopen(url, urllib.urlencode(request))

    ret = None
    if return_type == 'json':
        ret = json.load(req)
    else:
        ret = req.read()
    return ret

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

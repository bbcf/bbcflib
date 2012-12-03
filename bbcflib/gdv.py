"""
===================
Module: bbcflib.gdv
===================

Python API for the GDV genome data viewer.
"""

import json, urllib, urllib2
from bbcflib.common import normalize_url

default_url = 'http://gdv.epfl.ch/pygdv'

################################################################################
def _gdv_request(**kw):
    req_keys = ['mail','key','project_key','project_id','name','assembly',
                'url','fsys','trackname','force','extension','delfile']
    request = dict((k,kw[k]) for k in req_keys if kw.get(k))
    kw['serv_url'] = normalize_url(kw['serv_url'])
    url = "/".join([kw[k] for k in ['serv_url','obj','action','id'] if k in kw])
    req = urllib2.urlopen(url, urllib.urlencode(request))
    if kw.get('return_type','') == 'json': return json.load(req)
    return req.read()

def get_project(mail, key, project_key, serv_url=default_url):
    '''
    Retrieve project information from its key.

    :param mail: email to login via TEQUILA.
    :param key: GDV user key.
    :param project_key: the project key.
    :param serv_url: GDV's url.
    :rtype: JSON
    '''
    return _gdv_request(mail=mail, key=key, project_key=project_key, serv_url=serv_url,
                        obj='projects', action='get', return_type='json')

def new_project(mail, key, name, assembly_id, serv_url=default_url):
    '''
    Create a new project on GDV.

    :param mail: email to login via TEQUILA.
    :param key: GDV user key.
    :param name: project name.
    :param assembly_id: GenRep assembly identifier (must be BBCF_VALID).
    :rtype: JSON
    '''
    return _gdv_request(mail=mail, key=key, name=name, assembly=assembly_id,
                        serv_url=serv_url, obj='projects', action='create', return_type='json')


def delete_project(mail, key, project_id, serv_url=default_url):
    return _gdv_request(mail=mail, key=key, serv_url=serv_url, id=project_id, 
                        obj='projects', action='delete', return_type='json')


def single_track(mail, key, assembly_id=None, project_id=None,
                 url=None, name=None, extension=None,
                 delete_target=False, delete_source=False,
                 serv_url=default_url):
    '''
    Create a new track on GDV.

    :param mail: email to login via TEQUILA.
    :param key: GDV user key.
    :param assembly_id: GenRep assembly identifier (must be BBCF_VALID), optional if a project_id is specified.
    :param project_id: the project identifier to add the track to.
    :param url: an url or path pointing to a file.
    :param extension: file extension.
    :param name: the track name.
    :param delete_target: Force the file to be recomputed (boolean).
    :param delete_source: If true and file is a local parth, the original file will be removed after job success.
    :rtype: JSON
    '''
    _url = None
    _fsys = None
    if url.startswith(("http://","https://","ftp://")): _url = url
    else: _fsys = url
    return _gdv_request(mail=mail, key=key, 
                        assembly=assembly_id, project_id=project_id, 
                        url=_url, fsys=_fsys, trackname=name, extension=extension,
                        force=delete_target, delfile=delete_source, 
                        serv_url=serv_url, obj='tracks', action='create', return_type='json')

def multiple_tracks(mail, key, assembly_id=None, project_id=None,
                    urls=[], names=[], extensions=[],
                    force=False, serv_url=default_url):
    '''
    Add multiple tracks to GDV.

    :param urls: a list of urls.
    :param extensions: a list of extensions, in the same order as the urls.
    :param names: a list of track names, in the same order as the urls.

    For other params: see single_track
    '''
    tracks = []
    if isinstance(names,basestring): tracknames = tracknames.split()
    if isinstance(extensions,basestring): extensions = extensions.split()
    if isinstance(urls,basestring): urls = urls.split()
    ntracks = len(urls)
    tr = names+[None]*(ntracks-len(names))
    ex = extensions+[None]*(ntracks-len(extensions))
    tracks = [single_track(mail, key,
                           assembly_id=assembly_id, project_id=project_id,
                           url=u, name=tr[n], extension=ex[n], force=force, serv_url=serv_url)
              for n,u in enumerate(urls)]
    return tracks

def delete_track(mail, key, track_id, serv_url=default_url):
    return _gdv_request(mail=mail, key=key, serv_url=serv_url, id=track_id, 
                        obj='tracks', action='delete', return_type='json')



#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

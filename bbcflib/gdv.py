import urllib, urllib2

############ GDV requests ############
def create_gdv_project( gdv_key, gdv_email,
                        name, run_key, nr_assembly_id,
                        gdv_url="http://svitsrv25.epfl.ch/gdv/post" ):
    request = { "id": "gdv_post",
                "mail": gdv_email,
                "key": gdv_key,
                "command": "new_project",
                "type": "htschipseq",
                "obfuscated": run_key,
                "name": name,
                "seq_id": nr_assembly_id }
    project_id = urllib2.urlopen( gdv_url, urllib.urlencode(request) ).read()
    return project_id

def add_gdv_track( gdv_key, gdv_email,
                   files, project_id,
                   server_url="http://htsstation.vital-it.ch/chipseq/results/",
                   gdv_url="http://svitsrv25.epfl.ch/gdv/post" ):
    request = { "id": "gdv_post",
                "mail": gdv_email,
                "key": gdv_key,
                "command": "add_sqlite",
                "project_id": project_id,
                "url": '',
                "datatype": "quantitative" }
    for f in files:
        request['url'] = server_url+f
        urllib2.urlopen( gdv_url, urllib.urlencode(request) ).read()
    return 
############################################################

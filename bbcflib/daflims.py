"""
===============
bbcflib.daflims
===============

This module provides an interface to the DNA Array Facility's LIMS.
To use it, create a DAFLIMS object with your username and password to
the facility, then call the ``fetch`` method on the object to get a
particular file.

Files in the LIMS are identified by four fields: facility, machine, run, and lane.  For instance, lane 3 of run 91 sequenced on R2D2 at lgtf.  To fetch this file to 'myfile.fastq' in the current working directory, we would write::

    d = DAFLIMS(username=..., password=...)
    d.fetch('lgtf', 'R2D2', 91, 3, 'myfile.fastq')
"""
import urllib2
import tarfile
import re

class DAFLIMS(object):
    def __init__(self, username=None, password=None, config=None, section="daflims"):
        if (username==None or password==None) and config==None:
            raise TypeError("Must provide a username and password, or a ConfigParser.")
        elif config != None:
            if username == None: self.username = config.get(section, 'daflims_username')
            else:                self.username = username

            if password == None: self.password = config.get(section, 'daflims_password')
            else:                self.password = password
        else:
            self.username = username
            self.password = password

    def fetch(self, facility, machine, run, lane, to=None):
        if not(isinstance(facility, str)):
            raise ValueError("facility must be a string, found %s" % str(facility))
        if not(isinstance(machine, str)):
            raise ValueError("machine must be a string, found %s" % str(machine))
        if not(isinstance(run, int)):
            raise ValueError("run must be an integer, found %s" % str(run))
        if not(isinstance(lane, int)):
            raise ValueError("lane must be an integer, found %s" % str(lane))

        base_url="http://uhts-%s.vital-it.ch" % facility
        auth_handler = urllib2.HTTPDigestAuthHandler()
        auth_handler.add_password(realm="UHTS-LIMS-ws",
                                  uri=base_url,
                                  user=self.username,
                                  passwd=self.password)
        opener = urllib2.build_opener(auth_handler)
        urllib2.install_opener(opener)
        url = "/".join([base_url,"ws","symlinkname", machine, str(run), str(lane)])
        s = urllib2.urlopen(url).read()
        status = re.search(r'==(\w+)\s',s).groups()[0]
        links = re.search(r'\n(.*)\n',s).groups()[0].split("\t")
        if status == "DATA":
            link_name = links[0]
        else:
            raise ValueError("Request "+url+"\n"+str(links))
        url = "/".join([base_url,"ws","lanedesc"]+sample_descr[1:4])
        s = urllib2.urlopen(url).read()
        status = re.search(r'==(\w+)\s',s).groups()[0]
        lanedesc = re.search(r'\n(.*)\n',s).groups()[0].split("\t")
        if status == "DATA":
            lib_name = lanedesc[4]
        else:
            raise ValueError("Request "+url+"\n"+lanedesc)
        url = "/".join([base_url,"symlink",link_name])
        tar = tarfile.open(fileobj=urllib2.urlopen(url),mode="r|gz")
        if to == None:
            file_loc = os.getcwd() + unique_filename_in(os.getcwd())
        elif os.path.isdir(to):
            file_loc = to + unique_filename_in(to)
        else:
            file_loc = to
        tar.extractall(path=file_loc)
        fastqname = tar.getnames()[0]
        tar.close()
        return {lib_name: file_loc+"/"+fastqname}

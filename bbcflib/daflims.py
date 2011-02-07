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
    d.fetch('lgtf', 'R2D2', 91, 3, 'path/to/save/to')
"""
import urllib2
import tarfile
import re
import os
from urlparse import urlparse
from bein import unique_filename_in

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
        self.url_opener = None

    def open_url(self, url):
        """Same as run_method, but doesn't read the URL, just passes its handle back."""
        auth_handler = urllib2.HTTPDigestAuthHandler()
        parse_url = urlparse(url)
        base_url = parse_url.scheme + '://' + parse_url.netloc
        auth_handler.add_password(realm="UHTS-LIMS-ws",
                                  uri=base_url,
                                  user=self.username,
                                  passwd=self.password)
        opener = urllib2.build_opener(auth_handler)
        urllib2.install_opener(opener)
        return urllib2.urlopen(url)

    def run_method(self, method, facility, *args):
        """Runs *method* on the DAF LIMS with the rest of the arguments to the function.

        The rest of the arguments are converted to strings with str before being used
        in the URL.
        """
        base_url="http://uhts-%s.vital-it.ch" % facility
        url = "/".join([base_url,"ws",method] + [str(x) for x in args])
        return self.open_url(url).read()

    def check_description(self, facility, machine, run, lane):
        if not(isinstance(facility, str)):
            raise ValueError("facility must be a string, found %s" % str(facility))
        if not(isinstance(machine, str)):
            raise ValueError("machine must be a string, found %s" % str(machine))
        if not(isinstance(run, int)):
            raise ValueError("run must be an integer, found %s" % str(run))
        if not(isinstance(lane, int)):
            raise ValueError("lane must be an integer, found %s" % str(lane))

    def symlinkname(self, facility, machine, run, lane):
        """Fetch the names of the tar files matching the given description.

        The DAF LIMS stores three tar files for each combination of
        facility, machine, run, and lane, one for FASTQ, one for
        ELAND, and one for QSEQ.  symlinkname returns a list of
        dictionaries, each containing the keys 'fastq', 'eland', and
        'qseq'.
        """
        self.check_description(facility, machine, run, lane)
        response = self.run_method("symlinkname", facility, machine, run, lane).splitlines()
        if re.search('==DATA', response[0]) == None:
            raise ValueError(("symlinkname method failed on DAFLIMS (facility='%s', " + \
                              "machine='%s', run=%d, lane=%d): %s") % (facility, machine, run, lane,
                                                                     '\n'.join(response[1:])))
        if len(response) > 2:
            raise ValueError("symlinkname method returned multiple records: %s" % ('\n'.join(response[1:])))
        else:
            q = response[1].split('\t')
            return {'fastq': q[0], 'eland': q[1], 'qseq': q[2]}

    def lanedesc(self, facility, machine, run, lane):
        self.check_description(facility, machine, run, lane)
        response = self.run_method("lanedesc", facility, machine, run, lane).splitlines()
        if re.search('==DATA', response[0]) == None:
            raise ValueError(("lanedesc method failed on DAFLIMS (facility='%s', " + \
                              "machine='%s', run=%d, lane=%d): %s") % (facility, machine, run, lane,
                                                                     '\n'.join(response[1:])))
        if len(response) > 2:
            raise ValueError("lanedesc method returned multiple records: %s" % ('\n'.join(response[1:])))
        else:
            q = response[1].split('\t')
            return {'machine': q[0],
                    'run': int(q[1]),
                    'cycle': int(q[2]),
                    'lane': int(q[3]),
                    'library': q[4],
                    'quantity (/pM)': float(q[5]),
                    'protocol': q[6],
                    'organism': q[7],
                    'project': q[8],
                    'PI firstname': q[9],
                    'PI lastname': q[10],
                    'submitter firstname': q[11],
                    'submitter lastname': q[12],
                    'run type': q[13],
                    'NCBI ID': q[14]}


    def fetch_symlink(self, link_name, to=None):
        if to == None:
            target = os.path.join(os.getcwd(), unique_filename_in())
        elif os.path.isdir(to):
            target = os.path.join(to, unique_filename_in())
        else:
            target = to

        # We get a tar archive reading from the URL
        url = self.open_url(link_name)
        tar = tarfile.open(fileobj=url, mode='r|gz')

        tar_filename = tar.next()

        input_file = tar.extractfile(tar_filename)
        with open(target, 'w') as output_file:
            while True:
                chunk = input_file.read(4096000)
                if chunk == '':
                    break
                else:
                    output_file.write(chunk)
        output_file.close()
        input_file.close()
        tar.close()
        return target

    def fetch_structured(self, type, facility, machine, run, lane, to=None):
        """Fetches a FASTQ file to *to*, and returns all relevant information."""
        if type != 'fastq' and type != 'eland' and type != 'qseq':
            raise ValueError("type must be one of 'fastq', 'eland', or 'qseq'; found %s" % str(type))
        url = self.symlinkname(facility, machine, run, lane)[type]
        info = self.lanedesc(facility, machine, run, lane)
        filename = self.fetch_symlink(url, to)
        info['path'] = filename
        return info

    def fetch_fastq(self, facility, machine, run, lane, to=None):
        return self.fetch_structured('fastq', facility, machine, run, lane, to)

    def fetch_eland(self, facility, machine, run, lane, to=None):
        return self.fetch_structured('eland', facility, machine, run, lane, to)

    def fetch_qseq(self, facility, machine, run, lane, to=None):
        return self.fetch_structured('qseq', facility, machine, run, lane, to)
        

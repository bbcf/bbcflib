"""
=======================
Module: bbcflib.daflims
=======================

The DNA Array Facility at the University of Lausanne provides all
their data in a LIMS accessible from VITAL-IT.  This module hides the
details of fetching files from that LIMS.  Files are identified by
four fields:

  * The facility where they were generated (either ``'lgtf'`` for Lausanne,
    or ``'gva'`` at the University of Geneva).

  * The machine at the facility on which they were generated (``'R2D2'`` or
    ``'C3PO'`` in Lausanne, or ``'HWUSI-EAS691'`` in Geneva).

  * The run number on that machine (an integer).

  * The lane in that run (also an integer).

Connecting to the LIMS requires a username and password, and will
often only be possible from certain hosts.  A connection is
represented as a DAFLIMS object, which provides methods for fetching
the three kinds of files stored in the LIMS:

  * ``fetch_fastq``: The FASTQ file of all reads.  You probably want
    this.

  * ``fetch_eland``: The output of aligning all reads against some
    genome with Eland, a proprietary tool from Illumina.

  * ``fetch_qseq``: QSeq is the text format output directly by the
    sequencer.

For example, to fetch the FASTQ file of run 91, lane 3 sequenced on
R2D2 in Lausanne, write::

    d = DAFLIMS(username=..., password=...)
    d.fetch_fastq('lgtf', 'R2D2', 91, 3, '/path/to/save/to.fastq')

The last argument works much like a target given to the Unix command
``cp``: if it specifies a directory, the file is given a random name
in that directory.  If it is given a filename, the file is written to
that name.  If the last argument is omitted, the file is written to a
random name in the current working irectory.

  .. autoclass:: DAFLIMS
"""

# Built-in modules #
import os, re, tarfile, urllib2
from urlparse import urlparse

# Other modules #
from bein import unique_filename_in

################################################################################
class DAFLIMS(object):
    """Connect to the DNA Array Facility's LIMS.

    Either specify *username* and *password*, or provide a
    ``ConfigParser`` object to the *config* argument.  The
    ConfigParser should define ``daflims_username`` and
    ``daflims_password`` in the sections ``daflims`` (you can override
    the section with the *section* argument).  If you specify a
    *config* argument, and one or both of *username* and *password*,
    the latter override the ConfigParser.

      .. py:method:: fetch_fastq(facility, machine, run, lane, to=None)

      .. py:method:: fetch_eland(facility, machine, run, lane, to=None)

      .. automethod:: fetch_qseq
    """
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

    def _open_url(self, url):
        """Returns a file-like object connected to a *url*.

        The *url* is assumed to be part of the DNA Array Facility's
        LIMS, and logs in with the ``DAFLIMS`` object's username and
        password.
        """
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

    def _run_method(self, method, facility, *args):
        """Returns a string from running *method* on the LIMS.

        Methods on the LIMS are accessed with a uniform URL format.
        *method* is the method to run (for example ``"symlinkname"`` or
        ``"lanedesc"``).  *facility* is the facility to run against, and
        any further arguments to the function are appended as path
        components of the URL.
        """
        base_url="http://uhts-%s.vital-it.ch" % facility
        url = "/".join([base_url,"ws",method] + [str(x) for x in args])
        return self._open_url(url).read()

    def _check_description(self, facility, machine, run, lane):
        """Check if the identifier of a file is valid.

        *facility* and *machine* should be strings.  *run* and *lane*
        should be integers.  This is checked by a couple of
        functions, so we factor it out here.
        """
        if not(isinstance(facility, str)):
            raise ValueError("facility must be a string, found %s" % str(facility))
        if not(isinstance(machine, str)):
            raise ValueError("machine must be a string, found %s" % str(machine))
        if not(isinstance(run, int)):
            raise ValueError("run must be an integer, found %s" % str(run))
        if not(isinstance(lane, int)):
            raise ValueError("lane must be an integer, found %s" % str(lane))

    def _symlinkname(self, facility, machine, run, lane):
        """Fetch the URLs to access data in the LIMS.

        Returns a dictionary with the keys ``'fastq'``, ``'eland'``,
        and ``'qseq'``, each referring to a URL in the LIMS where that
        file is stored.
        """
        self._check_description(facility, machine, run, lane)
        response = self._run_method("symlinkname", facility, machine, run, lane).splitlines()
        if re.search('==DATA', response[0]) == None or len(response)<2:
            raise ValueError(("symlinkname method failed on DAFLIMS (facility='%s', " + \
                              "machine='%s', run=%d, lane=%d): %s") % (facility, machine, run, lane,
                                                                     '\n'.join(response[1:])))
        else:
            q = response[1].split('\t')
            return {'fastq': q[0], 'eland': q[1], 'qseq': q[2]}

    def _lanedesc(self, facility, machine, run, lane):
        """Fetch the metadata of particular data set in the LIMS.

        Returns a dictionary with the following keys (which refer to strings unless otherwise notes):

          * ``'machine'``
          * ``'run'`` (integer)
          * ``'lane'`` (integer)
          * ``'cycle'`` (integer)
          * ``'quantity (/pM)'`` (float, the concentration of DNA in picomolar)
          * ``'library'``
          * ``'project'``
          * ``'protocol'``
          * ``'run type'`` (``'ChipSeq'``, etc.)
          * ``'PI firstname'``
          * ``'PI lastname'``
          * ``'submitter firstname'``
          * ``'submitter lastname'``
          * ``'organism'``
          * ``'NCBI ID'`` (the genome Eland aligns against)
        """
        self._check_description(facility, machine, run, lane)
        response = self._run_method("lanedesc", facility, machine, run, lane).splitlines()

        # '==DATA' on the first line means success.  It should be
        # followed by one line of data.  Anything else is an error.
        if re.search('==DATA', response[0]) == None:
            raise ValueError(("lanedesc method failed on DAFLIMS (facility='%s', " + \
                              "machine='%s', run=%d, lane=%d): %s") % (facility, machine, run, lane,
                                                                     '\n'.join(response[1:])))
        if len(response) > 2:
            raise ValueError("lanedesc method returned multiple records: %s" % ('\n'.join(response[1:])))

        # If the response is valid, parse the fields into a dictionary.
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


    def _fetch_symlink(self, link_name, to=None):
        """Fetch the data from a file in the LIMS into *to*.

        *link_name* is a URL to a .tar.gz file in the LIMS.  These
        .tar.gz files all contain only one file, which we write to
        *to*.  If *to* is omitted, then the data is written to a
        randomly named file in the current working directory.  If
        *to* is a directory, the data is written to a randomly named
        file in that directory.  Otherwise *to* is taken as the full
        path to the file to write to.

        ``_fetch_symlink`` returns the path to the output file,
        including its filename.
        """
        if to == None:
            target = os.path.join(os.getcwd(), unique_filename_in())
        elif os.path.isdir(to):
            target = os.path.join(to, unique_filename_in())
        else:
            target = to

        # *link_name* is a URL to a .tar.gz file in the LIMS
        # containing exactly one file.  We stream from the HTTP
        # connection to a local file in 4kb chunks.
        url = self._open_url(link_name)
        tar = tarfile.open(fileobj=url, mode='r|gz')

        # Since the tar file contains exactly one file, calling
        # ``next()`` on the tar gives us the file we want.  We cannot
        # use ``getnames()[0]`` or similar methods, since they scan
        # all the way through the file, and we cannot rewind on HTTP
        # responses.
        tar_filename = tar.next()

        # extractfile returns a file-like object we can stream from.
        input_file = tar.extractfile(tar_filename)
        with open(target, 'w') as output_file:
            while True:
                chunk = input_file.read(4096)
                if chunk == '':
                    break
                else:
                    output_file.write(chunk)
        output_file.close()
        input_file.close()
        tar.close()
        return target

    def _fetch_structured(self, type, facility, machine, run, lane, to=None):
        """Fetch a file from its description.

        The return value is a dictionary containing the following keys
        (all of which refer to strings unless otherwise noted.

          * ``'path'`` (the path to the fetched file)
          * ``'machine'``
          * ``'run'`` (integer)
          * ``'lane'`` (integer)
          * ``'cycle'`` (integer)
          * ``'quantity (/pM)'`` (float, the concentration of DNA in picomolar)
          * ``'library'``
          * ``'project'``
          * ``'protocol'``
          * ``'run type'`` (``'ChipSeq'``, etc.)
          * ``'PI firstname'``
          * ``'PI lastname'``
          * ``'submitter firstname'``
          * ``'submitter lastname'``
          * ``'organism'``
          * ``'NCBI ID'`` (the genome Eland aligns against)

        *type* must be one of ``'fastq'``, ``'eland'``, or ``'qseq'``.
        """
        if type != 'fastq' and type != 'eland' and type != 'qseq':
            raise ValueError("type must be one of 'fastq', 'eland', or 'qseq'; found %s" % str(type))
        url = self._symlinkname(facility, machine, run, lane)[type]
        info = self._lanedesc(facility, machine, run, lane)
        filename = self._fetch_symlink(url, to)
        info['path'] = filename
        return info

    def fetch_fastq(self, facility, machine, run, lane, to=None):
        """Fetch a file from the LIMS to *to*.

        If *to* is omitted, then the data is written to a
        randomly named file in the current working directory.  If
        *to* is a directory, the data is written to a randomly named
        file in that directory.  Otherwise *to* is taken as the full
        path to the file to write to.

        The return value is a dictionary containing the following keys
        (all of which refer to strings unless otherwise noted.

          * ``'path'`` (the path to the fetched file)
          * ``'machine'``
          * ``'run'`` (integer)
          * ``'lane'`` (integer)
          * ``'cycle'`` (integer)
          * ``'quantity (/pM)'`` (float, the concentration of DNA in picomolar)
          * ``'library'``
          * ``'project'``
          * ``'protocol'``
          * ``'run type'`` (``'ChipSeq'``, etc.)
          * ``'PI firstname'``
          * ``'PI lastname'``
          * ``'submitter firstname'``
          * ``'submitter lastname'``
          * ``'organism'``
          * ``'NCBI ID'`` (the genome Eland aligns against)
        """
        return self._fetch_structured('fastq', facility, machine, run, lane, to)

    def fetch_eland(self, facility, machine, run, lane, to=None):
        """Fetch a file from the LIMS to *to*.

        If *to* is omitted, then the data is written to a
        randomly named file in the current working directory.  If
        *to* is a directory, the data is written to a randomly named
        file in that directory.  Otherwise *to* is taken as the full
        path to the file to write to.

        The return value is a dictionary containing the following keys
        (all of which refer to strings unless otherwise noted.

          * ``'path'`` (the path to the fetched file)
          * ``'machine'``
          * ``'run'`` (integer)
          * ``'lane'`` (integer)
          * ``'cycle'`` (integer)
          * ``'quantity (/pM)'`` (float, the concentration of DNA in picomolar)
          * ``'library'``
          * ``'project'``
          * ``'protocol'``
          * ``'run type'`` (``'ChipSeq'``, etc.)
          * ``'PI firstname'``
          * ``'PI lastname'``
          * ``'submitter firstname'``
          * ``'submitter lastname'``
          * ``'organism'``
          * ``'NCBI ID'`` (the genome Eland aligns against)
        """
        return self._fetch_structured('eland', facility, machine, run, lane, to)

    def fetch_qseq(self, facility, machine, run, lane, to=None):
        """Fetch a file from the LIMS to *to*.

        If *to* is omitted, then the data is written to a
        randomly named file in the current working directory.  If
        *to* is a directory, the data is written to a randomly named
        file in that directory.  Otherwise *to* is taken as the full
        path to the file to write to.

        The return value is a dictionary containing the following keys
        (all of which refer to strings unless otherwise noted.

          * ``'path'`` (the path to the fetched file)
          * ``'machine'``
          * ``'run'`` (integer)
          * ``'lane'`` (integer)
          * ``'cycle'`` (integer)
          * ``'quantity (/pM)'`` (float, the concentration of DNA in picomolar)
          * ``'library'``
          * ``'project'``
          * ``'protocol'``
          * ``'run type'`` (``'ChipSeq'``, etc.)
          * ``'PI firstname'``
          * ``'PI lastname'``
          * ``'submitter firstname'``
          * ``'submitter lastname'``
          * ``'organism'``
          * ``'NCBI ID'`` (the genome Eland aligns against)
        """
        return self._fetch_structured('qseq', facility, machine, run, lane, to)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

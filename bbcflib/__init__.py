"""
===========
bbcflib 0.1
===========

bbcflib is a set of Python modules for accessing facilities used by
the Bioinformatics and Biostatistics Core Facility (BBCF) at the EPFL.
At the moment it provides modules for BBCF's GenRep, the LIMS of the
DNA Array Facility at UNIL, the standard frontend code for workflows
written by Fabrice David, and sending job description emails.

All the functionality can be imported with::

    import bbcflib

which reexports the important functions from all the submodules.  It
reexports ``ConfigParser`` from ``ConfigParser``, and the following
modules from bbcflib:

* ``GenRep`` and ``Assembly`` from  :doc:`bbcflib_genrep`
* ``EmailReport`` from :doc:`bbcflib_email`
* ``DAFLIMS`` from :doc:`bbcflib_daflims`
* ``Frontend`` from :doc:`bbcflib_frontend`

"""

from ConfigParser import ConfigParser
from genrep import GenRep, Assembly
from email import EmailReport
from daflims import DAFLIMS
from frontend import Frontend


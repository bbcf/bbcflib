"""
===========
bbcflib 1.1
===========

bbcflib is a set of Python modules for accessing facilities used by
the Bioinformatics and Biostatistics Core Facility (BBCF) at the EPFL.
It provides modules for BBCF's 

* *GenRep* (genome data repository), 
* the LIMS of the DNA Array Facilities at UNIL and UNIGE, 
* the standard frontend code for workflows written by Fabrice David, 
* sending job completion emails.
* the *GDV* api for viewing genome tracks
* *mapseq* to map reads to reference genomes or transcriptomes
* *chipseq* to run ChIP-seq analyses

All the functionality can be imported with::

    import bbcflib

which reexports the important functions from all the submodules.  It
reexports ``ConfigParser`` from ``ConfigParser``, and the following
modules from bbcflib:

* ``GenRep`` and ``Assembly`` from  :doc:`bbcflib_genrep`
* ``EmailReport`` from :doc:`bbcflib_email`
* ``DAFLIMS`` from :doc:`bbcflib_daflims`
* ``Frontend`` from :doc:`bbcflib_frontend`
* ``mapseq`` from :doc:`bbcflib_mapseq`
* ``chipseq`` from :doc:`bbcflib_chipseq`
* ``gdv`` from :doc:`bbcflib_gdv`
"""

from ConfigParser import ConfigParser
from genrep import GenRep, Assembly
from email import EmailReport
from daflims import DAFLIMS
from frontend import Frontend, Job


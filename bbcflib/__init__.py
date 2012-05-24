"""
===========
Description
===========

bbcflib is a set of Python modules for accessing facilities used by
the Bioinformatics and Biostatistics Core Facility (BBCF) at the EPFL.
It provides modules for BBCF's

* *GenRep* (genome data repository),
* the LIMS of the DNA Array Facilities at UNIL and UNIGE,
* the standard frontend code for workflows written by Fabrice David,
* sending job completion emails.
* the *GDV* api for viewing genome tracks
* *mapseq* to map reads to reference genomes
* *chipseq* to run ChIP-seq analyses
* *rnaseq* to map reads to reference transcriptomes and compute statistics of differential expression
* *snp* to search for SNP with respect to a reference genome

All the functionality can be imported with::

    import bbcflib.module

where module is one of:

* ``GenRep`` and ``Assembly`` from  :doc:`bbcflib_genrep`
* ``EmailReport`` from :doc:`bbcflib_email`
* ``DAFLIMS`` from :doc:`bbcflib_daflims`
* ``Frontend`` from :doc:`bbcflib_frontend`
* ``mapseq`` from :doc:`bbcflib_mapseq`
* ``chipseq`` from :doc:`bbcflib_chipseq`
* ``rnaseq`` from :doc:`bbcflib_rnaseq`
* ``gdv`` from :doc:`bbcflib_gdv`
* ``snp`` from :doc:`bbcflib_snp`

============
Installation
============

bbcflib requires:

* Python >=2.6
* mailer >=0.6 (http://pypi.python.org/pypi/mailer)
* bein (http://bbcf.epfl.ch/bein/)
* numpy (http://numpy.scipy.org/)
* pysam (http://code.google.com/p/pysam/)

Latest source code is available from GitHub::

    http://github.com/bbcf/bbcflib

by clicking on "Downloads", or by cloning the git repository with::

    $ git clone https://github.com/bbcf/bbcflib.git


=======
License
=======

bbcflib is released under the GNU General Public License 3.0. A copy
of this license is in the LICENSE.txt file.
"""

b'This module needs Python 2.6 or later.'

__version__ = '2.0.0'

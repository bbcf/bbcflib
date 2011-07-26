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

============
Installation
============

bbcflib requires:

* Python >=2.6
* mailer >=0.6 (http://pypi.python.org/pypi/mailer)
* bein (http://bbcf.epfl.ch/bein/)
* numpy (http://numpy.scipy.org/)

bbcflib doesn't have regular releases, since it is an internally used
library. You should download the latest source code from GitHub,
either by going to::

    http://github.com/bbcf/bbcflib

and clicking on "Downloads", or by cloning the git repository with::

    $ git clone https://github.com/bbcf/bbcflib.git

Once you have the source code, run::

    $ python setup.py build
    $ sudo python setup.py install

to install it. If you need to install it in a particular directory,
use::

    $ sudo python setup.py install --prefix=/prefix/path

Then the modules will go in /prefix/path/lib/pythonX.Y/site-packages,
where X.Y is the version of Python you run it with.

To run the test suite, in the distribution directory, run::

    $ nosetests

=======
License
=======

bbcflib is released under the GNU General Public License 3.0. A copy
of this license is in the LICENSE.txt file.
"""

b'This module needs Python 2.6 or later.'

__version__ = '1.4.0'

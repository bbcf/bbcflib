bbcflib 1.3.0

Copyright 2010 Frederick Ross <bbcf at gmail dot com>,
2011 BBCF <webmaster dot bbcf at epfl dot ch>

bbcflib is a set of Python modules for accessing facilities used by
the Bioinformatics and Biostatistics Core Facility˘(BBCF) at the EPFL.

bbcflib is released under the GNU General Public License 3.0. A copy
of this license is in the LICENSE.txt file.

Installation
============

bbcflib requires:
* Python >=2.6
* mailer >=0.6 (http://pypi.python.org/pypi/mailer)
* bein (http://bbcf.epfl.ch/bein/)
* numpy (http://numpy.scipy.org/)

bbcflib doesn't have regular releases, since it is an internally used
library. You should download the latest source code from GitHub,
either by going to

    http://github.com/bbcf/bbcflib

and clicking on "Downloads", or by cloning the git repository with

    $ git clone https://github.com/bbcf/bbcflib.git

Once you have the source code, run

    $ python setup.py build
    $ sudo python setup.py install

to install it. If you need to install it in a particular directory,
use

    $ sudo python setup.py install --prefix=/prefix/path

Then the modules will go in /prefix/path/lib/pythonX.Y/site-packages,
where X.Y is the version of Python you run it with.

After installation process you need update your mime database by following command:
    $ sudo update-mime-database /usr/share/mime
if during installation you use --datadir options you need give to update-mime-database command the same path:
    $ python setup.py build
    $ sudo python setup.py install --datadir=/usr/local/share
    $ sudo update-mime-database /usr/local/share/mime

To run the test suite, in the distribution directory, run

    $ nosetests --with-doctest

Full documentation
==================

The full documentation can be found [on our website](http://bbcf.epfl.ch/bbcflib).

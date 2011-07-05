from distutils.core import setup

setup(
      name             = 'bbcflib',
      version          = '1.3.1',
      description      = 'Utility modules for deploying workflows in the BBCF',
      long_description = open('README.md').read(),
      license          = 'GNU General Public License 3.0',
      url              = 'http://bbcf.epfl.ch/bbcflib',
      author           = 'EPFL BBCF',
      author_email     = 'webmaster.bbcf@epfl.ch',
      install_requires = ['mailer', 'bein', 'rpy2', 'pysam', 'scipy', 'BeautifulSoup'],
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['bbcflib',
                          'bbcflib.tests',
                          'bbcflib.track',
                         ],
)

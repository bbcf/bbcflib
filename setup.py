from distutils.core import setup

setup(
      name             = 'bbcflib',
      version          = '1.3.0',
      description      = 'Utility modules for deploying workflows in the BBCF',
      long_description = open('README.md').read(),
      license          = 'GNU General Public License 3.0',
      url              = 'http://bbcf.epfl.ch/bbcflib',
      author           = 'EPFL BBCF, Fred Ross, Jacques Rougemont, Lucas Sinclair, Jonathan Mercier.',
      author_email     = 'webmaster.bbcf@epfl.ch',
      install_requires = ['mailer', 'bein', 'rpy2', 'BeautifulSoup'],
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['bbcflib',
                          'bbcflib.tests',
                          'bbcflib.track',
                         ],
)

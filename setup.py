from distutils.core import setup

setup(
        name            =   'bbcflib',
        version         =   '1.4.1',
        description     =   'Utility modules for deploying workflows in the BBCF',
        long_description=   open('README.md').read(),
        license         =   'GNU General Public License 3.0',
        url             =   'http://bbcf.epfl.ch/bbcflib',
        author          =   'EPFL BBCF, Jacques Rougemont, Fred Ross, Lucas Sinclair, Jonathan Mercier, Julien Delafontaine, Yohan Jarosz.',
        author_email    =   'webmaster.bbcf@epfl.ch',
        install_requires=   ['mailer', 'bein', 'rpy2', 'pysam', 'scipy', 'BeautifulSoup', 'cogent', 'mysql-python', 'sqlalchemy'],
        classifiers     =   ['Topic :: Scientific/Engineering :: Bio-Informatics'],
        packages        =   ['bbcflib',
                             'bbcflib.btrack',
                             'bbcflib.btrack.extras',
                             'bbcflib.btrack.formats',
                            ],
    )

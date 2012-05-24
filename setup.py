from distutils.core import setup

setup(
        name            =   'bbcflib',
        version         =   '2.0.0',
        description     =   'Utility modules for deploying workflows in the BBCF',
        long_description=   open('README.md').read(),
        license         =   'GNU General Public License 3.0',
        url             =   'http://bbcf.epfl.ch/bbcflib',
        author          =   'EPFL BBCF',
        author_email    =   'webmaster.bbcf@epfl.ch',
        install_requires=   ['mailer', 'bein', 'rpy2', 'pysam', 'scipy', 'numpy', 'sqlalchemy', 'sqlite3'],
        classifiers     =   ['Topic :: Scientific/Engineering :: Bio-Informatics'],
        packages        =   ['bbcflib',
                             'bbcflib.btrack',
                             'bbcflib.bFlatMajor',
                             'bbcflib.bFlatMajor.stream',
                             'bbcflib.bFlatMajor.numeric',
                             'bbcflib.bFlatMajor.figure' 
                            ],
    )

from distutils.core import setup, Command
from unittest       import TextTestRunner, TestLoader
from glob           import glob
from os.path        import splitext, basename, join as pjoin, walk
import os, pdb

class TestCommand(Command):
    user_options = [ ]

    def initialize_options(self):
        self._dir = os.getcwd()

    def finalize_options(self):
        pass

    def run(self):
        '''
        Finds all the tests modules in tests/, and runs them.
        '''
        testfiles = [ ]
        for t in glob(pjoin(self._dir, 'tests', '*.py')):
            if not t.endswith('__init__.py'):
                testfiles.append('.'.join(
                    ['tests', splitext(basename(t))[0]])
                )
        pdb.set_trace()
        tests = TestLoader().loadTestsFromNames(testfiles)
        t = TextTestRunner(verbosity = 1)
        t.run(tests)

class CleanCommand(Command):
    user_options = [ ]

    def initialize_options(self):
        self._clean_me = [ ]
        for root, dirs, files in os.walk('.'):
            for f in files:
                if f.endswith('.pyc'):
                    self._clean_me.append(pjoin(root, f))

    def finalize_options(self):
        pass

    def run(self):
        for clean_me in self._clean_me:
            try:
                os.unlink(clean_me)
            except:
                pass

setup(
        name            =   'bbcflib',
        version         =   '1.3.1',
        description     =   'Utility modules for deploying workflows in the BBCF',
        long_description=   open('README.md').read(),
        license         =   'GNU General Public License 3.0',
        url             =   'http://bbcf.epfl.ch/bbcflib',
        author          =   'EPFL BBCF, Fred Ross, Jacques Rougemont, Lucas Sinclair, Jonathan Mercier.',
        author_email    =   'webmaster.bbcf@epfl.ch',
        install_requires=   ['python-mailer', 'python-magic','bein', 'rpy2', 'pysam', 'scipy', 'BeautifulSoup'],
        classifiers     =   ['Topic :: Scientific/Engineering :: Bio-Informatics'],
        packages        =   ['bbcflib',
                             'bbcflib.track',
                             'bbcflib.track.extras',
                             'bbcflib.track.formats',
                             'bbcflib.track.magic'
                            ],
        package_dir     =   {'bbcflib.track.magic':'bbcflib/track/magic'},
        package_data    =   {'bbcflib.track.magic':['magic_data']},
        cmdclass        = { 'test': TestCommand, 'clean': CleanCommand }
    )

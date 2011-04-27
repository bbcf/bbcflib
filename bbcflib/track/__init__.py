"""
=========================
Subpackage: bbcflib.track
=========================

Provides easy read/write access to genomic tracks in a fashion that is independant from the underlying format.
Currently the following formats are implemented as read/write:

* Bio SQLite (http://bbcf.epfl.ch/twiki/bin/view/BBCF/SqLite)
* BED        (http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
* WIG        (http://genome.ucsc.edu/goldenPath/help/wiggle.html)

More formats can be added easily.
 
To get access to the information contained inside already existing tracks, you would do the following whatever the format of the track is::

    from bbcflib.track import Track
    with rpgenes as Track('tracks/rp_genes.sql'):
        rpgenes.read('chr3')

Optionally you can supply a name for every track you load, to help you keep track of your tracks::

    from bbcflib.track import Track
    with ribigenes as Track('tracks/ribi_genes.sql', name='Ribosome genesis from SGD'):
        ribigenes.read('chr7')

If your track is in a format that is missing chromosome information, you will need to supply an extra chromosome file::

    from bbcflib.track import Track
    with saccer as Track('tracks/yeast_genes.bed', chrfile='tracks/chrs/yeast.chr'):
        saccer.read('chr4')

For instance, the cumulative base coverage of features on chromosome two can be calculated like this::

    from bbcflib.track import Track
    with saccer as Track('tracks/yeast_genes.sql'):
        base_coverage = sum([f[1] - f[0] for f in saccer.read('chr2')])

To create a new track and then write to it, you would do the following::

    from bbcflib.track import new
    with mypeaks as new('tracks/inverted.sql', 'sql', name='Inverted strands'):
        mypeaks.write('chr1')

For instance, to make a new track from an old one, and invert the strand of every feature::

    from bbcflib.track import Track, new
    def invert_strands(data):
        for feature in data:
            yield (feature[0], feature[1], feature[2], feature[3], feature[4] == '+' and '-' or '+')
    with a as Track('tracks/orig.sql', name='Normal strands'):
        with b as new('tracks/inverted.sql', name='Inverted strands'):
            for chrom in a:
                b.write(chrom, invert_strands(a.read(chrom))) 
"""

import os, sys
from .. import common as com

#-----------------------------------------------------------------------------#   
def _determine_format(path):
    '''Try to guess the format of a track given its path. Returns a three letter extension'''
    return 'sql'

def _import_implementation(format):
    '''Try to import the implementation of a given format'''
    try:
        if not hasattr(sys.modules[__package__], format):
            __import__(__package__ + '.' + format)
        return sys.modules[__package__ + '.' + format]
    except ImportError, AttributeError:
        raise Exception("The format '" + format + "' is not supported at the moment")

###########################################################################   
class Track(object):
    '''The class used to access genomic data. It can load a track from disk, whatever the format is.

            * *path* is the path to track file.
            * *name* is an option string specifing the name of the track to load.
            * *chrfile* is the path to chromosome file. This is specified when the underlying format is missing chromosome length information. 
        
        Examples::

            with rpgenes as Track('tracks/rp_genes.sql'):
                rpgenes.read('chr3')

        Once a track is loaded you have access to the following attributes:

           * *path* is the file system path to the underlying file.
           * *type* is either ``qualitative`` or ``quantiative``.
           * *name* is given upon creation of the track.
           * *format* is either ``sql``, ``bed`` or ``wig``
           * *all_chrs* is a list of all available chromosome. For instance:
                ``['chr1, 'chr2', 'chr3']``
           * *meta_chr* is a dictionary of chromosome meta data (information like length, etc). For instance: 
                ``[{'name': 'chr1', 'length': 197195432}, {'name': 'chr2', 'length': 129993255}]``
           * *meta_track* is a dictionary of meta data associated to the track (information like the source, etc). For instance:
                 ``{'datatype': 'quantiative', 'source': 'SGD'}``
           * *chr_file* is the path to a chromosome file if one was given.
    '''

    def __new__(cls, path, name=None, chrfile=None):
        '''Internal factory-like method that is called before creating a new instance of Track.
           This function determines the format of the file that is provided and returns an
           instance of the approriate child class.'''
        if cls is Track:
            format         = _determine_format(path)
            implementation = _import_implementation(format)
            instance       = implementation.GenomicFormat(path, name, chrfile)
        else:
            instance = super(Track, cls).__new__(cls)
        return instance

    def __iter__(self):
        ''' Called when trying to iterate the class'''
        return iter(self.all_chrs)
    
    def __enter__(self):
        ''' Called when entering a <with> statement'''
        pass

    def __exit__(self, type, value, traceback):
        '''Called when exiting a <with> statement'''
        pass

    def __init__(self, path, name=None, chrfile=None):
        # Set attributes #
        self.path     = path
        self.name     = name
        self.chr_file = chrfile
        # Check existance #
        if not os.path.exists(path):
            raise Exception("The file '" + path + "' cannot be found")
        if os.path.isdir(path):
            raise Exception("The location '" + path + "' is a directory")
        # Call child function #
        self.load()
        # Sort chromosomes #
        self.all_chrs.sort(key=com.natural_sort)
        # Test variables #
        if self.type not in ['quantitative', 'qualitative']:
            raise Exception("The type of the track is invalid: " + type + ".")

    #-----------------------------------------------------------------------------#   
    def read(self, selection=None, fields=None):
        '''Read data from the genomic file.

        * *selection* can be the name of a chromosome, in which case all the data on that chromosome will be returned. It can also be a dictionary specifying a region in which case only features contained in that region will be returned. When left empty, the data from all chromosome is returned 

        * *fields* is a list of fields which will influence the order in which the information is returned. The default for quantitative tracks is ['start', 'stop', 'name', 'strand', 'score'] and ['start', 'stop', 'score'] for quantitative tracks.

        Examples::

            t.read()
            t.read('chr2')
            t.read({chr:'chr1', 'start':10000, 'stop':15000})
            t.read('chr3', ['name', strand'])
            t.read({chr:'chr5', 'start':0, 'stop':200}, ['start', 'score'])

        ``read`` returns a generator object yielding tuples.
        '''
        pass

    def write(self, chrom, data, fields):
        '''Writes data to a genomic file.

        * *chrom* is the name of the chromosome on which one wants to write. For instance, if one is using the BED format this will become the first column, while if one is using the SQL format this will become the name of the table to be created.

        * *data* must be an iterable object that yields features. As an example, the ``read`` function of this class produces such objects. 

        * *fields* is a list of fields which will influence the number of columns for every feature in the file. This should not change within the same track.

        Examples::

            t.write('chr1', [('10','20','A','0','+'), ('40','50','B','0','-')])
            def example_generator():
                for i in xrange(5):
                    yield ('10','20','X',i,'+') 
            t.write('chr2', example_generator())

        ``write`` returns nothing.
        '''
        pass

    def convert(self, path, format='sql'):
        '''Convert a track to a given format.
       
           * *path* is the file path where the new track will be created

           * *format* is the format into which the track will be converted.

           Examples::
            
               with t as Track('tracks/rp_genes.bed'):
                   t.convert('tracks/rp_genes.sql', 'sql')
               with t as Track('tracks/ribi_genes.sql'):
                   t.convert('tracks/rp_genes.bed', 'bed')
        
            ``convert`` returns nothing.
        '''
        pass

    #-----------------------------------------------------------------------------#
    @property
    def meta_chr(self): 
        print 'called getter'
        return self._x

    @meta_chr.setter
    def meta_chr(self, value): 
        print 'called setter'
        self._x = value

    @property
    def meta_track(self): 
        print 'called getter'
        return self._x

    @meta_track.setter
    def meta_track(self, value): 
        print 'called setter'
        self._x = value

    #-----------------------------------------------------------------------------#
    def load(self):
        raise NotImplementedError

    def new(self):
        raise NotImplementedError

###########################################################################   
def new(path, format='sql', type='qualitative', name='Unamed'):
    '''Create a new empty track in preparation for writing to it.

        * *path* is the file path where the new track will be created

        * *format* is the format into which the *track* will be converted.

        Examples::
            
            t = new('tmp/track.sql'):
            t = new('tracks/peaks.sql', 'sql', name='High affinity peaks'):
            t = new('tracks/scores.sql', 'sql', type='quantitative', name='Signal along the genome'):
        
        ``new`` returns a Track instance.

    '''
    pass

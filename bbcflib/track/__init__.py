"""
=========================
Subpackage: bbcflib.track
=========================

Provides easy read/write access to genomic tracks in a fashion that is independent from the underlying format.
Currently the following formats are implemented as read/write:

* Bio SQLite (http://bbcf.epfl.ch/twiki/bin/view/BBCF/SqLite)
* BED        (http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
* WIG        (http://genome.ucsc.edu/goldenPath/help/wiggle.html)

More formats can be added easily.
 
To get access to the information contained inside already existing tracks, you would do the following whatever the format of the track is::

    from bbcflib.track import Track
    with Track('tracks/rp_genes.sql') as rpgenes:
        data = rpgenes.read('chr3')

Optionally you can supply a name for every track you load, to help you keep track of your tracks::

    from bbcflib.track import Track
    with Track('tracks/ribi_genes.sql', name='Ribosome genesis from SGD') as ribigenes:
        data = ribigenes.read('chr7')

If your track is in a format that is missing chromosome information, you will need to supply an extra chromosome file::

    from bbcflib.track import Track
    with Track('tracks/yeast_genes.bed', chrfile='tracks/chrs/yeast.chr') as saccer:
        data = saccer.read('chr4')

For instance, the cumulative base coverage of features on chromosome two can be calculated like this::

    from bbcflib.track import Track
    with Track('tracks/yeast_genes.sql') as saccer:
        base_coverage = sum([f[1] - f[0] for f in saccer.read('chr2')])

To create a new track and then write to it, you would do the following::

    from bbcflib.track import new
    with new('tracks/rap1_peaks.sql', 'sql', name='Rap1 Peaks') as mypeaks:
        mypeaks.write('chr1', [('10','20','A','0','+')])

For instance, to make a new track from an old one, and invert the strand of every feature::

    from bbcflib.track import Track, new
    def invert_strands(data):
        for feature in data:
            yield (feature[0], feature[1], feature[2], feature[3], feature[4] == '+' and '-' or '+')
    with Track('tracks/orig.sql', name='Normal strands') as a:
        with new('tracks/inverted.sql', name='Inverted strands') as b:
            for chrom in a:
                b.write(chrom, invert_strands(a.read(chrom)))

To convert a track from a slow format (e.g. BED) to a database format (e.g. SQL) you first create a track object and call the convert method on it::

    from bbcflib.track import Track
    with Track('tracks/rp_genes.bed') as rpgenes:
        rpgenes.convert('tracks/rp_genes.sql', 'sql')
"""

import os, sys
from .. import common as com

#-----------------------------------------------------------------------------#   
def _determine_format(path):
    '''Try to guess the format of a track given its path. Returns a three letter extension'''
    # Get extension #
    extension = os.path.splitext(path)[1][1:]
    # If no extension then try magic #
    known_format_extensions = {
        'SQLite 3.x database':                       'sql',
        'SQLite database (Version 3)':               'sql',
        'Hierarchical Data Format (version 5) data': 'hdf5',
    }
    if not extension: 
        try:
            import magic
        except ImportError:
            raise Exception("The format of the track '" + path + "' cannot be determined")
        # Let the user customize magic #
        if os.path.exists('magic'):
            m = magic.Magic('magic')
        else:
            m = magic.Magic()
        # Check the result of magic #
        type = m.from_file(path)
        try:
            extension = known_format_extensions[type]
        except KeyError as err:
            raise Exception("The format of the track '" + path + "' resolves to " + type + " which is not supported at the moment.")
    # Synonyms #
    if extension == 'db': extension = 'sql'
    # Return the format #
    return extension

def _import_implementation(format):
    '''Try to import the implementation of a given format'''
    try:
        if not hasattr(sys.modules[__package__], format):
            __import__(__package__ + '.' + format)
        return sys.modules[__package__ + '.' + format]
    except ImportError, AttributeError:
        raise Exception("The format '" + format + "' is not supported at the moment")

def join_read_queries(track, selections, fields):
    ''' Join read results when selection is a list'''
    def _add_chromsome(sel, data):
        if type(sel) == str: chrom = (sel,)
        else:                chrom = (sel['chr'],)
        for f in data: yield chrom + f
    for sel in selections:
        for f in _add_chromsome(sel, track.read(sel, fields)): yield f

###########################################################################   
class Track(object):
    '''The class used to access genomic data. It can load a track from disk, whatever the format is.

            * *path* is the path to track file.
            * *name* is an option string specifying the name of the track to load.
            * *chrfile* is the path to chromosome file. This is specified when the underlying format is missing chromosome length information. 
        
        Examples::

            with Track('tracks/rp_genes.sql') as rpgenes:
                data = rpgenes.read('chr3')

        Once a track is loaded you have access to the following attributes:

           * *path* is the file system path to the underlying file.
           * *type* is either ``qualitative`` or ``quantitative``.
           * *name* is given upon creation of the track.
           * *format* is either ``sql``, ``bed`` or ``wig``
           * *all_chrs* is a list of all available chromosome. For instance:
                ``['chr1, 'chr2', 'chr3']``
           * *meta_chr* is a dictionary of chromosome meta data (information like length, etc). For instance: 
                ``[{'name': 'chr1', 'length': 197195432}, {'name': 'chr2', 'length': 129993255}]``
           * *meta_track* is a dictionary of meta data associated to the track (information like the source, etc). For instance:
                 ``{'datatype': 'quantitative', 'source': 'SGD'}``
           * *chr_file* is the path to a chromosome file if one was given.
    '''

    qualitative_fields = ['start', 'end', 'name', 'score', 'strand'] 
    quantiative_fields = ['start', 'end', 'score'] 
    field_types = {
        'start':        'integer',
        'end':          'integer',
        'score':        'real',
        'strand':       'integer',
        'name':         'text',
        'thick_start':  'integer',
        'thick_end':    'integer',
        'item_rgb':     'text',
        'block_count':  'integer',
        'block_sizes':  'text',
        'block_starts': 'text',
    }
    
    #-----------------------------------------------------------------------------#   
    def __new__(cls, path, name=None, chrfile=None):
        '''Internal factory-like method that is called before creating a new instance of Track.
           This function determines the format of the file that is provided and returns an
           instance of the appropriate child class.'''
        if cls is Track:
            format         = _determine_format(path)
            implementation = _import_implementation(format)
            instance       = implementation.GenomicFormat(path, name, chrfile)
        else:
            instance = super(Track, cls).__new__(cls)
        return instance

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

    def __iter__(self):
        ''' Called when trying to iterate the class'''
        return iter(self.all_chrs)
    
    def __enter__(self):
        ''' Called when entering a <with> statement'''
        return self

    def __exit__(self, type, value, traceback):
        '''Called when exiting a <with> statement'''
        self.unload(type, value, traceback)

    #-----------------------------------------------------------------------------#   
    def read(self, selection=None, fields=None):
        '''Read data from the genomic file.

        * *selection* can be the name of a chromosome, in which case all the data on that chromosome will be returned. It can also be a dictionary specifying a region in which case only features contained in that region will be returned. To combine multiple selections you can specify a list including chromosome names and region dictionaries. As exepected, if such is the case, the joined data from those selections will be returned with an added 'chr' field in front since the results may span several chromosomes. When *selection* is left empty, the data from all chromosome is returned.

        * *fields* is a list of fields which will influence the order in which the information is returned. The default for quantitative tracks is ``['start', 'end', 'name', 'score', 'strand']`` and ``['start', 'end', 'score']`` for quantitative tracks.

        Adding the parameter ``'inclusion':'strict'`` to a region dictionary will return only features exactly contained inside the interval instead of features simply included in the interval.

        Examples::

            data = t.read()
            data = t.read('chr2')
            data = t.read(['chr1','chr2','chr3'])
            data = t.read({'chr':'chr1', 'start':10000, 'stop':15000})
            data = t.read({'chr':'chr1', 'start':10000, 'stop':15000, 'inclusion':'strict'})
            data = t.read('chr3', ['name', 'strand'])
            data = t.read({'chr':'chr5', 'start':0, 'stop':200}, ['strand', 'start', 'score'])

        ``read`` returns a generator object yielding tuples.
        '''
        raise NotImplementedError

    def write(self, chrom, data, fields):
        '''Writes data to a genomic file.

        * *chrom* is the name of the chromosome on which one wants to write. For instance, if one is using the BED format this will become the first column, while if one is using the SQL format this will become the name of the table to be created.

        * *data* must be an iterable object that yields tuples of the correct length. As an example, the ``read`` function of this class produces such objects. 

        * *fields* is a list of fields which will influence the number of columns for every feature in the file and hence the length of the tuples to be generated. The value of *fields* should not change within the same track.

        Examples::

            t.write('chr1', [('10','20','A','0','+'), ('40','50','B','0','-')])
            def example_generator():
                for i in xrange(5):
                    yield ('10','20','X',i,'+') 
            t.write('chr2', example_generator())

        ``write`` returns nothing.
        '''
        raise NotImplementedError

    def convert(self, path, format='sql', type=None, mean_scores=False):
        '''Convert a track to a given format.
       
           * *path* is the file path where the new track will be created

           * *format* is the format into which the track will be converted.
           
           * *type* is the datatype into which the track will be converted. It is either ``qualitative`` or ``quantitative``. Certain formats have default datatypes. For instance a BED file will become a qualitative SQL track, while a WIG file will become a quantitative SQL track. However, this can be overwritten. In the case of a quantitative BED track, the information contained in the names and strand columns will be discarded. In the case of qualitative WIG track, all features will be missing names.

           * *mean_scores* is a boolean value that defaults to False. When converting a track to a non-natural datatype, in particular a BED file to a quantitative track, you might want to ignore the errors produced by overlapping features and resolve the conflicts by meaning the score in that region. This is done by setting *mean_scores* to True.

           Examples::
            
               with Track('tracks/rp_genes.bed') as t:
                   t.convert('tracks/rp_genes.sql', 'sql')
               with Track('tracks/ribi_genes.sql') as t:
                   t.convert('tracks/rp_genes.bed', 'bed')
               with Track('tracks/yeast_genes.bed') as t:
                   t.convert('tracks/yeast_genes.sql', 'sql', 'quantitative')
               with Track('tracks/rap1_peaks.bed') as t:
                   t.convert('tracks/rap1_peaks.sql', 'sql', 'quantitative', True)
        
           ``convert`` returns nothing.
        '''
        raise NotImplementedError

    #-----------------------------------------------------------------------------#
    @property
    def type(self):
        raise NotImplementedError

    @property
    def fields(self):
        raise NotImplementedError
    
    @property
    def meta_chr(self): 
        raise NotImplementedError

    @meta_chr.setter
    def meta_chr(self, value): 
        raise NotImplementedError

    @property
    def meta_track(self): 
        raise NotImplementedError

    @meta_track.setter
    def meta_track(self, value): 
        raise NotImplementedError

    #-----------------------------------------------------------------------------#
    def load(self):
        pass
   
    def unload(self):
        pass
 
    @property
    def default_fields(self):
        return getattr(Track, self.type + '_fields')

###########################################################################   
def new(path, format='sql', type='qualitative', name='Unnamed'):
    '''Create a new empty track in preparation for writing to it.

        * *path* is the file path where the new track will be created

        * *format* is the format into which the *track* will be created.
        
        * *name* is an optional name for the track.

        Examples::
            
            t = new('tmp/track.sql')
            t = new('tracks/peaks.sql', 'sql', name='High affinity peaks')
            t = new('tracks/scores.sql', 'sql', type='quantitative', name='Signal along the genome')
        
        ``new`` returns a Track instance.
    '''
    if os.path.exists(path):
        raise Exception("The location '" + path + "' is already taken")
    implementation = _import_implementation(format)
    implementation.create(path, type, name)
    return Track(path, name=name)

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#

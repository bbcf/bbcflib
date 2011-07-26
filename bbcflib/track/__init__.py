"""
Provides easy read/write access to genomic tracks in a fashion that is independent from the underlying format.
Currently the following formats are implemented as read/write:

* Bio SQLite (http://bbcf.epfl.ch/twiki/bin/view/BBCF/SqLite)
* BED        (http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
* WIG        (http://genome.ucsc.edu/goldenPath/help/wiggle.html)
* bedGraph   (http://genome.ucsc.edu/goldenPath/help/bedgraph.html)
* bigWig     (http://genome.ucsc.edu/goldenPath/help/bigWig.html)

More formats can be added easily.

To get access to the information contained inside already existing tracks, you would do the following whatever the format of the track is::

    from bbcflib import track
    with track.load('tracks/rp_genes.sql') as rpgenes:
        data = rpgenes.read('chr3')

Optionally you can supply a name for every track you load, to help you keep track of your tracks::

    from bbcflib import track
    with track.load('tracks/ribi_genes.sql', name='Ribosome genesis from SGD') as ribigenes:
        data = ribigenes.read('chr7')

If your track is in a format that is missing chromosome information (such as the length of every chromosome), you can supply an assembly name or a chromosome file::

    from bbcflib import track
    with track.load('tracks/yeast_genes.bed', chrmeta='sacCer2') as saccer:
        data = saccer.read('chr4')
    with track.load('tracks/yeast_genes.bed', chrmeta='tracks/chrs/yeast.chr') as saccer:
        data = saccer.read('chr4')

For instance, the cumulative base coverage of features on chromosome two can be calculated like this::

    from bbcflib import track
    with track.load('tracks/yeast_genes.sql') as saccer:
        base_coverage = sum([f[1] - f[0] for f in saccer.read('chr2')])

To create a new track and then write to it, you would do the following::

    from bbcflib.track import new
    with new('tracks/rap1_peaks.sql', 'sql', name='Rap1 Peaks') as mypeaks:
        mypeaks.write('chr1', [(10, 20, 'A', 0.0, 1)])

For instance, to make a new track from an old one, and invert the strand of every feature::

    from bbcflib import track, new
    def invert_strands(data):
        for feature in data:
            yield (feature[0], feature[1], feature[2], feature[3], feature[4] == 1 and -1 or 1)
    with track.load('tracks/orig.sql', name='Normal strands') as a:
        with new('tracks/inverted.sql', name='Inverted strands') as b:
            for chrom in a:
                b.write(chrom, invert_strands(a.read(chrom)))

To convert a track from a format (e.g. BED) to an other format (e.g. SQL) you first load the track and call the export method on it::

    from bbcflib import track
    with track.load('tracks/rp_genes.bed') as rpgenes:
        rpgenes.export('tracks/rp_genes.sql', 'sql')

To set the chromosome metadata or the track metadata you simply asign to that attribute::

    from bbcflib import track
    with track.load('tracks/scores.sql') as t:
        t.chrmeta    = ``{'chr1': {'length': 197195432}, 'chr2': {'length': 129993255}}``
        t.attributes = {'datatype': 'quantitative', 'source': 'UCSC'}

It is important to note that the general numbering convention of features on a chromosome varies depending on the source of the data. For instance, UCSC and Ensembl differ in this point such that an interval labeled `(start=4,end=8)` will span four base pairs according to UCSC but will span five base pairs according to Ensembl. The representation that the this packages sticks to is explained `here <http://bbcf.epfl.ch/twiki/bin/view/BBCF/NumberingConvention>`_.
"""

__all__ = ['load', 'new']

# Built-in modules #
import os

###########################################################################
def load(path, format=None, name=None, chrmeta=None, datatype=None, readonly=False):
    '''Loads a track from disk, whatever the format is.

            * *path* is the path to track file to load.
            * *format* is an optional string specifying the format of the track to load when it cannot be guessed from the file extension.
            * *name* is an optional string specifying the name of the track to load. This can help with keep track of multiple different tracks.
            * *chrmeta* is an optional parameter useful in cases where the track to be loaded is missing extra meta data such chromosome length information. One can include this at the load phase by specifying one of three things: *chrmeta* can be the name of an assembly such as ``sacCer2``. *chrmeta* can be a dictionary where each key is chromosome names that specifies the ``length`` as integers. *chrmeta* can be the path to chromosome file. The chromosome file is structured as tab-separated text file containing two columns: the first specifies a chromosomes name and the second its length as an integer.
            * *datatype* is an optional variable that can take the value of either ``qualitative`` or ``quantitative``. It is only useful when loading a track that is ambiguous towards its datatype, as can be certain text files. For instance, In the case of WIG track becoming qualitative, all features will be missing names, but overlapping features will suddenly be authorized.
            * *readonly* is an optional boolean variable that defaults to ``False``. When set to ``True``, any operation attempting to write to the track will silently be ignored.

        Examples::

            from bbcflib import track
            with track.load('tracks/rp_genes.sql') as rpgenes:
                data = rpgenes.read()
            with track.load('tracks/yeast', 'sql', 'S. cer. genes') as yeast:
                data = yeast.read()
            with track.load('tracks/peaks.bed', 'bed', chrmeta='hg19') as peaks:
                data = peaks.read()
            with track.load('tracks/scores.wig', 'wig', chrmeta='tracks/cser.chr', datatype='qualitative') as scores:
                data = scores.read()
            with track.load('tracks/repeats.sql', readonly=True) as rpgenes:
                data = rpgenes.read()

        ``load`` returns a Track instance.
    '''
    return Track(path, format, name, chrmeta, datatype, readonly)

def new(path, format=None, name=None, chrmeta=None, datatype=None):
    '''Creates a new empty track in preparation for writing to it.

            * *path* is the path to track file to create.
            * *format* is an optional string specifying the format of the track to load when it cannot be guessed from the file extension.
            * *name* is an optional string specifying the name of the track to load. This can help with keep track of multiple different tracks.
            * *chrmeta* is an optional parameter useful in cases where the track to be loaded is missing extra meta data such chromosome length information. One can include this at the creation phase by specifying one of three things: *chrmeta* can be the name of an assembly such as ``sacCer2``. *chrmeta* can be a dictionary where each key is chromosome names that specifies the ``length`` as integers. *chrmeta* can be the path to chromosome file. The chromosome file is structured as tab-separated text file containing two columns: the first specifies a chromosomes name and the second its length as an integer.
            * *datatype* is an optional variable that can take the value of either ``qualitative`` or ``quantitative``. It is only useful when creating a track can can support both types, such as the BioSQLite format.

        Examples::

            from bbcflib import track
            with track.new('tmp/track.sql') as t:
                t.write('chr1', [(10, 20, 'Gene A', 0.0, 1)])
            with track.new('tracks/peaks.sql', 'sql', name='High affinity peaks') as t:
                t.write('chr5', [(500, 1200, 'Peak1', 11.3, 0)])
            with track.new('tracks/scores.sql', 'sql', chrmeta='sacCer2' datatype='quantitative',) as t:
                t.write('chr1', [(10, 20, 500.0)])

        ``new`` returns a Track instance.
    '''
    if os.path.exists(path): raise Exception("The location '" + path + "' is already taken")
    if not format: format = os.path.splitext(path)[1][1:]
    cls = import_implementation(format).TrackFormat
    cls.create(path)
    return Track(path, format, name, chrmeta, datatype, readonly=False, empty=True)

###########################################################################
class Track(object):
    '''Once a track is loaded you have access to the following attributes:

           * *path* is the file system path to the underlying file.
           * *datatype* is either ``qualitative`` or ``quantitative``.
           * *name* is given upon creation of the track.
           * *format* is always ``sql`` as, when loading a ``bed`` track for instance, a transparent underlying Bio SQLite track is sliently created.
           * *all_chrs* is a list of all available chromosome. For instance:
                ``['chr1, 'chr2', 'chr3']``
           * *chrmeta* is a  dictionary of meta data associated to each chromosome (information like length, etc). For instance:
                ``{'chr1': {'length': 197195432}, 'chr2': {'length': 129993255}}``
           * *attributes* is a dictionary of meta data associated to the track (information like the source, etc). For instance:
                 ``{'datatype': 'quantitative', 'source': 'SGD'}``

        The track object itself is iterable and will yield the name of all chromosomes.

        Examples::

            from bbcflib import track
            with track.load('tracks/rp_genes.sql') as rpgenes:
                for chrom in rpgenes: print chrom
                if 'chrY' in rpgenes: print 'Male'
                if len(rpgenes) != 23: print 'Aneuploidy'
    '''

    #-----------------------------------------------------------------------------#
    def read(self, selection=None, fields=None, order='start,end', cursor=False):
        '''Reads data from the genomic file.

        * *selection* can be several things.

        *selection* can be the name of a chromosome, in which case all the data on that chromosome will be returned.

        *selection* can also be a dictionary specifying: regions, score intervals or strands. Indeed, you can specify either region in which case only features contained in that region will be returned or a dictionary specifying a score interval in which case only features contained in that score boundaries will be returned. You can also specify a strand. The dictionary can specify one or several of these arguemts. See examples for more details

        Adding the parameter ``'inclusion':'strict'`` to a region dictionary will return only features exactly contained inside the interval instead of features simply included in the interval.

       To combine multiple selections you can specify a list including chromosome names and region dictionaries. As expected, if such is the case, the joined data from those selections will be returned with an added 'chr' field in front since the results may span several chromosomes.

        When *selection* is left empty, the data from all chromosome is returned.

        * *fields* is a list of fields which will influence the length of the tuples returned and the way in which the information is returned. The default for quantitative tracks is ``['start', 'end', 'name', 'score', 'strand']`` and ``['start', 'end', 'score']`` for quantitative tracks.

        * *order* is a sublist of *fields* which will influence the order in which the tuples are yielded. By default results are sorted by ``start`` and, secondly, by ``end``.

        * *cursor* is a boolean which should be set true if you are performing several operations on the same track at the same time. This is the case, for instance when you are chaining a read operation to a write operation.

        Examples::

            from bbcflib import track
            with track.load('tracks/example.sql') as t:
                data = t.read()
                data = t.read('chr2')
                data = t.read('chr3', ['name', 'strand'])
                data = t.read(['chr1','chr2','chr3'])
                data = t.read({'chr':'chr1', 'start':100})
                data = t.read({'chr':'chr1', 'start':10000, 'end':15000})
                data = t.read({'chr':'chr1', 'start':10000, 'end':15000, 'inclusion':'strict'})
                data = t.read({'chr':'chr1', 'strand':1})
                data = t.read({'chr':'chr1', 'score':(10,100)})
                data = t.read({'chr':'chr1', 'start':10000, 'end':15000 'strand':-1 'score':(10,100)})
                data = t.read({'chr':'chr5', 'start':0, 'end':200}, ['strand', 'start', 'score'])
            # Duplicate a chromosome
            with track.load('tracks/copychrs.sql') as t:
                t.write('chrY', t.read('chrX', cursor=True))

        ``read`` returns a generator object yielding tuples.
        '''
        raise NotImplementedError

    def write(self, chrom, data, fields=None):
        '''Writes data to a genomic file.

        * *chrom* is the name of the chromosome on which one wants to write. For instance, if one is using the BED format this will become the first column, while if one is using the SQL format this will become the name of the table to be created.

        * *data* must be an iterable object that yields tuples of the correct length. As an example, the ``read`` function of this class produces such objects.

        * *fields* is a list of fields which will influence the number of columns for every feature in the file and hence the length of the tuples to be generated. The value of *fields* should not change within the same track.

        Examples::

            from bbcflib import track
            with track.load('tracks/example.sql') as t:
                t.write('chr1', [(10, 20, 'A', 0.0, 1), (40, 50, 'B', 0.0, -1)])
                def example_generator():
                    for i in xrange(5):
                        yield (10, 20, 'X', i, 1)
                t.write('chr2', example_generator())

        ``write`` returns nothing.
        '''
        raise NotImplementedError

    def remove(self, chrom=None):
        '''Removes data from a given chromosome.

        * *chrom* is the name of the chromosome that one wishes to delete or a list of chromosomes to delete.

        Called with no arguments, will remove every chromosome.

        Examples::

            from bbcflib import track
            with track.load('tracks/example.sql') as t:
                t.remove('chr1')
            with track.load('tracks/example.sql') as t:
                t.remove(['chr1', 'chr2', 'chr3'])
            with track.load('tracks/example.sql') as t:
                t.remove()

        ``remove`` returns nothing.
        '''
        raise NotImplementedError

    def count(self, selection=None):
        '''Counts the number of features or entries in a given selection.

        * *selection* is the name of a chromosome, a list of chromosomes, a particular span or a list of spans. In other words, a value similar to the *selection* parameter of the *read* method.

        Called with no arguments, will count every feature in a track.

        Examples::

            from bbcflib import track
            with track.load('tracks/example.sql') as t:
                num = t.count('chr1')
            with track.load('tracks/example.sql') as t:
                num = t.count(['chr1','chr2','chr3'])
            with track.load('tracks/example.sql') as t:
                num = t.count({'chr':'chr1', 'start':10000, 'end':15000})

        ``count`` returns an integer.
        '''
        raise NotImplementedError

    def export(self, path, format=None):
        '''Exports a track to a given location and format.

           * *path* is the file path where the new track will be created

           * *format* is the format into which the track will be converted.

           Examples::

               from bbcflib import track
               with track.load('tracks/rp_genes.sql') as t:
                   t.export('tracks/rp_genes.bed')
               with track.load('tracks/rp_genes.bed') as t:
                   t.write('chr1', generate_data())

           ``export`` returns nothing but a new file is created at the specified *path* while the current track object is left untouched.
        '''
        if os.path.exists(path): raise Exception("The location '" + path + "' is already taken")
        if not format: format = os.path.splitext(path)[1][1:]
        with new(path, format) as t:
            t.chrmeta    = self.chrmeta
            t.attributes = self.attributes
            fields = self.fields
            for chrom in self.all_chrs: t.write(chrom, self.read(chrom), fields)

    def convert(self, path, format=None):
        '''Converts a track to a given format dynamically.

           * *path* is the file path where the new track will be created

           * *format* is the format into which the track will be converted.

           Examples::

               from bbcflib import track
               with track.load('tracks/rp_genes.bed') as t:
                   t.convert('tracks/rp_genes.sql', 'sql')
                   t.write('chr1', generate_data())

           ``convert`` returns nothing but the track object is mutated and changes format dynamically. You can thus continue using your track after the conversion.
        '''
        if os.path.exists(path): raise Exception("The location '" + path + "' is already taken")
        if not format: format = os.path.splitext(path)[1][1:]
        if format == self.format:
            raise Exception("The track '" + path + "' cannot be converted to the " + format + " format because it is already in that format.")
        self.mutate_source(path, format)

    def ucsc_to_ensembl(self):
        '''Converts all entries of a track from the UCSC standard to the Ensembl standard effectively adding one to every start position.

           Examples::

               from bbcflib import track
               with track.load('tracks/example.sql') as t:
                   t.ucsc_to_ensembl()

           ``ucsc_to_ensembl`` returns nothing.
        '''
        raise NotImplementedError

    def ensembl_to_ucsc(self):
        '''Converts all entries of a track from the Ensembl standard to the UCSC standard effectively subtracting one from every start position.

           Examples::

               from bbcflib import track
               with track.load('tracks/rp_genes.bed') as t:
                   t.ensembl_to_ucsc()

           ``ensembl_to_ucsc`` returns nothing.
        '''
        raise NotImplementedError

    #-----------------------------------------------------------------------------#
    def mutate_source(self, path, format):
        '''Change a <Track> instance to a given format dynamically'''
        self.format = format
        cls = import_implementation(format).TrackFormat
        cls.create(path)
        cls.mutate_destination(self, path, format)

    @classmethod
    def mutate_destination(cls, self, path, format):
        '''Until other formats are added, this should be never called.'''
        raise NotImplementedError

    #-----------------------------------------------------------------------------#
    qualitative_fields  = ['start', 'end', 'name', 'score', 'strand']
    quantitative_fields = ['start', 'end', 'score']
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
    def __new__(cls, path, format=None, name=None, chrmeta=None, datatype=None, readonly=False, empty=False):
        '''Internal factory-like method that is called before creating a new instance of Track.
           This function determines the format of the file that is provided and returns an
           instance of the appropriate child class.'''
        if cls is Track:
            # Check existence #
            if not os.path.exists(path): raise Exception("The file '" + path + "' cannot be found")
            elif os.path.isdir(path): raise Exception("The location '" + path + "' is a directory")
            # Get child class #
            if not format: format = determine_format(path)
            implementation = import_implementation(format)
            instance       = super(Track, cls).__new__(implementation.TrackFormat)
            # Set attributes #
            instance.format   = format
            instance.modified = False
            instance.readonly = readonly
            instance._chrmeta    = ChromMetaData()
            instance._attributes = TrackMetaData()
        else:
            instance = super(Track, cls).__new__(cls)
        return instance

    def __init__(self, path, format=None, name=None, chrmeta=None, datatype=None, readonly=False, empty=False):
       pass

    def __enter__(self):
        return self

    def __exit__(self, errtype, value, traceback):
        self.unload(errtype, value, traceback)

    def __iter__(self):
        return iter(self.all_chrs)

    def __len__(self):
        return len(self.all_chrs)

    def __contains__(self, key):
        return key in self.all_chrs

# Internal modules #
from .track_util import determine_format, import_implementation
from .track_meta import ChromMetaData, TrackMetaData

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#

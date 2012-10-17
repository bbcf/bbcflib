Genomic tracks
==============

Here is a short tutorial showing how to manage track files with the Python library **btrack** from the **bbcflib** package.

What is it useful for?
----------------------

Bioinformaticians have to deal with large files in multiple formats.
This involves tedious conversions, manipulations, and hundreds of very similar scripts that each of us rewrites constantly.
Also, shell commands have their limits, and most of the time it is simpler/mandatory to use a real programming language instead (here we choose Python).
The purpose of **btrack** is to provide an immediate access to the file's content, neither having to think about formats' specificities, nor convert or decode binary formats.

What formats are supported?
---------------------------

All kinds of raw text files organized in columns can be read (if column fields are specified),
e.g. **csv**, **sam**, or tab-delimited files.
The following formats are automatically recognized and decoded:
**bed**, **wig**, **bedGraph**, **bigWig**, **BAM**, **sqlite**, **sga**, **gff** (**gtf**).

URLs pointing to such files (ex.: http://genome.ucsc.edu/goldenPath/help/examples/bedExample2.bed)
and gzipped files are handled automatically.

Glossary
--------

* **track**: any kind of file that records at least a position in the genome, usually many of them, either intervals (e.g. gene positions) or scores (e.g. density of reads in a give region). By extension, any object containing the same information.
* **stream**: an iterator, similar to a cursor pointing to lines of the track file one by one, never reading back. It can be read over only once.
* **field**: a name describing the content of a column of a track file, for instance the most usual 'chr', 'start', 'end', 'score', 'name' or 'strand'.

How does is work?
-----------------

The library is made of two main functions: ``track`` and ``convert``, and two main classes: ``Track`` and ``FeatureStream``.

When one calls the ``track`` function on a file name, it creates a ``Track`` instance that knows the format of the file and its field names, records genomic information about the species (if specified). It is the interface that gives access to the file's content, similarly to a ``file`` object: it does not itself contain the data, but one can ``read`` and ``write`` it.

Reading a Track object returns a ``FeatureStream`` instance, which is iterable over the data, line by line.
On purpose, it does **not** load all the data in memory as a list would do.

A ``FeatureStream`` is basically an iterator that yields tuples of the type ('chr1',12,14,0.5). Once it has been manipulated at will, it can be written to another ``Track`` using the latter's ``write`` method.

The ``convert`` method can translate a file from a given format to another

How do I use it?
----------------

1. Create (open) a track::

    >>> from bbcflib.btrack import track
    >>> t = track("myfile.bed")

    >>> t
    <bbcflib.btrack.text.BedTrack object at 0x1004a1710>
    >>> t.fields
    ['chr', 'start', 'end', 'name', 'score', 'strand']
    >>> t.format
    'bed'

2. Read the content of a track::

    >>> s = t.read()

    >>> s
    <bbcflib.btrack.FeatureStream object at 0x101244350>
    >>> s.next()
    ('chrII', 45975, 46367, 'YBL092W', 6.0, 1)
    >>> s.next()
    ('chrII', 59818, 60735, 'YBL087C', 8.0, -1)
    >>> for x in s: print x
    ('chrII', 88521, 89123, 'YBL072C', 1.2, -1)
    ('chrII', 168426, 169379, 'YBL027W', 9.0, 1)
    ('chrII', 300166, 301254, 'YBR031W', 11.0, 1)
    ('chrII', 332829, 333810, 'YBR048W', 4.0, 1)
    ...

    >>> s.next()
    StopIteration     # Already read entirely
    >>> s = t.read()  # Reset the stream
    >>> s.next()      # Read again from the beginning
    ('chrII', 59818, 60735, 'YBL087C', 8.0, -1)

3. Write a stream to a new empty track::

    >>> out = track("newfile.wig")
    >>> out.write(s)  # 's' is a FeatureStream
    >>> out.close()

    >>> out
    <bbcflib.btrack.bin.WigTrack object at 0x101769e90>

    # Note: "newfile.wig" is created only at the time you write to it
    # Note: a `Track` may not be written entirely until you `close` it!

4. Convert a track::

    >>> from bbcflib.btrack import convert
    >>> convert("myfile.bed", "myfile.wig")

5. Add genomic information to a Track (from GenRep)::

    >>> t = track("myfile.bed", chrmeta='mm9')  # Mouse assembly name
    >>> t.chmeta
    {'chrY': {'length': 15902555, 'ac': '2752_NC_000087.6'},
     'chrX': {'length': 166650296, 'ac': '2751_NC_000086.6'},
     'chr13': {'length': 120284312, 'ac': '2744_NC_000079.5'},
    ...
    >>> t.assembly
    <bbcflib.genrep.Assembly object at 0x10179b310>
    >>> t.assembly.name
    u'mm9'

See :func:`bbcflib.genrep.Assembly` for more on genomic meta info.

6. Make a selection from a track::

    t = track("myfile.bed")

    # Read only one chromosome:
    s = t.read('chr7')

    # Read only some fields:
    s = t.read(fields=['start','score'])

    # Read only features which either are on chr1 and start within 1000 bp
    # from the beginning of the chromosome, or are on chr2 and end between
    # 3907400 and 4302000:
    sel = [{'chr':'chr1','start':(1,1000)},
           {'chr':'chr2','end':(3907400,4302000)}]
    s = t.read(selection=sel)

7. Read a custom text file::

    t = track("myfile", format='txt', separator='\t',
                        fields=['seq','name','start','info'])

8. Loop on chromosomes::

    t = track("myfile.bed", chrmeta='mm9')
    for chrom in t.chrmeta:
        s = t.read(chrom)
        ...

Advanced features
-----------------

* Streams can be created programmatically, without reference to a track file, either using a list, or an iterator::

    from bbcflib.btrack import FeatureStream
    s = FeatureStream([('chr1',12,13,'a'),('chr1',23,28,'b')],
                      fields=['chr','start','end','name'])

    def generator():
        for x in [10,20,30]:
            yield ('chr1',x,x+5)

    s = FeatureStream(generator(), fields=['chr','start','end'])


* Items are converted to a specific type upon reading and writing, depending on the field name. 
  The conversion function are given in a dictionary called ``intypes`` (converting from text to Python object) and ``outtypes``
  (converting from Python to a text format). For example, the default type for a 'score' field is *float*. 
  If your file contains scores like "NA" which are not convertible with *float()*, then you can specify::

        >>> t = track("myfile.bedgraph",intypes={'score':str})
        >>> t.read().next()
        ('chr1', 1, 101, 'NA')

   Similarly you can convert when writing to file::

        >>> t = track("myfile.bedgraph",outtypes={'score': lambda x=0: "%s" %int(x+.5)})
        >>> t.write([('chr1',10,14,23.56)])
        "chr1    10      14      24"

* To convert from Ensembl-formatted file to the UCSC interval convention::

        >>> t = track("myfile.bedgraph")
        >>> ensembl_to_ucsc(t.read()).next()
        ('chr1', 0, 101, 1.0)
	>>> stream = FeatureStream([('chr1',10,14,23.56)],fields=t.fields)
        >>> t.write(ucsc_to_ensembl(stream),mode='append')
        "chr1    11      14      23.56"



bFlatMajor: data manipulations
------------------------------

**btrack** basically parses track files but does not transform the original data.
To manipulate your data, the **bbcflib** library provides powerful tools to concatenate, intersect, annotate, etc.
It will always take ``FeatureStream`` objects as input, so first open the track using ``btrack.track``,
then ``read`` it and provide the ouput stream to one of **bFlatMajor**'s functions.
Most of them will also return streams, so that you can pass it to another function,
and write the final result to a new ``Track``.

For more info, see **bFlatMajor**'s :doc:`developer documentation <bbcflib_bFlatMajor>` .

Miscellaneous notes
-------------------

* Handling BAM files requires `samtools <http://samtools.sourceforge.net/>`_ .
* Handling bigWig files requires UCSC's *bigWigToBedGraph* (for reading) and *bedGraphToBigWig* (for writing) - look `here <http://genome.ucsc.edu/goldenPath/help/bigWig.html>`_.
* Looping on chromosomes is necessary for several manipulations (see :doc:`bbcflib.bFlatMajor <bbcflib_bFlatMajor>`).
* The ``Track`` class is the parent of multiple subclasses, one for each type of track file (such as :func:`bbcflib.btrack.text.BedTrack` or :func:`bbcflib.btrack.sql.SqlTrack`).
* Look at the :doc:`developer documentation <bbcflib_btrack>` for more details.




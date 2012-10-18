Track manipulations
===================

Here is a short tutorial showing how to manipulate streams of data with the Python library
**bFlatMajor** from the **bbcflib** package.

What is it useful for?
----------------------

The **btrack** library from the **bbcflib** package reads track-like files and returns iterators
over the data (see the documentation `here <http://bbcf.epfl.ch/bbcflib/tutorial_btrack.html>`_).

Usually, one wants to modify this data through a sequence of manipulations before writing it back.
**bFlatMajor** provides a collection of useful functions that take streams as input and perform
common manipulations such as concatenate, intersect or reorder.

Before starting
---------------

This tutorial assumes that you already went through
`btrack's tutorial <http://bbcf.epfl.ch/bbcflib/tutorial_btrack.html>`_,
and thus a *stream* in this context is a short name for a ``btrack.FeatureStream`` object.

A ``FeatureStream`` is an iterator, thus **a stream can be read only once**.
To read the data a second time, one must recreate the stream (read the track again).

How do I use it?
----------------

First, create a stream from scratch or from a track file::

    >>> from bbcflib.btrack import FeatureStream
    >>> s = FeatureStream([('chr1',12,13),('chr1',24,28)], fields=['chr','start','end'])

    # or

    >>> from bbcflib.btrack import track
    >>> t = track("myfile.bed")
    >>> s = t.read(fields=['chr','start','end'])

Then import the function(s) that you need and give your stream as input::

    >>> from bbcflib.bFlatMajor.common import duplicate
    >>> s1 = duplicate(s,'chr','newfield') # copies the 'chr' field to a new one

    >>> s1
    <bbcflib.btrack.FeatureStream object at 0x534248907>
    >>> s1.next()
    ('chr1',12,13,'chr1')
    >>> s1.fields
    ['chr','start','end','chr']

Of course, one can chain functions as long as they return streams::

    >>> from bbcflib.bFlatMajor.common import apply
    >>> s2 = apply(s1,'newfield', lambda x:"aaa") # rename all 'newfield' entries to 'aaa'
    >>> s3 = apply(s2,'end', lambda x:x+12)       # add 12 to all 'end' entries
    >>> s3.next()
    ('chr1',24,40,'aaa')

    # or equivalently:

    >>> s3 = apply(apply(duplicate(s,...) ,...) ,...)

Finally, write the result to a new file using **btrack**::

    >>> from bbcflib.btrack import track
    >>> t = track("newfile.bed", fields=s1.fields)
    >>> t.write(s1)
    >>> t.close()

For many of **bFlatMajor**'s functions,

1. The track must be sorted w.r.t. chromosome, start, end (in this priority order).
   This can be done with a shell ``sort``, but we propose to use the inner
   ``bbcflib.bFlatMajor.common.sorted_stream`` on the stream itself for this purpose
   (will not modify the original file).

2. The function must be applied chromosome by chromosome. Typically::

    from bbcflib.bFlatMajor.common import fusion
    from bbcflib.btrack import track
    t = track("byfile.bed", chrmeta='mm9')
    out = track("newfile.bed")
    for chr in t.chrmeta:
        s = t.read(chr)
        s1 = fusion(s)
        out.write(s1)

    # Running ``fusion`` on the whole genome would paste together
    # regions from different chromosomes.

Both concern every function that has to compare two regions' coordinates.
In general, we advise to always loop on the chomosomes list.

How do I find the function I need?
----------------------------------

**bFlatMajor**'s functions are classified in four submodules:

* **common**: low-level, usual manipulations, usually called implicitly inside of other functions.
* **stream**: functions that return streams.
* **numeric**: functions that return vectors of matrices (*numpy* arrays).
* **figure**: functions that create plots (using a binding to R).

Here are brief descriptions of the main functions (subject to changes):

bFlatMajor.common functions:
----------------------------

``from bbcflib.bFlatMajor.common import *``

* ``copy``: return n independant copies of the input stream.
* ``select``: keeps only the specified fields.
* ``reorder``: change the fields' order.
* ``apply``: apply a custom function to all entries of the specified field(s).
* ``duplicate``: copies one of the fields and its entries to a new one.
* ``concat_fields``: concatenates two fields and their respective entries.
* ``split_field``: when possible, splits a field and its entries into two distinct ones.
* ``map_chromosomes``: translates chromosome names to GenRep standard (e.g. 'chr1').
* ``score_threshold``: filters scores with respect to a threshold.
* ``unroll``: returns one score per genomic position.
* ``sorted_stream``: sorts the stream, by default w.r.t chr, start and end.
* ``shuffled``: returns a stream of randomly located features similar to the original stream.
* ``fusion``: fuses every two overlapping regions A,B into a single one A|B.
* ``cobble``: breaks every two overlapping regions A,B into three: A - A|B - B.

bFlatMajor.stream functions:
----------------------------

``from bbcflib.bFlatMajor.stream import *``

* ``getNearestFeature``: find the nearest gene to each of the input's features.
* ``concatenate``: makes a single stream from the union of several ones.
* ``selection``: filters elements of a stream w.r.t. some given criteria.
* ``neighborhood``: enlarges each of the input's regions.
* ``intersect``: returns the intersection of several streams.
* ``merge_scores``: returns a stream with per-base average (or sum) of several score tracks.
* ``filter_scores``: keeps only scores belonging to a given set of regions.
* ``score_by_feature``: attribute to each given region the sum or average of (independantly) given scores that span the region.
* ``window_smoothing``: applies to the scores a smoothing filter along the sequence.

bFlatMajor.numeric functions:
----------------------------

``from bbcflib.bFlatMajor.stream import *``

* ``score_array``: returns a vector of scores, one for each unique name in the stream.
* ``correlation``: calculates the auto-correlation.
* ``feature_matrix``: returns an array with names as rows and scores as columns, one column for each input score stream.
* ``summed_feature_matrix``: returns an array with for each input score stream, the average score over all features.

bFlatMajor.figure functions:
----------------------------

``from bbcflib.bFlatMajor.figure import *``

* ``scatterplot``: scatter plot (2-d points).
* ``lineplot``: same, but points are bounded by lines.
* ``boxplot``: box plot (quantile plot).
* ``heatmap``: heat map (2-d colored matrix).
* ``pairs``: a scatter plot of each pair of variables one against the other.

Common errors
-------------

* **StopIteration**: The stream is empty, but one tries to read its next element.
* **IndexError**: Most of the time, this is due to an incoherence with the number of fields,
  or a required field that was not found.
* **TypeError**: Common fields, such as 'chr','start','end','frame','strand','score', have
  specific types (resp. str,int,int,int,int,float). Ensure that if you give such a name to a field,
  its entries have the right type.
* **ValueError**: Can have a lot of different causes, but often due to conversion issues
  (see **TypeError**). Ensure that numeric entries are not surrounded by quotes.

Advanced features
-----------------

* Fields order (how to access a field's entry - s.fields.index('field_name'))
* How to use ``combine``
* Build your own function
* The need for ``fusion`` or ``cobble``

More documentation
------------------

* For more details on how each individual function works,
  look at the :doc:`developer documentation <bbcflib_bFlatMajor>`.



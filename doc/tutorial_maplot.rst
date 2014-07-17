M-A plots
=========

Here is a short tutorial showing how to use `maplot.py <https://github.com/delafont/maplot>`_.

Introduction
------------

This script permits to create an `MA-plot` to compare transcription levels of a set of genes
(or other genomic features) in two different conditions C1 ad C2 - typically for the analysis
of RNA-seq data.

Example ref: `Robinson and Oshlack, Genome Biology, 2010 <http://genomebiology.com/2010/11/3/R25>`_

MA-plots show the mean expression of both samples on the horizontal axis (in log10 scale), and the
expression ratio between the two conditions on the vertical axis (in log2 scale). Each gene is a point
on the graph. In average the expression is weak on the right of the graph, and strong on the right.
A gene is over-expressed in condition C1 if the corresponding point is high over the horizontal axis,
and conversely under-expressed in C2 if it is under the axis.

Quantile splines are also drawn to help detecting significantly differentially expressed genes.
Typically one is interested in points over the 95%/99% percentile, or below the 5%/1%.

Usage
-----

Suppose you have a text file that contains unique names in the first column, and scores for two different
conditions in columns 2 and 3 (default, see below). Here is the basic command::

    maplot.py data.txt

This will produce a png image called *maplot_xxxx.png* (see below for annotation).
But more interestingly, you can get an interactive figure where you can click on a gene to get its name::

    maplot.py -m interactive

It may very well happen that your file has its columns in a different order. In this case you can use
the -c (for scores) and -n (for names) options::

    maplot.py -c 3,4 -l 6 data.txt

If you enter several values for names (e.g. -n 1,2,3), values from these columns will be concatenated.
This may be useful if no column contains unique identifiers, but a combination of them yields unique ones.

Also the column separator character may not be a tabulation, or it may not be automatically correctly
detected. You can specify it with the -s option (for tab, try '``C^V``' or '\\t')::

    maplot.py -s ';' data.txt

Different kinds of normalization can generate quite different shapes and scales. In particular, it is
supposed by default that expression levels are given in raw read counts. But RPKM scores will be drawn more
nicely if you specify it with the -f option::

    maplot.py -f rpkm data.txt

If the figure still is not automatically framed as you wish, you can choose bounds on the x and y axes::

    maplot.py --xmin -1 --xmax 2 --ymin -5 --ymax 5 data.txt

Splines (the blue lines) can also be cut on the left and on the right if they unwantedly prolungate or
badly fit in the extremities::

    maplot.py --smin 0 --smax 3 data.txt

Depending on your data, quantile splines may interpolate badly. There are several options to correct
this: -d to change the degree of the spline (default: 4), -b to change the number of bins you split the
data points into to calculate the quantiles (default: 30), -q to simply not display the splines::

    maplot.py -d 3 -b 10 data.txt

Quantiles lines are only approximations, going through points representing given percentiles in each bin.
To get the list of genes that have scores above or below a certain percentile, use the -e option.
A text file will be produced, named *extreme_ratios_xxxx.txt*. For example, for the points beyond
the **5** %-95% interval::

    maplot.py -e 5 data.txt

One can give more than one data file (up to six). Each data set will be plotted in a different color::

    maplot.py data1.txt data2.txt

In 'normal' mode (static .png figure), you can choose to annotate or not one of the datasets.
For a list of three datasets d1,d2,d3, if you want to annotate only d3, use --annotate 001.
It is advised to annotate only very small secondary datasets, otherwise everything risks to be covered by
the labels. A good usage of this is to use the -e option first, and re-enter the file produced as a
second dataset, annotating only the latter::

    maplot.py --annotate 0,1 data.txt extreme_ratios_1234.txt

Finally, you can add a title to your figure::

    maplot.py -t myTitle data.txt

.. image:: [Figure]

Contact
-------

Julien Delafontaine,
EPFL, Bioinformatics and Biostatistics Core Facility (BBCF),
2011-2012,
julien.delafontaine@yandex.com

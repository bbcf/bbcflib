RNA-seq
=======

Here is a short tutorial showing how to launch a RNA-seq analysis from the interface http://htsstation.vital-it.ch/rnaseq/.


New Job
-------

A RNA-seq analysis works from aligned data, given as BAM file(s) through the BAM URL field (there is one BAM file per run). Each BAM file represents a sample ("run"); several samples that were produced in the same conditions (replicates) form a "group".

.. warning:: Alignments have to be performed on the **exonome**.

.. note::  If the alignment was performed on something else than the exonome, you can still retrieve the read counts on the features you aligned to in the output file that is supposed to contain the counts on exons. Make sure in this case to check only the Exons box of the Count level option (see below).


The BAM URLs can be given directly as a http:// or ftp:// address accessible from outside. You can add manually as many groups and as many runs per group you want by using the links `Add group of runs` and `Add run in this group`. Each sample will then be labeled GroupName.RunId in the output files. Make sure to use short group names, without spaces (prefer "_" character to separate words) and without any special character in it (e.g. "%&?!" ).

If you used the HTSstation `mapping module <http://htsstation.vital-it.ch/mapseq/>`_ to do the mapping, you can copy the 20-random characters keys obtained as a result into the `Mapping key` field, and validate using the link `Add data from Mapping`. In such case, all relevant fields will be automatically filled in (see tutorial of our `mapping module <http://htsstation.vital-it.ch/mapseq/>`_ for more details about those fields). To add several samples, successively enter the correponding keys and click on `Add data from Mapping`.

.. image:: images/RNAseq_newjob.png


Select an assembly from the list. Make sure you are selecting the one that was used for the mapping! If your assembly is not listed, please send us an `email <mailto:webmaster.bbcf@epfl.ch>`_.

Name your analysis in the `Run description` field. Please, use short names, without spaces (prefer "_" character to separate words) and without any special characters (e.g., "%&?!" ... ).
Submit your e-mail in order to receive an message upon completion of the pipeline.

.. image:: images/RNAseq_generals.png

Finally, click on the `Create` button and confirm to launch the job.

.. image:: images/RNAseq_create.png


Results
_______

When the job finishes successfully, you will receive an e-mail with a link to the page where you can download the results. Results include mapping reports for each sample, and tab-delimited files containing counts and rpkm for genes, exons and transcripts for all samples at once.

From there you can also create an interactive MA-plot to look for differential expression: using the checkboxes, select the two samples you want to compare, and click on the `Compute MA-plot` button.

.. image:: images/RNAseq_create_maplot.png

On the graph's page, click on the points you are interested in to display its name in the column on the right. Click on it again to remove it from the list. Click on the name to get information about the selected feature from Ensembl. Note that the graph may take a long time to load and react if there are a lot of features to draw.

.. image:: images/RNAseq_maplot.png

RNA-seq Tutorial
================

Here is a short tutorial showing how to launch a RNA-seq analysis from the interface http://htsstation.vital-it.ch/rnaseq/.


New Job
-------
A RNA-seq analysis works from aligned data, given as BAM file(s) through the BAM URL field (there is one BAM file per run).

.. warning:: The alignment has to be performed on the EXONOME!

The URL can be given directly (as a http:// or ftp:// address accessible from outside) or retrieved by using the `Mapping key` obtained when running the `mapping module <http://htsstation.vital-it.ch/mapseq/>`_. In such case, fields related to the sequencing facility (`#Run, Facility, Machine, Run number and Lane number`) might be automatically filled in if relevant (see tutorial of mapseq for more details about those fields).
The image below shows how such key can be used. In this example, data were not coming directly from a sequencing facility, so the corresponding fields are left empty.

.. image:: images/RNAseq_newJob1.png

Each BAM file represents a sample; several samples that were produced in the same conditions (replicates) can be grouped. You can add as many groups and as many runs per group you want by using the links `Add group of runs` and `Add run in this group` (see tutorial xxx for a more detailed explanation of what are groups and runs). Each sample will then be labeled GroupName.RunId in the output files. Make sure to use short group names, without spaces (prefer "_" character to separate words) and without any special character in it (e.g. "%&?!" ).

Finally, you will have to select an assembly from the list. Make sure you are selecting the one that was used for the mapping! If your assembly is not listed, please send us an `email <mailto:webmaster.bbcf@epfl.ch>`_.

.. note::  If the alignment was performed on something else than the exonome, you can still retrieve the read counts on the features you aligned to in the output file that is supposed to contain the counts on exons. Make sure in this case to check only the Exons box of the Pileup level option.


General Parameters
------------------

Run description: name your analysis. Please, use short names, without spaces (prefer "_" character to separate words) and without any special characters (e.g., "%&?!" ... ).
Submit your e-mail in order to receive an message upon completion of the pipeline.

.. image:: images/RNAseq_newJob2.png


RNA-Seq options
---------------

Check Genes/Exons/Transcripts pileup levels to get the counts/rpkm respectively on genes, exons and transcripts. Each of the three will be a distinct tab-delimited output file. The only reason not to check all boxes is to save time if you know you are interested in only one of them.

.. image:: images/RNAseq_newJob3.png


Then click on Create and confirm to launch the job.


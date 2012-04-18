SNP
===

Here is a short tutorial showing how to search for SNPs from the interface http://htsstation.vital-it.ch/snp/.


New Job
-------

An RNA-seq analysis works from aligned data, given as BAM file(s) through the BAM URL field (there is one BAM file per run). Each BAM file represents a sample ("run"); several samples that were produced in the same conditions (replicates) form a "group".

.. warning:: Alignments have to be performed on the **exonome**.

.. note:: If the mapping was performed independantly, modify the BAM header so that the identifier of each exon is of the form `@SQ  SN:exonID|geneID|start|end|strand  LN:length`.

The BAM URLs can be given directly as an `http://` or `ftp://` address accessible from outside. You can add manually as many groups and as many runs per group you want by using the links `Add group of runs` and `Add run in this group`. Each sample will then be labeled *group_name.run_index* in the output files. Make sure to use short group names, without spaces (prefer "_" character to separate words) and without any special character in it (e.g. "%&?!" ).

If you used the HTSstation `mapping module <http://htsstation.vital-it.ch/mapseq/>`_ to do the mapping, you can copy the 20-random characters keys obtained as a result into the `Mapping key` field, and validate using the link `Add data from Mapping`. In such case, all relevant fields will be automatically filled in (see tutorial of our `mapping module <http://htsstation.vital-it.ch/mapseq/>`_ for more details about those fields). To add samples from other independant mappings, successively enter the correponding keys and click on `Add data from Mapping`.

.. image:: images/RNAseq_newjob.png

Then select an assembly from the list. Make sure you are selecting the one that was used for the mapping. If your assembly is not listed, please send us an `email <mailto:webmaster.bbcf@epfl.ch>`_.

Name your analysis in the `Run description` field. Please, use short names, without spaces (prefer "_" character to separate words) and without any special characters (e.g., "%&?!" ... ).
Submit your e-mail in order to receive an message upon completion of the pipeline.

.. image:: images/RNAseq_generals.png

Finally, click on the `Create` button and confirm to launch the job.

.. image:: images/RNAseq_create.png


Results
-------

When the job finishes successfully, you will receive an e-mail with a link to the page where you can download the results. Results consist in tab-delimited files containing counts and rpkm for genes, exons and transcripts, a differential expression analysis for each pair of groups in the experiment, and mapping reports of reads that did not map to the exonome first.


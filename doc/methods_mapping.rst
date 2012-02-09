Mapping
=======

This module is a wrapper to `Bowtie <http://bowtie-bio.sourceforge.net/manual.shtml>`_. We always allow multiple matches per read, and document this in the `bam` output by adding the optional *NH* flag, see :py:func:`bbcflib.mapseq.add_nh_flag`.

We use two `C` programs to analyze `bam` files:
 * `bamstat <http://github.com/bbcf/bbcfutils/blob/master/C/bamstat.cc>`_
 * `bam2wig <http://github.com/bbcf/bbcfutils/blob/master/C/bam2wig.cc>`_

The first generates statistics on reads mapping. In particular, the number of reads per strand, if unbalanced, can indicated a sequencing bias (see the Bowtie documentation).

The second creates genomic density profiles per strand or merged. We have implemented several optional features described in `Leleu et al. <http://www.ncbi.nlm.nih.gov/pubmed/20861161>`_:
 * Correct for PCR artifacts: only at most **n** reads per strand-specific genomic position are kept, where **n** is computed as the 95% percentile of a Poisson distribution with the same mean as the average genome coverage. This is done by parsing the position-ordered output `bam` and counting the successive lines with the same position, strand and library, see :py:func:`bbcflib.mapseq.remove_duplicate_reads`.
 * Read density smoothing: the read length hass no biological meaning, in principle only the first maping position should be counted in the genomic density. However extending the coverage to *n* bases downstream makes a smoother density which is visually more appealing. The value of this extension is arbitrary and can less or more than the actual read length.
 * Down-weight multiple hits: using the `NH` flag described above, we count *1/NH* read at each alignement, so that the sum of all scores in the genome is always equal to number of reads times the read extension.
 * Normalization: use total numer of mapped reads times *1e-7* as a normalization factor.
 * Strand shifting and merging: reads mapping to the reverse and forward strand correspond to the ends of the DNA fragment making the sequencing library, it is generally advised to shift the mapping position downstream by half the difference between fragment size and read extension.

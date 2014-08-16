Using rnacounter to count reads in genomic intervals
====================================================

`rnacounter` estimates abundances of genes and their different transcripts
from read densities. Exons and introns can also be quantified.
It requires a genome-level BAM file and a
GTF/GFF file describing the exon structure, such as those provided by Ensembl or GenRep.
The GTF is assumed to be sorted at least w.r.t. chromosome name,
and the chromosome identifiers in the GTF must be the same as the BAM references.
The method used is described in [<ref>].

Most basic usage::

   rnacounter test.bam test.gtf > counts_table.txt

The command above will create a tab-delimited text file of gene counts, their RPKM
and annotation information sur as the genomic location.
Many options are then available, listed below.

Use `rnacounter join` to merge several output files produced using **the same GTF**,
to create a single table with counts from all samples::

   rnacounter join tab1.txt tab2.txt ... > tab_all.txt


Options:
--------

* `-h, --help` and `-v, --version`::

    Display information about the program usage / the version currently installed.

* `-s, --stranded`::

    If the protocol was strand-specific and this option is provided,
    sense and antisense counts are both reported in two consecutive lines
    with a different tag in the last column.
    They can be split afterwards by piping the result as for instance with
    `... | grep 'antisense'`.
    Using the `--threshold` option together with `--stranded`
    will exclude only elements with both sense and antisense counts under the threshold.

* `-n, --normalize`::

    RPKM are automatically calculated together with raw read counts. RPKM are counts
    divided by the length of the transcript as well as by a sample-specific
    normalization constant, usually the total number of aligned reads in the sample (default).
    This value can be changed to a user-defined integer.
    Typically, if you want to compare the same gene in several samples,
    the normalization will cancel out anyway
    and giving `-n 1` will speed up the process since it will skip counting the alignments.
    Some stats programs also require raw counts anyway and do their own normalization.
    To get FPKM instead, see `--fraglength`.

* `-f, --fraglength`::

    Since in a transcript of length L there are only L-F+1 different positions where
    one can cut a fragment of length F, one should correct for this bias before RPKM
    calculation (then usually called FPKM). Typical fragment lengths are around 350;
    default value is 1 (no correction).

* `--nh`::

    A flag "NH" can be added to BAM files to indicate the number of times the read
    could be mapped to different locations in the genome. Adding this option
    will take this number into account by adding 1/NH instead of 1 to an exon read count.

* `--noheader`::

    By default the program adds one line with column descriptors on top of the output file.
    For easier piping the result in some other program, such as `cut`, one can choose to
    not add the header by adding this option.

* `--threshold`::

    Features with counts inferior or equal to the given threshold (positive number)
    will not be reported in the ouput. By default everything is reported
    - even with zero counts.

* `--gtf_type`::

    Usually one uses standard (Ensembl etc.) GTF files to count reads in
    exons/genes/transcripts. The only lines of interest are then the ones with
    value "exon" (default) in the 3rd column. If you are counting something else
    or provided your own, differently formatted GTF, with this option you can specify
    the 3rd column value of the lines to consider.

* `--format`::

    One can also give an annotation file in BED format, in which case each line
    is considered as an independant, disjoint interval with no splicing structure.
    Default is "gtf", can be changed to "bed".

* `-t, --type`::

    The type of feature you want to count reads in. Can be "genes" (default),
    "transcripts", "exons" or "introns".
    One can give multiple comma-separated values, in which case all
    the different features will be mixed in the output but can easily be split
    using the last column tag, as for instance with `... | grep 'exon'`.
    Then if `--method` is specified it must have the same number of values as `type`,
    also as a comma-separated list, or a single one that is applied to all types.

* `-c, --chromosomes`::

    Consider only a subset of all chromosomes by providing a comma-separated list
    of chromosome names (that must match those of the GTF and BAM).

* `-o, --output`::

    The output is `stdout` by default (output directly to screen), which permits
    redirection to a file. Alternatively one can redirect the standard output to
    a file using this option. If the file name already exists, it will be overwritten.

* `-m, --method`::

    Feature counts are inferred from the counts on (slices of) exons
    with the chosen `--method`: "raw" (default, HTSeq-like) or
    "nnls" (non-negative least squares, see [<ref>]).


Miscellaneous notes:
--------------------

* Overlapping regions:
  In "raw" counting mode, regions spanned by two or more genes, together with the
  alignements inside these regions, are ignored - as in HTSeq's "union" mode.

* Exons and introns:
  Because annotated exons often overlap a lot, in "raw" mode, "exon" counts are actually
  that of their disjoint slices, and their name in the output table is formatted as
  "exon1|exon2" if a slice is spanned by exon1 and exon2. In "nnls" mode, exon counts
  are inferred from disjoint slices as for genes.

  Intronic regions also annotated as exons in some alternative transcripts are
  ignored whatever the chosen method is. Because they don't have official IDs,
  introns slices are given names following this pattern:
  "T-i7" means that is is the 7th intron of transcript T,
  "T1-i1|T2-i3|..." means it is part of several transcripts.

* Non-integer counts:
  The fact that some reads cross exon boundaries as well as considering the NH flag
  make the reported number not be integers. They still represent count data and can
  be rounded afterwards if necessary.

* Custom input:
  If your GTF does not represent exons but custom genomic intervals to simply count
  reads in, provide at least a unique `exon_id` in the attributes as a feature name,
  and the type field (column 3) must be set to 'exon' or specified with the
  `--gtf_ftype` option. If not specified, `gene_id`, `transcript_id` and `exon_id`
  will all get the value of `exon_id`.

* Paired-end support:
  At the moment alignments of paired-end reads are not treated specially, i.e.
  all reads are considered as single-end.


Examples:
---------

* Probably the best way to get isoforms counts::

    rnacounter -t transcripts -m nnls --nh -f 350 sample.bam mouse.gtf > transcript_counts.txt

* Compare gene counts between two conditions, HTSeq-like::

    rnacounter -n 1 group1.bam mouse.gtf > gene_counts1.txt
    rnacounter -n 1 group2.bam mouse.gtf > gene_counts2.txt
    rnacounter join gene_counts1.txt gene_counts2.txt > gene_counts.txt

  Then send it to DESeq/EdgeR/whatever other stats program that asks for such a table.


Reference:
----------

<?>


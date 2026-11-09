.. _commands/profile:

profile
=======

Calculate percentage of reads mapped to each bin.

Produces a table indicating the percentage of reads **mapped** to an assembly which are assigned to each
genome bin or assigned to unbinned contigs. This information is also used to determine the percentage of each
genome bin relative to all genome bins under consideration. This is a useful indication of the relative
percentage of different populations within a community when the majority of populations are represented by a
genome bin and the majority of reads map to the assembly. If the majority of reads do not map to the
assembly, these results will not be a reliable indication of the relative proportion of different
populations. This command requires a file indicating the coverage profile of all sequences within the genome
bins. This file can be created with the :ref:`coverage <commands/coverage>` command.

The columns of the output table are (the last four columns are given for each BAM file):

* **Bin Id**: unique identifier of bin
* **Bin size (Mbp)**: size of bin in Mbp
* **mapped reads**: number of reads mapped to contigs comprising the bin
* **% mapped reads**: (reads mapped to bin) / (total number of reads mapped to assembly)
* **% binned populations**: estimates the proportion of a bin relative to all recovered bins. This is
  determined from the percentage of reads mapped to a bin normalized for the size of the bin and considering
  only reads mapped to bins. Specifically, it is calculated as [(% mapped reads) / (bin size)] * (1 / C),
  where C is the sum of size-adjusted bin coverages over all bins (i.e., the sum of
  (% mapped reads) / (bin size) over all bins).
* **% community**: estimates the proportion of a bin relative to the number of reads mapped to assembled
  contigs and adjusted for the size of the bin. Assuming the majority of reads map to an assembled contig,
  this is an estimate of the relative proportion of a bin (i.e., population) within the community. It is
  calculated as (% binned populations) * (100 - (percentage of reads assigned to unbinned contigs)).

Arguments
---------

.. literalinclude:: help/profile.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm coverage ./bins coverage.tsv example_1.bam
    checkm profile coverage.tsv

.. _commands/unbinned:

unbinned
========

Identify unbinned sequences.

Given a set of genome bins along with a FASTA file of sequences, determines which of these sequences are not
currently contained in a genome bin. This is useful for determining which sequences from an assembly were not
binned by an automated binning algorithm. The unbinned sequences are written to a FASTA file and their
length and GC to a tab-separated statistics file.

Arguments
---------

.. literalinclude:: help/unbinned.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm unbinned ./bins seqs.fna unbinned.fna unbinned_stats.tsv

.. _commands/ssu_finder:

ssu_finder
==========

Identify SSU (16S/18S) rRNAs in sequences.

Identifies SSU (16S and 18S) rRNA genes residing on sequences, if these sequences are contained within a
genome bin. A summary of the identified genes is written to ``<output_dir>/ssu_summary.tsv`` and their
sequences to ``<output_dir>/ssu.fna``.

Arguments
---------

.. literalinclude:: help/ssu_finder.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm ssu_finder seqs.fna ./bins ./ssu_finder

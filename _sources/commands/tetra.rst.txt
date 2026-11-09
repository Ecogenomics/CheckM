.. _commands/tetra:

tetra
=====

Calculate tetranucleotide signature of sequences.

Produces tetranucleotide signatures for all sequences within a FASTA file. Tetranucleotide signatures are
required by the :ref:`tetra_plot <commands/tetra_plot>`, :ref:`dist_plot <commands/dist_plot>`, and
:ref:`outliers <commands/outliers>` commands.

Arguments
---------

.. literalinclude:: help/tetra.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm tetra seqs.fna tetra.tsv

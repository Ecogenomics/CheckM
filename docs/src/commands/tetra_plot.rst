.. _commands/tetra_plot:

tetra_plot
==========

Create tetranucleotide distance (TD) histogram and delta-TD plot.

Provides a plot analogous to the :ref:`gc_plot <commands/gc_plot>` suitable for assessing the
tetranucleotide signatures of sequences within a genome bin. The Manhattan distance is used to determine the
difference between each sequence's tetranucleotide signature and the tetranucleotide signature of the entire
genome bin.

This plot requires a file indicating the tetranucleotide signature of all sequences within the genome bins.
This file can be created with the :ref:`tetra <commands/tetra>` command.
This command requires the output folder of the :ref:`analyze <commands/analyze>`,
:ref:`lineage_wf <commands/lineage_wf>`, or :ref:`taxonomy_wf <commands/taxonomy_wf>` command
(``results_dir``).

A plot is written for each bin to ``<output_dir>/<bin id>.tetra_dist_plots.<image_type>``.

Arguments
---------

.. literalinclude:: help/tetra_plot.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm tetra seqs.fna tetra.tsv
    checkm tetra_plot ./output ./bins ./plots tetra.tsv 95

.. _commands/dist_plot:

dist_plot
=========

Create image with GC, CD, and TD distribution plots together.

Produces a single figure combining the plots produced by :ref:`gc_plot <commands/gc_plot>`,
:ref:`coding_plot <commands/coding_plot>`, and :ref:`tetra_plot <commands/tetra_plot>`.

This plot requires a file indicating the tetranucleotide signature of all sequences within the genome bins.
This file can be created with the :ref:`tetra <commands/tetra>` command.
This command requires the output folder of the :ref:`analyze <commands/analyze>`,
:ref:`lineage_wf <commands/lineage_wf>`, or :ref:`taxonomy_wf <commands/taxonomy_wf>` command
(``results_dir``).

A plot is written for each bin to ``<output_dir>/<bin id>.ref_dist_plots.<image_type>``.

Arguments
---------

.. literalinclude:: help/dist_plot.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm tetra seqs.fna tetra.tsv
    checkm dist_plot ./output ./bins ./plots tetra.tsv 95

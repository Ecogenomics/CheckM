.. _commands/marker_plot:

marker_plot
===========

Plot position of marker genes on sequences.

Plots the position of marker genes on sequences within a genome bin. This provides information regarding the
extent to which marker genes are collocated. The number of marker genes within each fixed-size window along a
sequence is indicated with different colours. Sequences without any marker genes are not shown.

This command requires the output folder of the :ref:`analyze <commands/analyze>`,
:ref:`lineage_wf <commands/lineage_wf>`, or :ref:`taxonomy_wf <commands/taxonomy_wf>` command
(``results_dir``), and the :ref:`qa <commands/qa>` command must also have been run on this folder (the
workflows run ``qa``).

A plot is written for each bin with marker genes to ``<output_dir>/<bin id>.marker_pos_plot.<image_type>``.

Arguments
---------

.. literalinclude:: help/marker_plot.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm marker_plot ./output ./bins ./plots

.. _commands/coding_plot:

coding_plot
===========

Create coding density (CD) histogram and delta-CD plot.

Provides a plot analogous to the :ref:`gc_plot <commands/gc_plot>` suitable for assessing the coding density
of sequences within a genome bin.

This command requires the output folder of the :ref:`analyze <commands/analyze>`,
:ref:`lineage_wf <commands/lineage_wf>`, or :ref:`taxonomy_wf <commands/taxonomy_wf>` command
(``results_dir``).

A plot is written for each bin to ``<output_dir>/<bin id>.coding_density_plots.<image_type>``.

Arguments
---------

.. literalinclude:: help/coding_plot.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm coding_plot ./output ./bins ./plots 95

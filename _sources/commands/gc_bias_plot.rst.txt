.. _commands/gc_bias_plot:

gc_bias_plot
============

Plot bin coverage as a function of GC.

Plots the coverage of windows within each genome bin as a function of their GC content, using reads mapped
to the sequences in a BAM file. The BAM file must be sorted and indexed (i.e. a ``.bai`` index file must be
present; see :ref:`coverage <commands/coverage>`).

A plot is written for each bin to ``<output_dir>/<bin id>.gc_bias_plot.<image_type>``.

Arguments
---------

.. literalinclude:: help/gc_bias_plot.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm gc_bias_plot ./bins ./plots example.bam

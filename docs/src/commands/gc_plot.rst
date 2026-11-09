.. _commands/gc_plot:

gc_plot
=======

Create GC histogram and delta-GC plot.

Provides a two-pane plot suitable for assessing the GC distribution of sequences within a genome bin. The
first pane is a histogram of the percentage of non-overlapping windows (5 kbp by default; see ``-w``) with a
given percent GC. A typical genome will produce a unimodal distribution. The second pane plots each sequence
in the genome bin as a function of its deviation from the average GC of the entire genome (x-axis) and
sequence length (y-axis). The dashed red lines indicate the expected deviation from the mean GC as a function
of length. This expected deviation is pre-calculated from a set of trusted reference genomes and the
percentile plotted is provided as an argument to this command (``dist_value``). A good default value to use
for this distribution parameter is 95.

A plot is written for each bin to ``<output_dir>/<bin id>.gc_plots.<image_type>``.

Arguments
---------

.. literalinclude:: help/gc_plot.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm gc_plot ./bins ./plots 95

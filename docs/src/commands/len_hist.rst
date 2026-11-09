.. _commands/len_hist:

len_hist
========

Sequence length histogram.

Produces a histogram of the number of sequences within a genome bin at different sequence length intervals.
This provides additional information regarding the quality of an assembled genome.

A plot is written for each bin to ``<output_dir>/<bin id>.len_hist.<image_type>``.

Arguments
---------

.. literalinclude:: help/len_hist.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm len_hist ./bins ./plots

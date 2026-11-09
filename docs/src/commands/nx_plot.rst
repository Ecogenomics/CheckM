.. _commands/nx_plot:

nx_plot
=======

Create Nx-plots.

Produces a plot indicating the Nx value of a genome bin for all values of x. This provides a more
comprehensive view of the quality of an assembly than simply considering
`N50 <https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics>`_.

A plot is written for each bin to ``<output_dir>/<bin id>.nx_plot.<image_type>``.

Arguments
---------

.. literalinclude:: help/nx_plot.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm nx_plot ./bins ./plots

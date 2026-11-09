.. _commands:

Commands
========

CheckM is executed from the command line and consists of a series of commands in order to support a number of
different analyses and workflows. These commands are organized into several related groups, as listed by
``checkm -h``. Help for each command can be shown with ``checkm <command> -h``.

The two most common workflows are to assess genomes using either *lineage-specific* or *taxonomic-specific*
marker sets (see :ref:`workflows`).


Lineage-specific marker set
---------------------------

.. toctree::
   :maxdepth: 1

   tree
   tree_qa
   lineage_set


Taxonomic-specific marker set
-----------------------------

.. toctree::
   :maxdepth: 1

   taxon_list
   taxon_set


Apply marker set to genome bins
-------------------------------

.. toctree::
   :maxdepth: 1

   analyze
   qa


Common workflows
----------------

These workflows combine the above commands.

.. toctree::
   :maxdepth: 1

   lineage_wf
   taxonomy_wf


.. _commands/plots:

Plots
-----

CheckM can produce a number of plots for assessing the quality of genome bins.

Reference distribution plots:

.. toctree::
   :maxdepth: 1

   gc_plot
   coding_plot
   tetra_plot
   dist_plot

General plots:

.. toctree::
   :maxdepth: 1

   nx_plot
   len_hist
   marker_plot
   gc_bias_plot


Bin exploration and modification
--------------------------------

.. toctree::
   :maxdepth: 1

   unique
   merge
   outliers
   modify


Utility functions
-----------------

CheckM also provides a number of additional commands that may be useful for exploring genome bins. Some of
these commands also produce summary statistics of sequences required by specific plotting commands.

.. toctree::
   :maxdepth: 1

   unbinned
   coverage
   tetra
   profile
   ssu_finder


Reference data and testing
--------------------------

.. toctree::
   :maxdepth: 1

   data
   test


Removed commands
----------------

Commands removed in CheckM v1.1.x (e.g., ``bin_qa_plot``, ``par_plot``, ``cov_pca``, ``tetra_pca``,
``len_plot``, ``bin_union``, ``bin_compare``, ``join_tables``) are listed in :ref:`python3/removed`.


checkm -h
---------

.. literalinclude:: help/checkm.txt
   :language: text

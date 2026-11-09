.. _commands/tree:

tree
====

Place bins in the reference genome tree.

This is the first step of the :ref:`lineage-specific workflow <workflows/lineage>`. Genes are called in each
bin with Prodigal (unless ``--genes`` is used), phylogenetically informative marker genes are identified with
HMMER, and each bin is placed in the reference genome tree with pplacer.

* ``--reduced_tree``: use the reduced reference genome tree, which requires < 16 GB of memory
  (see :ref:`installing/system_requirements`).
* ``--ali``: generate a HMMER alignment file for each bin in
  ``<output folder>/bins/<bin id>/hmmer.tree.ali.txt``.
* ``--nt``: generate nucleotide gene sequences for each bin in ``<output folder>/bins/<bin id>/genes.fna``.
* ``--genes``: bins contain genes as amino acids instead of nucleotide contigs.
* ``--pplacer_threads``: number of threads used by pplacer. By default, the number of threads given by
  ``-t`` is used, to a maximum of 64. pplacer's memory usage increases only slightly with additional
  threads.

Called genes are written to ``<output folder>/bins/<bin id>/genes.*``. Bins smaller than 100 kbp are
processed with Prodigal's metagenomic mode (``-p meta``).

.. note::

   The following heuristic is used to establish the translation table used by Prodigal: use table 11 unless
   the coding density using table 4 is more than 5 percentage points higher than when using table 11 and the
   coding density under table 4 is > 70%. Distinguishing between tables 4 and 25 is challenging so CheckM
   does not attempt to distinguish between these two tables. If you know the correct translation table for
   your genomes, it is recommended that you call genes outside of CheckM and provide CheckM with the protein
   sequences (see ``--genes``).

Arguments
---------

.. literalinclude:: help/tree.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm tree -t 8 -x fa ./bins ./output

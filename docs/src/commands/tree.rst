.. _commands/tree:

tree
====

Place bins in the reference genome tree.

This is the first step of the :ref:`lineage-specific workflow <workflows/lineage>`. Genes are called in each
bin with Prodigal (unless ``--genes`` is used), phylogenetically informative marker genes are identified with
HMMER, and each bin is placed in the reference genome tree with `gPlace <https://github.com/cmc-aau/gPlace>`_.

* ``--ali``: generate a HMMER alignment file for each bin in
  ``<output folder>/bins/<bin id>/hmmer.tree.ali.txt.gz``.
* ``--nt``: generate nucleotide gene sequences for each bin in ``<output folder>/bins/<bin id>/genes.fna.gz``.
* ``--genes``: bins contain genes as amino acids instead of nucleotide contigs.
* ``-t``: number of threads, which is also the number of threads used by gPlace.
* ``--gplace_max_memory``: memory gPlace may use to place bins, e.g. ``64G`` or ``500M`` (a plain number is
  GiB). If the placement needs more memory, gPlace places bins tile by tile within this limit, which gives
  identical results. By default, gPlace sets this from the available memory. If the limit is too small,
  CheckM reports gPlace's error message, which gives the smallest limit that can be used (see
  :ref:`installing/system_requirements`).

Called genes are written to ``<output folder>/bins/<bin id>/genes.*.gz`` (see :ref:`outputs/files`). Bins
smaller than 100 kbp are processed with Prodigal's metagenomic mode (``-p meta``).

.. note::

   CheckM <= 1.2.10 placed bins with pplacer and provided the ``--reduced_tree`` and ``--pplacer_threads``
   options, which were removed in CheckM 1.3.0.

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
    checkm tree -t 8 -x fa --gplace_max_memory 8G ./bins ./output

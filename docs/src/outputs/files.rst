.. _outputs/files:

Output files
============

The output folder of the :ref:`lineage_wf <commands/lineage_wf>` and :ref:`taxonomy_wf <commands/taxonomy_wf>`
workflows contains:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - File
     - Description
   * - ``quality_report.tsv``
     - Completeness, contamination, strain heterogeneity and genome statistics of each bin: the
       tab-separated table produced by ``qa -o 2 --tab_table`` (see :ref:`outputs/qa`).
   * - ``lineage.ms.gz`` or ``<taxon>.ms.gz``
     - Marker file used to assess the bins (see :ref:`commands/lineage_set` and :ref:`commands/taxon_set`).
   * - ``checkm.log``
     - Log of the commands run on the output folder.
   * - ``bins/<bin id>/``
     - Called genes (``genes.faa.gz``, ``genes.gff.gz``, and ``genes.fna.gz`` with ``--nt``) and HMMER
       results (``hmmer.tree.txt.gz``, ``hmmer.analyze.txt.gz``, and the ``*.ali.txt.gz`` alignments with
       ``--ali``) of each bin.
   * - ``storage/``
     - Genome statistics (``bin_stats.tree.tsv.gz``, ``bin_stats.analyze.tsv.gz``, ``bin_stats_ext.tsv.gz``),
       marker gene positions (``marker_gene_stats.tsv.gz``), alignments of marker genes found multiple times
       in a bin (``aai_qa.tar.gz``, used to calculate strain heterogeneity), and HMM information
       (``*.pkl.gz``).
   * - ``storage/tree/``
     - Placement of the bins in the reference genome tree (:ref:`lineage_wf <commands/lineage_wf>` only):
       marker gene alignments (``<marker>.masked.faa.gz``), the concatenated alignment
       (``concatenated.fasta.gz``), gPlace's placements (``concatenated.pplacer.json.gz``) and tree
       (``concatenated.tre.gz``), and gPlace's logs (``pplacer.out``, ``gplace.log``,
       ``gplace.warnings.log``).

The ``bins`` and ``storage`` folders are needed to re-run :ref:`qa <commands/qa>` and the
:ref:`plotting commands <commands/plots>` on the output folder. Use ``--rm_tmp_files`` to remove them once a
workflow has finished if only ``quality_report.tsv`` is needed.


Compressed files
----------------

From CheckM 1.3.0, most output files are gzip-compressed, which makes the output folder ~3.5 times smaller
(e.g., 1.2 GB instead of 4.4 GB for 1,000 genomes). They can be read with standard tools, e.g.:

.. code-block:: bash

    zcat ./output/storage/bin_stats_ext.tsv.gz | head
    tar -tzf ./output/storage/aai_qa.tar.gz

CheckM commands read both compressed outputs and the uncompressed outputs of CheckM <= 1.2.10, and accept
file names without the ``.gz`` extension (e.g., ``checkm qa ./output/lineage.ms ./output``). If both a
compressed and an uncompressed marker file exist (e.g., ``lineage.ms`` and ``lineage.ms.gz`` in an output
folder reused from an earlier version of CheckM), the most recently modified file is used.

.. note::

   CheckM <= 1.2.10 wrote these files uncompressed (e.g., ``lineage.ms``, ``bins/<bin id>/genes.faa`` and
   ``storage/bin_stats_ext.tsv``) and stored the alignments used to calculate strain heterogeneity in the
   ``storage/aai_qa`` folder. Scripts that read these files directly need to be updated for CheckM >= 1.3.0.

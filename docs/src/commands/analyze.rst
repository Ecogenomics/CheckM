.. _commands/analyze:

analyze
=======

Identify marker genes in bins and calculate genome statistics.

Identifies the marker genes in a marker file within each bin. The marker file can be a lineage-specific
marker file produced by :ref:`lineage_set <commands/lineage_set>`, a taxonomic-specific marker file produced
by :ref:`taxon_set <commands/taxon_set>`, or a HMM file of :ref:`custom marker genes <workflows/custom>`.
Completeness and contamination are then estimated with the :ref:`qa <commands/qa>` command.

* ``--ali``: generate a HMMER alignment file for each bin in
  ``<output folder>/bins/<bin id>/hmmer.analyze.ali.txt.gz``.
* ``--nt``: generate nucleotide gene sequences for each bin in ``<output folder>/bins/<bin id>/genes.fna.gz``.
* ``--genes``: bins contain genes as amino acids instead of nucleotide contigs.

Outputs:

* called gene files are placed in ``<output folder>/bins/<bin id>/genes.*.gz`` (see :ref:`outputs/files`)
* HMMER result files are placed in ``<output folder>/bins/<bin id>/hmmer.*``

Genes are called with Prodigal as described for the :ref:`tree <commands/tree>` command.

Arguments
---------

.. literalinclude:: help/analyze.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm analyze ./output/lineage.ms ./bins ./output
    checkm analyze -t 8 -x fa bacteria.ms ./bins ./output

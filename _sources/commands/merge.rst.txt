.. _commands/merge:

merge
=====

Identify bins with complementary sets of marker genes.

Identifies genome bins with complementary sets of marker genes. Merging such bins will result in a notable
increase in completeness with only a marginal, or no, increase in contamination. Caution must be exercised
before merging two bins.

To identify complementary sets of marker genes, a common set of markers must be sought in each bin, so this
command requires a taxonomic-specific marker set (see :ref:`taxon_set <commands/taxon_set>`) or a HMM file of
:ref:`custom marker genes <workflows/custom>`; lineage-specific marker files are not supported. We generally
consider both the bacterial and archaeal sets produced by the ``taxon_set`` command.

It is entirely possible that two bins will have complementary sets of marker genes, but should not be
merged. We have observed this situation many times. Additional information should be used to confirm the
merging of bins. We only merge bins after verifying that they have similar genomic characteristics (e.g.,
GC, coverage) and are placed in similar locations within a reference genome tree. This information is
available using the :ref:`qa <commands/qa>` and :ref:`tree_qa <commands/tree_qa>` commands, respectively.

Pairs of bins meeting the ``--delta_comp``, ``--delta_cont``, ``--merged_comp``, and ``--merged_cont``
criteria are written to ``<output_dir>/merger.tsv``.

Arguments
---------

.. literalinclude:: help/merge.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm taxon_set domain Bacteria bacteria.ms
    checkm merge bacteria.ms ./bins ./output

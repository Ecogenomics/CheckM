.. _commands/lineage_set:

lineage_set
===========

Infer lineage-specific marker sets for each bin.

This step of the :ref:`lineage-specific workflow <workflows/lineage>` uses the placement of each bin in the
reference genome tree, as determined by the :ref:`tree <commands/tree>` command, to infer a lineage-specific
marker set for each bin. The marker sets are written to the specified marker file, which is used by the
:ref:`analyze <commands/analyze>` and :ref:`qa <commands/qa>` commands.

* ``-u, --unique``: minimum number of unique phylogenetic markers required to use a lineage-specific marker
  set, otherwise a domain-level marker set is used.
* ``-m, --multi``: maximum number of multi-copy phylogenetic markers before defaulting to a domain-level
  marker set.
* ``--force_domain``: use domain-level marker sets for all bins.
* ``--no_refinement``: do not perform lineage-specific marker set refinement.

Arguments
---------

.. literalinclude:: help/lineage_set.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm lineage_set ./output ./output/lineage.ms

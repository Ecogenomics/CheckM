.. _commands/taxon_list:

taxon_list
==========

List available taxonomic-specific marker sets.

Produces a table indicating all taxa for which a :ref:`taxonomic-specific marker set <workflows/taxonomy>`
can be produced with the :ref:`taxon_set <commands/taxon_set>` command. Use ``--rank`` to restrict the list
to a taxonomic rank (e.g., life, domain, phylum, ...).

Arguments
---------

.. literalinclude:: help/taxon_list.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm taxon_list
    checkm taxon_list --rank phylum

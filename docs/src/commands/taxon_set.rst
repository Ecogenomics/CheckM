.. _commands/taxon_set:

taxon_set
=========

Generate a taxonomic-specific marker set.

Writes the marker set for the specified taxon to a marker file, which is used by the
:ref:`analyze <commands/analyze>` and :ref:`qa <commands/qa>` commands
(see :ref:`workflows/taxonomy`). Use :ref:`taxon_list <commands/taxon_list>` to see the available taxa.
The marker file is gzip-compressed, and ``.gz`` is appended to its name if not already present (e.g.,
``bacteria.ms`` is written as ``bacteria.ms.gz``); later commands accept the name with or without ``.gz``.

* ``rank``: taxonomic rank of the desired taxonomic-specific marker set (e.g., domain)
* ``taxon``: taxon of interest (e.g., Bacteria)
* ``marker_file``: output file describing the taxonomic-specific marker set

Arguments
---------

.. literalinclude:: help/taxon_set.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm taxon_set domain Bacteria bacteria.ms
    checkm taxon_set phylum Cyanobacteria cyanobacteria.ms

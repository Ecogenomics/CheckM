.. _commands/taxonomy_wf:

taxonomy_wf
===========

Runs taxon_set, analyze, qa.

Runs the :ref:`taxonomic-specific workflow <workflows/taxonomy>`, which consists of the
:ref:`taxon_set <commands/taxon_set>`, :ref:`analyze <commands/analyze>`, and :ref:`qa <commands/qa>`
commands. All bins are assessed with the marker set of the specified taxon. Use
:ref:`taxon_list <commands/taxon_list>` to see the available taxa. See the individual commands for a
description of their options.

The marker file ``<taxon>.ms`` (e.g. ``Bacteria.ms``) is written to the output folder, and the ``qa`` summary
of bin completeness and contamination (``-o 1``) is printed to the console or to the file given with ``-f``.

Arguments
---------

.. literalinclude:: help/taxonomy_wf.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm taxonomy_wf domain Bacteria ./bins ./output
    checkm taxonomy_wf phylum Cyanobacteria -t 8 -x fa ./bins ./output

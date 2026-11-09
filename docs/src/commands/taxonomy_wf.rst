.. _commands/taxonomy_wf:

taxonomy_wf
===========

Runs taxon_set, analyze, qa.

Runs the :ref:`taxonomic-specific workflow <workflows/taxonomy>`, which consists of the
:ref:`taxon_set <commands/taxon_set>`, :ref:`analyze <commands/analyze>`, and :ref:`qa <commands/qa>`
commands. All bins are assessed with the marker set of the specified taxon. Use
:ref:`taxon_list <commands/taxon_list>` to see the available taxa. See the individual commands for a
description of their options.

The results are written to ``quality_report.tsv`` in the output folder: the tab-separated table produced by
``qa -o 2 --tab_table`` (see :ref:`outputs/qa`). The marker file ``<taxon>.ms.gz`` (e.g. ``Bacteria.ms.gz``)
is also written to the output folder, so :ref:`qa <commands/qa>` can be re-run to produce other tables.

* ``--rm_tmp_files``: remove the ``bins`` and ``storage`` folders once ``quality_report.tsv`` has been
  written. ``qa`` and the plotting commands can then no longer be run on the output folder.

.. note::

   CheckM <= 1.2.10 printed the ``qa -o 1`` table to the console or to the file given with ``-f``. The
   ``-f`` and ``--tab_table`` options were removed in CheckM 1.3.0.

Arguments
---------

.. literalinclude:: help/taxonomy_wf.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm taxonomy_wf domain Bacteria ./bins ./output
    checkm taxonomy_wf phylum Cyanobacteria -t 8 -x fa ./bins ./output

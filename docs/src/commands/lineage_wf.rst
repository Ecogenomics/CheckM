.. _commands/lineage_wf:

lineage_wf
==========

Runs tree, lineage_set, analyze, qa.

Runs the :ref:`lineage-specific workflow <workflows/lineage>`, which consists of the
:ref:`tree <commands/tree>`, :ref:`lineage_set <commands/lineage_set>`, :ref:`analyze <commands/analyze>`,
and :ref:`qa <commands/qa>` commands. This is the recommended workflow for assessing the completeness and
contamination of genome bins. See the individual commands for a description of their options.

The results are written to ``quality_report.tsv`` in the output folder: the tab-separated table produced by
``qa -o 2 --tab_table`` (see :ref:`outputs/qa`). The marker file ``lineage.ms.gz`` is also written to the
output folder, and the :ref:`qa <commands/qa>` command can be re-run on the output folder to produce other
tables:

.. code-block:: bash

    checkm qa -o 1 ./output/lineage.ms ./output

* ``--rm_tmp_files``: remove the ``bins`` and ``storage`` folders once ``quality_report.tsv`` has been
  written. This saves disk space, but ``qa`` and the plotting commands can then no longer be run on the
  output folder. A folder containing the input genomes is never removed.
* ``--gplace_max_memory``: limit the memory used to place bins in the reference genome tree (see
  :ref:`commands/tree` and :ref:`installing/system_requirements`).

.. note::

   CheckM <= 1.2.10 printed the ``qa -o 1`` table to the console (or to the file given with ``-f``, as a
   tab-separated table with ``--tab_table``). These options, and ``--reduced_tree`` and
   ``--pplacer_threads``, were removed in CheckM 1.3.0.

Arguments
---------

.. literalinclude:: help/lineage_wf.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm lineage_wf -t 8 -x fa ./bins ./output
    checkm lineage_wf -t 8 -x fa --rm_tmp_files ./bins ./output
    checkm lineage_wf -t 8 -x fa --gplace_max_memory 8G ./bins ./output

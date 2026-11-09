.. _commands/lineage_wf:

lineage_wf
==========

Runs tree, lineage_set, analyze, qa.

Runs the :ref:`lineage-specific workflow <workflows/lineage>`, which consists of the
:ref:`tree <commands/tree>`, :ref:`lineage_set <commands/lineage_set>`, :ref:`analyze <commands/analyze>`,
and :ref:`qa <commands/qa>` commands. This is the recommended workflow for assessing the completeness and
contamination of genome bins. See the individual commands for a description of their options.

The marker file ``lineage.ms`` is written to the output folder, and the ``qa`` summary of bin completeness and
contamination (``-o 1``) is printed to the console or to the file given with ``-f``. The
:ref:`qa <commands/qa>` command can be re-run on the output folder to produce other tables:

.. code-block:: bash

    checkm qa -o 2 ./output/lineage.ms ./output

If you are on a machine with < 40 GB of memory, use ``--reduced_tree`` (see
:ref:`installing/system_requirements`).

Arguments
---------

.. literalinclude:: help/lineage_wf.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm lineage_wf -t 8 -x fa ./bins ./output
    checkm lineage_wf -t 8 --reduced_tree --tab_table -f checkm_results.tsv ./bins ./output

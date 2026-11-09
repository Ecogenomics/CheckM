.. _commands/qa:

qa
==

Assess bins for contamination and completeness.

Estimates the completeness, contamination, and strain heterogeneity of each bin from the marker genes
identified by the :ref:`analyze <commands/analyze>` command, and reports the results in one of several
tables. The ``qa`` command must be given the marker file and output folder used with ``analyze``. The
:ref:`lineage_wf <commands/lineage_wf>` and :ref:`taxonomy_wf <commands/taxonomy_wf>` commands write their
marker file to the output folder (``lineage.ms`` or ``<taxon>.ms``), so ``qa`` can be re-run on the output
of a workflow to produce other tables.

The ``-o, --out_format`` option specifies the desired output (see :ref:`outputs/qa_formats`):

1. summary of bin completeness, contamination, and strain heterogeneity
2. extended summary of bin quality (includes GC, genome size, coding density, ...)
3. summary of bin quality for increasingly basal lineage-specific marker sets
4. list of marker genes for each bin along with the number of times each marker was identified
5. list of bin ID, marker gene ID, and called gene ID for each identified marker gene
6. list of marker genes present multiple times in a bin
7. list of marker genes present multiple times on the same scaffold
8. list indicating the position of each marker gene within a bin
9. FASTA file of marker genes identified in each bin

The statistics reported in these tables are described in :ref:`outputs/qa`.

Other options:

* ``--exclude_markers``: file specifying markers to exclude from marker sets. Each marker to exclude should
  be listed on a separate line of the file.
* ``--individual_markers``: treat markers as independent (i.e., ignore the co-located set structure).
* ``--skip_adj_correction``: do not combine hits to adjacent genes (see note below).
* ``--skip_pseudogene_correction``: do not filter hits that are likely to be pseudogenes. By default, hits
  with an alignment shorter than 30% of the length of the marker HMM are ignored.
* ``--aai_strain``: amino acid identity (AAI) threshold used to identify strain heterogeneity.
* ``-a, --alignment_file``: produce a file showing the alignment of multi-copy genes and their AAI, which can
  be used to further assess strain heterogeneity.
* ``--ignore_thresholds``: ignore model-specific score thresholds.
* ``-e, --e_value``: e-value cut-off (not used if model-specific thresholds are specified).
* ``-l, --length``: minimum alignment length as a fraction of the marker HMM length (not used if
  model-specific thresholds are specified).
* ``-c, --coverage_file``: file containing the coverage of each sequence; coverage information is appended
  to table type 2 when this file is provided (see :ref:`coverage <commands/coverage>`).
* ``-f, --file``: print results to file instead of the console.
* ``--tab_table``: for tabular outputs, print a tab-separated values table instead of a table formatted for
  console output.

.. note::

   Adjacent called genes matching the same marker gene may indicate a true duplication event, a gene calling
   error, or an assembly error. Unless ``--skip_adj_correction`` is used, hits to the same marker gene on a
   run of consecutive genes of a scaffold are combined into a single hit, so they are not treated as
   contamination. The IDs of the combined genes are joined with a pair of ampersands (``&&``). Output
   format 7 lists marker genes present multiple times on the same scaffold, which is useful for identifying
   true gene duplication events, gene calling errors, or assembly errors.

Arguments
---------

.. literalinclude:: help/qa.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm qa ./output/lineage.ms ./output
    checkm qa -o 2 --tab_table -f qa_extended.tsv ./output/lineage.ms ./output

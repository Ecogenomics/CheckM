.. _commands/outliers:

outliers
========

[Experimental] Identify outliers in bins relative to reference distributions.

Produces a table indicating sequences within genome bins that are outliers in either GC, tetranucleotide
distance (TD), or coding density (CD) space relative to the expected distribution of these genomic
statistics. The expected distribution was pre-calculated from reference genomes and the percentile used to
identify outliers can be set with the ``-d`` flag (default = 95). Use ``-r`` to report sequences that are
outliers in ``any`` (default) or ``all`` reference distributions.

This command requires the output folder of the :ref:`analyze <commands/analyze>`,
:ref:`lineage_wf <commands/lineage_wf>`, or :ref:`taxonomy_wf <commands/taxonomy_wf>` command
(``results_dir``) and a file indicating the tetranucleotide signature of all sequences within the genome
bins. This file can be created with the :ref:`tetra <commands/tetra>` command. The output file can be passed
to the :ref:`modify <commands/modify>` command to remove the outlying sequences from a bin.

Arguments
---------

.. literalinclude:: help/outliers.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm tetra seqs.fna tetra.tsv
    checkm outliers ./output ./bins tetra.tsv outliers.tsv

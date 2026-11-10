.. _quick_start:

Quick start
===========

Typical workflow
----------------

CheckM works on a directory of genome bins in FASTA format. By default, CheckM assumes genomes consist of
contigs/scaffolds in nucleotide space and that the files to process end with the extension ``fna``. You can
specify a different extension with the ``-x`` flag. CheckM calls genes internally using
`Prodigal <https://github.com/hyattpd/Prodigal>`_, taking care to identify genes with recoded stop codons.
You can call genes externally and provide CheckM with FASTA files containing genes in amino acid space. To
specify this, use the ``--genes`` flag. Again, you may need to change the extension CheckM looks for (e.g.,
``-x faa``).

Instead of a directory, genomes can also be specified with a tab-separated file listing a genome ID and the
path to the genome FASTA file on each line (and, with ``--genes``, the path to the protein FASTA file as a
third column). Bin IDs are always taken from the FASTA file names with the extension removed.

CheckM consists of a series of commands in order to support a number of different analyses and workflows
(see :ref:`commands`). If you are in a rush to get started, the standard workflow for CheckM is:

.. code-block:: bash

    checkm lineage_wf <bin folder> <output folder>

For a full list of options, run ``checkm lineage_wf -h`` (see :ref:`commands/lineage_wf`). To speed up
processing, use the ``-t`` flag to specify the desired number of threads. Placing genomes in the reference
genome tree requires approximately 18 GB of memory with 64 threads; on a machine with less memory, use
``--gplace_max_memory`` (see :ref:`installing/system_requirements`).

The results are written to the tab-separated file ``quality_report.tsv`` in the output folder (see
:ref:`outputs/qa`). Use ``--rm_tmp_files`` to remove the intermediate ``bins`` and ``storage`` folders once the
workflow has finished, which saves disk space but means ``qa`` and the plotting commands can no longer be run
on the output folder (see :ref:`outputs/files`).

After performing a lineage (``checkm lineage_wf``) or taxonomy (``checkm taxonomy_wf``) workflow, you can
re-run the :ref:`qa <commands/qa>` command to produce a number of different output tables. The ``qa``
command requires a marker file to run. This file is produced during the workflow and is stored in the
CheckM output directory specified (``lineage.ms.gz`` for ``lineage_wf`` or ``<taxon>.ms.gz`` for
``taxonomy_wf``; the name can be given without ``.gz``). Plots can be produced with the
:ref:`plotting commands <commands/plots>`.


Example usage
-------------

Assume you have putative genomes in the directory ``/example_dir/bins`` with ``fa`` as the file extension
and want to store the CheckM results in ``/example_dir/checkm``. To process these genomes with 8 threads,
simply run:

.. code-block:: bash

    checkm lineage_wf -t 8 -x fa /example_dir/bins /example_dir/checkm

Or, to process files of called genes in amino acid space which have the extension ``faa``, use:

.. code-block:: bash

    checkm lineage_wf --genes -t 8 -x faa <bin folder> <output folder>

The results are written to ``/example_dir/checkm/quality_report.tsv``. The completeness, contamination, and
other statistics reported by CheckM are described in :ref:`outputs/qa`.

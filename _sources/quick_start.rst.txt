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
processing, use the ``-t`` flag to specify the desired number of threads. If you are on a machine with
< 40 GB of memory, the ``--reduced_tree`` flag can be used, which reduces the memory requirements to less
than 16 GB (see :ref:`installing/system_requirements`).

After performing a lineage (``checkm lineage_wf``) or taxonomy (``checkm taxonomy_wf``) workflow, you can
re-run the :ref:`qa <commands/qa>` command to produce a number of different output tables. The ``qa``
command requires a marker file to run. This file is produced during the workflow and is stored in the
CheckM output directory specified (``lineage.ms`` for ``lineage_wf`` or ``<taxon>.ms`` for
``taxonomy_wf``). Plots can be produced with the :ref:`plotting commands <commands/plots>`.


Example usage
-------------

Assume you have putative genomes in the directory ``/home/donovan/bins`` with ``fa`` as the file extension
and want to store the CheckM results in ``/home/donovan/checkm``. To process these genomes with 8 threads,
simply run:

.. code-block:: bash

    checkm lineage_wf -t 8 -x fa /home/donovan/bins /home/donovan/checkm

Or, to process files of called genes in amino acid space which have the extension ``faa``, use:

.. code-block:: bash

    checkm lineage_wf --genes -t 8 -x faa <bin folder> <output folder>

To write the results to a tab-separated file instead of printing a table to the console, add
``--tab_table -f checkm_results.tsv``. The completeness, contamination, and other statistics reported by
CheckM are described in :ref:`outputs/qa`.

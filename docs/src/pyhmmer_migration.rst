.. _pyhmmer_migration:

Migration to pyhmmer
====================

CheckM >= 1.4.0 runs `HMMER <http://hmmer.org/>`_ through `pyhmmer <https://pyhmmer.readthedocs.io>`_, a
Python library that includes HMMER 3.4. CheckM <= 1.3.2 ran the HMMER programs (``hmmsearch``, ``hmmalign``,
``hmmfetch`` and ``nhmmer``), which had to be installed separately.

Completeness, contamination and strain heterogeneity estimates are unchanged. The search for marker genes
takes less time, but needs more memory.


What you need to do
-------------------

* **Installing:** nothing. pyhmmer (>= 0.12.0) is installed with CheckM by pip, Bioconda and Docker (see
  :ref:`installing`). It is available for Linux and macOS. The HMMER programs are no longer used, so they do
  not need to be on your system path.
* **gPlace:** CheckM 1.4.0 also requires gPlace >= 0.1.2 (CheckM 1.3.x required >= 0.1.0). Update it with
  ``pip install --upgrade gtdb-gplace`` in the environment where gPlace is installed.
* **Commands:** no arguments were added, removed or changed.
* **Scripts that read CheckM's output:** check the changes to the output files below. The names and formats
  of the files are unchanged.


Changes to results and output files
-----------------------------------

CheckM 1.4.0 was compared with CheckM 1.3.2 on the two sets of test genomes described in
:ref:`performance`:

* ``quality_report.tsv`` of ``lineage_wf`` is identical for the 128 genomes and the 1,000 genomes.
* ``quality_report.tsv`` of ``taxonomy_wf`` is identical for the 128 genomes.
* ``ssu_summary.tsv`` and ``ssu.fna`` of ``ssu_finder`` are identical for the 128 genomes.

pyhmmer was also compared with the programs of HMMER 3.4 on the 128 genomes. The output of ``hmmfetch`` and
``hmmalign`` is identical. The following files differ:

* **HMMER tables** (e.g., ``bins/<bin id>/hmmer.analyze.txt.gz``): E-values that ``hmmsearch`` writes as
  ``0`` can be written as very small numbers (e.g., ``5.2e-313``). This affected 177 of the 1,125,861 hits
  to the 128 genomes. All other values are identical.
* **Reports written with** ``--ali`` (e.g., ``bins/<bin id>/hmmer.analyze.ali.txt.gz``): pyhmmer does not
  write this report, so CheckM writes it in the layout of the output of ``hmmsearch``. It has the same hits
  and alignments, but not the header and search statistics of ``hmmsearch``. The spacing of some column
  headers differs, and E-values of ``0`` can be very small numbers, as above.
* **Hit tables of** ``ssu_finder``: CheckM writes these tables. The hits are those of ``nhmmer`` run with a
  single thread (``--cpu 0``).


Performance
-----------

``lineage_wf`` was run with 64 threads, as described in :ref:`performance`. CheckM 1.4.0 used pyhmmer 0.12.3
and CheckM 1.3.2 used the programs of HMMER 3.4.

.. list-table:: Wall time
   :header-rows: 1
   :widths: 46 12 18 18

   * - Step
     - Genomes
     - CheckM 1.3.2
     - CheckM 1.4.0
   * - Calling genes and searching for the 43 marker genes of the reference genome tree (``tree``)
     - 128
     - 75 s
     - 75 s
   * - Searching for the lineage-specific marker genes (``analyze``)
     - 128
     - 197 s
     - 154 s
   * - ``lineage_wf``
     - 128
     - 299 s
     - 254 s
   * - Calling genes and searching for the 43 marker genes of the reference genome tree (``tree``)
     - 1,000
     - 254 s
     - 239 s
   * - Searching for the lineage-specific marker genes (``analyze``)
     - 1,000
     - 780 s
     - 693 s
   * - ``lineage_wf``
     - 1,000
     - 1,120 s
     - 1,018 s

* **Time:** searching for the lineage-specific marker genes is the longest step of ``lineage_wf``, and takes
  22% (128 genomes) and 11% (1,000 genomes) less time with pyhmmer. Most of the time of the ``tree`` step is
  spent calling genes with Prodigal, so it changes little. ``lineage_wf`` takes 15% and 9% less time.
* **Memory:** with 1,000 genomes, searching for the lineage-specific marker genes needed 5.9 GiB of memory
  with pyhmmer and 3.1 GiB with ``hmmsearch``. This memory increases with the number of threads. It is less
  than the memory used to place genomes in the reference genome tree (13.7 GiB), so the peak memory of
  ``lineage_wf`` is not changed by pyhmmer. It matters if the memory of gPlace is limited with
  ``--gplace_max_memory``: ``lineage_wf`` with ``--gplace_max_memory 6G`` then needed 6.3 GiB with CheckM
  1.4.0 and 5.3 GiB with CheckM 1.3.2.
* ``ssu_finder``: took 15 s for the 128 genomes instead of 10 s. pyhmmer searches a sequence file with one
  thread for each HMM, so CheckM searches parts of the file in separate processes.

The lower peak memory of ``lineage_wf`` in CheckM 1.4.0 without ``--gplace_max_memory`` (13.7 GiB instead of
about 17 GiB for 1,000 genomes) is due to gPlace 0.1.2, not pyhmmer.


Reference
---------

Larralde M, Zeller G. 2023. PyHMMER: a Python library binding to HMMER for efficient sequence analysis.
Bioinformatics 39: btad214.

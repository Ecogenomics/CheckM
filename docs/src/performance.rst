.. _performance:

Performance
===========

The run time, memory and disk space used by ``lineage_wf`` with 64 threads were measured for CheckM 1.2.7,
1.2.10, and >=1.3.0 on two sets of genomes randomly selected from the non-isolate
(metagenome-assembled and single-cell) genomes in `GlobDB <https://globdb.org/>`_ release 232:

* 128 genomes (112 bacteria and 16 archaea)
* 1,000 genomes (875 bacteria and 125 archaea)

CheckM >= 1.3.0 places genomes in the reference genome tree with `gPlace <https://github.com/cmc-aau/gPlace>`_
instead of pplacer. CheckM >= 1.4.0 runs HMMER through `pyhmmer <https://pyhmmer.readthedocs.io>`_ instead of
the HMMER programs. CheckM 1.2.10 and all versions >= 1.3.0 gave identical marker lineages, completeness,
contamination and strain heterogeneity for all genomes. CheckM 1.2.7 differed for 1 of the 128 genomes and
6 of the 1,000 genomes due to bug fixes in CheckM 1.2.8 (see :ref:`quality_estimate_changes`).


Results
-------

All runs used 64 threads.

**128 genomes**

.. list-table::
   :header-rows: 1
   :widths: 16 12 18 18 14 16

   * - CheckM version
     - Wall time
     - Peak memory
     - Placement memory
     - Output size
     - Release date
   * - 1.2.7
     - 9.7 min
     - 43.1 GiB
     - 43.0 GiB
     - 0.62 GB
     - 2026-11-08
   * - 1.2.10
     - 8.5 min
     - 43.2 GiB
     - 43.1 GiB
     - 0.62 GB
     - 2026-11-09
   * - 1.3.0
     - 5.0 min
     - 16.7 ± 1.2 GiB
     - 16.6 ± 1.2 GiB
     - 0.18 GB
     - 2026-11-10
   * - 1.3.1
     - 4.9 min
     - 16.7 ± 1.2 GiB
     - 16.6 ± 1.2 GiB
     - 0.17 GB
     - 2026-11-19
   * - 1.3.2
     - 5.0 min
     - 16.7 ± 1.2 GiB
     - 16.6 ± 1.2 GiB
     - 0.17 GB
     - 2026-11-19
   * - 1.4.0
     - 4.2 min
     - 13.6 GiB
     - 13.5 GiB
     - 0.17 GB
     - 2026-11-23

**1,000 genomes**

.. list-table::
   :header-rows: 1
   :widths: 16 12 18 18 14 16

   * - CheckM version
     - Wall time
     - Peak memory
     - Placement memory
     - Output size
     - Release date
   * - 1.2.7
     - 41.2 min
     - 46.8 GiB
     - 46.6 GiB
     - 4.46 GB
     - 2026-11-08
   * - 1.2.10
     - 24.3 min
     - 46.7 GiB
     - 46.6 GiB
     - 4.43 GB
     - 2026-11-09
   * - 1.3.0
     - 19.6 min
     - 17.1 ± 1.4 GiB
     - 17.0 ± 1.4 GiB
     - 1.23 GB
     - 2026-11-10
   * - 1.3.1
     - 19.0 min
     - 17.1 ± 1.4 GiB
     - 17.0 ± 1.4 GiB
     - 1.21 GB
     - 2026-11-19
   * - 1.3.2
     - 18.7 min
     - 17.1 ± 1.4 GiB
     - 17.0 ± 1.4 GiB
     - 1.21 GB
     - 2026-11-19
   * - 1.4.0
     - 17.0 min
     - 13.7 GiB
     - 13.6 GiB
     - 1.21 GB
     - 2026-11-23

**1,000 genomes with** ``--gplace_max_memory 6G``

.. list-table::
   :header-rows: 1
   :widths: 16 12 18 18 14 16

   * - CheckM version
     - Wall time
     - Peak memory
     - Placement memory
     - Output size
     - Release date
   * - 1.3.1
     - 19.0 min
     - 8.3 ± 0.0 GiB
     - 4.9 ± 1.1 GiB
     - 1.21 GB
     - 2026-11-19
   * - 1.3.2
     - 19.3 min
     - 5.7 ± 0.5 GiB
     - 4.9 ± 1.1 GiB
     - 1.21 GB
     - 2026-11-19
   * - 1.4.0
     - 17.1 min
     - 6.3 GiB
     - 4.2 GiB
     - 1.21 GB
     - 2026-11-23

* **Peak memory** is the peak memory of the whole workflow.
* **Placement memory** is the peak memory of the step that places genomes in the reference genome tree
  (pplacer for CheckM <= 1.2.10, gPlace for CheckM >= 1.3.0). This step sets the peak memory of the workflow,
  except when gPlace's memory is limited with ``--gplace_max_memory``: the peak is then usually reached while
  calling genes with 64 threads at the start of the workflow (8.3 GiB for CheckM 1.3.1). CheckM 1.3.2 limits
  the memory used to call genes (see :ref:`changelog`), which lowers this to 5.3 GiB. In the run of CheckM
  1.4.0, the peak of 6.3 GiB was also reached while calling genes.
* **Mean ± standard deviation:** the memory used by gPlace 0.1.0 differs between runs on the same genomes,
  although the placements are identical (e.g., 15.6 to 18.4 GiB in 5 runs with 1,000 genomes). For CheckM
  1.3.0, 1.3.1 and 1.3.2, which were measured with gPlace 0.1.0, the memory of gPlace is therefore given as
  the mean ± standard deviation of 5 runs, and the peak memory is calculated for each of these runs (see
  Methods). These versions give gPlace the same input and so have the same placement memory. CheckM 1.4.0
  requires gPlace >= 0.1.2, whose memory differed by less than 0.05 GiB between 5 runs, so its values are
  from a single run of the workflow. Values for pplacer are also single runs.
* **Output size** is the size of the output folder. Most output files are gzip-compressed in CheckM >= 1.3.0
  (see :ref:`outputs/files`).
* **Release date** is the date the version was released (see :ref:`changelog`).

Placing genomes in the tree took 3.7 min (128 genomes) and 5.5 min (1,000 genomes) with pplacer, and 12 to
14 s and 17 to 19 s with gPlace. Limiting gPlace to 6 GiB of memory gave identical results. gPlace then
places genomes tile by tile, which took about 1 min for 1,000 genomes with gPlace 0.1.0 and 19 s with
gPlace 0.1.2. gPlace 0.1.0 could exceed this limit (up to 6.2 GiB); gPlace 0.1.2 used 4.2 GiB.

CheckM 1.4.0 searches for marker genes with pyhmmer, which takes less time than the HMMER programs but needs
more memory. With 1,000 genomes, identifying marker genes in the ``analyze`` step took 11.6 min instead of
13.1 min in CheckM 1.3.2, and needed 5.9 GiB of memory instead of 3.1 GiB. This is less than the memory of
gPlace, so it only matters when gPlace's memory is limited with ``--gplace_max_memory``. The lower peak
memory of CheckM 1.4.0 without this limit is due to gPlace 0.1.2. See :ref:`pyhmmer_migration` for details.

With 1,000 genomes, CheckM 1.4.0 takes 30% less time than CheckM 1.2.10 and needs 71% less memory and 73% less
disk space. The memory used by gPlace depends little on the number of genomes, as most of it holds the
reference genome tree. The memory used while calling genes, searching for marker genes and placing genomes
increases with the number of threads.


Methods
-------

* **Machine:** 2x AMD EPYC 7402 24-core processors (48 cores, 96 hardware threads) and 503 GiB of memory,
  running Linux. CheckM was run with 64 threads.
* **Commands:** ``checkm lineage_wf -t 64 -x fna <genome folder> <output folder>``, with
  ``--pplacer_threads 64`` for CheckM <= 1.2.10. Runs were made one at a time.
* **Software:** Prodigal 2.6.3, HMMER 3.4 (CheckM <= 1.3.2), pyhmmer 0.12.3 (CheckM 1.4.0), pplacer
  1.1.alpha20 (CheckM <= 1.2.10), gPlace 0.1.0 (CheckM 1.3.0, 1.3.1 and 1.3.2) and gPlace 0.1.2
  (CheckM 1.4.0).
* **Memory:** the memory of all processes started by CheckM was sampled once a second from the memory cgroup
  of the run. This counts memory shared between processes once, which matters as pplacer and gPlace use
  worker processes that share the reference genome tree; the largest single process (as reported by
  ``/usr/bin/time -v``) used 34.9 GiB with pplacer and 12.0 GiB with gPlace. Memory used by the kernel to
  cache files is not included. The memory in use on the whole machine agreed to within 0.2 GiB.
* **Memory of gPlace 0.1.0:** gPlace was run 5 times for each set of genomes and memory limit, with the command
  and input used by CheckM (``gplace place -j 64``, with ``--max_memory`` when limited), and measured as
  above. For each of these runs, the peak memory of the workflow is the larger of the memory of gPlace plus
  the memory held by CheckM during placement (0.10 GiB for 128 genomes, 0.15 GiB for 1,000 genomes), and the
  peak memory of the other steps of the workflow (a single run of each CheckM version). Wall time and output
  size are from a single run of each workflow. gPlace 0.1.2 was run 5 times in the same way and used 13.5 GiB
  (128 genomes), 13.6 GiB (1,000 genomes) and 4.2 GiB (1,000 genomes, limit of 6 GiB) in every run.

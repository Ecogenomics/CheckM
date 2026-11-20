.. _performance:

Performance
===========

The run time, memory and disk space used by ``lineage_wf`` with 64 threads were measured for CheckM 1.2.7,
1.2.10, 1.3.0, 1.3.1 and 1.3.2 on two sets of genomes randomly selected from the non-isolate
(metagenome-assembled and single-cell) genomes in `GlobDB <https://globdb.org/>`_ release 232:

* 128 genomes (112 bacteria and 16 archaea)
* 1,000 genomes (875 bacteria and 125 archaea)

CheckM >= 1.3.0 places genomes in the reference genome tree with `gPlace <https://github.com/cmc-aau/gPlace>`_
instead of pplacer. CheckM 1.2.10 and all versions >= 1.3.0 gave identical marker lineages, completeness,
contamination and strain heterogeneity for all genomes. CheckM 1.2.7 differed for 1 of the 128 genomes and
6 of the 1,000 genomes due to bug fixes in CheckM 1.2.8 (see :ref:`quality_estimate_changes`).


Results
-------

All runs used 64 threads.

.. list-table::
   :header-rows: 1
   :widths: 24 10 12 18 18 14

   * - CheckM version
     - Genomes
     - Wall time
     - Peak memory
     - Placement memory
     - Output size
   * - 1.2.7
     - 128
     - 9.7 min
     - 43.1 GiB
     - 43.0 GiB
     - 0.62 GB
   * - 1.2.10
     - 128
     - 8.5 min
     - 43.2 GiB
     - 43.1 GiB
     - 0.62 GB
   * - 1.3.0
     - 128
     - 5.0 min
     - 16.7 ± 1.2 GiB
     - 16.6 ± 1.2 GiB
     - 0.18 GB
   * - 1.3.1
     - 128
     - 4.9 min
     - 16.7 ± 1.2 GiB
     - 16.6 ± 1.2 GiB
     - 0.17 GB
   * - 1.3.2
     - 128
     - 5.0 min
     - 16.7 ± 1.2 GiB
     - 16.6 ± 1.2 GiB
     - 0.17 GB
   * - 1.2.7
     - 1,000
     - 41.2 min
     - 46.8 GiB
     - 46.6 GiB
     - 4.46 GB
   * - 1.2.10
     - 1,000
     - 24.3 min
     - 46.7 GiB
     - 46.6 GiB
     - 4.43 GB
   * - 1.3.0
     - 1,000
     - 19.6 min
     - 17.1 ± 1.4 GiB
     - 17.0 ± 1.4 GiB
     - 1.23 GB
   * - 1.3.1
     - 1,000
     - 19.0 min
     - 17.1 ± 1.4 GiB
     - 17.0 ± 1.4 GiB
     - 1.21 GB
   * - 1.3.2
     - 1,000
     - 18.7 min
     - 17.1 ± 1.4 GiB
     - 17.0 ± 1.4 GiB
     - 1.21 GB
   * - 1.3.0 (``--gplace_max_memory 8G``)
     - 1,000
     - 20.3 min
     - 8.5 ± 0.1 GiB
     - 6.4 ± 1.2 GiB
     - 1.23 GB
   * - 1.3.1 (``--gplace_max_memory 6G``)
     - 1,000
     - 19.0 min
     - 8.3 ± 0.0 GiB
     - 4.9 ± 1.1 GiB
     - 1.21 GB
   * - 1.3.2 (``--gplace_max_memory 6G``)
     - 1,000
     - 19.3 min
     - 5.7 ± 0.5 GiB
     - 4.9 ± 1.1 GiB
     - 1.21 GB

* **Peak memory** is the peak memory of the whole workflow.
* **Placement memory** is the peak memory of the step that places genomes in the reference genome tree
  (pplacer for CheckM <= 1.2.10, gPlace for CheckM >= 1.3.0). This step sets the peak memory of the workflow,
  except when gPlace's memory is limited with ``--gplace_max_memory``: the peak is then usually reached while
  calling genes with 64 threads at the start of the workflow (8.4 GiB for CheckM 1.3.0, 8.3 GiB for 1.3.1).
  CheckM 1.3.2 limits the memory used to call genes (see :ref:`changelog`), which lowers this to 5.3 GiB.
* **Mean ± standard deviation:** the memory used by gPlace differs between runs on the same genomes, although
  the placements are identical (e.g., 15.6 to 18.4 GiB in 5 runs with 1,000 genomes). The memory of gPlace is
  therefore given as the mean ± standard deviation of 5 runs, and the peak memory is calculated for each of
  these runs (see Methods). CheckM 1.3.0, 1.3.1 and 1.3.2 give gPlace the same input and so have the same
  placement memory. Values for pplacer are single runs.
* **Output size** is the size of the output folder. Most output files are gzip-compressed in CheckM >= 1.3.0
  (see :ref:`outputs/files`).

Placing genomes in the tree took 3.7 min (128 genomes) and 5.5 min (1,000 genomes) with pplacer, and 13 to
14 s and 17 to 19 s with gPlace. Limiting gPlace to 8 or 6 GiB of memory increased this to 47 to 61 s for
1,000 genomes, as gPlace then places genomes tile by tile, and gave identical results. gPlace can exceed this
limit: it used up to 8.5 GiB with a limit of 8 GiB and up to 6.2 GiB with a limit of 6 GiB.

With 1,000 genomes, CheckM 1.3.2 takes 23% less time than CheckM 1.2.10 and needs 63% less memory and 73% less
disk space. The memory used by gPlace depends little on the number of genomes, as most of it holds the
reference genome tree. The memory used while calling genes and placing genomes increases with the number of
threads.


Methods
-------

* **Machine:** 2x AMD EPYC 7402 24-core processors (48 cores, 96 hardware threads) and 503 GiB of memory,
  running Linux. CheckM was run with 64 threads.
* **Commands:** ``checkm lineage_wf -t 64 -x fna <genome folder> <output folder>``, with
  ``--pplacer_threads 64`` for CheckM <= 1.2.10. Runs were made one at a time.
* **Software:** HMMER 3.4, Prodigal 2.6.3, pplacer 1.1.alpha20 (CheckM <= 1.2.10) and gPlace 0.1.0
  (CheckM >= 1.3.0).
* **Memory:** the memory of all processes started by CheckM was sampled once a second from the memory cgroup
  of the run. This counts memory shared between processes once, which matters as pplacer and gPlace use
  worker processes that share the reference genome tree; the largest single process (as reported by
  ``/usr/bin/time -v``) used 34.9 GiB with pplacer and 12.0 GiB with gPlace. Memory used by the kernel to
  cache files is not included. The memory in use on the whole machine agreed to within 0.2 GiB.
* **Memory of gPlace:** gPlace was run 5 times for each set of genomes and memory limit, with the command
  and input used by CheckM (``gplace place -j 64``, with ``--max_memory`` when limited), and measured as
  above. For each of these runs, the peak memory of the workflow is the larger of the memory of gPlace plus
  the memory held by CheckM during placement (0.10 GiB for 128 genomes, 0.15 GiB for 1,000 genomes), and the
  peak memory of the other steps of the workflow (a single run of each CheckM version). Wall time and output
  size are from a single run of each workflow.

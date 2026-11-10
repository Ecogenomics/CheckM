.. _performance:

Performance
===========

The run time, memory and disk space used by ``lineage_wf`` were measured for CheckM 1.2.7, 1.2.10 and 1.3.0
on two sets of genomes randomly selected from the non-isolate (metagenome-assembled and single-cell) genomes in
`GlobDB <https://globdb.org/>`_ release 232:

* 128 genomes (112 bacteria and 16 archaea)
* 1,000 genomes (875 bacteria and 125 archaea)

CheckM 1.3.0 places genomes in the reference genome tree with `gPlace <https://github.com/cmc-aau/gPlace>`_
instead of pplacer. CheckM 1.3.0 and 1.2.10 gave identical marker lineages, completeness, contamination and
strain heterogeneity for all genomes. CheckM 1.2.7 differed for 1 of the 128 genomes and 6 of the 1,000 genomes
due to bug fixes in CheckM 1.2.8 (see :ref:`quality_estimate_changes`).


Results
-------

.. list-table::
   :header-rows: 1
   :widths: 22 12 14 16 16 20

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
     - 16.8 GiB
     - 16.7 GiB
     - 0.18 GB
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
     - 18.3 GiB
     - 18.2 GiB
     - 1.23 GB
   * - 1.3.0 (``--gplace_max_memory 8G``)
     - 1,000
     - 20.3 min
     - 8.4 GiB
     - 6.2 GiB
     - 1.23 GB

* **Peak memory** is the peak memory of the whole workflow.
* **Placement memory** is the peak memory of the step that places genomes in the reference genome tree
  (pplacer for CheckM <= 1.2.10, gPlace for CheckM 1.3.0). This step sets the peak memory of the workflow,
  except when gPlace's memory is limited with ``--gplace_max_memory``: the peak (8.4 GiB) is then reached
  while calling genes with 64 threads at the start of the workflow.
* **Output size** is the size of the output folder. Most output files are gzip-compressed from CheckM 1.3.0
  (see :ref:`outputs/files`).

Placing genomes in the tree took 3.7 min (128 genomes) and 5.5 min (1,000 genomes) with pplacer, and 14 s and
19 s with gPlace. Limiting gPlace to 8 GiB of memory increased this to 53 s for 1,000 genomes, as gPlace then
places genomes tile by tile, and gave identical results.

With 1,000 genomes, CheckM 1.3.0 takes 19% less time than CheckM 1.2.10 and needs 61% less memory and 72% less
disk space. The memory used by gPlace depends little on the number of genomes, as most of it holds the
reference genome tree. The memory used while calling genes and placing genomes increases with the number of
threads.


Methods
-------

* **Machine:** 2x AMD EPYC 7402 24-core processors (96 threads) and 502 GB of memory, running Linux.
* **Commands:** ``checkm lineage_wf -t 64 -x fna <genome folder> <output folder>``, with
  ``--pplacer_threads 64`` for CheckM <= 1.2.10. Runs were made one at a time.
* **Software:** HMMER 3.4, Prodigal 2.6.3, pplacer 1.1.alpha20 (CheckM <= 1.2.10) and gPlace 0.1.0
  (CheckM 1.3.0).
* **Memory:** the memory of all processes started by CheckM was sampled once a second from the memory cgroup
  of the run. This counts memory shared between processes once, which matters as pplacer and gPlace use
  worker processes that share the reference genome tree; the largest single process (as reported by
  ``/usr/bin/time -v``) used 34.9 GiB with pplacer and 12.0 GiB with gPlace. Memory used by the kernel to
  cache files is not included. The memory in use on the whole machine agreed to within 0.2 GiB.

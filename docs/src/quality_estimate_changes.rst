.. _quality_estimate_changes:

Quality estimate changes
========================

CheckM 1.2.8 fixed two bugs that can change the completeness, contamination and strain heterogeneity
estimates of a genome. The changes are small and affect few genomes. CheckM 1.2.8, 1.2.9, 1.2.10 and 1.3.0
give identical estimates, so the changes below apply to any CheckM >= 1.2.8 compared with CheckM <= 1.2.7.


Bug fixes that change estimates
-------------------------------

* **Pfam same-clan filter.** When two Pfam markers hit the same gene, CheckM keeps only the best hit if the
  two Pfams belong to the same Pfam clan. CheckM <= 1.2.7 treated two Pfams without a clan as members of the
  same clan, so a hit to an unrelated marker could be discarded. With the fix, such a marker can be found in
  a genome (higher completeness) or found in an additional copy (higher contamination).
* **Gene fragments.** Hits to the same marker on a run of consecutive genes (e.g., a gene split into
  several fragments) are now combined into a single hit, as these are most likely one gene split by a gene
  calling or assembly error. CheckM <= 1.2.7 combined such hits in pairs, in an order that depended on hits
  to other markers, so a gene split into three or more fragments could be counted as more than one copy of
  the marker. With the fix, contamination and strain heterogeneity can be lower.

The other fixes in CheckM 1.2.8 (see :ref:`changelog`) did not change the estimates of any test genome.


Test
----

``lineage_wf`` was run with CheckM 1.2.7 and 1.3.0 on two sets of genomes randomly selected from the
non-isolate (metagenome-assembled and single-cell) genomes in `GlobDB <https://globdb.org/>`_ release 232:

* 128 genomes (112 bacteria and 16 archaea)
* 1,000 genomes (875 bacteria and 125 archaea)

These are the genomes used for the :ref:`performance` results. The marker lineage, completeness,
contamination and strain heterogeneity of each genome were compared between the two versions.


Results
-------

The estimates changed for 1 of the 128 genomes and 6 of the 1,000 genomes (0.6%). The marker lineage did not
change for any genome, and no change was larger than 0.18 percentage points for completeness or
contamination. No genome crossed a MIMAG quality threshold (completeness of 50% or 90%, contamination of 5%
or 10%).

.. list-table::
   :header-rows: 1
   :widths: 13 11 24 13 13 14 12

   * - Genome
     - Genome set
     - Marker lineage
     - Completeness
     - Contamination
     - Strain heterogeneity
     - Cause
   * - GCA_035327525
     - 128
     - k__Bacteria (UID203)
     - 75.99 → 76.15
     -
     -
     - Clan filter
   * - GCA_018701055
     - 1,000
     - c__Gammaproteobacteria (UID4267)
     - 74.87 → 74.90
     -
     -
     - Clan filter
   * - GCA_030455605
     - 1,000
     - g__Lactobacillus (UID377)
     -
     - 0.41 → 0.51
     -
     - Clan filter
   * - GCA_045561055
     - 1,000
     - c__Gammaproteobacteria (UID4267)
     -
     - 1.58 → 1.62
     - 16.67 → 14.29
     - Clan filter
   * - GCA_964384205
     - 1,000
     - k__Bacteria (UID2495)
     -
     - 5.54 → 5.64
     -
     - Clan filter
   * - GCA_963783575
     - 1,000
     - k__Bacteria (UID203)
     -
     - 3.93 → 3.77
     -
     - Gene fragments
   * - GCA_965202795
     - 1,000
     - k__Bacteria (UID3187)
     -
     - 0.18 → 0.00
     -
     - Gene fragments

Values are CheckM 1.2.7 → CheckM 1.3.0; empty cells are unchanged.

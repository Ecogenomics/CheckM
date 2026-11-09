.. _outputs/qa:

Reported statistics
===================

The following statistics are reported by the :ref:`qa <commands/qa>` command. Output format 1
(``-o 1``, the default and the output of the workflows) reports the bin ID, marker set information, marker
counts, completeness, contamination, and strain heterogeneity. Output format 2 (``-o 2``) additionally reports
the genome statistics. See :ref:`outputs/qa_formats` for the other output formats.


Marker set and quality estimates
--------------------------------

* **Bin Id**: unique identifier of genome bin (derived from the input FASTA file name).
* **Marker lineage**: indicates the taxonomic rank of the lineage-specific marker set used to estimate genome
  completeness, contamination, and strain heterogeneity. More detailed information about the placement of a
  genome within the reference genome tree can be obtained with the :ref:`tree_qa <commands/tree_qa>`
  command. The UID indicates the branch within the reference tree used to infer the marker set applied to
  estimate the bin's quality. For :ref:`custom marker genes <workflows/custom>`, the marker lineage is
  ``N/A``.
* **# genomes**: number of reference genomes used to infer the lineage-specific marker set (``-1`` for custom
  marker genes).
* **# markers**: number of marker genes within the inferred lineage-specific marker set.
* **# marker sets**: number of co-located marker sets within the inferred lineage-specific marker set.
* **0, 1, 2, 3, 4, 5+**: number of marker genes identified 0, 1, 2, 3, 4, or 5 or more times.
* **Completeness**: estimated completeness of genome as determined from the presence/absence of marker genes
  and the expected collocalization of these genes (see Methods in the
  `CheckM manuscript <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/>`_ for details).
* **Contamination**: estimated contamination of genome as determined by the presence of multi-copy marker
  genes and the expected collocalization of these genes (see Methods in the
  `CheckM manuscript <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/>`_ for details).
* **Strain heterogeneity**: estimated strain heterogeneity as determined from the number of multi-copy marker
  pairs which exceed a specified amino acid identity threshold (default = 90%; see ``--aai_strain``). High
  strain heterogeneity suggests the majority of reported contamination is from one or more closely related
  organisms (i.e. potentially the same species), while low strain heterogeneity suggests the majority of
  contamination is from more phylogenetically diverse sources (see Methods in the
  `CheckM manuscript <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/>`_ for more details).


Genome statistics
-----------------

The following statistics are also reported in output format 2:

* **Genome size (bp)**: number of nucleotides (including unknowns specified by N's) in the genome.
* **# ambiguous bases**: number of ambiguous (N's) bases in the genome.
* **# scaffolds**: number of scaffolds within the genome.
* **# contigs**: number of contigs within the genome as determined by splitting scaffolds at any position
  consisting of more than 10 consecutive ambiguous bases.
* **N50 (scaffolds)**: N50 statistic as calculated over all scaffolds.
* **N50 (contigs)**: N50 statistic as calculated over all contigs.
* **Mean scaffold length (bp)**: mean length of scaffolds within the genome.
* **Mean contig length (bp)**: mean length of contigs within the genome.
* **Longest scaffold (bp)**: the longest scaffold within the genome.
* **Longest contig (bp)**: the longest contig within the genome.
* **GC**: number of G/C nucleotides relative to all A, C, G, and T nucleotides in the genome.
* **GC std (scaffolds > 1kbp)**: standard deviation of the GC content of scaffolds longer than 1 kbp.
* **Coding density**: the number of nucleotides within a coding sequence (CDS) relative to all nucleotides
  in the genome.
* **Translation table**: indicates which
  `genetic code <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`_ was used to translate
  nucleotides into amino acids (see :ref:`tree <commands/tree>`).
* **# predicted genes**: number of predicted coding sequences (CDS) within the genome as determined using
  `Prodigal <https://github.com/hyattpd/Prodigal>`_.
* **Coverage (<BAM>)**, **Coverage std (<BAM>)**: mean and standard deviation of the coverage of the bin in
  each BAM file. Only reported when a coverage file produced by the :ref:`coverage <commands/coverage>`
  command is given with ``-c``.

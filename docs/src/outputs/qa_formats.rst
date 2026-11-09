.. _outputs/qa_formats:

qa output formats
=================

The :ref:`qa <commands/qa>` command can produce the following outputs, selected with ``-o, --out_format``.
Results are printed to the console, or to a file with ``-f``. Formats 1, 2, and 3 are printed as a table
formatted for console output unless ``--tab_table`` is given, in which case a tab-separated values table is
produced. Formats 4 to 8 are always printed as tab-separated text, and format 9 is a FASTA file. The statistics are described
in :ref:`outputs/qa`.

1. **Summary of bin completeness, contamination, and strain heterogeneity**: Bin Id, Marker lineage,
   # genomes, # markers, # marker sets, 0-5+, Completeness, Contamination, and Strain heterogeneity.
2. **Extended summary of bin quality**: the statistics of format 1 along with genome statistics (e.g., GC,
   genome size, coding density). Coverage information is appended when a coverage file is given with ``-c``.
3. **Summary of bin quality for increasingly basal lineage-specific marker sets**: as format 1, with each
   lineage-specific marker set considered for the bin on a separate row.

   * Node Id: unique ID of the internal node in the genome tree from which lineage-specific markers were
     inferred

4. **List of marker genes for each bin along with the number of times each marker was identified**:

   * Node Id: unique ID of the internal node in the genome tree from which lineage-specific markers were
     inferred
   * Marker lineage: indicates the lineage used for inferring the marker set
   * *Useful for identifying lineage-specific gene loss or duplication.*

5. **List of bin ID, marker gene ID, and called gene ID for each identified marker gene.**
6. **List of marker genes present multiple times in a bin.**
7. **List of marker genes present multiple times on the same scaffold.**

   * *Useful for identifying true gene duplication events, gene calling errors, or assembly errors. See the
     note on adjacent marker genes in* :ref:`qa <commands/qa>`.

8. **List indicating the position of each marker gene within a bin.**
9. **FASTA file of marker genes identified in each bin.** The header of each sequence gives the bin ID, the
   contig, the gene (``geneId``, ``start``, ``end``, ``strand``, ``protlen``), and the marker hit
   (``marker``, ``mstart``, ``mend``). Do not use ``--tab_table`` with this format, as only the FASTA
   output is supported.

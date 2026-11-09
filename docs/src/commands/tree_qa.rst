.. _commands/tree_qa:

tree_qa
=======

Assess phylogenetic markers found in each bin.

This optional step of the :ref:`lineage-specific workflow <workflows/lineage>` reports the placement of
each bin in the reference genome tree. It must be run on the output folder of the
:ref:`tree <commands/tree>` command. The ``-o, --out_format`` option specifies the desired output:

1. brief summary of genome tree placement indicating the number of unique phylogenetically informative
   markers found, the number of markers found multiple times, and a taxon string indicating the placement
   of each bin within the genome tree
2. detailed summary of genome tree placement giving a more detailed indication of where each bin is within
   the genome tree, general characteristics about each bin (e.g., GC, genome size, coding density), and
   general characteristics about all reference genomes descendant from the parental node of each bin
   (e.g., mean and standard deviation of GC)
3. genome tree in Newick format decorated with IMG genome IDs which can be used to examine the
   phylogenetic neighbours of each bin
4. genome tree in Newick format decorated with taxonomy strings which can be used to examine the
   phylogenetic neighbours of each bin
5. multiple sequence alignment of reference genomes and bins which can be used to infer a *de novo* genome
   tree

Use ``-f`` to write results to a file instead of the console, and ``--tab_table`` to write tabular outputs
as a tab-separated values table instead of a table formatted for console output.

Arguments
---------

.. literalinclude:: help/tree_qa.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm tree_qa ./output
    checkm tree_qa -o 2 --tab_table -f tree_qa.tsv ./output

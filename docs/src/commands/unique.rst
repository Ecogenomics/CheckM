.. _commands/unique:

unique
======

Ensure no sequences are assigned to multiple bins.

Checks each putative genome and ensures no sequence has been assigned to multiple genomes. For most automated
binning methods, the assignment of a sequence to multiple putative genomes would indicate a serious binning
error. In practice, this command verifies that each sequence name is unique across all putative genomes.

Arguments
---------

.. literalinclude:: help/unique.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm unique ./bins

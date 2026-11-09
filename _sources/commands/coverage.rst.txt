.. _commands/coverage:

coverage
========

Calculate coverage of sequences.

Produces coverage profiles for all sequences within a set of genome bins. This command requires sorted and
indexed BAM files (i.e. a ``.bai`` index file must be present for each BAM file) produced with a tool such
as `BWA <https://github.com/lh3/bwa>`_. Coverage profiles are required by the
:ref:`profile <commands/profile>` command and can be added to the extended ``qa`` table (``qa -o 2 -c``; see
:ref:`qa <commands/qa>`).

Arguments
---------

.. literalinclude:: help/coverage.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm coverage ./bins coverage.tsv example_1.bam example_2.bam

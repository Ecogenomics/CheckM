.. _commands/modify:

modify
======

[Experimental] Modify sequences in a bin.

This is an experimental command that allows sequences to be added to (``-a``) or removed from (``-r``) a
genome bin. It is also compatible with the :ref:`outliers <commands/outliers>` command: all sequences within a
putative genome identified as outliers can be removed by passing the output of ``outliers`` with ``-o``
(this cannot be combined with ``-a`` or ``-r``). Caution is warranted here as we have not explored the rate
at which sequences identified as outliers do not belong in a genome bin.

Arguments
---------

.. literalinclude:: help/modify.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm modify -r seq_id1 -r seq_id2 seqs.fna bin.fna new_bin.fna
    checkm modify -o outliers.tsv seqs.fna bin.fna new_bin.fna

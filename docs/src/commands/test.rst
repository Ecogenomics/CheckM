.. _commands/test:

test
====

Run CheckM on a test genome (E. coli) to verify operation.

Runs the *E. coli* K-12 W3110 genome, which is included with the CheckM reference data, through the standard
CheckM pipeline and verifies the resulting output files. This verifies that CheckM, its reference data, and
all third-party dependencies are installed correctly (see :ref:`installing/testing`). The output directory
can be removed once the test has run.

Arguments
---------

.. literalinclude:: help/test.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm test ~/checkm_test_results

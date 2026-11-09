.. _installing/testing:

Testing the installation
========================

A general test of CheckM which verifies all third-party dependencies can be run using the
:ref:`test <commands/test>` command:

.. code-block:: bash

    checkm test ~/checkm_test_results

This runs the *E. coli* K-12 W3110 genome through the standard CheckM pipeline and verifies the resulting
output files. The output directory can be removed once the test has run.

Unit tests
----------

Unit tests are also provided in the ``checkm/test`` directory of the source code. These are designed to aid
in development and are not needed to verify an installation. They use Python's ``unittest`` module and can
be run from the root of a CheckM source checkout. ``CHECKM_DATA_PATH`` must be set to the
:ref:`reference data <installing/reference_data>` before running them:

.. code-block:: bash

    export CHECKM_DATA_PATH=/path/to/my_checkm_data
    python3 -m unittest discover -s checkm/test -t .

.. _installing/reference_data:

Reference data
==============

CheckM relies on a number of precalculated data files which can be downloaded from either:

* `<https://data.ace.uq.edu.au/public/CheckM_databases>`_
* `<https://zenodo.org/record/7401545>`_

The reference data (``checkm_data_2015_01_16.tar.gz``) must be decompressed into a directory. For example:

.. code-block:: bash

    mkdir -p ~/.checkm   # or wherever you want to put it
    cd ~/.checkm
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar -xzf checkm_data_2015_01_16.tar.gz


Setting the location of the reference data
------------------------------------------

The recommended way to tell CheckM where the reference data is located is the ``CHECKM_DATA_PATH``
environment variable:

.. code-block:: bash

    export CHECKM_DATA_PATH=/path/to/my_checkm_data

Alternatively, the location can be stored in CheckM's ``DATA_CONFIG`` file with the
:ref:`data <commands/data>` command:

.. code-block:: bash

    checkm data setRoot <checkm_data_dir>

``checkm data setRoot`` must be able to write to the ``DATA_CONFIG`` file (located in the installed
``checkm`` Python package, unless a different file is given as a second argument) and to the reference data
directory, as it creates a manifest file (``.dmanifest``) in this directory.

CheckM determines the location of the reference data in the following order:

1. the ``CHECKM_DATA_PATH`` environment variable, if it is set (``checkm data setRoot`` has no effect while
   it is set)
2. the directory stored in ``DATA_CONFIG`` by ``checkm data setRoot``
3. ``~/.checkm`` (after printing a message about it)

The location of the reference data used by a run is written to the CheckM log.

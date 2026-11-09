.. _installing/bioconda:

Bioconda
========

A Conda package for CheckM is available on `Bioconda <https://anaconda.org/bioconda/checkm-genome>`_ as
``checkm-genome``. It installs CheckM along with its :ref:`dependencies <installing/dependencies>`:

.. code-block:: bash

    conda create -n checkm -c conda-forge -c bioconda checkm-genome
    conda activate checkm

The Bioconda package has been generously put together and maintained by community members (if this is you,
please let us know so we can acknowledge you here!). It may lag behind the latest release on PyPI.

CheckM >= 1.3.0 requires `gPlace <https://github.com/cmc-aau/gPlace>`_ instead of pplacer. If gPlace is not
installed with the Bioconda package, install it from PyPI (this requires Python >= 3.12):

.. code-block:: bash

    pip3 install gtdb-gplace

Finally, download and configure the :ref:`reference data <installing/reference_data>`.

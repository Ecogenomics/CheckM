.. _installing/bioconda:

Bioconda
========

A Conda package for CheckM is available on `Bioconda <https://anaconda.org/bioconda/checkm-genome>`_ as
``checkm-genome``. It installs CheckM along with HMMER, Prodigal, and pplacer:

.. code-block:: bash

    conda create -n checkm -c conda-forge -c bioconda checkm-genome
    conda activate checkm

The Bioconda package has been generously put together and maintained by community members (if this is you,
please let us know so we can acknowledge you here!). It may lag behind the latest release on PyPI.

Finally, download and configure the :ref:`reference data <installing/reference_data>`.

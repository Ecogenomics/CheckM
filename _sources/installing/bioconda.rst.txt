.. _installing/bioconda:

Bioconda
========

A Conda package for CheckM is available on `Bioconda <https://anaconda.org/bioconda/checkm-genome>`_ as
``checkm-genome``. It installs CheckM along with all of its :ref:`dependencies <installing/dependencies>`,
including gPlace, HMMER and Prodigal.

The package can be installed with either `mamba <https://mamba.readthedocs.io>`_ or
`conda <https://docs.conda.io>`_. We recommend mamba, which resolves dependencies faster than conda.
mamba is included in `Miniforge <https://github.com/conda-forge/miniforge>`_.

With mamba (recommended):

.. code-block:: bash

    mamba create -n checkm -c conda-forge -c bioconda checkm-genome
    mamba activate checkm

With conda:

.. code-block:: bash

    conda create -n checkm -c conda-forge -c bioconda checkm-genome
    conda activate checkm

The Bioconda package has been generously put together and maintained by community members (if this is you,
please let us know so we can acknowledge you here!). It may lag behind the latest release on PyPI.

Finally, download and configure the :ref:`reference data <installing/reference_data>`.

.. _installing/pip:

pip
===

CheckM is a Python 3 program and is available on `PyPI <https://pypi.org/project/checkm-genome/>`_ as
``checkm-genome``:

.. code-block:: bash

    pip3 install checkm-genome

This will install CheckM and all other required Python libraries
(see :ref:`installing/python_dependencies`). If ``pip`` complains about breaking system packages, install
CheckM into a virtual environment:

.. code-block:: bash

    python3 -m venv .venv       # or: uv venv, if you prefer https://pypi.org/project/uv/
    source .venv/bin/activate   # or: whatever the venv command tells you to run
    pip3 install checkm-genome  # or: uv pip install checkm-genome

The :ref:`bioinformatic tool dependencies <installing/dependencies>` (HMMER, Prodigal, and pplacer) need to
be installed separately and placed on your system path. If you are not too concerned about having the
latest versions, your package manager may provide them:

* Debian, Ubuntu, etc.: ``apt install hmmer prodigal pplacer``
* Fedora: ``dnf install hmmer`` (Prodigal and pplacer need to be obtained separately, e.g. from Bioconda)

Alternatively, a Conda environment can be set up with the dependencies from Bioconda and CheckM installed
with ``pip``:

.. code-block:: bash

    conda create -n checkm python=3.9
    conda activate checkm
    conda install -c bioconda numpy matplotlib pysam
    conda install -c bioconda hmmer prodigal pplacer
    pip3 install checkm-genome

Finally, download and configure the :ref:`reference data <installing/reference_data>`.

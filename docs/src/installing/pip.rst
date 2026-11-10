.. _installing/pip:

pip
===

CheckM is a Python program (CheckM >= 1.3.0 requires Python >= 3.11) and is available on `PyPI <https://pypi.org/project/checkm-genome/>`_ as
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

The :ref:`bioinformatic tool dependencies <installing/dependencies>` (HMMER, Prodigal, and gPlace) need to
be installed separately and placed on your system path. gPlace is installed from PyPI and requires
Python >= 3.12:

.. code-block:: bash

    pip3 install gtdb-gplace

If you are not too concerned about having the latest versions, your package manager may provide HMMER and
Prodigal:

* Debian, Ubuntu, etc.: ``apt install hmmer prodigal``
* Fedora: ``dnf install hmmer`` (Prodigal needs to be obtained separately, e.g. from Bioconda)

Alternatively, a Conda environment can be set up with the dependencies from Bioconda and CheckM installed
with ``pip``. This can be done with mamba (recommended) or conda (replace ``mamba`` with ``conda``):

.. code-block:: bash

    mamba create -n checkm -c conda-forge -c bioconda python=3.12 hmmer prodigal
    mamba activate checkm
    pip3 install checkm-genome gtdb-gplace

Finally, download and configure the :ref:`reference data <installing/reference_data>`.

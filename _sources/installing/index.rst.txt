.. _installing:

Installing CheckM
=================

The latest release of CheckM is |release|. CheckM >= 1.3.0 requires Python >= 3.11 (CheckM 1.1.0 to 1.2.10
require Python 3).

CheckM can be installed with :ref:`pip <installing/pip>` or :ref:`Bioconda <installing/bioconda>`, or run
from a :ref:`Docker <installing/docker>` image built with the Dockerfile in the CheckM repository. All
methods require the :ref:`CheckM reference data <installing/reference_data>`. We do not recommend
installing CheckM from the master branch, which may be unstable; please install an official release.

.. toctree::
   :maxdepth: 1

   pip
   bioconda
   docker
   reference_data
   testing


.. _installing/system_requirements:

System requirements
-------------------

CheckM is designed to run on Linux. The limiting requirement for CheckM is memory. Inference of
:ref:`lineage-specific marker sets <workflows/lineage>` places genomes in the reference genome tree with
`gPlace <https://github.com/cmc-aau/gPlace>`_, which requires approximately **12 GB** of memory, most of it for
the reference tree, so it hardly depends on the number of genomes (12.0 GB for 128 genomes and 12.1 GB for
1,000 genomes with 64 threads). On machines with less memory, use ``--gplace_max_memory`` (e.g.,
``--gplace_max_memory 6G``) to have gPlace place genomes tile by tile within this limit, which gives identical
results: 128 genomes can be placed with 6 GB of memory.

CheckM <= 1.2.10 placed genomes with pplacer, which required approximately 40 GB of memory with the full
reference genome tree, or less than 16 GB with the reduced tree (``--reduced_tree``, removed in CheckM 1.3.0).

System requirements are **far more modest** if you plan to make use of
:ref:`taxonomic-specific marker sets <workflows/taxonomy>` or your own
:ref:`custom marker genes <workflows/custom>`, as this bypasses the need to place genomes in the reference
genome tree.

If you plan to process a large number of genomes, see :ref:`faq/large_numbers`.


.. _installing/dependencies:

Bioinformatic tool dependencies
-------------------------------

CheckM requires the following programs to be
`added to your system path <https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path>`_:

* `HMMER <http://hmmer.org/>`_ (>= 3.1b1)
* `Prodigal <https://github.com/hyattpd/Prodigal>`_ (2.60 or >= 2.6.1)

  * the executable must be named ``prodigal`` and not ``prodigal.linux``

* `gPlace <https://github.com/cmc-aau/gPlace>`_ (>= 0.1.0)

  * gPlace requires Python >= 3.12 and can be installed with ``pip install gtdb-gplace``. It can be
    installed in a separate environment as long as the ``gplace`` command is on your system path.
  * gPlace replaced pplacer in CheckM 1.3.0 and gives the same placements as pplacer, bit for bit, using
    less memory and time. CheckM <= 1.2.10 requires `pplacer <https://matsen.fhcrc.org/pplacer/>`_ (>= 1.1)
    and ``guppy`` instead.

These programs are installed automatically in the :ref:`Docker <installing/docker>` image, but must be
installed separately when CheckM is installed with :ref:`pip <installing/pip>`. See
:ref:`installing/bioconda` for the Bioconda package.


.. _installing/python_dependencies:

Python dependencies
-------------------

The following Python libraries are required and are installed automatically by ``pip``:

* numpy >= 1.21.3
* scipy >= 1.7.3
* matplotlib >= 3.5.1
* pysam >= 0.19.0
* dendropy >= 4.5.2
* tqdm >= 4.62
* xopen >= 2.0 (reads and writes the gzip-compressed :ref:`output files <outputs/files>`)
* setuptools >= 81


Running CheckM
--------------

Once CheckM and its reference data are installed, CheckM is ready to run. For a list of CheckM commands,
type:

.. code-block:: bash

    checkm

See :ref:`quick_start` to get started.


.. _installing/upgrade:

Upgrading CheckM
----------------

You can upgrade CheckM through `pip <https://pypi.org/project/checkm-genome/>`_:

.. code-block:: bash

    pip3 install checkm-genome --upgrade --no-deps

Upgrading CheckM v1 does not require new reference data.


.. _installing/kbase:

Web version through KBase
-------------------------

The `CheckM lineage workflow <https://narrative.kbase.us/#catalog/apps/kb_Msuite/run_checkM_lineage_wf>`_ has
been made available through `KBase <https://kbase.us/>`_ for those looking for a web-based solution.

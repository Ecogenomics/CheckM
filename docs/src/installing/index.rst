.. _installing:

Installing CheckM
=================

The latest release of CheckM is |release|. CheckM >= 1.1.0 requires Python 3.

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
:ref:`lineage-specific marker sets <workflows/lineage>` using the full reference genome tree requires
approximately **40 GB** of memory. However, a reduced genome tree (``--reduced_tree``) can also be used to
infer lineage-specific marker sets, which is suitable for machines with as little as **16 GB** of memory. We
recommend using the full tree if possible, though our results suggest that the same lineage-specific marker
set will be selected for the vast majority of genomes regardless of the underlying reference tree.

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

* `pplacer <https://matsen.fhcrc.org/pplacer/>`_ (>= 1.1)

  * ``guppy``, which is part of the pplacer package, must also be on your system path
  * pplacer binaries can be found on the
    `pplacer GitHub <https://github.com/matsen/pplacer/releases/tag/v1.1.alpha17>`_ page

These programs are installed automatically with the :ref:`Bioconda <installing/bioconda>` package and in the
:ref:`Docker <installing/docker>` image, but must be installed separately when CheckM is installed with
:ref:`pip <installing/pip>`.


.. _installing/python_dependencies:

Python dependencies
-------------------

The following Python libraries are required and are installed automatically by ``pip``:

* numpy >= 1.21.3
* scipy >= 1.7.3
* matplotlib >= 3.5.1
* pysam >= 0.19.0
* dendropy >= 4.5.2
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

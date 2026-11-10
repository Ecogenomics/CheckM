******
CheckM
******

CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells,
or metagenomes. It provides robust estimates of genome completeness and contamination by using collocated
sets of genes that are ubiquitous and single-copy within a phylogenetic lineage. Assessment of genome quality
can also be examined using plots depicting key genomic characteristics (e.g., GC, coding density) which
highlight sequences outside the expected distributions of a typical genome. CheckM also provides tools for
identifying genome bins that are likely candidates for merging based on marker set compatibility, similarity
in genomic characteristics, and proximity within a reference genome tree.

CheckM is open source and released under the
`GNU General Public License (Version 3) <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.
The source code is available on `GitHub <https://github.com/Ecogenomics/CheckM>`_.


CheckM v2 has been released
===========================

.. important::

   CheckM v2 has been released and can be found at `<https://github.com/chklovski/CheckM2>`_.
   CheckM v1 (documented here) estimates genome completeness and contamination based on the presence or
   absence of marker genes, i.e., genes that are typically ubiquitous and single copy. In contrast,
   `CheckM v2 <https://github.com/chklovski/CheckM2>`_ uses machine learning models to estimate
   completeness and contamination. Our `assessment <https://www.nature.com/articles/s41592-023-01940-w>`_
   indicates CheckM v2 will be more accurate in general. However, you may consider running both CheckM v1
   and v2 in some cases since these are independent methodologies for estimating genome quality. It can also
   be insightful to explore the marker genes that are missing in your genomes, especially if you have
   multiple related genomes and are interested in systematic loss of an otherwise ubiquitous gene.

   Outside of minor updates to fix compatibility issues with external dependencies, no maintenance is
   intended for CheckM v1.


Latest release
==============

The latest release of CheckM is |release|. See the :ref:`changelog` for the changes in each release.

CheckM >= 1.3.0 requires Python >= 3.11. We do not recommend installing CheckM from the master branch,
which may be unstable. Please install an official release of CheckM (see :ref:`installing`).


Running CheckM
==============

1. Install CheckM and its reference data (:ref:`installing`).

2. Run the lineage-specific workflow on a directory of genomes (:ref:`quick_start`):

   .. code-block:: bash

       checkm lineage_wf -t 8 -x fa <bin folder> <output folder>

3. See :ref:`workflows` for the available workflows, :ref:`commands` for all commands, and
   :ref:`outputs` for a description of the reported statistics. A list of commands is also shown by
   running ``checkm`` (or ``checkm -h``), and help for an individual command by running
   ``checkm <command> -h``.

Information about obtaining improved quality estimates for CPR (Patescibacteria) genomes can be found in
:ref:`workflows/cpr`.

Please report bugs through the `GitHub issues <https://github.com/Ecogenomics/CheckM/issues>`_ system
(see :ref:`faq`).


Citing CheckM
=============

If you use CheckM in your research, please cite CheckM and the software it relies on, as described in
:ref:`references`.


.. toctree::
   :caption: Getting started
   :maxdepth: 1

   installing/index
   quick_start
   faq


.. toctree::
   :caption: Running CheckM
   :maxdepth: 1

   workflows
   commands/index
   outputs/index


.. toctree::
   :caption: About
   :maxdepth: 1

   changelog
   references


.. toctree::
   :caption: Additional information
   :maxdepth: 1

   performance
   quality_estimate_changes
   python3

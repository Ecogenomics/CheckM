.. _workflows:

Workflows
=========

CheckM is executed from the command line and consists of a series of :ref:`commands <commands>` in order to
support a number of different analyses and workflows. The two most common workflows are to assess genomes
using either *lineage-specific* or *taxonomic-specific* marker sets. Genome quality can also be assessed using
custom marker genes, such as the marker set for CPR genomes provided with CheckM.


.. _workflows/lineage:

Lineage-specific workflow
-------------------------

The recommended workflow for assessing the completeness and contamination of genome bins is to use
lineage-specific marker sets. This workflow consists of 4 mandatory (M) steps and 1 recommended (R) step:

.. code-block:: text

    (M) checkm tree <bin folder> <output folder>
    (R) checkm tree_qa <output folder>
    (M) checkm lineage_set <output folder> <marker file>
    (M) checkm analyze <marker file> <bin folder> <output folder>
    (M) checkm qa <marker file> <output folder>

The :ref:`tree <commands/tree>` command places genome bins into a reference genome tree. All genomes to be
analyzed must reside in a single bin directory (or be listed in a tab-separated file; see
:ref:`quick_start`). CheckM assumes genome bins are in FASTA format with the extension ``fna``, though this
can be changed with the ``-x`` flag.

The ``tree`` command can optionally be followed by the :ref:`tree_qa <commands/tree_qa>` command, which will
indicate the number of phylogenetically informative marker genes found in each genome bin along with a
taxonomic string indicating its approximate placement in the tree. If desired, genome bins with few
phylogenetically informative marker genes may be removed in order to reduce the computational requirements
of the following commands. Alternatively, if only genomes from a particular taxonomic group are of interest,
these can be moved to a new directory and analyzed separately.

The :ref:`lineage_set <commands/lineage_set>` command creates a marker file indicating lineage-specific
marker sets suitable for evaluating each genome. This marker file is passed to the
:ref:`analyze <commands/analyze>` command in order to identify marker genes and calculate genome statistics.
Finally, the :ref:`qa <commands/qa>` command estimates the completeness and contamination of each genome bin
and can be used to produce different tables summarizing the quality of each genome bin.

For convenience, the 4 mandatory steps can be executed using the :ref:`lineage_wf <commands/lineage_wf>`
command:

.. code-block:: bash

    checkm lineage_wf <bin folder> <output folder>

This workflow writes the marker file ``lineage.ms`` to the output folder.


.. _workflows/taxonomy:

Taxonomic-specific workflow
---------------------------

In some cases it is convenient to analyze all genome bins with the same marker set. A common example would
be a set of genomes from the same taxonomic group. The workflow for using a taxonomic-specific marker set
consists of 3 mandatory (M) steps and 1 recommended (R) step:

.. code-block:: text

    (R) checkm taxon_list
    (M) checkm taxon_set <rank> <taxon> <marker file>
    (M) checkm analyze <marker file> <bin folder> <output folder>
    (M) checkm qa <marker file> <output folder>

The :ref:`taxon_list <commands/taxon_list>` command produces a table indicating all taxa for which a marker
set can be produced. All supported taxa at a given taxonomic rank can be listed by passing ``taxon_list``
the ``--rank`` flag:

.. code-block:: bash

    checkm taxon_list --rank phylum

The :ref:`taxon_set <commands/taxon_set>` command is used to produce marker sets for a specific taxon:

.. code-block:: bash

    checkm taxon_set phylum Cyanobacteria cyanobacteria.ms

The marker file produced by the ``taxon_set`` command is passed to the :ref:`analyze <commands/analyze>`
command in order to identify marker genes within each genome bin. All putative genomes to be analyzed must
reside in a single bin directory (or be listed in a tab-separated file). CheckM assumes genomes are in FASTA
format with the extension ``fna``, though this can be changed with the ``-x`` flag. Finally, the
:ref:`qa <commands/qa>` command estimates completeness and contamination and can be used to produce
different tables summarizing the quality of each genome bin.

For convenience, the above workflow can be executed in a single step with the
:ref:`taxonomy_wf <commands/taxonomy_wf>` command:

.. code-block:: bash

    checkm taxonomy_wf <rank> <taxon> <bin folder> <output folder>

This workflow writes the marker file ``<taxon>.ms`` (e.g. ``Bacteria.ms``) to the output folder.


.. _workflows/custom:

Using custom marker genes
-------------------------

CheckM supports using custom marker genes for assessing genome completeness and contamination. The desired
marker genes must be specified as hidden Markov models (HMMs) constructed with `HMMER <http://hmmer.org/>`_
(HMMER3 format). Genome quality can be assessed using these marker genes as follows:

.. code-block:: bash

    checkm analyze <custom HMM file> <bin folder> <output folder>
    checkm qa <custom HMM file> <output folder>

The HMM file is passed to the :ref:`analyze <commands/analyze>` command in order to identify marker genes in
each genome bin. Finally, the :ref:`qa <commands/qa>` command estimates the completeness and contamination
of each genome bin and can be used to produce different tables summarizing the quality of each genome bin.

When using custom marker genes, all markers are treated as a single marker set. In the output of ``qa``, the
marker lineage is ``N/A``, ``# genomes`` is ``-1`` (indicating it is unknown how many genomes were used to
create the custom marker set), and ``# marker sets`` is ``1``.


.. _workflows/cpr:

Using the CPR marker set
------------------------

CheckM tends to underestimate the completeness of CPR (Patescibacteria) genomes as they are often missing
many of the "universal" set of 104 bacterial marker genes used by CheckM. The 43 marker genes proposed by
`Brown et al., 2015 <https://www.nature.com/articles/nature14486>`_ likely provide improved estimates of the
quality of CPR genomes. This marker set has been provided by Adarsh Singh and
`Nate Cira <https://www2.rowland.harvard.edu/people/nate-cira>`_ for use in CheckM and is available in the
`custom_marker_sets <https://github.com/Ecogenomics/CheckM/tree/master/custom_marker_sets>`_ directory of
the CheckM repository as ``cpr_43_markers.hmm`` (the markers are described in the Excel spreadsheet
``cpr_43_markers.xlxs``). It can be applied to CPR genomes as a :ref:`custom marker set <workflows/custom>`:

.. code-block:: bash

    checkm analyze cpr_43_markers.hmm <bin folder> <output folder>
    checkm qa cpr_43_markers.hmm <output folder>

As an example, this marker set was applied to two CPR genomes:

* `Candidatus Roizmanbacteria bacterium GCA_001790095.1 <https://www.ncbi.nlm.nih.gov/assembly/GCA_001790095.1/>`_,
  which has an estimated completeness of 95.4% with the CPR marker genes and 72.1% when using the 104
  bacterial marker genes
* `Candidatus Paceibacter normanii GCA_000398025.1 <https://www.ncbi.nlm.nih.gov/assembly/GCF_000398025.1/>`_,
  which has an estimated completeness of 69.8% with the CPR marker genes and 58.6% when using the 104
  bacterial marker genes

As for any custom marker set, the marker lineage is reported as ``N/A``, ``# genomes`` as ``-1``, and
``# marker sets`` as ``1``.

.. _outputs:

Outputs
=======

The genome quality estimates and statistics reported by CheckM are produced by the :ref:`qa <commands/qa>`
command, which is also run as the final step of the :ref:`lineage_wf <commands/lineage_wf>` and
:ref:`taxonomy_wf <commands/taxonomy_wf>` workflows. The workflows write the results to
``quality_report.tsv`` in the output folder. The files in the output folder are described in
:ref:`outputs/files`.

.. toctree::
   :maxdepth: 1

   files
   qa
   qa_formats

Each command that writes to an output directory also writes a log file, ``checkm.log``, to this directory.

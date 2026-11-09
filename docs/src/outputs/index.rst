.. _outputs:

Outputs
=======

The genome quality estimates and statistics reported by CheckM are produced by the :ref:`qa <commands/qa>`
command, which is also run as the final step of the :ref:`lineage_wf <commands/lineage_wf>` and
:ref:`taxonomy_wf <commands/taxonomy_wf>` workflows.

.. toctree::
   :maxdepth: 1

   qa
   qa_formats

Each command that writes to an output directory also writes a log file, ``checkm.log``, to this directory.

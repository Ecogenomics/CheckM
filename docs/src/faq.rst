.. _faq:

FAQ and support
===============


.. _faq/large_numbers:

Processing large numbers of genomes
-----------------------------------

If you plan to process a large number of genomes, you may wish to break these into smaller batches.
With CheckM 1.3.0, running 1,000 genomes at a time with 64 threads takes ~20 minutes and 18 GB of memory
(see :ref:`performance`).
Exceeding the available memory of your machine will cause CheckM to use swap space (as per any
program), which will substantially increase the time to process genomes; use ``--gplace_max_memory`` to limit
the memory used to place genomes in the reference genome tree (see :ref:`installing/system_requirements`).
CheckM's output files are gzip-compressed (see :ref:`outputs/files`), and ``--rm_tmp_files`` can be used to
remove the intermediate files of a workflow once it has finished.

It has also been reported that processing large numbers of genomes at once can cause a Python recursion
error. As such, splitting large runs into batches of 1,000 genomes is recommended.


How do I improve the quality estimates for CPR genomes?
-------------------------------------------------------

Use the CPR marker set provided with CheckM, as described in :ref:`workflows/cpr`.


Where is the CheckM reference data?
-----------------------------------

See :ref:`installing/reference_data` for how to download the reference data and tell CheckM where it is
located.


.. _faq/bug_report:

Reporting a bug
---------------

When submitting a bug, please provide as much information as possible and, if feasible, a simplified version
of the data set which is causing the issue. The more information you provide, the quicker we will be able to
resolve the issue. Please report bugs through the
`GitHub issues <https://github.com/Ecogenomics/CheckM/issues>`_ system. The log file ``checkm.log`` written to
the output directory is useful to include with a bug report.

Known issues are tracked in the `GitHub issues <https://github.com/Ecogenomics/CheckM/issues>`_.


Feature requests
----------------

CheckM v1 is no longer actively developed: outside of minor updates to fix compatibility issues with
external dependencies, no maintenance is intended. New features are not planned.
`CheckM2 <https://github.com/chklovski/CheckM2>`_ is recommended for new analyses.


Contact information
-------------------

Inquiries can be sent to Donovan Parks (donovan.parks [at] gmail.com).


Contributors
------------

The following people have directly submitted enhancements to the CheckM code base or other substantial
contributions:

* Hunter Cameron: `<https://github.com/Ecogenomics/CheckM/pull/31>`_

Massive thanks to `baudrly <https://github.com/baudrly>`_, `Vini Salazar <https://github.com/vinisalazar>`_,
and `Asaf Peer <https://github.com/asafpr>`_ for the initial Python 2 to 3 porting of CheckM. Other
contributions are acknowledged in the :ref:`changelog`.

.. _python3:

Migration to Python 3
=====================

CheckM has been ported to Python 3 to accommodate Python 2 reaching end of life on January 1, 2020.
CheckM >= 1.1.0 requires Python 3, and Python 2 is no longer supported.

Many thanks to `baudrly <https://github.com/baudrly>`_, `Vini Salazar <https://github.com/vinisalazar>`_, and
`Asaf Peer <https://github.com/asafpr>`_ for the initial Python 2 to 3 port.


.. _python3/validation:

Python 2 to 3 validation
------------------------

The port of CheckM to Python 3 was validated on a set of 1,000 genomes randomly selected from the
`GTDB <https://gtdb.ecogenomic.org/>`_ R89 representative genomes. Results were compared to those generated
with CheckM v1.0.18, the last Python 2 version of CheckM. Identical results were obtained for the
``lineage_wf``, ``taxonomy_wf``, and ``ssu_finder`` commands across this set of test genomes. Other CheckM
commands were run on a small set of 3 genomes to verify they run to completion under Python 3.


.. _python3/removed:

Removed functionality
---------------------

The following commands were removed in CheckM v1.1.x in order to simplify the code base and focus CheckM and
support requests on critical functionality:

* ``bin_qa_plot``: non-critical, rarely used plot which does not scale to the large numbers of MAGs now being
  recovered
* ``par_plot``: non-critical plot and the same information is better presented in the reference
  distribution plots
* ``cov_pca``, ``tetra_pca``: alternatives to these static plots exist in tools such as
  `Anvi'o <http://merenlab.org/software/anvio/>`_
* ``len_plot``: rarely used plot which is largely redundant with the ``len_hist`` and ``nx_plot`` plots
* ``bin_union``, ``bin_compare``: feature-rich alternatives now exist, such as
  `DAS Tool <https://github.com/cmks/DAS_Tool>`_ and `UniteM <https://github.com/dparks1134/UniteM>`_

The ``join_tables`` command is also no longer available.

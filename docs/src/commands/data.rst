.. _commands/data:

data
====

Set the location of the CheckM reference data.

The ``data setRoot`` action stores the location of the CheckM reference data in CheckM's ``DATA_CONFIG``
file. A different ``DATA_CONFIG`` file can be given as an optional second argument. Setting the
``CHECKM_DATA_PATH`` environment variable is an alternative to this command and takes precedence over it.
See :ref:`installing/reference_data` for details.

Arguments
---------

.. literalinclude:: help/data.txt
   :language: text

Example
-------

.. code-block:: bash

    checkm data setRoot /path/to/my_checkm_data

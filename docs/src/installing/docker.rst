.. _installing/docker:

Docker
======

A `Dockerfile <https://github.com/Ecogenomics/CheckM/blob/master/docker/Dockerfile>`_ provided by a
community member is available in the ``docker`` directory of the CheckM repository. The image installs
HMMER, Prodigal, and pplacer from the Ubuntu package repository, installs the latest release of CheckM
from PyPI, and downloads the :ref:`CheckM reference data <installing/reference_data>` into the image
(``CHECKM_DATA_PATH`` is set accordingly).

No pre-built image is provided. To build the image:

.. code-block:: bash

    git clone https://github.com/Ecogenomics/CheckM.git
    cd CheckM/docker
    docker build -t checkm .

CheckM can then be run by mounting the directory containing your genomes into the container, e.g.:

.. code-block:: bash

    docker run --rm -v /path/to/my_data:/data checkm \
        checkm lineage_wf -t 8 -x fa /data/bins /data/checkm_output

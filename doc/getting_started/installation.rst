.. _installation:

Installation
============

A CADET installation consists of two parts: The CADET core simulator and a frontend.

Install CADET core simulator
----------------------------

The core simulator can be compiled from source, or you can download pre-built binaries.
If you want to extend or modify CADET (e.g., add a custom binding model), you will need to build CADET from source.

Install pre-built binaries
^^^^^^^^^^^^^^^^^^^^^^^^^^
CADET can be installed via conda from the ``conda-forge channel``.

``conda install -c conda-forge cadet``

This requires a working `conda installation <https://docs.anaconda.com/anaconda/install/index.html>`_.

Optionally, use `mamba <https://github.com/mamba-org/mamba>`_ which uses a faster dependency solver than ``conda``.

``mamba install -c conda-forge cadet``

Install from source
^^^^^^^^^^^^^^^^^^^

- :ref:`build_linux`
- :ref:`build_windows`
- :ref:`build_osx`

.. _cadet_python:

Install a frontend
------------------

As of now, only a Python frontend is provided.

Install CADET-Python
^^^^^^^^^^^^^^^^^^^^

The easiest way to create CADET simulations is to use the `CADET-Python <https://github.com/modsim/CADET-python>`_ frontend.
For this purpose, we recommend installing `Anaconda <https://www.anaconda.com/>`_.
Anaconda is a high-performance scientific distribution of Python that includes many common packages needed for scientific and engineering work.
Download the installer from their `website <https://www.anaconda.com/>`_ and run it for the local user.

To install ``CADET-Python``, open an `Anaconda Shell` and execute:

.. code-block:: bash

    pip install cadet-python

If you would also like to use `CADET-Match <https://github.com/modsim/CADET-Match>`_ for parameter estimation, run:

.. code-block:: bash

    pip install cadetmatch

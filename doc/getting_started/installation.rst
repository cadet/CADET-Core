.. _installation:

Installation
============

A CADET installation consists of two parts: The CADET core simulator and a frontend.

Install CADET core simulator
----------------------------

The core simulator can be compiled from source, or you can download pre-built binaries.
At the moment, only pre-built binaries for MS Windows are provided.
If you want to extend or modify CADET (e.g., add a custom binding model), you will need to build CADET from source.

Install pre-built binaries
^^^^^^^^^^^^^^^^^^^^^^^^^^
- Download the latest `release <https://github.com/modsim/CADET/releases>`_ for your platform (currently only binaries for windows are provided).
- Unzip the archive to your destination directory


Install from source
^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    build_linux
    build_windows
    build_osx

.. _cadet_python:

Install a frontend
------------------

As of now, a MATLAB and a Python frontend are provided.
In general, we recommend to use the Python frontend.
Note that the MATLAB interface is no longer actively developed and will be deprecated in a later version.

Install CADET-Python
^^^^^^^^^^^^^^^^^^^^

The easiest way to create CADET simulations is to use the `CADET-Python <https://github.com/modsim/CADET-python>`_ frontend.
For this purpose, we recommend installing `Anaconda <https://www.anaconda.com/>`_.
Anaconda is a high-performance scientific distribution of Python that includes many common packages needed for scientific and engineering work.
Download the installer from their `website <https://www.anaconda.com/>`_ and run it for the local user.

Before installing ``CADET-Python``, open an `Anaconda Shell` and and install these additional packages:

.. code-block:: bash

    conda install numpy scipy matplotlib gitpython jupyterlab ipywidgets
    conda install -c conda-forge jupyter_contrib_nbextensions jupyter_nbextensions_configurator

Moreover, we need to allow some additional channels for installing ``CADET-Python``:

.. code-block:: bash

    conda config --add channels anaconda-fusion
    conda config --add channels conda-forge

Then, to install ``CADET-Python`` run:

.. code-block:: bash

    conda install -c immudzen cadet-python

If you would also like to use `CADET-Match <https://github.com/modsim/CADET-Match>`_ for parameter estimation, run:

.. code-block:: bash

    conda install -c immudzen cadetmatch


Install CADET-MI
^^^^^^^^^^^^^^^^

.. warning::
    Note that the MATLAB interface is no longer actively developed and will be deprecated in a later version.

The MATLAB frontend is distributed with the pre-built binaries.
Run MATLAB and call ``installCADET()`` in the command window.

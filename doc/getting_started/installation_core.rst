.. _installation_core:

CADET-Core Installation
=======================

The core simulator can be compiled from source, or you can download pre-built binaries.
If you want to extend or modify CADET-Core (e.g., add a custom binding model), you will need to build CADET-Core from source.

.. note::

   On macOS ARM64 systems, CADET-Core must be built from source for now; see :ref:`build_osx` for instructions.

Install pre-built binaries
^^^^^^^^^^^^^^^^^^^^^^^^^^
CADET-Core can be installed via `conda <https://github.com/conda-forge/miniforge>`_ from the ``conda-forge channel``.

``conda install -c conda-forge cadet``

Install from source
^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   build_linux
   build_windows
   build_osx

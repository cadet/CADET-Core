CADET Introduction
==================

Performing a forward simulation comprises several steps: 
  * Setting up the model including all parameters
  * Defining connectivity and dynamic events
  * Setting up the simulator and actually running the simulation 
  * Evaluating results (e.g., plotting)

For this purpose, we recommend using `CADET-Process <https://cadet-process.readthedocs.io/>`_, an object oriented Python frontend for CADET.
CADET still must be downloaded (or built from source) as explained in the :ref:`installation guide <installation_core>`.

.. toctree::
   :maxdepth: 1

   tutorials/breakthrough


CADET-Core developers who might want to test their extensions of the simulator should use CADET-Python, a plain file based API for CADET.

.. toctree::
   :maxdepth: 1

   /../developer_guide/cadet_python

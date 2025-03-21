.. _installation_frontend:

Frontend Installation
=====================

CADET provides a Python API, called ``CADET-Python``, which can be used to configure the model according to :ref:`file_format`.
Setting up the model using ``CADET-Python`` can become very tedious, especially for systems, and is outdated now that an actual frontend is available:

We recommend using the ``CADET-Process`` frontend, which facilitates modeling processes using an object oriented model builder.
This interface layer provides convenient access to all model parameters in the system.
It automatically checks validity of the parameter values and sets reasonable default values where possible.

Install CADET-Process (recommended for users)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install ``CADET-Process``, open an `conda shell` and execute:

.. code-block:: bash

    pip install CADET-Process

Install CADET-Python (recommended for CADET-Core developers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install ``CADET-Python``, open an `conda shell` and execute:

.. code-block:: bash

    pip install CADET-Python 


Additional Resources
^^^^^^^^^^^^^^^^^^^^

- `CADET-Process documentation <https://cadet-process.readthedocs.io/>`_
- `CADET-Process Workshop <https://github.com/cadet/CADET-Workshop>`_
- `CADET-Python Tutorial <https://github.com/cadet/CADET-Python-Tutorial>`_
.. _getting_started:

Getting started
===============

This section details the steps to install CADET, as pre-built binaries or building from scratch, and also provides the series of tutorials that will help in building complete models from scratch in CADET. 

A CADET installation consists of two parts: The CADET-Core simulator and a frontend.

.. toctree::
   :maxdepth: 1

   installation_core

After the core simulator is installed, you can verify the installation via

.. code-block:: bash

    createLWE
    cadet-cli LWE.h5

The first command calls the ``createLWE``, which generates an ``h5`` configuration file for a fully configured `Load-Wash-Elute` process (some parameters can be changed, run ``createLWE --help`` for input options). The model is then simulated via the second command.

Next, in order to build your own model from scratch, we recommend installing a frontend.

.. toctree::
   :maxdepth: 2

   installation_frontend


Finally, the following introductory tutorials will guide you through the process of building a complete model from scratch.

.. toctree::
   :maxdepth: 1

   introduction

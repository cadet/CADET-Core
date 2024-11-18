CADET-Core
======

.. image:: https://img.shields.io/github/release/cadet/cadet-core.svg
   :target: https://github.com/cadet/cadet-core/releases

.. image:: https://github.com/cadet/cadet-core/actions/workflows/ci.yml/badge.svg?branch=master
   :target: https://github.com/cadet/cadet-core/actions/workflows/ci.yml?query=branch%3Amaster

.. image:: https://anaconda.org/conda-forge/cadet/badges/downloads.svg
   :target: https://anaconda.org/conda-forge/cadet

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8179015.svg
   :target: https://doi.org/10.5281/zenodo.8179015

.. image:: https://img.shields.io/badge/JuRSE_Code_Pick-Oct_2024-blue.svg
   :target: https://www.fz-juelich.de/en/rse/community-initiatives/jurse-code-of-the-month/october-2024

- **Website (including documentation):** https://cadet.github.io
- **Forum:** https://forum.cadet-web.de
- **Source:** https://github.com/cadet/cadet-core
- **Bug reports:** https://github.com/cadet/cadet-core/issues
- **Demo:** https://www.cadet-web.de 
- **Newsletter:** https://cadet-web.de/newsletter/

Installation
------------
CADET-Core can be installed via conda from the ``conda-forge`` channel.

``conda install -c conda-forge cadet``

This requires a working `conda installation <https://docs.anaconda.com/anaconda/install/index.html>`_.

Optionally, use `mamba <https://github.com/mamba-org/mamba>`_ which uses a faster dependency solver than ``conda``.

``mamba install -c conda-forge cadet``

`Additional information <https://cadet.github.io/master/getting_started/installation>`_ and a `tutorial <https://cadet.github.io/master/getting_started/tutorials/breakthrough>`_ are available to guide you through the installation and the first steps of using CADET.

Citing
------------
The development of CADET has been a collaborative effort, with multiple dedicated individuals contributing their expertise to create a powerful and versatile open-source software tool.
Countless hours of hard work have been invested to provide the scientific community with a valuable resource.
As an open-source project, CADET relies on the support and recognition from users and researchers to thrive.
Therefore, we kindly ask that any publications or projects leveraging the capabilities of CADET acknowledge its creators and their contributions by citing an adequate selection of our publications.

**General:**

- Leweke, S.; von Lieres, E.: `Chromatography Analysis and Design Toolkit (CADET) <https://doi.org/10.1016/j.compchemeng.2018.02.025>`_, Computers and Chemical Engineering **113** (2018), 274–294.

- von Lieres, E.; Andersson, J.: `A fast and accurate solver for the general rate model of column liquid chromatography <https://doi.org/10.1016/j.compchemeng.2010.03.008>`_, Computers and Chemical Engineering **34,8** (2010), 1180–1191.

**Numerics:**

- Breuer, J. M.; Leweke, S.; Schmölder, J.; Gassner, G.; von Lieres, E.: `Spatial discontinuous Galerkin spectral element method for a family of chromatography models in CADET <https://doi.org/10.1016/j.compchemeng.2023.108340>`_, Computers and Chemical Engineering **177** (2023), 108340.

- Leweke, S.; von Lieres, E.: `Fast arbitrary order moments and arbitrary precision solution of the general rate model of column liquid chromatography with linear isotherm <http://dx.doi.org/10.1016/j.compchemeng.2015.09.009>`_, Computers and Chemical Engineering **84** (2016), 350–362.

- Püttmann, A.; Schnittert, S.; Naumann, U.; von Lieres, E.: `Fast and accurate parameter sensitivities for the general rate model of column liquid chromatography <http://dx.doi.org/10.1016/j.compchemeng.2013.04.021>`_, Computers and Chemical Engineering **56** (2013), 46–57.

Additionally, to ensure reproducibility of your work, we recommend citing the zenodo doi corresponding to the specific CADET release that you used.

Selected applications that demonstrate the capabilities and use-cases of CADET-Core are highlighted in the `documentation <https://cadet.github.io>`_.

Ongoing Development
-------------------

We do our best to provide you with a stable API. However, CADET is actively developed and breaking changes can sometimes be unavoidable. For non-developers, it is recommended to upgrade from release to release instead of always working with the most recent commit.

Bugs
----

Please report any bugs that you find `here <https://github.com/cadet/cadet-core/issues>`_. Or, even better, fork the repository on `GitHub <https://github.com/cadet/cadet-core>`_ and create a pull request (PR) with the fix. 

Donations
---------

`Donations <https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=FCQ2M89558ZAG>`_ for helping to host, maintain, and further develop the CADET project are highly appreciated.


License
-------

CADET-Core is licensed under the GPL v3.0 License (see `LICENSE.txt <https://github.com/cadet/cadet-core/blob/master/LICENSE.txt>`_).

Copyright (C) 2004-present CADET-Core Authors (see `AUTHORS.md <https://github.com/cadet/cadet-core/blob/master/AUTHORS.md>`_).

Acknowledgments
---------------

Please refer to the `list of contributors <https://github.com/cadet/cadet-core/blob/master/CONTRIBUTING.md>`_ who helped building and funding this project.


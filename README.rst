CADET-Core
==========

.. image:: https://img.shields.io/github/release/cadet/cadet-core.svg
   :target: https://github.com/cadet/cadet-core/releases

.. image:: https://github.com/cadet/cadet-core/actions/workflows/ci.yml/badge.svg?branch=master
   :target: https://github.com/cadet/cadet-core/actions/workflows/ci.yml?query=branch%3Amaster

.. image:: https://codecov.io/gh/cadet/CADET-Core/graph/badge.svg?token=NKLJL03PA5
   :target: https://codecov.io/gh/cadet/CADET-Core/graph

.. image:: https://joss.theoj.org/papers/282a9ec56a0680e51ed4f1fa8fda3650/status.svg
   :target: https://joss.theoj.org/papers/282a9ec56a0680e51ed4f1fa8fda3650

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

To stay up to date on CADET news, development updates, and workshops, please subscribe to our `Newsletter <http://eepurl.com/hzN4EP>`_.

Installation
------------
CADET-Core can be installed via conda from the ``conda-forge`` channel.

``conda install -c conda-forge cadet``

This requires a working `conda installation <https://github.com/conda-forge/miniforge>`_.

`Additional information <https://cadet.github.io/master/getting_started/installation>`_ and a `tutorial <https://cadet.github.io/master/getting_started/tutorials/breakthrough>`_ are available to guide you through the installation and the first steps of using CADET.

Citing
------------
The development of CADET-Core has been a collaborative effort, with multiple dedicated individuals contributing their expertise to create a powerful and versatile open-source software tool.
Countless hours of hard work have been invested to provide the scientific community with a valuable resource.
As an open-source project, CADET-Core relies on the support and recognition from users and researchers to thrive.
Therefore, we kindly ask that any publications or projects leveraging the capabilities of CADET-Core acknowledge its creators and their contributions by citing an adequate selection of our publications.

**General:**

- Leweke, S.; Breuer, J.; Schmölder, J.; Jäpel, R.; Lanzrath, H.; Rao, J.; Hassan, J.; Zhang, W.; Berger, A.; Heymann, W.; von Lieres, E.: `CADET-Core 5.0: High-Performance Solver for Advanced Biotechnology Process Modeling <https://doi.org/10.21105/joss.07881>`_, Journal of Open Source Software **10** (2025), 7881.

- Leweke, S.; von Lieres, E.: `Chromatography Analysis and Design Toolkit (CADET) <https://doi.org/10.1016/j.compchemeng.2018.02.025>`_, Computers and Chemical Engineering **113** (2018), 274–294.

- von Lieres, E.; Andersson, J.: `A fast and accurate solver for the general rate model of column liquid chromatography <https://doi.org/10.1016/j.compchemeng.2010.03.008>`_, Computers and Chemical Engineering **34,8** (2010), 1180–1191.

**Major extensions:**

- Breuer, J. M.; Leweke, S.; Schmölder, J.; Gassner, G.; von Lieres, E.: `Spatial discontinuous Galerkin spectral element method for a family of chromatography models in CADET <https://doi.org/10.1016/j.compchemeng.2023.108340>`_, Computers and Chemical Engineering **177** (2023), 108340.

- Zhang, W.; Przybycien T., Breuer J. M. , Leweke S. , von Lieres E.: `Solving crystallization/precipitation population balance models in CADET, part II: Size-based Smoluchowski coagulation and fragmentation equations in batch and continuous modes <https://doi.org/10.1016/j.compchemeng.2024.108860>`_, Computers and Chemical Engineering **192** (2025), 108860.

- Zhang, W.; Przybycien T., Schmölder J. , Leweke S. , von Lieres E.: `Solving crystallization/precipitation population balance models in CADET, part I: Nucleation growth and growth rate dispersion in batch and continuous modes on nonuniform grids <https://doi.org/10.1016/j.compchemeng.2024.108612>`_, Computers and Chemical Engineering **183** (2024), 108612.

- Püttmann, A.; Schnittert, S.; Naumann, U.; von Lieres, E.: `Fast and accurate parameter sensitivities for the general rate model of column liquid chromatography <http://dx.doi.org/10.1016/j.compchemeng.2013.04.021>`_, Computers and Chemical Engineering **56** (2013), 46–57.

Additionally, to ensure reproducibility of your work, we recommend citing the `zenodo doi <https://doi.org/10.5281/zenodo.8179015>`_ corresponding to the specific CADET-Core release that you used.

For a comprehensive list and guidance on citing CADET-Core publications, please refer to the publications section of the `documentation <https://cadet.github.io/master/publications.html>`_.

Ongoing Development
-------------------

We do our best to provide you with a stable API. However, CADET-Core is actively developed and breaking changes can sometimes be unavoidable. For non-developers, it is recommended to upgrade from release to release instead of always working with the most recent commit.

Bugs
----

Please report any bugs that you find `here <https://github.com/cadet/cadet-core/issues>`_. Or, even better, fork the repository on `GitHub <https://github.com/cadet/cadet-core>`_ and create a pull request (PR) with the fix. 

Donations
---------

`Donations <https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=FCQ2M89558ZAG>`_ for helping to host, maintain, and further develop the CADET-Core project are highly appreciated.


Copyright and License Notice
----------------------------

Copyright (C) 2008-present: The CADET-Core Authors (see `AUTHORS.md <https://github.com/cadet/cadet-core/blob/master/AUTHORS.md>`_).

As of 24th June 2024, newly added files are licensed under the GNU AGPL v3.0-or-later.

Files existing in the repository before 24th June 2024 remain licensed under the GNU GPL v3.0-or-later unless their file header states otherwise.

When a file contains an SPDX license identifier, that identifier governs the license for that file.

This repository contains files licensed under both the GNU General Public
License (GPL) and the GNU Affero General Public License (AGPL), as described
above.

You may redistribute and/or modify each file under the terms of the license
that applies to that file, i.e. either the GNU GPL or the GNU AGPL as
published by the Free Software Foundation, either version 3 of the applicable
license, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the applicable license for more details.

You should have received a copy of the applicable license along with this program, the full texts are available in `LICENSE.txt <LICENSE.txt>`_ (AGPL) and `LICENSE-GPL.txt <LICENSE-GPL.txt>`_ (GPL).
If not, see <https://www.gnu.org/licenses/>.

Except as contained in this notice, the name of a copyright holder shall not be used in advertising
or otherwise to promote the sale, use, or other dealings in this Software without prior written
authorization of the copyright holder.


Acknowledgments
---------------

Please refer to the `list of contributors <https://github.com/cadet/cadet-core/blob/master/AUTHORS.md>`_ who helped building and funding this project.


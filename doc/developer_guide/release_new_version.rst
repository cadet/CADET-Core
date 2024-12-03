.. _release_new_version:

CADET-Core version release
==========================

CADET-Core releases follow the semantic versioning system, which is documented `here <https://semver.org/>`_.

Release checklist
-----------------

- Run the release tests:

  - The release tests contain extensive testing that is not included in our CI, such as EOC tests.
    Running these tests might take a while and this should be done on the server.
  - Some tests are implemented in CADET-Core with a `ReleaseCI` flag, which should be run by executing the ci workflow.
  - More tests are implemented in Python, the code can be found in `CADET-Verification <https://github.com/cadet/CADET-Verification>`_
  Compare the results with the previous run.
  The release process can only be continued if the results are reasonable.

- Run performance benchmarks:

  - If numerical algorithms were refactored or if performance-critical infrastructure was changed, you should run performance benchmarks to compare the latest release with the planned new one.
    To this end, you can refer to the performance benchmark templates in CADET-Reference, e.g. the `benchmark for the modified Newton method <https://jugit.fz-juelich.de/IBG-1/ModSim/cadet/cadet-reference/-/tree/benchmark_modified_newton?ref_type=heads>`_

- Create a version bump commit, which will be the target commit for the release.
  The bump commit contains:
  
  - Update of the version number in the `version.txt`, `zenodo.json`, `cadet.hpp` and cadet.doxyfile, compare to last `bump version` commit
  - Update of the authors list if needed: CONTRIBUTING.md and zenodo.json
  - Update of the copyright (years)
  - Update of the file format if needed

- Create the release on github `here <https://github.com/cadet/CADET-Core/releases/new>`_.

  - Add the version number according to the semantic versioning system as the tag and set the master branch as target.
  - Add release notes with these categories:

    - Added: New features, enhancements, or functionalities introduced in this release.
    - Fixed: Bug fixes and corrections made to resolve issues from previous versions.
    - Changed: Modifications to existing features and breaking changes for major releases including changes in the interface.
    - Updated: Improvements to documentation, minor tweaks, or other updates that don’t fit into the other categories.

- Check success of zenodo archiving:

  - Upon release, Zenodo automatically archives the release, generating a version-specific DOI (Digital Object Identifier) for it and storing a copy of the source code, along with any associated files.
    The `concept DOI <https://doi.org/10.5281/zenodo.8179015>`_, which is also given in the repository README, does not change but represents the repository as a whole and always points to the latest version.

Release of binaries on conda-forge
----------------------------------

To ensure CADET-Core is accessible to a broad community, it is available as a Python package on conda-forge.
Other software, such as our frontend, `CADET-Process`, and our Python interface, `CADET-Python`, import this package.

- go to github.com/conda-forge/cadet-feedstock
- create a new branch and open a PR:
- change the recipe/meta.yaml file:

  - generate sha key via ``curl -sL https://github.com/cadet/{{ name }}/archive/refs/tags/v{{ version }}.tar.gz | openssl sha256``, which requires open ssl und curl
  - update version number and sha key
  - set build number to zero (build: number: 0)

- Upon opening the PR, a todo checklist is automatically generated. After solving the todos, comment `@conda-forge-admin, please rerender` in the PR conversation. Automatic checks will be run and the bot will tell us the changes are fine. Then we can merge the PR, triggering the release on conda-forge.











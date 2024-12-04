CADET-Core release
==================

CADET-Core releases follow the *semantic versioning system*, which is documented [here](https://semver.org/).

**CADET-Core is released following the so-called *GitLab flow*:**
In GitLab flow, feature and hotfix branches contain work for new features and bug fixes which will be merged back into the master branch immediately when they’re finished, reviewed, and approved. The master branch is ready to be deployed, but not necessarily the source of truth for a new release: Each release has an associated release branch that is based off the master branch. Changes that are intended to be released are cherry-picked from master into the release branch. The release is made on the release branch. If desired, some release branches can be kept and maintained, which is currently not planned for the CADET-Core project.


**To make a release, the tasks in the following checklist have to be executed sequentially.**


Release on github
-----------------

- [ ] Run the release tests:

  - [ ] The release tests contain extensive testing that is not included in our CI, such as EOC tests.
    Running these tests might take a while and this should be done on the server.
  - [ ] Some tests are implemented in CADET-Core with a `ReleaseCI` flag. These tests need to pass and can be run by executing the `cd.yml` workflow.
  - [ ] More tests are implemented in Python, the code can be found in [CADET-Verification](https://github.com/cadet/CADET-Verification)
  Compare the results with the previous run.
  The release process can only be continued if the results are reasonable.

- [ ] Run performance benchmarks:

  - [ ] If numerical algorithms were refactored or if performance-critical infrastructure was changed, you should run performance benchmarks to compare the latest release with the planned new one.
    To this end, you can refer to the performance benchmark templates in CADET-Reference, e.g. the [benchmark for the modified Newton method](https://jugit.fz-juelich.de/IBG-1/ModSim/cadet/cadet-reference/-/tree/benchmark_modified_newton?ref_type=heads)

- [ ] Open a PR and create a version bump commit `Bump version to vX.X.X`:
  
  - [ ] Update of the version number in the `version.txt`, `zenodo.json`, `cadet.hpp` and cadet.doxyfile, compare to last `bump version` commit
  - [ ] Update of the authors list if needed: CONTRIBUTING.md and zenodo.json
  - [ ] Update of the copyright (years)
  - [ ] Update of the file format if needed

- [ ] Merge the PR into master.

- [ ] Create a Release branch off of master and cherry-pick the desired commits onto it.

- [ ] Create the release on github [here](https://github.com/cadet/CADET-Core/releases/new):

  - [ ] Set the release branch as target and manually specify a tag as vX.X.X with the version number according to the semantic versioning system.
  - [ ] Add release notes with these categories:

    - [ ] Added: New features, enhancements, or functionalities introduced in this release.
    - [ ] Fixed: Bug fixes and corrections made to resolve issues from previous versions.
    - [ ] Changed: Modifications to existing features and breaking changes for major releases including changes in the interface.
    - [ ] Updated: Improvements to documentation, minor tweaks, or other updates that don’t fit into the other categories.

- [ ] Check success of zenodo archiving: Upon release, Zenodo automatically archives the release, generating a version-specific DOI (Digital Object Identifier) for it and storing a copy of the source code, along with any associated files. The [concept DOI](https://doi.org/10.5281/zenodo.8179015), which is also given in the repository README, does not change but represents the repository as a whole and always points to the latest version.

Release of binaries on conda-forge
----------------------------------

To ensure CADET-Core is accessible to a broad community, it is available as a Python package on conda-forge.
Other software, such as our frontend, `CADET-Process`, and our Python interface, `CADET-Python`, import this package.

- [ ] Go to github.com/conda-forge/cadet-feedstock
- [ ] Create a new branch and open a PR:
- [ ] Change the recipe/meta.yaml file:

  - [ ] Generate sha key via ``curl -sL https://github.com/cadet/{{ name }}/archive/refs/tags/v{{ version }}.tar.gz | openssl sha256``, which requires open ssl und curl
  - [ ] Update version number and sha key
  - [ ] Set build number to zero (build: number: 0)

- [ ] Upon opening the PR, a todo checklist is automatically generated. After solving the todos, comment `@conda-forge-admin, please rerender` in the PR conversation. Automatic checks will be run and the bot will tell us the changes are fine. Then we can merge the PR, triggering the release on conda-forge.


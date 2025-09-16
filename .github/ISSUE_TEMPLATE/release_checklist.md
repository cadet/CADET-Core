---
name: Release checklist
about: Checklist to make a release
title: ''
labels: 'maintenance'
assignees: ''

---

# CADET-Core Release Checklist

CADET-Core follows the semantic versioning system described at [semver.org](https://semver.org/).

CADET-Core is released using the GitLab flow. In GitLab flow, feature and hotfix branches contain work for new features and bug fixes. These branches are merged back into the master branch immediately after they are finished, reviewed, and approved. The master branch is ready to be deployed but is not necessarily the source of truth for a new release. Each release has an associated release branch that is based on the master branch. Changes that are intended to be released are cherry-picked from master onto the release branch. The release is made on the release branch. Some release branches can be maintained, for example branch `v5.0.X` to support the last major release.

The following checklist describes the steps to execute sequentially for creating a new release.

---

## Preparation

- [ ] Create a release branch `vX.X.X` (typically from `master`) and cherry-pick the commits to include
- [ ] Create a version bump commit `Bump version to vX.X.X`:
  - Update the version number in `version.txt`, `zenodo.json` (two places), `cadet.hpp` and `cadet.doxyfile`
  - Update the authors list if needed in `CONTRIBUTING.md` and `zenodo.json`
  - Update the file format in `driver.hpp` if required
  - If the release contains all commits from `master`, merge the bump commit into `master`

---

## Pre-release build and Docker image

- [ ] Manually dispatch the `build_docker_containers` workflow ([cd_docker.yml](https://github.com/cadet/CADET-Core/blob/master/.github/workflows/cd_docker.yml)) and input the commit hash of the last commit on the release branch.  
  This triggers the Docker workflow to build an image containing:
  - The CADET-Core version at this commit
  - The latest released version of CADET-Python
  - The latest released version of CADET-Process

---

## Continuous Deployment tests

- [ ] Run the CD tests as described in the [CADET-Verification README](https://github.com/cadet/CADET-Verification).
- [ ] In the [CADET-Verification Output repository](https://github.com/cadet/CADET-Verification-Output):
  - Open the automatically created branch for this release and compare the results with the previous release. To this end, there are comparison utility functions as described in the [CADET-Verification README](https://github.com/cadet/CADET-Verification)
  - If the numerical engine or code structure changed, check the required simulation times
- [ ] If tests fail or show unexpected behaviour:
  - Fix the issues on the release branch
  - Repeat the tests until results are as expected
  - Merge fixes into `master` and resolve conflicts if needed
- [ ] Add the new run to the [meta-issue](https://github.com/cadet/CADET-Verification-Output/issues/1) and upload the log generated from the convergence data comparison
- [ ] Link the output branch in the release issue

---

## Creating the release on GitHub

- [ ] Go to [GitHub Releases](https://github.com/cadet/CADET-Core/releases/new):
  - Set the release branch as the target
  - Specify the tag `vX.X.X` according to semantic versioning
  - Add release notes with sections for Added, Fixed, Changed, and Updated
  - Publish the release.
- [ ] Check creation of docker image with the new release
- [ ] Verify Zenodo archiving:
  - Confirm that a version-specific DOI was created
  - Ensure that the source code and associated files are archived
  - Note that the [concept DOI](https://doi.org/10.5281/zenodo.8179015) remains constant

---

## Release of binaries on conda-forge

To ensure CADET-Core is accessible to a broad community, it is available as a Python package on conda-forge.  
Other software, such as CADET-Process and CADET-Python, import this package.

- [ ] Go to your fork of [cadet-feedstock](https://github.com/conda-forge/cadet-feedstock) or create one if it does not exist
- [ ] Create a new branch on your fork and change the file `recipe/meta.yaml`:
  - install openSSL
  - Generate the SHA256 key (replace `{{ version }}` with the semantic version number):  
    ```bash
    curl -sL https://github.com/cadet/CADET-Core/archive/refs/tags/v{{ version }}.tar.gz | openssl sha256
    ```
  - Update the version number and SHA256 key
  - Set the build number to zero (`build: number: 0`)
- [ ] Open a PR onto the main branch of the conda-forge repo, and complete the automatically generated checklist
  - Note: to check if the license file is included, check https://github.com/cadet/CADET-Core/archive/refs/tags/v{{ version }}.tar.gz if the License file is in the location specified in the meta.yml in the variable `license_file`
- [ ] Wait for the automatic checks to pass
- [ ] Merge the pull request to trigger the conda-forge release
- [ ] Double check if the new version of cadet-core is on conda-forge

## Follow-up
- [ ] Create or update a CADET forum post announcing the release, including release notes
- [ ] If this release checklist was updated, add these changes to the corresponding issue template


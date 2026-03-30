.. _release_strategy:

Code development and release strategy
=====================================

Agile code development
^^^^^^^^^^^^^^^^^^^^^^

PRs are merged into master as soon as they are complete, i.e. no half-baked solutions, tested, documented, updated C-API (if required).
This happens independent of when we want to release the code.

Due to our limited resources and since we’re not paid by our users, we need to avoid technical debt arising from trying to satisfy everybody at once.
Specifically, we do not implement multiple interfaces on master since they all need to be tested and later deprecated.

Release strategy and Code distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Patch releases** latest minor/major release is supported by patch releases, i.e. bug fixes and minor improvements are made on the corresponding release branch **and** on master.

**Minor releases**: New minor releases are branched off of master with the desired code changes to be cherry-picked.
If an interface breaking change is desired to be released in a minor version release, the code must be reworked to implement the old interface.
Added features implementing a new interface can be added to a minor version release.
Only major version relases may contain breaking changes to the interface.

**Major version pre-releases**: The above rule prevents the direct release (i.e. release without additional work) of code that breaks the interface until the next major version is published.
To distribute the code e.g. for publications, collaborations, early-access, we can create a pre-release of the major version, requiring users to adapt to the new interface as it currently is.
In this case, we can distribute binaries or the docker image and the corresponding documentation to collaborators.
This is a pre-release that we also want to have on zenodo, in contrast to pre-releases made for testing/verification.
It should be made clear to users that this is not a formal release and that the interface may still evolve.

**Major release**: Major releases may contain breaking changes to the interface.
Foreseeable breaking changes to established feature should thus be bundled in a major release.

This strategy is meant to reduce development time and maintenance overhead and shifts tedious/painful rework to where it pays off - collaboration, publication.

Making a new release
^^^^^^^^^^^^^^^^^^^^

Every release of CADET-Core follows the `release checklist
<https://github.com/CADET/CADET-Core/blob/master/.github/ISSUE_TEMPLATE/release_checklist.md>`_,
ensuring proper testing, documentation and deployment.

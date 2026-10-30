.. _CAPI:

Changing the CAPI
==================

The CADET-Core C API is designed to provide a stable and user-friendly interface for users who want to interact with CADET-Core from other programming languages.
One example for this is the implementation of the `CADET-Python software <https://github.com/cadet/CADET-python>`_, which implements a Python interface to CADET-Core using the C API.


**An example** for a change to the C API is given by the `CADET-Core PR #366 <https://github.com/cadet/CADET-Core/pull/366>`_ and corresponding `CADET-Python companion PR #71 <https://github.com/cadet/CADET-Python/pull/71>`_.


General Requirements
~~~~~~~~~~~~~~~~~~~~

- Changes to the CAPI should be made with care, as they can potentially break the interface for users of the C API and thus must be made **backwards compatible** or require a **major CADET-Core version** bump.
- For backwards compatibility, **multiple versions** of the C API can be **maintained in parallel**, see e.g. the various minor, patch or pre-release versions in ``CAPIv1.cpp``.
- Every change to the CAPI must **bump the C API version**, which is defined in ``cadet/LibVersionInfo.hpp``.
- Part of updating the C API is to **implement the corresponding support in CADET-Python**.
- Before a change to the C API is merged, it must be **tested** by implementing a **test case in CADET-Python** that uses the new C API functionality and ensures that it works as expected.


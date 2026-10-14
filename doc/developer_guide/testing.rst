.. _testing:

CADET-Core testing
==================

Run the tests
-------------

If you want to run tests in CADET-Core you need to ensure that in the `CMakeSettings.json` file:

1. The ``cmakeCommandArgs`` contain:
     1. ``-DENABLE_TESTS=ON`` to enable building the test runner
     2. ``-DENABLE_STATIC_LINK_LAPACK=ON -DENABLE_STATIC_LINK_DEPS=ON`` to create statically linked dependencies
2. The `variables` field contains:
    .. code-block:: json

        {
          {
            "name": "HDF5_USE_STATIC_LIBRARIES",
            "value": "1",
            "type": "STRING"
          },
          {
            "name": "BUILD_SHARED_LIBS",
            "value": "0",
            "type": "STRING"
          }
        }

After building, you can find the testRunner.exe in ``CADET-root/build/test/Release`` or, when build in debug mode, in ``CADET-root/build/test/Debug``.

To debug specific tests (with flag [testHere]) from the Visual Studio IDE, you can add the following configuration to the launch.vs.json file mentioned in the :ref:`debugging` section:

.. code-block:: json

    {
      "type": "default",
      "project": "CMakeLists.txt",
      "projectTarget": "testRunner.exe (test/Debug/testRunner.exe)",
      "name": "testRunner.exe (test/Debug/testRunner.exe)",
      "args": [
        "[testHere]"
      ]
    }

Select the testRunner.exe in the startup item dropdown menu and you can start debugging tests with the specified flag.
To see the options, run `testrunner --help`. To get the available test flags for instance, execute the testrunner with the corresponding option, i.e. `testrunner -t`.
You can provide multiple arguments like `"[testHere][testHere2]"` to run all test that have both flags, or with comma separation `"[testHere],[testHere2]"` to run all tests with the first and/or the second flag.

Adding tests for your model
---------------------------

Every model and numerical method implemented in CADET-Core has to be tested adequatly.
Here, we distinguish between unit tests and numerical reference tests.

Unit tests
^^^^^^^^^^

Unit tests are designed to verify the correctness of individual units or components of code in isolation.
Unit tests in CADET-Core typically encompass tests for the Jacobian implementation and consistent initialization of the model, which can be adapted from corresponding existing tests.
Testing the Jacobian typically involves comparing the analytical Jacobian to the AD Jacobian to verify that the residual implementation is consistent to the analytical Jacobian.
To this end, it might be necessary to increase the maximum number of AD directions for your test case, which can be done via the cmake argument `NUM_MAX_AD_DIRS`, as described in the :ref:`build_options`.

Numerical reference tests
^^^^^^^^^^^^^^^^^^^^^^^^^

Every major feature of the model or numerical method should be tested in a separate test case, preferably recreating examples from literature, including your own publication on the corresponding CADET extension.
This section provides a guide on creating reference tests, including convergence tables and research data management (RDM).
We strongly emphasize that the procedure described here should be taken into consideration not just as part of the software testing but also for method/model validation, especially when a paper publication is planned.

We utilize `CADET-Verification <https://github.com/cadet/CADET-Verification>`_ and `CADET-RDM <https://github.com/cadet/CADET-RDM>`_ to ensure reproducibility of the tests and thus maintainability and stability of the added modules.

**1. Define your model setup in CADET-Verification:**
There are many benchmark cases already implemented in CADET-Verification, please refer to them and adopt style, naming conventions and integration into the codebase.

**2. Add tests to CADET-Verification script:**
Add (potentially reduced/minimallized) test cases, e.g. case studies from a paper publication to `CADET-Verification <https://github.com/cadet/CADET-Verification>`_.
These tests should be added to the main scipt of CADET-Verification (refer to the README file), so that they are executed as part of the CADET-Core CD pipeline.
For numerical additions to the code, please add tests that verify these methods, e.g. experimental order-of-convergence (EOC) tests.

**3. Add numerical reference tests to CADET-Core:**
Use the code added to CADET-Verification and generate reference solutions to the benchmark cases.
Then, add numerical reference tests (refer to those already implemented), comparing a newly computed solution to the stored reference solution.
Use low discretization criteria for this, so that the test executes quickly.
Add these tests to the CI workflow by adding a `[CI]` flag.

Manufactured solution
^^^^^^^^^^^^^^^^^^^^^
If no analytical solution is available for your model, you can use a so-called "manufactured solution" to verify your implementation.
A manufactured solution was implemented e.g. for the radial flow model, see test/testRadialKernel.cpp.

General CADET-Core testing procedure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Create a new test file ``test\NewModelMethod.cpp``, an easy way to begin with is to copy one of the files that implements tests which are the closest to the ones that you are planning on, e.g. LumpedRateModelWithPores to test a new unit operation.
Add your tests to the testrunner executable by adding ``NewModelMethod.cpp`` to the list in the command ``add_executable(testrunner`` within the ``test\CmakeList.txt`` file.
Note that every test needs an unique name, which is specified for each test by ``TEST_CASE("My first test", "[FLAG1],[FLAG2]")``.
Flags are used as options for the testrunner.exe and are specified within square brackets.
Reuse existing flags and add new ones for your extension.
The ``[CI]`` flag is used for tests that shall be rerun as part of our github continuous integration (CI) pipeline.

Maintenance of the tests
------------------------

Some changes will break the tests without them being necessarily wrong. A change in the numerics for instance, will most likely break some tests.
This can be fixed by carefully adapting the absolute and relative tolerances for the broken tests. These changes should not change the magnitude of the tolerances, except if this is within an acceptable and expected new tolerance).

Test coverage
-------------

``Codecov`` is used to analyze test coverage through an automated workflow (see ``coverage.yml``), which runs on changes to the ``test/coverage`` branch and can also be triggered manually (workflow dispatch).
Coverage reports are automatically uploaded to ``Codecov``, providing detailed insights into tested code areas, uncovered lines, and coverage trends over time.
Tests marked with the ``[CI]`` flag are automatically included in the coverage analysis; otherwise, they can be added manually in the workflow file.
The goal is to maintain complete test coverage by ensuring every new functionality or code change includes corresponding tests executed with Codecov.
The file ``.codecov.yml`` defines a target coverage percentage. If the coverage falls below that, the workflow fails.
This percentage should be updated as we continuously expand our test coverage.
Further, the file defines ignored paths, such as the `ThirdParty` directory as we dont consider test coverage of third party software.

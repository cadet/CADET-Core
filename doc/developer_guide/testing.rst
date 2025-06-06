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

We utilize `CADET-Database <https://jugit.fz-juelich.de/IBG-1/ModSim/cadet/cadet-database>`_ and `CADET-Verification <https://github.com/cadet/CADET-Verification>`_ to ensure reproducibility of the tests.

We recommend using this procedure not only for the tests but also for the publication to make good standards regarding RDM.

The general procedure how to add a reference test is as follows:

**1. Define your model setup in CADET-Database:**
The easiest way to do this, is to clone the `CADET-Database <https://jugit.fz-juelich.de/IBG-1/ModSim/cadet/cadet-database>`_, checkout the core_tests branch and store your model setup file as a json file in the cadet-config folder.
You can create this json config file by translating the standard CADET h5 file to a json using python, see the utility/h5ToJson.py script for reference.
Make sure to give your model setup a meaningful and unique name (follow the naming logic of other setups in that folder).
Note that the return data should only specify the output required for the reference test, i.e. usually the outlet of a single unit.

**2. Add a CADET-Verification script:**
This python script reads the model setup from CADET-Database and generates all the reference data that will be used in your tests.
This includes a specific low resolution reference solution as well as convergence (EOC) tables computed using a high precision solution.
If analytical solutions are available for the considered model, you should use that as the high precision solution.
Otherwise, you additionally need to generate a very high resolution numerical reference solution, preferably with an accuracy of at least :math:`\mathcal{L}_\text{inf} \approx 1e-8`.

**3. Add numerical reference tests to CADET-Core:**
These tests should read the model setups previously defined in CADET-Database and run them with the same numerical specification as used to compute the reference solutions computed by the CADET-Verification script.
The resulta are compared to the reference solution generated by CADET-Verification.
This type of tests ensures that the model is still functional and that the numerics for this model have not changed.
Hence, every major feature of the model should be tested in a separete tests.
This way, we make sure that ongoing CADET-Core development doesnt break the model and these tests should be included in the CI pipeline by adding the [CI] flag as described in the implementation procedure section.

**4. Add EOC tests to CADET-Core (optional):**
These tests should be part of the paper publication which introduces the new model implemented in CADET-Core and can also be included in the CADET-Core tests.
Verifying the experimental order of convergence (EOC) is widely considered the most rigorous and best scientific practice in model and method validation, which is why we recommend including the EOC tables in your publication.
The convergence tests should not be added to the standard CI but only be rerun on release, i.e. by adding the [ReleaseCI] flag.
Details on how to compute EOC tables can be found elsewhere, please also refer to the already implemented EOC tests in CADET-Verification.

Manufactured solution
^^^^^^^^^^^^^^^^^^^^^
If no analytical solution is available for your model, you can use a so-called "manufactured solution" to verify your implementation.
A manufactured solution was implemented e.g. for the radial flow model, see test/testRadialKernel.cpp.

Implementation procedure
^^^^^^^^^^^^^^^^^^^^^^^^
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

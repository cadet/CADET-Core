.. _build_options

CADET build options
===================

The cmakeSettings.json provides standard build configurations for a ``DEBUG`` build, a ``RELEASE`` build and a ``RELEASE_with_Debug_Info`` build.
The ``DEBUG`` build compiles code suited for debuggin, i.e. the code is compiled such that every line of code is individually executed.
This way, its possible to go through the code by setting breakpoints in an IDE such as visual studio.
Additionally, the ``DEBUG`` build will print lot of additional information, such as the current step size and time point of the simulation.

The ``RELEASE`` build compiles optimized code and should be used to actually use CADET for simulations.

The ``RELEASE_with_Debug_Info`` build compiles optimized code but the debug information will still be printed.

CMake command arguments are specified for every build configuration in the cmakeSettings.json.

todo: describe functionality of every one.
DVCPKG_TARGET_TRIPLET=x64-windows-static -DENABLE_STATIC_LINK_LAPACK=ON -DENABLE_STATIC_LINK_DEPS=ON -DENABLE_TESTS=OFF -DENABLE_ANALYTIC_JACOBIAN_CHECK

todo: there are some additional arguments, that are not used in the current cmakeSettings.json, these should also be described
DENABLE_THREADING
-> Parallelized code will be compiled, using the TBB library. Note that the non-parallelized code is faster compared to the parallelized code when only one thread is being used. The number of threads is specified in the filed ``N_THREADS``.
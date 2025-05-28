# Build CADET-Core on MacOS

## Prerequisites

* CMake (>= 3.12.0)
* GCC >= 7.0, Clang >= 3.9, or Intel C++ 18.0
* Optional: Git

Assumed directory structure:

<pre>
|- CADET-Core
|    - src
|    - include
|    - [...]
|    - build
|    - install
</pre>

Note that the version numbers of the files and packages below are subject to change and will not always reflect the most recent version.

Also note that you have to use the same compiler for all packages. This is especially important if some of the packages are installed via a package manager such as [Homebrew](http://brew.sh/) which uses the system compiler (Clang).

## Build dependencies

```
brew update > /dev/null || true
brew install cmake --without-docs
brew install hdf5
brew install tbb
brew install eigen
```

### LAPACK

You can either use the native LAPACK implementation provided by Mac OS X (vecLib, Accelerate)
or install the freely available [Intel MKL](https://software.intel.com/sites/campaigns/nest/) which is very fast and probably faster than Accelerate.

## Build CADET-Core

- Clone the CADET-Core source code `git clone https://github.com/cadet/cadet-core.git CADET-Core`
- Create the directories `CADET-Core/build` and `CADET-Core/install`

* If using Intel MKL, execute `export MKLROOT=/opt/intel/mkl`
* Using standard LAPACK: Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>/install" ../`

    Using MKL (sequential): Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>/install" -DBLA_VENDOR=Intel10_64lp_seq ../`

    Using MKL (parallel): Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>/install" -DBLA_VENDOR=Intel10_64lp ../`

There are further compile flags that can be passed to the above cmake command, please refer to the build options described in the developer guide.

* Execute `make`
* Execute `make install`

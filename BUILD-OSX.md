## Prerequisites

* CMake (>= 3.12.0)
* GCC >= 7.0, Clang >= 3.9, or Intel C++ 18.0
* Optional: MATLAB R2009a or greater
* Optional: Git

Assumed directory structure:

<pre>
&lt;ROOT&gt;
   |-code
   |-cadet
   |-build
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

## Build CADET

* Download release of CADET or checkout from git
* Place the source in `<ROOT>/code` and create the directory `<ROOT>/build`
* Open a terminal and change to `<ROOT>/build`

* If using Intel MKL, execute `export MKLROOT=/opt/intel/mkl`
* Using standard LAPACK: Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>/cadet" ../code/`
 
    Using MKL (sequential): Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>/cadet" -DBLA_VENDOR=Intel10_64lp_seq ../code/`
 
    Using MKL (parallel): Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>/cadet" -DBLA_VENDOR=Intel10_64lp ../code/`
* Execute `make`
* Execute `make install`


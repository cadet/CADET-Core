## Prerequisites

* CMake (>= 3.12.0)
* GCC >= 7.0, Clang >= 3.9, or Intel C++ 18.0
* Optional: Git

Assumed directory structure:

<pre>
&lt;ROOT&gt;
   |-code
   |-cadet
   |-build
</pre>

Note that the version numbers of the files and packages below are subject to change and will not always reflect the most
recent version.

## Install dependencies
```
sudo apt-get update
sudo apt -y install build-essential libhdf5-dev liblapack-dev libblas-dev libtbb-dev libsuperlu-dev libeigen3-dev
```

### LAPACK

You can either use a LAPACK implementation provided by your distribution or install the freely available 
[Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=linux&distributions=online) which is very fast and probably faster than your 
distribution's implementation.

Obtain LAPACK from your distribution:
* Install the packages (LAPACK and BLAS) of your distribution (e.g., `liblapack3`, `liblapack-dev`, `libblas3`, `libblas-dev` for Ubuntu and Debian). Note that some packages only provide reference (i.e., slow) implementations and others (e.g., ATLAS, GOTO) perform much faster.


## Build CADET

- Download or git clone the CADET source code
- Place the source in `<ROOT>\code`
- Create the directories `<ROOT>\build` and `<ROOT>\cadet`


- Open a terminal and change to `<ROOT>/build`
- If using MKL, execute `export MKLROOT=/opt/intel/mkl`
- Using standard LAPACK: Execute `cmake -DCMAKE_INSTALL_PREFIX="../cadet" ../code/`
 
   - Using MKL (sequential): Execute `cmake -DCMAKE_INSTALL_PREFIX="../cadet" -DBLA_VENDOR=Intel10_64lp_seq ../code/`
 
   - Using MKL (parallel): Execute `cmake -DCMAKE_INSTALL_PREFIX="../cadet" -DBLA_VENDOR=Intel10_64lp ../code/`

- Execute `make`
- Execute `make install`

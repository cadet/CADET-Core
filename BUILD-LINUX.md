# Build CADET-Core on Linux

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

Note that the version numbers of the files and packages below are subject to change and will not always reflect the most
recent version.

## Install dependencies

### Debian/Ubuntu/...

```
sudo apt-get update
sudo apt -y install build-essential cmake libhdf5-dev libsuperlu-dev libeigen3-dev
```
#### LAPACK

You can either use a LAPACK implementation provided by your distribution or install the freely available Intel MKL.

For Intel run:
```
sudo apt -y install intel-mkl
```

For distro defaults run:
```
sudo apt -y install liblapack3 liblapack-dev libblas3 libblas-dev
```

### Fedora

Required packages:
```
sudo dnf install cmake gcc-c++ git eigen3-devel hdf5-devel lapack-devel
```

Recommended packages:
```
sudo dnf install SuperLU-devel suitesparse-devel
```

## Build CADET-Core

- Clone the CADET-Core source code `git clone https://github.com/cadet/cadet-core.git CADET-Core`
- Create the directories `CADET-Core/build` and `CADET-Core/install`

- Open a terminal and change to `CADET-Core/build`
- If using MKL, execute `export MKLROOT=/opt/intel/mkl`
- To compile:

	- Using standard LAPACK: Execute `cmake -DCMAKE_INSTALL_PREFIX="../install" ../`

	- Using MKL (sequential): Execute `cmake -DCMAKE_INSTALL_PREFIX="../install" -DBLA_VENDOR=Intel10_64lp_seq ../`

	- Using MKL (parallel): Execute `cmake -DCMAKE_INSTALL_PREFIX="../install" -DBLA_VENDOR=Intel10_64lp ../`

There are further compile flags that can be passed to the above cmake command, please refer to the build options described in the developer guide.

- To build:
	- Execute `make`
- To install:
	- Execute `make install`

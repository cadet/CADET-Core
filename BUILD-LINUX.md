# Build CADET-Core on Linux

## Prerequisites

* CMake (>= 3.12.0)

* GCC >= 7.0, Clang >= 3.9, or Intel C++ 18.0
* Optional: Git

Assumed directory structure:

<pre>
|- CADET
|    - src
|    - include
|    - [...]
|    - build
|    - install
</pre>

Note that the version numbers of the files and packages below are subject to change and will not always reflect the most
recent version.

## Install dependencies

```
sudo apt-get update
sudo apt -y install build-essential cmake libhdf5-dev libsuperlu-dev libeigen3-dev
```

### LAPACK

You can either use a LAPACK implementation provided by your distribution or install the freely available Intel MKL

for Intel run

```
sudo apt -y install intel-mkl
```

for distro defaults run

```
sudo apt -y install liblapack3 liblapack-dev libblas3 libblas-dev
```

## Build CADET

- Clone the CADET source code `git clone https://github.com/cadet/cadet-core.git CADET`
- Create the directories `CADET/build` and `CADET/install`

- Open a terminal and change to `CADET/build`
- If using MKL, execute `export MKLROOT=/opt/intel/mkl`
- To compile:
	- Using standard LAPACK: Execute `cmake -DCMAKE_INSTALL_PREFIX="../install" ../`

	- Using MKL (sequential): Execute `cmake -DCMAKE_INSTALL_PREFIX="../install" -DBLA_VENDOR=Intel10_64lp_seq ../`

	- Using MKL (parallel): Execute `cmake -DCMAKE_INSTALL_PREFIX="../install" -DBLA_VENDOR=Intel10_64lp ../`

- To build:
	- Execute `make`
- To install:
	- Execute `make install`

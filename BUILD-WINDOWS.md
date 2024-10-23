## Prerequisites

* Microsoft Visual Studio **2022** (Community Edition, or better)
* Intel OneAPI MKL
* Optional: Git
* Optional but not generally recommended*: Intel OneAPI TBB

*For most use-cases it is more efficient to parallelize by running multiple CADET simulations instead 
of parallelizing within one CADET simulation. Including the parallelization code in CADET can lead to performance
losses, even if parallelization within CADET is not used. 
Therefore, we recommend not including the parallelization library TBB 
unless you know your simulations are large enough to benefit from it.

Assumed directory structure:

<pre>
|- CADET
|    - src
|    - include
|    - [...]
|    - build
</pre>

Note that the version numbers of the files and packages below are subject to change and will not always reflect the most
recent version.


## Visual Studio:
We are using Visual Studio because it is the easiest way to install all required tools, compilers etc. in one step.

- Download [MS Visual Studio](https://visualstudio.microsoft.com/de/downloads/)
- Install Visual Studio with the workload "Desktop Developments with C++"

## Intel oneAPI Base Toolkit

- Download
  the [Intel oneAPI Base Toolkit online installer](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=window&distributions=webdownload&options=online).
- During the installation, select "Custom Installation" and deselect everything except for:
    - Intel® oneAPI Math Kernel Library (MKL)
    - optional*: Intel® oneAPI Threading Building Blocks (TBB)

## Optional:  clink

- Clink provides text completion, history, and line-editing to Windows Command Prompt
- Download the [clink installer exe](https://github.com/chrisant996/clink) and install clink.

## Prepare CADET code

- Clone the CADET source code into a `CADET` folder: 
  - `git clone https://github.com/cadet/cadet-core.git CADET`
- Create the directories `CADET\build` and `CADET\install`

## Build CADET in Visual Studio

- Open Visual Studio and open the `CADET` folder
- Navigate to "Tools" - "Command Line" and open either a "Developer Command Prompt" or "Developer PowerShell"
- Execute the following command:
    - `vcpkg integrate install` (this only needs to be run _once_ per PC and will require admin privileges)
- At the top, where it says `DEBUG`, select `aRELEASE` instead
- Wait for `vcpkg` to install all the dependencies. The first time this is done on your PC it can take ~15-30 minutes
- Wait for `cmake generation` to finish (see `output` window)
- From the status bar at the top select `Build`, `Build all`
- Once that finishes, select `Build`, `Install CadetFramework`
- The binaries will be located in `CADET\out\install\aRELEASE\bin`


## Build CADET from the command line

- Open Visual Studio and select "continue without code"
- Navigate to "Tools" - "Command Line" and open either a "Developer Command Prompt" or "Developer PowerShell"
- Execute the following commands:
- For Command Prompt:
    - `cd CADET\build`
    - `vcpkg integrate install` (this only needs to be run _once_ and will require admin privileges)
    - `set MKLROOT="C:/Program Files (x86)/Intel/oneAPI/mkl/latest"`
    - `cmake -DCMAKE_INSTALL_PREFIX=..\out\install\aRELEASE -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="%VCPKG_ROOT%/scripts/buildsystems/vcpkg.cmake" -DVCPKG_TARGET_TRIPLET=x64-windows-static -DENABLE_STATIC_LINK_LAPACK=ON -DENABLE_STATIC_LINK_DEPS=ON -DBLA_VENDOR=Intel10_64lp_seq --fresh ../`
        - If you want to use parallelization and have installed TBB, instead
          execute `set TBBROOT="C:/Program Files (x86)/Intel/oneAPI/tbb/latest"` and
        - `cmake -DCMAKE_INSTALL_PREFIX=..\out\install\aRELEASE -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="%VCPKG_ROOT%/scripts/buildsystems/vcpkg.cmake" -DVCPKG_TARGET_TRIPLET=x64-windows-static -DENABLE_STATIC_LINK_LAPACK=ON -DENABLE_STATIC_LINK_DEPS=ON -DBLA_VENDOR=Intel10_64lp --fresh ../`
    - `msbuild.exe INSTALL.vcxproj /p:Configuration=Release;Platform=x64`
  
- For PowerShell:
    - `cd CADET\build`
    - `vcpkg integrate install` (this only needs to be run _once_ and will require admin privileges)
    - `$ENV:MKLROOT = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"`
    - `cmake -DCMAKE_INSTALL_PREFIX="..\out\install\aRELEASE" -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="$ENV:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" -DVCPKG_TARGET_TRIPLET=x64-windows-static -DENABLE_STATIC_LINK_LAPACK=ON -DENABLE_STATIC_LINK_DEPS=ON -DBLA_VENDOR=Intel10_64lp_seq "../" --fresh`
      - If you want to use parallelization and have installed TBB, instead
           execute `$ENV:TBBROOT = "C:\Program Files (x86)\Intel\oneAPI\tbb\latest"`
      and 
      - `cmake -DCMAKE_INSTALL_PREFIX="..\out\install\aRELEASE" -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="$ENV:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" -DVCPKG_TARGET_TRIPLET=x64-windows-static -DENABLE_STATIC_LINK_LAPACK=ON -DENABLE_STATIC_LINK_DEPS=ON -DBLA_VENDOR=Intel10_64lp "../" --fresh`
    - `msbuild.exe INSTALL.vcxproj /p:Configuration="Release;Platform=x64"`
- The binaries will be located in `CADET\out\install\aRELEASE\bin`

## Test build results
- Navigate to the install location `cd CADET\out\install\aRELEASE\bin`
- Run:
  - `cadet-cli.exe --version`
  - `createLWE.exe`
  - `cadet-cli.exe LWE.h5`
- And confirm the output of the LWE.h5 by opening it in HDF5view or loading it in CADET-Process.
- If you get no printed return from the first command, run cadet-cli.exe by double-clicking it in the file explorer.
This raises error messages that are not raised from a cmd or PowerShell window.
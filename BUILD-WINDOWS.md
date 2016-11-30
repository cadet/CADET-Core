# Prerequisites

* CMake (>= 3.1.0)
* Microsoft Visual Studio 2015 (Community Edition, or better)
* Optional: MATLAB R2009a or greater
* Optional: Git

Assumed directory structure:

<pre>
&lt;ROOT&gt;
   |-Libs
   |   |- sundials
   |   |- hdf5
   |   |- lapack
   |-code
   |-cadet
   |-build
</pre>

Note that the version numbers of the files and packages below are subject to change and will not always reflect the most recent version.

# Build dependencies

## HDF5

* Make sure that no HDF5 libraries are installed (remove already existing HDF5 installations)
* Download CMake-enabled source from https://support.hdfgroup.org/HDF5/release/cmakebuild.html or https://support.hdfgroup.org/HDF5/release/cmakebuild5110.html
* Unzip and make sure that the directory path does not contain blank spaces
* Open VS2015x64 Command Prompt and change to the unzipped directory
* Execute `ctest -S HDF5config.cmake,BUILD_GENERATOR=VS201564,INSTALLDIR="<ROOT>\Libs\hdf5" -C Release -V`
* Extract the created `HDF5-1.8.17-win64.zip` file to `<ROOT>/Libs/hdf5` such that you have `<ROOT>/Libs/hdf5/lib`
* We need to remove the `lib` prefix of the `.lib` files in `<ROOT>/Libs/hdf5/lib` (i.e., rename `libhdf5.lib` to `hdf5.lib`). Change to `<ROOT>/Libs/hdf5/lib` and execute `rename "lib*" "///*"`

## SUNDIALS

* Download SUNDIALS source from http://computation.llnl.gov/projects/sundials/sundials-software
* Unzip
* Open VS2015x64 Command Prompt and change to the parent directory of the unzipped directory
* Create a new folder `sundialsbuild` and change to it
* Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>\Libs\sundials" -DEXAMPLES_ENABLE=OFF -DOPENMP_ENABLE=ON -DBUILD_SHARED_LIBS=OFF -G "Visual Studio 14 Win64" -DCMAKE_C_FLAGS="/GL" -DCMAKE_STATIC_LINKER_FLAGS="/LTCG" ..\sundials-2.7.0\`
* Execute `msbuild.exe INSTALL.vcxproj /p:Configuration=Release;Platform=x64`

## LAPACK

In the following, CLAPACK is built and used. You can also install the freely available [Intel MKL](https://software.intel.com/sites/campaigns/nest/) which is very fast and certainly much faster than CLAPACK.

* Download `clapack-3.2.1-CMAKE.tgz` from https://icl.cs.utk.edu/lapack-for-windows/clapack/index.html#build and unzip
* Open VS2015x64 Command Prompt and change to the parent directory of the unzipped directory
* Create a new folder `clapackbuild` and change to it
* Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>\Libs\clapack" -G "Visual Studio 14 Win64" -DCMAKE_C_FLAGS="/GL" -DCMAKE_CXX_FLAGS="/GL" -DCMAKE_STATIC_LINKER_FLAGS="/LTCG" ..\clapack-3.2.1-`
* Execute `msbuild.exe CLAPACK.sln /p:Configuration=Release;Platform=x64`
* Execute `msbuild.exe INSTALL.vcxproj /p:Configuration=Release;Platform=x64`
* Rename `<ROOT>\Libs\clapack\lib\libf2c.lib` to `<ROOT>\Libs\clapack\lib\f2c.lib`

# Build CADET

* Download release of CADET or checkout from git
* Place the source in `<ROOT>\code` and create the directory `<ROOT>\build`
* Open VS2015x64 Command Prompt and change to `<ROOT>\build`
* Execute `set SUNDIALS_ROOT=<ROOT>\Libs\sundials`
* Execute `set HDF5_ROOT=<ROOT>\Libs\hdf5`
* Using CLAPACK: Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>\cadet" -DBUILD_TESTS=OFF -DCMAKE_LIBRARY_PATH="<ROOT>\Libs\clapack\lib" -G "Visual Studio 14 Win64" -DBLA_VENDOR=CLAPACK ..\code\`
 
    Using MKL (sequential): Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>\cadet" -DBUILD_TESTS=OFF -DCMAKE_LIBRARY_PATH="<MKL_ROOT>\compilers_and_libraries_2016\windows\mkl\lib\intel64_win" -G "Visual Studio 14 Win64" -DBLA_VENDOR=Intel10_64lp_seq ..\code\`
 
    Using MKL (parallel): Execute `cmake -DCMAKE_INSTALL_PREFIX="<ROOT>\cadet" -DBUILD_TESTS=OFF -DCMAKE_LIBRARY_PATH="<MKL_ROOT>\compilers_and_libraries_2016\windows\mkl\lib\intel64_win" -G "Visual Studio 14 Win64" -DBLA_VENDOR=Intel10_64lp ..\code\`
* Execute `msbuild.exe CadetFramework.sln /p:Configuration=Release;Platform=x64`
* Execute `msbuild.exe INSTALL.vcxproj /p:Configuration=Release;Platform=x64`


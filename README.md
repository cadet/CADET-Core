![CADET Logo](doc/logo/CADET-GitHub.png "Chromatography Analysis and Design Toolkit")

# Chromatography Analysis and Design Toolkit

The Chromatography Analysis and Design Toolkit (CADET) is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. The core of the CADET software is a fast and accurate solver for the General Rate Model (GRM) of packed bed liquid chromatography. The CADET solver covers a wide range of GRM variants, combining different transport and binding models with state-of-the-art mathematical algorithms and scientific computing techniques. 

The CADET framework also comprises a MATLAB interface and standard routines for parameter estimation, process optimization and experimental design. CADET is freely distributed (under the terms of the GPLv3) as a contribution to the scientific community. If you find it useful for your own work, we would appreciate acknowledgements of the CADET software and citations of our papers (see the [BibTeX citations](https://github.com/modsim/cadet/wiki/Referencing-CADET)):

* von Lieres, E.; Andersson, J.: [A fast and accurate solver for the general rate model of column liquid chromatography](http://dx.doi.org/10.1016/j.compchemeng.2010.03.008), Computers and Chemical Engineering 34,8 (2010), 1180–1191.
* Püttmann, A.; Schnittert, S.; Naumann, U.; von Lieres, E.: [Fast and accurate parameter sensitivities for the general rate model of column liquid chromatography](http://dx.doi.org/10.1016/j.compchemeng.2013.04.021), Computers and Chemical Engineering 56,13 (2013), 46-57.
* Püttmann, A.; Schnittert, S.; Leweke, S.; von Lieres, E.: [Utilizing algorithmic differentiation to efficiently compute chromatograms and parameter sensitivities](http://dx.doi.org/10.1016/j.ces.2015.08.050), Chemical Engineering Science, 139 (2016), 152–162.

## Features

* Fast and accurate solution of the partial differential algebraic equations (PDAE) based on the finite volume method and the WENO scheme
* Computation of derivatives with respect to model parameters by a forward sensitivity approach combined with algorithmic differentiation (AD)
* Shared memory parallelization using Intel TBB
* Native MATLAB interface (MEX) with many examples for rapid development
* Modular design allows for easy experimenting with the code (e.g., adding a new binding model or unit operation)
* Supports XML and HDF5 as data formats
* Multi-platform: Works on Windows, Linux, and Mac OS X

## Get CADET

Download the [latest release](https://github.com/modsim/cadet/releases) for your platform.
Check the [tutorials](https://github.com/modsim/cadet/wiki/tutorials) on how to install CADET.

## Requirements and Dependencies

CADET has been successfully built and run on the following platforms:

* MS Windows 7
* Linux (Ubuntu 12.04, Ubuntu 16.04, SLES11 SP3, SLES12)
* Mac OS X (10.6.8 Snow Leopard, 10.9 Mavericks).

What you need to build CADET:

* A C++11 capable compiler (e.g., GCC >= 4.7, Clang >= 3.3, MS Visual Studio >= 2015)
* [CMake](http://cmake.org/) >= 3.1.0
* IDAS of the [SUNDIALS](http://computation.llnl.gov/casc/sundials/main.html) package
* A LAPACK implementation (e.g., [LAPACK](http://www.netlib.org/lapack/index.html), [CLAPACK](https://icl.cs.utk.edu/lapack-for-windows/clapack/index.html#build), [Intel MKL](https://software.intel.com/sites/campaigns/nest/), or a native one provided by your operating system)
* [HDF5](http://www.hdfgroup.org/HDF5/)
* MATLAB (R2010b or higher, optional but highly recommended for interfacing the simulator) 

## Tutorials and Instructions

Please find instructions on how to build, install, or use CADET in the [Wiki](https://github.com/modsim/cadet/wiki).

For example, there are tutorials on

* [How to install CADET](https://github.com/modsim/cadet/wiki/Install-CADET)
* [How to build CADET](https://github.com/modsim/cadet/wiki/Build-CADET)

## Online Simulator

Try CADET in the browser using the [online simulator web service](http://www.cadet-web.de).

## Commercial software using CADET

For commercial software built on top of CADET and designed for experimentalists, [KBI](http://www.kbibiopharma.com/) has built [![Reveal Chromatography](https://www.cadet-web.de/static/simulation/RevealChromLogo.png)](http://bit.ly/reveal_prod).

## Ongoing Development

Since CADET is actively developed, **do not expect a stable API**. Breaking changes and extensive restructuring may occur in any commit and release.
For non-developers it is recommended to upgrade from release to release instead of always working with the most recent commit.

## Donations

[Donations](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=FCQ2M89558ZAG) for helping to host, maintain and further develop the CADET project are highly appreciated.

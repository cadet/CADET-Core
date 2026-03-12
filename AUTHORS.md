# Contributors

The following individuals have made significant contributions to the CADET-Core software and are listed in *chronological order*. Together they form the "CADET-Core authors" as mentioned in the copyright statements and [LICENSE.md](LICENSE.txt) file.

Major additions or modifications to CADET-Core are explicitly documented here.

## Special thanks to everyone who has helped with this project

* [Eric von Lieres](https://github.com/lieres)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-0309-8408)[^1]: Supervision, concepts, first MATLAB implementation, user interface design
* Joel Andersson[^1]: Domain decomposition (GRM), linear solver (GRM), first C implementation
* Sebastian Schnittert[^1]: First C++ implementation, first Matlab interface, first file format
* Andreas Püttmann[^1]: Parameter sensitivities, Jacobians via AD, AD techniques (block- & band-compression)
* [Samuel Leweke](https://github.com/sleweke)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-9471-4511)[^1]: Complete rewrite of v3.0 from scratch, build system, code & framework design and architecture, Matlab interface, file format, frontends (MEX, CLI), domain decomposition (multi-unit systems), reactions, 2D general rate model (paid by GE Healthcare), multiple particle types (paid by GE Healthcare)
* [William Heymann](https://github.com/immudzen)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-5093-0797)[^1]: Networks of unit operations, compiler flag optimization
* [Salah Azzouzi](https://github.com/azzouzis)[^1]
* [Johannes Schmölder](https://github.com/schmoelder)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-0446-7209)[^1]
* [Jayghosh Rao](https://github.com/jayghoshter)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-1216-9394)[^1]
* [Jazib Hassan](https://github.com/jazib-hassan-juelich)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-2741-5460)[^1]
* [Jan Michael Breuer](https://github.com/jbreue16)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-1999-2439)[^2][^1]: Spatial DG discretization (1D, 2D, particles), Crystallization module, refactored/generalized chromatography units
* [Ronald Jäpel](https://github.com/ronald-jaepel)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-4912-5176)[^1]
* [Hannah Lanzrath](https://github.com/hannahlanzrath)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-2675-9002)[^1]
* [Antonia Berger](https://github.com/AntoniaBerger)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0009-0002-0207-9042)[^1]: Refactored/generalized reactions system
* [Wendi Zhang](https://github.com/WFlynnZ)[![ORCID iD](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-6127-1408)[^3]: Crystallization module

[^1]: Forschungszentrum Juelich GmbH, IBG-1: Biotechnology, Juelich, Germany
[^2]: University of Cologne, Department of Mathematics and Computer Science, Cologne, Germany
[^3]: Rensselaer Polytechnic Institute, Chemical and Biological Engineering, Troy, New York, USA

## Funding Acknowledgement

* GE Healthcare sponsored the development of core-shell beads, multiple particle types, and 2D general rate model.
* Bayer AG contributed to bug fixes, initial C API, initial Python frontend.
* Johannes Schmölder has received support from the IMI2/ EU/EFPIA joint undertaking Inno4Vac (grant no. 101007799).
* Jan Michael Breuer has received funding by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – 548805630.

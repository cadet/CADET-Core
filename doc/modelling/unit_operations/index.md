(unit-operation-models)=

# Unit operation models

A comprehensive collection of chromatographic unit-operation models is provided through the Python-based tool [CADET-Equations](https://github.com/cadet/CADET-Equations)
The tool provides a graphical interface and allows users to export the corresponding governing equations in LaTeX or PDF format.
The application is [currently hosted on streamlit](https://cadet-equations-74ko8eryoxmsbqspggxrj2.streamlit.app/) and is planned to be integrated into the CADET website in the future.

A short comparison of the most prominent unit operation model features
is given in [](#table-features-unit-operations).

(table-features-unit-operations)=

| Unit operation model | Radial dispersion | Pore diffusion | Film diffusion | Particle geometries | Multiple particle types |
| --- | --- | --- | --- | --- | --- |
| [General Rate Model](#general-rate-model-model) | × | ✓ | ✓ | ✓ | ✓ |
| [Lumped Rate Model with Pores](#lumped-rate-model-with-pores-model) | × | × | ✓ | ✓ | ✓ |
| [Lumped Rate Model without Pores](#lumped-rate-model-without-pores-model) | × | × | × | × | × |
| [2D General Rate Model](#d-general-rate-model-model) | ✓ | ✓ | ✓ | ✓ | ✓ |
| [CSTR](#cstr-model) | × | × | × | × | ✓ |
| [Multi-Channel Transport Model](#multi-channel-transport-model-model) | × | × | × | × | × |
| [Crystallization](#ffcrystallization) | × | × | × | × | × |

Moreover, the pseudo unit operations [](#inlet-model), and [](#outlet-model) act as sources and sinks for the system.
We further note that radial flow model variants are available for the LRM, LRMP and GRM.


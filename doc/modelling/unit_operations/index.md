(unit-operation-models)=

# Unit operation models

A comprehensive collection of chromatographic unit-operation models is provided through the Python-based tool [CADET-Equations](https://github.com/cadet/CADET-Equations).
The tool provides a graphical interface and allows users to export the corresponding governing equations in LaTeX or PDF format.
The application is [hosted as a web application](https://equations.cadet-web.de/) and will be integrated into the CADET modeling documentation in the future.

A short comparison of the most prominent unit operation model features
is given in {numref}`table-features-unit-operations`.

(table-features-unit-operations)=

```{list-table} Supported features of the different unit operations models
:header-rows: 1

* - Unit operation model
  - Radial dispersion
  - Pore diffusion
  - Film diffusion
  - Particle geometries
  - Multiple particle types
* - {ref}`general-rate-model-model`
  - ×
  - ✓
  - ✓
  - ✓
  - ✓
* - {ref}`lumped-rate-model-with-pores-model`
  - ×
  - ×
  - ✓
  - ✓
  - ✓
* - {ref}`lumped-rate-model-without-pores-model`
  - ×
  - ×
  - ×
  - ×
  - ×
* - {ref}`d-general-rate-model-model`
  - ✓
  - ✓
  - ✓
  - ✓
  - ✓
* - {ref}`cstr-model`
  - ×
  - ✓
  - ✓
  - ✓
  - ✓
* - {ref}`multi-channel-transport-model-model`
  - ×
  - ×
  - ×
  - ×
  - ×
* - {ref}`ffcrystallization`
  - ×
  - ×
  - ×
  - ×
  - ×
```

Moreover, the pseudo unit operations {ref}`inlet-model`, and {ref}`outlet-model` act as sources and sinks for the system.
We further note that column geometry variants, specifically axial-flow, radial-flow and concical frustum, are available for the LRM, LRMP and GRM.

```{toctree}
:hidden:
:glob:

general_rate_model
lumped_rate_model_without_pores
lumped_rate_model_with_pores
2d_general_rate_model
multi_channel_transport_model
cstr
inlet
outlet
```


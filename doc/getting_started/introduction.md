# CADET Introduction

Performing a forward simulation comprises several steps:
: - Setting up the model including all parameters
  - Defining connectivity and dynamic events
  - Setting up the simulator and actually running the simulation
  - Evaluating results (e.g., plotting)

For this purpose, we recommend using [CADET-Process](https://cadet-process.readthedocs.io/), an object oriented Python frontend for CADET.
CADET still must be downloaded (or built from source) as explained in the [installation guide](#installation-core).

CADET-Core developers who might want to test their extensions of the simulator should use CADET-Python, a plain file based API for CADET.


.. _gaussian_process_regression_config:

Gaussian Process Regression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = GAUSSIAN_PROCESS_REGRESSION**

For information on model equations, refer to :ref:`gaussian_process_regression`.

``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =======================
**Type:** int        **Range:** {0,1}           **Length:** 1/NTOTALBND
===================  =========================  =======================

``CP_NDIM``
   Number of input dimensions per training point used in ``CP_VALS``.
   Must equal number of components, user sanity check.

===================  =========================  =======================
**Type:** int        **Range:** NCOMP           **Length:** 1
===================  =========================  =======================

``CS_NDIM``
   Number of output dimensions per training point used in ``CS_VALS``.
   Must equal total number of bound states, user sanity check.

===================  =========================  =======================
**Type:** int        **Range:** NTOTALBOUND     **Length:** 1
===================  =========================  =======================

``CP_VALS``
   Flattened pore-phase concentration training inputs used by the Gaussian
   process regression model. The values represent the training input points
   :math:`X` used to evaluate the kernel function. The array is interpreted
   as an ``NTRAIN`` major vector.

**Unit:** :math:`mol~m_{MP}^{-3}`

===================  =========================  ============================
**Type:** double     **Range:** unrestricted    **Length:** NTRAIN * CP_NDIM
===================  =========================  ============================

``CS_VALS``
   Solid-phase training targets corresponding to ``CP_VALS``.
   These values form the training output vector used to compute the GPR
   coefficient vector :math:`\alpha = (K + \sigma_n^2 I)^{-1} y`.

**Unit:** :math:`mol~m_{SP}^{-3}`

===================  =========================  ============================
**Type:** double     **Range:** unrestricted    **Length:** NTRAIN * CS_NDIM
===================  =========================  ============================

``GPR_KKIN``
   Linear-driving-force coefficients in component-major ordering.

**Unit:** :math:`s^{-1}`

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** NCOMP
===================  =========================  =======================

``KERNEL``
   Selects the kernel function used by the Gaussian process regression
   model. Supported values are ``MLP``, ``RBF``, ``RBF_Linear``, and
   ``MLP_Linear``.

===================  ================================================  ====================
**Type:** string     **Range:** {MLP, RBF, RBF_Linear, MLP_Linear}     **Length:** 1
===================  ================================================  ====================


Kernel-Specific Hyperparameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each bound state requires its own kernel configuration specified using the
suffix ``_BNDXXX``, where ``XXX`` is the zero-padded bound state index
(e.g., ``_BND000``, ``_BND001``, ``_BND002``).

Bound state indices are assigned sequentially across components: if component 0
has 2 bound states and component 1 has 1 bound state, the indices are:

- Bound state 0: component 0, state 0
- Bound state 1: component 0, state 1
- Bound state 2: component 1, state 0

The following table summarizes which parameters are required for each kernel:

+--------------------------------+------------+-----------+------------------+------------------+
| **Kernel**                     | **MLP**    | **RBF**   | **RBF_Linear**   | **MLP_Linear**   |
+================================+============+===========+==================+==================+
| MLP_WEIGHT_VARIANCE_BND_XXX    | ✓          | —         | —                | ✓                |
+--------------------------------+------------+-----------+------------------+------------------+
| MLP_BIAS_VARIANCE_BND_XXX      | ✓          | —         | —                | ✓                |
+--------------------------------+------------+-----------+------------------+------------------+
| MLP_VARIANCE_BND_XXX           | ✓          | —         | —                | ✓                |
+--------------------------------+------------+-----------+------------------+------------------+
| RBF_VARIANCE_BND_XXX           | —          | ✓         | ✓                | —                |
+--------------------------------+------------+-----------+------------------+------------------+
| RBF_LENGTHSCALE_BND_XXX        | —          | ✓         | ✓                | —                |
+--------------------------------+------------+-----------+------------------+------------------+
| LINEAR_VARIANCE_BND_XXX        | —          | —         | ✓                | ✓                |
+--------------------------------+------------+-----------+------------------+------------------+
| GAUSSIAN_NOISE_VARIANCE_BND_XXX| ✓          | ✓         | ✓                | ✓                |
+--------------------------------+------------+-----------+------------------+------------------+

``MLP_WEIGHT_VARIANCE_BND_XXX``
   Weight variance hyperparameter for the MLP (arc-cosine) kernel.
   Required for ``MLP`` and ``MLP_Linear`` kernels.

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** 1
===================  =========================  =======================

``MLP_BIAS_VARIANCE_BND_XXX``
   Bias variance hyperparameter for the MLP (arc-cosine) kernel.
   Required for ``MLP`` and ``MLP_Linear`` kernels.

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** 1
===================  =========================  =======================

``MLP_VARIANCE_BND_XXX``
   Output variance hyperparameter for the MLP (arc-cosine) kernel.
   Required for ``MLP`` and ``MLP_Linear`` kernels.

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** 1
===================  =========================  =======================

``RBF_VARIANCE_BND_XXX``
   Output variance hyperparameter for the RBF (squared exponential) kernel.
   Required for ``RBF`` and ``RBF_Linear`` kernels.

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** 1
===================  =========================  =======================

``RBF_LENGTHSCALE_BND_XXX``
   Lengthscale hyperparameter for the RBF kernel, expected as :math:`\ell^2`
   (squared lengthscale as exported from scikit-learn).
   Required for ``RBF`` and ``RBF_Linear`` kernels.

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** 1
===================  =========================  =======================

``LINEAR_VARIANCE_BND_XXX``
   Variance hyperparameter for the linear kernel.
   Required for ``RBF_Linear`` and ``MLP_Linear`` kernels.

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** 1
===================  =========================  =======================

``GAUSSIAN_NOISE_VARIANCE_BND_XXX``
   Noise variance hyperparameter :math:`\sigma_n^2` added to the kernel
   diagonal for numerical stability and observation noise modeling.
   Required for all kernel types.

===================  =========================  ====================
**Type:** double     **Range:** :math:`> 0`     **Length:** 1
===================  =========================  ====================

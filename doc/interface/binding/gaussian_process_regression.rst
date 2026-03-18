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

``CP_VALS``
   Flattened pore-phase concentration training inputs used by the Gaussian
   process regression model. The values represent the training input points
   :math:`X` used to evaluate the kernel function. The array is interpreted
   according to ``NDIM``.

**Unit:** :math:`mol~m_{MP}^{-3}`

===================  =========================  =======================
**Type:** double     **Range:** unrestricted    **Length:** NTRAIN * NDIM
===================  =========================  =======================

``CS_VALS``
   Solid-phase training targets corresponding to ``CP_VALS``.
   These values form the training output vector used to compute the GPR
   coefficient vector :math:`\alpha = (K + \sigma_n^2 I)^{-1} y`.

**Unit:** :math:`mol~m_{SP}^{-3}`

===================  =========================  =======================
**Type:** double     **Range:** unrestricted    **Length:** NTRAIN
===================  =========================  =======================

``TRAINED_PARAMS``
   Trained kernel hyperparameters of the Gaussian process regression model.
   The parameters are expected in the following order:

   - MLP weight variance
   - MLP bias variance
   - MLP variance
   - Linear variance
   - RBF variance
   - RBF lengthscale
   - Gaussian noise variance

   All entries must be provided, regardless of the selected kernel.

===================  =========================  =======================
**Type:** double     **Range:** kernel-dependent **Length:** 7
===================  =========================  =======================

``KERNEL``
   Selects the kernel function used by the Gaussian process regression
   model. Supported values are ``MLP``, ``RBF``, ``RBF_Linear``, and
   ``MLP_Linear``.

===================  ================================================  =========
**Type:** string     **Range:** {MLP, RBF, RBF_Linear, MLP_Linear}     **Length:** 1
===================  ================================================  =========

``NDIM``
   Number of input dimensions per training point used in
   ``CP_VALS``. Must be a positive integer.

===================  =========================  =======================
**Type:** int        **Range:** :math:`\geq 1`  **Length:** 1
===================  =========================  =======================

``GPR_KKIN``
   Linear-driving-force coefficients in component-major ordering.

**Unit:** :math:`s^{-1}`

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** NTOTALBND
===================  =========================  =======================

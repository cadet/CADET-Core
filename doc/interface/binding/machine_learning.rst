.. _machine_learning_config:

Machine Learning Binding
~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption - ADSORPTION_MODEL = MACHINE_LEARNING**

For information on model equations, refer to :ref:`machine_learning`.

``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =======================
**Type:** int        **Range:** {0,1}           **Length:** 1/NTOTALBND
===================  =========================  =======================

``ML_KKIN``
   Linear-driving-force coefficients in component-major ordering.
   Controls the rate at which the solid-phase loading approaches the
   neural network-predicted equilibrium.

**Unit:** :math:`s^{-1}`

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** NCOMP
===================  =========================  =======================

``LAYERS``
   Number of hidden layers in the neural network architecture.
   Currently supports 1 or 2 hidden layers.

===================  =========================  =======================
**Type:** int        **Range:** {1, 2}          **Length:** 1
===================  =========================  =======================

``NODES``
   Number of nodes per hidden layer. All hidden layers have the same
   number of nodes.

===================  =========================  =======================
**Type:** int        **Range:** :math:`\geq 1`  **Length:** 1
===================  =========================  =======================

``NORM_FACTOR``
   Normalization factors applied element-wise to the pore-phase
   concentration before feeding into the neural network. Each entry
   corresponds to one input dimension (component).

===================  =========================  =======================
**Type:** double     **Range:** :math:`> 0`     **Length:** NCOMP
===================  =========================  =======================

``POROSITY_FACTOR``
   Scaling factor applied to the neural network prediction. This can be
   used to account for porosity differences or unit conversions between
   training data and simulation conditions.

===================  =========================  =======================
**Type:** double     **Range:** :math:`> 0`     **Length:** 1
===================  =========================  =======================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/model_weights**

Neural network weights and biases are organized hierarchically by layer. All weight matrices must be stored in column-major (Fortran) order, matching the export format from common machine learning frameworks.

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/model_weights/layer_0**

First hidden layer parameters.

``KERNEL``
   Weight matrix :math:`W_1` connecting input to first hidden layer.
   Shape: (NODES x NCOMP). Stored in column-major order.

===================  =========================  =======================
**Type:** double     **Range:** unrestricted    **Length:** NODES * NCOMP
===================  =========================  =======================

``BIAS``
   Bias vector :math:`b_1` for first hidden layer.

===================  =========================  =======================
**Type:** double     **Range:** unrestricted    **Length:** NODES
===================  =========================  =======================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/model_weights/layer_1**

Second hidden layer parameters (for LAYERS=2) or output layer parameters (for LAYERS=1).

``KERNEL``
   Weight matrix :math:`W_2`.

   - For LAYERS=1: Shape (1 x NODES), connects hidden to output.
   - For LAYERS=2: Shape (NODES x NODES), connects first to second hidden layer.

   Stored in column-major order.

===================  =========================  ============================
**Type:** double     **Range:** unrestricted    **Length:** NODES (LAYERS=1)
                                                            or NODES*NODES (LAYERS=2)
===================  =========================  ============================

``BIAS``
   Bias vector :math:`b_2`.

   - For LAYERS=1: Output bias (length 1).
   - For LAYERS=2: Second hidden layer bias (length NODES).

===================  =========================  =======================
**Type:** double     **Range:** unrestricted    **Length:** 1 (LAYERS=1)
                                                            or NODES (LAYERS=2)
===================  =========================  =======================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/model_weights/layer_2**

Output layer parameters (only for LAYERS=2).

``KERNEL``
   Weight matrix :math:`W_3` connecting second hidden layer to output.
   Shape: (1 x NODES). Stored in column-major order.

===================  =========================  =======================
**Type:** double     **Range:** unrestricted    **Length:** NODES
===================  =========================  =======================

``BIAS``
   Bias scalar :math:`b_3` for output layer.

===================  =========================  =======================
**Type:** double     **Range:** unrestricted    **Length:** 1
===================  =========================  =======================
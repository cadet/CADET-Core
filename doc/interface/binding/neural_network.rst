.. _neural_network_config:

Neural Network Binding
~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption - ADSORPTION_MODEL = NEURAL_NETWORK**

For information on model equations, refer to :ref:`neural_network`.

``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =======================
**Type:** int        **Range:** {0,1}           **Length:** 1/NTOTALBND
===================  =========================  =======================

``NN_KKIN``
   Linear-driving-force coefficients in component-major ordering.
   Controls the rate at which the solid-phase loading approaches the
   neural network-predicted equilibrium.

**Unit:** :math:`s^{-1}`

===================  =========================  =======================
**Type:** double     **Range:** :math:`\geq 0`  **Length:** NCOMP
===================  =========================  =======================

``NLAYERS``
   Number of hidden layers in the neural network architecture.
   Currently supports 1 or 2 hidden layers.

===================  =========================  =======================
**Type:** int        **Range:** {1, 2}          **Length:** 1
===================  =========================  =======================

``NNODES``
   Number of nodes per hidden layer. All hidden layers have the same
   number of nodes.

===================  =========================  =======================
**Type:** int        **Range:** :math:`\geq 1`  **Length:** 1
===================  =========================  =======================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_YYY**

``NORM_FACTOR``
   Normalization factors applied pore-phase concentration before feeding into the neural network.

===================  =========================  =======================
**Type:** double     **Range:** :math:`> 0`     **Length:** 1
===================  =========================  =======================

``POROSITY_FACTOR``
   Scaling factor applied to the neural network prediction. This can be
   used to account for porosity differences or unit conversions between
   training data and simulation conditions.

===================  =========================  =======================
**Type:** double     **Range:** :math:`> 0`     **Length:** 1
===================  =========================  =======================

Neural network weights and biases are organized hierarchically by layer.
All weight matrices must be stored in column-major (Fortran) order.

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_YYY/layer_0**

First hidden layer parameters.

``KERNEL``
   Weight matrix :math:`W_1` connecting input to first hidden layer.
   Shape: (NNODES x NCOMP). Stored in column-major order.

===================  =============================  ===========================
**Type:** double     **Range:** :math:`\mathbb{R}`    **Length:** NNODES * NCOMP
===================  =============================  ===========================

``BIAS``
   Bias vector :math:`b_1` for first hidden layer.

===================  =============================  =======================
**Type:** double     **Range:** :math:`\mathbb{R}`    **Length:** NNODES
===================  =============================  =======================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_YYY/layer_1**

Second hidden layer parameters (for NLAYERS=2) or output layer parameters (for NLAYERS=1).

``KERNEL``
   Weight matrix :math:`W_2`.

   - For NLAYERS=1: Shape (1 x NNODES), connects hidden to output.
   - For NLAYERS=2: Shape (NNODES x NNODES), connects first to second hidden layer.

   Stored in column-major order.

===================  =============================  ================================
**Type:** double     **Range:** :math:`\mathbb{R}`    **Length:** NNODES (NLAYERS=1)
                                                        or NNODES*NNODES (NLAYERS=2)
===================  =============================  ================================

``BIAS``
   Bias vector :math:`b_2`.

   - For NLAYERS=1: Output bias (length 1).
   - For NLAYERS=2: Second hidden layer bias (length NNODES).

===================  =============================  ===========================
**Type:** double     **Range:** :math:`\mathbb{R}`    **Length:** 1 (NLAYERS=1)
                                                          or NNODES (NLAYERS=2)
===================  =============================  ===========================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_YYY/layer_2**

Output layer parameters (only for NLAYERS=2).

``KERNEL``
   Weight matrix :math:`W_3` connecting second hidden layer to output.
   Shape: (1 x NNODES). Stored in column-major order.

===================  =============================  =======================
**Type:** double     **Range:** :math:`\mathbb{R}`    **Length:** NNODES
===================  =============================  =======================

``BIAS``
   Bias scalar :math:`b_3` for output layer.

===================  =============================  ================
**Type:** double     **Range:** :math:`\mathbb{R}`    **Length:** 1
===================  =============================  ================
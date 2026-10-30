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
   Number of hidden layers :math:`N` in the neural network architecture.

===================  =========================  =======================
**Type:** int        **Range:** :math:`\geq 1`  **Length:** 1
===================  =========================  =======================

``NNODES``
   Number of nodes per hidden layer.
   All hidden layers share the same width.

===================  =========================  =======================
**Type:** int        **Range:** :math:`\geq 1`  **Length:** 1
===================  =========================  =======================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_XXX/layer_0**

``NORM_FACTOR``
   Normalization factors applied to the pore-phase
   concentration before feeding into the neural network.

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

Neural network weights and biases are organized hierarchically by bound state and layer.
For each bound state, there are ``NLAYERS + 1`` layer groups:
``layer_0`` through ``layer_{NLAYERS-1}`` are the hidden layers and
``layer_{NLAYERS}`` is the output layer. All weight matrices must be stored in
column-major (Fortran) order.

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_XXX/layer_0**

First hidden layer parameters. Always present.

``KERNEL``
   Weight matrix :math:`W_1` connecting the input to the first hidden layer.
   Shape: (NNODES x NCOMP). Stored in column-major order.

===================  =============================  ===========================
**Type:** double     **Range:** :math:`\mathbb{R}`  **Length:** NNODES * NCOMP
===================  =============================  ===========================

``BIAS``
   Bias vector :math:`b_1` for the first hidden layer.

===================  =============================  =======================
**Type:** double     **Range:** :math:`\mathbb{R}`  **Length:** NNODES
===================  =============================  =======================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_XXX/layer_l**
**(for l = 1, ..., NLAYERS-1)**

Intermediate hidden layer parameters. Present only when ``NLAYERS >= 2``.
Layer ``layer_l`` connects hidden layer :math:`l` to hidden layer :math:`l+1`.

``KERNEL``
   Weight matrix :math:`W_{l+1}` of shape (NNODES x NNODES). Stored in column-major order.

===================  =============================  ========================
**Type:** double     **Range:** :math:`\mathbb{R}`  **Length:** NNODES*NNODES
===================  =============================  ========================

``BIAS``
   Bias vector :math:`b_{l+1}` for hidden layer :math:`l+1`.

===================  =============================  =======================
**Type:** double     **Range:** :math:`\mathbb{R}`  **Length:** NNODES
===================  =============================  =======================

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_XXX/layer_NLAYERS**

Output layer parameters. Always present. Index equals ``NLAYERS``.

``KERNEL``
   Weight matrix :math:`W_{N+1}` connecting the last hidden layer to the scalar output.
   Shape: (1 x NNODES). Stored in column-major order.

===================  =============================  =======================
**Type:** double     **Range:** :math:`\mathbb{R}`  **Length:** NNODES
===================  =============================  =======================

``BIAS``
   Bias scalar :math:`b_{N+1}` for the output layer.

===================  =============================  ================
**Type:** double     **Range:** :math:`\mathbb{R}`  **Length:** 1
===================  =============================  ================
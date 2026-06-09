.. _neural_network:

Neural Network Binding
~~~~~~~~~~~~~~~~~~~~~~

The Neural Network (NN) binding model is a non-parametric, data-driven binding model that represents the equilibrium relation between pore-phase and solid-phase concentration using a trained feedforward neural network. The model predicts an equilibrium loading from a neural network trained on experimental or simulated adsorption data.

For each component :math:`i` and bound state :math:`m`, an equilibrium loading

.. math::

    c^{s,\ast}_{i,m} = f_{i,m}(c^p)

is constructed from user-provided trained neural network weights and biases. In the current implementation, the neural network is a feedforward architecture with exponential linear unit (ELU) activations, supporting either one or two hidden layers.

Network Architecture
********************

The neural network predictor uses the following architectures:

**Single hidden layer:**

.. math::

    f(x) = W_2 \cdot \text{ELU}(W_1 \cdot x + b_1) + b_2

**Two hidden layers:**

.. math::

    f(x) = W_3 \cdot \text{ELU}(W_2 \cdot \text{ELU}(W_1 \cdot x + b_1) + b_2) + b_3

where :math:`W_i` denote weight matrices, :math:`b_i` denote bias vectors, and the exponential linear unit (ELU) activation function is defined as:

.. math::

    \text{ELU}(z) = 
    \begin{cases}
        z & \text{if } z \geq 0 \\
        e^z - 1 & \text{if } z < 0
    \end{cases}

The input to the neural network is the normalized pore-phase concentration vector:

.. math::

    x = c^p \odot \alpha_{\text{norm}}

where :math:`\alpha_{\text{norm}}` is the normalization factor specified by ``NORM_FACTOR``, and :math:`\odot` denotes element-wise multiplication.

The network output is scaled by a porosity factor and shifted by an offset:

.. math::

    c^{s,\ast}(c^p) = \beta_{\text{poros}} \cdot \left(f(c^p \odot \alpha_{\text{norm}}) - f(0)\right)

where :math:`\beta_{\text{poros}}` is specified by ``POROSITY_FACTOR``.
The offset :math:`f(0)` is computed once during configuration to ensure that the predicted loading is zero when the pore-phase concentration is zero (enforcing physically valid boundary conditions).

Normalization
*************

The normalization factor :math:`\alpha_{\text{norm}}` scales the liquid-phase concentrations before they are passed to the ANN.
The ANN is trained and evaluated on these normalized concentrations rather than on the raw concentrations.
Normalization factors should be chosen such that the normalized ANN inputs are typically of order 1.
This generally improves training stability, prediction accuracy, and robustness of the nonlinear solver.
Values that are excessively large or small should be avoided.

**With** mechanistic knowledge, it may be practical to choose normalization factors that reflect characteristic concentration scales of the respective bound state (e.g., equilibrium constants, affinity parameters, or other physically meaningful scaling factors).
**Without** mechanistic knowledge, it may be practical to choose normalization factors purely for numerical conditioning. A common choice is to scale each input component to values of order unity, for example :math:`\alpha_{\text{norm}} = max(c^p)`

Porosity factor
***************

The porosity scaling factor :math:`\beta_{\text{poros}}` scales this prediction to the physical solid-phase loading used by the binding model.
Since equilibrium loadings are often reported on different reference volumes or masses, :math:`\beta_{\text{poros}}` is intended to be used as a simple correction factor to account for differences in porosity or unit conversions between the training data and the simulation conditions.

Kinetic Form
************

The model is used in a kinetic linear-driving-force form. For each component :math:`i` and bound state :math:`m`, the exchange term is based on the deviation of the current solid-phase loading :math:`c^s_{i,m}` from the neural network-predicted equilibrium loading :math:`c^{s,\ast}_{i,m}`:

.. math::

    \frac{\partial c^s_{i,m}}{\partial t}
    =
    k^{\mathrm{kin}}_{i,m}\left(c^{s,\ast}_{i,m}(c^p) - c^s_{i,m}\right).

Equivalently, in residual form the implementation evaluates

.. math::

    r_{i,m}
    =
    k^{\mathrm{kin}}_{i,m}\left(c^s_{i,m} - c^{s,\ast}_{i,m}(c^p)\right).

Thus, the neural network provides the equilibrium target, while the kinetic constant :math:`k^{\mathrm{kin}}_{i,m}` controls how fast the equilibrium state is approached. The kinetic parameter is configured through ``NN_KKIN``.

Jacobian Computation
********************

An analytical Jacobian is implemented, where the gradient with respect to the pore-phase concentration is computed as:

**Single hidden layer:**

.. math::

    \frac{\partial c^{s,\ast}}{\partial c^p} = W_1^\top \cdot \text{diag}(\text{ELU}'(z_1)) \cdot W_2^\top

**Two hidden layers:**

.. math::

    \frac{\partial c^{s,\ast}}{\partial c^p} = W_1^\top \cdot \text{diag}(\text{ELU}'(z_1)) \cdot W_2^\top \cdot \text{diag}(\text{ELU}'(z_2)) \cdot W_3^\top

where :math:`z_i` denotes the pre-activation at layer :math:`i`, and the ELU derivative is:

.. math::

    \text{ELU}'(z) = 
    \begin{cases}
        1 & \text{if } z \geq 0 \\
        e^z & \text{if } z < 0
    \end{cases}

The chain rule accounts for input normalization and output scaling when computing the full Jacobian contribution to the binding residual.

For more information on model parameters required to define in CADET file format, see :ref:`neural_network_config`.
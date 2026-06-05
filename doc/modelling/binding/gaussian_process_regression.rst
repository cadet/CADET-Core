.. _gaussian_process_regression:

Gaussian Process Regression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Gaussian process regression (GPR) model is a non-parametric, data-driven binding model that represents the equilibrium relation between pore-phase and solid-phase concentration by Gaussian process regression trained on tabulated data.
The model predicts equilibrium loadings from training inputs in pore-phase concentration space and corresponding solid-phase loadings.

Model Structure
***************

The GPR binding model creates **one independent GPR model per bound state**.
For a system with multiple components where each component may have multiple bound states, each bound state :math:`q_{i,m}` (component :math:`i`, bound state :math:`m`) is predicted by its own GPR model with its own trained hyperparameters.

For each bound state at global index :math:`\text{j}`, an equilibrium loading

.. math::

    c^{s,\ast}_{\text{j}} = f_{\text{j}}(c^p)

is constructed from user-provided training data:

- Pore-phase inputs (shared across all bound states) and solid-phase targets (one column per bound state)
- Bound-state-specific trained kernel hyperparameters

GPR Prediction
**************

The GPR predictor evaluates the equilibrium loading in the standard kernel form

.. math::

    c^{s,\ast}_{\text{j}}(c^p) = k(c^p, X)^\top \alpha_{\text{j}},

where:

- :math:`X` denotes the training inputs (shared pore-phase concentrations from ``CP_VALS``)
- :math:`k(c^p, X)` is the vector of kernel evaluations between the current pore-phase concentration and all training samples, using the kernel type and hyperparameters specific to bound state ``bndIdx``
- :math:`\alpha_{\text{j}}` is the coefficient vector obtained from the linear system for this bound state:

.. math::

    \alpha_{\text{j}} = \left(K_{\text{j}}(X,X) + \sigma_{n,\text{j}}^2 I\right)^{-1} y_{\text{j}}.

Here:

- :math:`K_{\text{j}}(X,X)` is the kernel matrix assembled from the training inputs using the kernel specific to bound state ``bndIdx``
- :math:`y_{\text{j}}` is the vector of solid-phase training values for this bound state (extracted from the corresponding column of ``CS_VALS``)
- :math:`\sigma_{n,\text{j}}^2` is the Gaussian noise variance for this bound state (``GAUSSIAN_NOISE_VARIANCE_BNDXXX``)

Kernel Functions
****************

Depending on the selected kernel, the covariance function is given by one of the following forms.

For the radial basis function kernel,

.. math::

    k_{\mathrm{RBF}}(x,x')
    =
    \sigma_{\mathrm{RBF}}^2
    \exp\!\left(
        -\frac{\lVert x-x' \rVert^2}{2\,\ell_{\mathrm{RBF}}}
    \right).

For the linear kernel,

.. math::

    k_{\mathrm{Lin}}(x,x')
    =
    \sigma_{\mathrm{Lin}}^2 \, x^\top x'.

For the multilayer perceptron kernel,

.. math::

    k_{\mathrm{MLP}}(x,x')
    =
    \sigma_{\mathrm{MLP}}^2 \frac{2}{\pi}
    \arcsin\!\left(
        \frac{\sigma_w^2 x^\top x' + \sigma_b^2}
        {\sqrt{\sigma_w^2 x^\top x + \sigma_b^2 + 1}
         \sqrt{\sigma_w^2 x'^\top x' + \sigma_b^2 + 1}}
    \right).

In addition, the implementation supports additive mixed kernels

.. math::

    k_{\mathrm{RBF+Lin}}(x,x') = k_{\mathrm{RBF}}(x,x') + k_{\mathrm{Lin}}(x,x'),

.. math::

    k_{\mathrm{MLP+Lin}}(x,x') = k_{\mathrm{MLP}}(x,x') + k_{\mathrm{Lin}}(x,x').

These kernel definitions are used both for prediction and for evaluation of the Jacobian contribution.

Offset Correction
*****************

An offset is computed once during configuration for each bound state as the GPR prediction at zero input and is subtracted from subsequent predictions:

.. math::

    c^{s,\ast}_{\text{j}}(c^p) = f_{\text{j}}(c^p) - f_{\text{j}}(0).

This shifts the model response such that the predicted loading is zero at the origin. This correction is applied independently for each bound state.

Kinetic Form
************

The model is used in a kinetic linear-driving-force form. For each bound state at global index ``bndIdx`` belonging to component :math:`i`, the exchange term is based on the deviation of the current solid-phase loading :math:`c^s_{\text{j}}` from the GPR-predicted equilibrium loading :math:`c^{s,\ast}_{\text{j}}`:

.. math::

    \frac{\partial c^s_{\text{j}}}{\partial t}
    =
    k^{\mathrm{kin}}_{i}\left(c^{s,\ast}_{\text{j}}(c^p) - c^s_{\text{j}}\right).

Equivalently, in residual form the implementation evaluates

.. math::

    r_{\text{j}}
    =
    k^{\mathrm{kin}}_{i}\left(c^s_{\text{j}} - c^{s,\ast}_{\text{j}}(c^p)\right).

Thus, the Gaussian process regression model provides the equilibrium target for each bound state, while the kinetic constant :math:`k^{\mathrm{kin}}_{i}` (shared across all bound states of component :math:`i`) controls how fast the equilibrium state is approached. The kinetic parameter is configured through ``GPR_KKIN`` in component-major ordering.


For more information on model parameters required to define in CADET file format, see :ref:`gaussian_process_regression_config`.
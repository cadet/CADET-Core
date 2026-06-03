.. _gaussian_process_regression:

Gaussian Process Regression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The Gaussian process regression (GPR) model is a non-parametric, data-driven binding model that represents the equilibrium relation between pore-phase and solid-phase concentration by Gaussian process regression trained on tabulated data. The model predicts an equilibrium loading from training inputs in pore-phase concentration space and corresponding solid-phase loadings.

For each component :math:`i` and bound state :math:`m`, an equilibrium loading

.. math::

    c^{s,\ast}_{i,m} = f_{i,m}(c^p)

is constructed from user-provided training data and trained kernel parameters. In the current implementation, the prediction is evaluated from pore-phase training inputs ``CP_VALS``, solid-phase targets ``CS_VALS``, and trained kernel hyperparameters ``TRAINED_PARAMS`` provided in the ``training_data`` scope. The input dimension is specified by ``NDIM`` and the kernel type by ``KERNEL``. Supported kernels are ``MLP``, ``RBF``, ``RBF_Linear``, and ``MLP_Linear``.

The GPR predictor evaluates the equilibrium loading in the standard kernel form

.. math::

    c^{s,\ast}(c^p) = k(c^p, X)^\top \alpha,

where :math:`X` denotes the training inputs, :math:`k(c^p, X)` is the vector of kernel evaluations between the current pore-phase concentration and all training samples, and :math:`\alpha` is obtained from the linear system

.. math::

    \alpha = \left(K(X,X) + \sigma_n^2 I\right)^{-1} y.

Here, :math:`K(X,X)` is the kernel matrix assembled from the training inputs, :math:`y` is the vector of solid-phase training values, and :math:`\sigma_n^2` is the Gaussian noise variance added to the diagonal before Cholesky factorization.

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

An offset is computed once during configuration as the GPR prediction at zero input and is subtracted from subsequent predictions. Thus, the equilibrium loading used by the binding model is

.. math::

    c^{s,\ast}(c^p) = f(c^p) - f(0).

This shifts the model response such that the predicted loading is zero at the chosen reference point.

Kinetic form
************

The model is used in a kinetic linear-driving-force form. For each component :math:`i` and bound state :math:`m`, the exchange term is based on the deviation of the current solid-phase loading :math:`c^s_{i,m}` from the GPR-predicted equilibrium loading :math:`c^{s,\ast}_{i,m}`:

.. math::

    \frac{\partial c^s_{i,m}}{\partial t}
    =
    k^{\mathrm{kin}}_{i,m}\left(c^{s,\ast}_{i,m}(c^p) - c^s_{i,m}\right).

Equivalently, in residual form the implementation evaluates

.. math::

    r_{i,m}
    =
    k^{\mathrm{kin}}_{i,m}\left(c^s_{i,m} - c^{s,\ast}_{i,m}(c^p)\right).

Thus, the Gaussian process regression model provides the equilibrium target, while the kinetic constant :math:`k^{\mathrm{kin}}_{i,m}` controls how fast the equilibrium state is approached.
The kinetic parameter is configured through ``GPR_KKIN``.


For more information on model parameters required to define in CADET file format, see :ref:`gaussian_process_regression_config`.
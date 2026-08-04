.. _parameter_dependencies:

Parameter Dependencies
======================

Many physical models require parameters that are not constant but depend on the operating conditions, spatial position, or the current solution state.
Parameter dependencies provide a general mechanism to express such relationships while preserving the underlying model formulation.
Instead of replacing governing equations, they modify selected model parameters by evaluating a prescribed functional relation.

A parameter dependency of a model parameter :math:`p` can generally be viewed as

.. math::

    p = f(\cdot),

where the function :math:`f` may depend on spatial coordinates, model parameters, solution variables, external functions, or combinations thereof.

Parameter dependencies enable the representation of empirical correlations and constitutive relations without introducing additional governing equations.


Current Implementation
----------------------

The dependency framework is designed as a generic extension mechanism for model parameters.
The currently implemented parameter dependencies represent selected examples that demonstrate the general concept and provide reusable implementation patterns for further extensions.
Only a limited subset of model parameters and dependency functions is currently supported.
Complete coverage of all model parameters is not intended.
The set of currently available parameter dependencies, supported functional forms, and corresponding configuration options is described in the interface documentation (see :ref:`parameter_dependencies_config`).
If a desired parameter dependency is not implemented, users are encouraged to `open an issue to CADET-Core <https://github.com/cadet/CADET-Core/issues>`_ or, even better, contribute a pull request.


Types of Dependencies
---------------------

Several categories of parameter dependencies can be distinguished based on the quantity on which a parameter depends.

External Function Dependence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A parameter may depend on an external function or profile :math:`F\colon \left[ 0, T_{\text{end}}\right] \times [0, L] \to \mathbb{R}`, where :math:`L` denotes the physical length of the unit operation, or :math:`T\colon \left[0, T_{\text{end}}\right] \to \mathbb{R}` if the unit operation model has no axial length.
The current implementation represents this dependence by a polynomial of third degree.
For example, the adsorption rate :math:`k_a` can be defined as

.. math::

    \begin{aligned}
        k_a(T) &= k_{a,3} T^3 + k_{a,2} T^2 + k_{a,1} T + k_{a,0}.
    \end{aligned}

This type of dependency is implemented in CADET for binding model parameters.
For more details, refer to :ref:`dependence-on-external-function_bind`.


Parameter-Parameter Dependence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A parameter may depend on another model parameter, such as film diffusion depending on the velocity in an axial flow column with constant porosity.
The parameter dependence is recomputed during every residual evaluation to account for changes of the independent parameter over time.

The independent parameter may itself vary spatially.
More generally, a parameter may depend explicitly on the spatial position or on another spatially varying parameter.
Examples include spatially varying porosity distributions, or parameters depending on the local velocity in a radial flow column.

From a modeling perspective, the dependent parameter then becomes a spatial field rather than a single scalar quantity.
Consequently, it must be represented on the computational mesh and handled consistently by the spatial discretization.
Evaluating the dependency only from (averaged) state values generally introduces an additional approximation and may reduce the accuracy of the numerical solution.


Parameter-State Dependence
^^^^^^^^^^^^^^^^^^^^^^^^^^

A parameter may depend on a solution variable of the model.
Typical examples include transport coefficients such as surface diffusion depending on salt concentration.
In this case, the parameter becomes part of the nonlinear model and must be re-evaluated whenever the solution changes.
State-dependent parameters may introduce additional nonlinearities, increasing the computational effort and the number of nonlinear iterations required for convergence during the nonlinear solution procedure.

If the state variable varies spatially, the dependent parameter likewise becomes a spatial field and must be represented on the computational mesh and handled consistently by the spatial discretization.
Evaluating the dependency only from (averaged) state values generally introduces an additional approximation and may reduce the accuracy of the numerical solution.

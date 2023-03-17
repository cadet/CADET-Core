.. _parameter_dependencies:

Parameter Dependencies
======================

Some parameters depend on other parameters (parameter-parameter dependency) or the solution variables (parameter-state dependency).
Parameter dependencies are defined in the unit operation scope.

Parameter-Parameter Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Group /input/model/unit_XXX
---------------------------

``COL_DISPERSION_DEP``

   Parameter dependence of column dispersion on the interstitial velocity. Available for the LRM, LRMP and GRM units (with FV discretization only at the moment)
   
   ================  =====================================  =============
   **Type:** string  **Range:** :math:`\texttt{POWER_LAW}`  **Length:** 1
   ================  =====================================  =============

``FILM_DIFFUSION_DEP``

   Parameter dependence of film diffusion on the interstitial velocity. Available for the LRMP unit (with FV discretization only at the moment)
   
   ================  =====================================  =============
   **Type:** string  **Range:** :math:`\texttt{POWER_LAW}`  **Length:** 1
   ================  =====================================  =============


**Correlations**
""""""""""""""""

Different types of parameter correlations are can be applied.
The following correlations can be used for all parameter-parameter dependencies, but we specify the required input fields only for ``COL_DISPERSION_DEP``, for the sake of conciseness.

**Power Law**

.. math::

    \begin{aligned}
        p_{dep} &= p_{dep} \cdot b \ |p_{on}^x|
    \end{aligned}

Here, :math:`p_{dep}` is the dependent parameter and :math:`p_{on}` is the parameter it depends on.

``COL_DISPERSION_DEP_BASE``

   Base :math:`b` of the power law parameter dependence. Optional, defaults to :math:`1.0`
   
   ================  =============================  =============
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  =============

``COL_DISPERSION_DEP_EXPONENT``

   Exponent :math:`x` of the power law parameter dependence
   
   ================  =============================  =============
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  =============

``COL_DISPERSION_DEP_ABS``

   Specifies whether or not the absolute value should be computed. Optional, defaults to :math:`1`
   
   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============


Parameter-State Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Group /input/model/unit_XXX
---------------------------

Parameter-State Dependencies are not fully implemented yet.

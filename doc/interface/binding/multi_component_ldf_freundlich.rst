.. _multi_component_ldf_freundlich_config:

Multi Component Linear Driving Force Freundlich
===============================================

**Group /input/model/unit_XXX/adsorption – ADSORPTION_MODEL = MULTI_COMPONENT_LDF_FREUNDLICH**

For information on model equations, refer to :ref:`multi_component_ldf_freundlich_model`.

``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =========================================
**Type:** int        **Range:** {0,1}  		    **Length:** 1/NTOTALBND
===================  =========================  =========================================

``MCLDFFRL_KLDF``
   Rate constants in linear driving force approach

**Unit:** :math:`s^{-1}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  =========================================

``MCLDFFRL_KF``
   Proportionality constants

**Unit:** :math:`m_{MP}^{3/n}~m_{SP}^{-3}~mol^{1-1/n}`

===================  =========================  ==================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** NCOMP
===================  =========================  ================================== 

``MCLDFFRL_EXP``
   Freundlich exponent

**Unit:** [-]

===================  =========================  ==================================
**Type:** double     **Range:** :math:`\gt 0`   **Length:** NCOMP
===================  =========================  ================================== 

``MCLDFFRL_A``
   Component influences in row-major ordering

**Unit:** [-]

===================  =========================  ==================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** :math:`\text{NCOMP}^2`
===================  =========================  ================================== 

``MCLDFFRL_TAU``
   Small constant that ensures numerical stability

**Unit:** :math:`mol~m_{MP}^{-3}`

===================  =========================  =========================================
**Type:** double     **Range:** :math:`\gt 0`   **Length:** 1
===================  =========================  =========================================

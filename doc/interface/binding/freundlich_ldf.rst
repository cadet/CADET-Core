.. _freundlich_ldf_config:

Freundlich LDF
~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/adsorption – ADSORPTION_MODEL = FREUNDLICH_LDF**


``IS_KINETIC``
   Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
   quasi-stationary. If a single value is given, the mode is set for all
   bound states. Otherwise, the adsorption mode is set for each bound
   state separately.

===================  =========================  =========================================
**Type:** int        **Range:** {0,1}           **Length:** 1/NTOTALBND
===================  =========================  =========================================

``FLDF_KKIN``
   Kinetic constants for each component


**Unit:** :math:`..`

===================  =========================  ==================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1/NTOTALBND
===================  =========================  ==================================  


``FLDF_KF``
   Freundlich coefficient for each component

**Unit:** :math:`mol^{1-n}g^{1}L^{-n}`

===================  =========================  ==================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1/NTOTALBND
===================  =========================  ==================================  


``FLDF_N``
   Freundlich exponent for each component

**Unit:** :math:`[-]`

===================  =========================  ==================================
**Type:** double     **Range:** :math:`\ge 0`   **Length:** 1/NTOTALBND
===================  =========================  ==================================  

     

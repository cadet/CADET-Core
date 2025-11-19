.. _activated_sludge_model:

Activated Sludge Model (ASM3h)
===============================

The Activated Sludge Model 3 extended (ASM3h) is a comprehensive reaction model for simulating biological 
wastewater treatment processes. It belongs to the Activated Sludge Models (ASM) family, a set of standardized 
models developed by the International Water Association (IWA). The ASM3h model describes the fate of 13 
biochemical components in wastewater treatment systems.

Model Components
~~~~~~~~~~~~~~~~

The ASM3h model tracks 13 biochemical components representing species and compounds in wastewater:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Abbreviation
     - Component
   * - :math:`SO`
     - Dissolved oxygen
   * - :math:`SS`
     - Readily biodegradable substrate
   * - :math:`SNH`
     - Ammonium
   * - :math:`SNO`
     - Nitrite and nitrate
   * - :math:`SN2`
     - Dinitrogen, released by denitrification
   * - :math:`SALK`
     - Alkalinity, bicarbonate
   * - :math:`SI`
     - Soluble inert organic
   * - :math:`XI`
     - Inert Particulate organic
   * - :math:`XS`
     - Slowly biodegradable substrate
   * - :math:`XH`
     - Heterotrophic biomass
   * - :math:`XSTO`
     - Organics stored by heterotrophs
   * - :math:`XA`
     - Autotrophic, nitrifying biomass
   * - :math:`XMI`
     - Mineral particulate matter from biomass

The net flux for each component is computed using the stoichiometric matrix :math:`S \in \mathbb{R}^{13 \times 13}` 
and a reaction rate vector :math:`\varphi_j(c)`, where subscript :math:`j` denotes the reaction index. 
The ASM3h model comprises 13 biochemical reactions detailed below.

Biochemical Reaction Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
    :header-rows: 1
    :widths: 15 35 50

    * - Reaction Number
      - Reaction Description
      - Rate Equation
    * - :math:`r_1`
      - Hydrolysis of organic structures
      - :math:`k_{h,20} \cdot f_T^{0.04} \cdot \frac{XS/XH}{XS/XH + K_X} \cdot XH`
    * - :math:`r_2`
      - Aerobic storage of SS
      - :math:`k_{STO} \cdot \frac{SO}{SO + K_{h,O2}} \cdot \frac{SS}{SS + K_{h,SS}} \cdot XH`
    * - :math:`r_3`
      - Anoxic storage of SS
      - :math:`k_{STO} \cdot \eta_{h,NO3} \cdot \frac{K_{h,O2}}{SO + K_{h,O2}} \cdot \frac{SNO}{SNO + K_{h,NO3}} \cdot \frac{SS}{SS + K_{h,SS}} \cdot XH`
    * - :math:`r_4`
      - Aerobic growth of XH
      - :math:`\mu_H \cdot \frac{SO}{SO + K_{h,O2}} \cdot \frac{XSTO/XH}{XSTO/XH + K_{h,STO}} \cdot \frac{SNH}{SNH + K_{h,NH4}} \cdot \frac{SALK}{SALK + K_{h,ALK}} \cdot XH`
    * - :math:`r_5`
      - Anoxic growth of XH (denitrification)
      - :math:`\mu_H \cdot \eta_{h,NO3} \cdot \frac{K_{h,O2}}{SO + K_{h,O2}} \cdot \frac{SNO}{SNO + K_{h,NO3}} \cdot \frac{XSTO/XH}{XSTO/XH + K_{h,STO}} \cdot \frac{SNH}{SNH + K_{h,NH4}} \cdot \frac{SALK}{SALK + K_{h,ALK}} \cdot XH`
    * - :math:`r_6`
      - Aerobic endogenous respiration of XH
      - :math:`b_H \cdot \frac{SO}{SO + K_{h,O2}} \cdot XH`
    * - :math:`r_7`
      - Anoxic endogenous respiration of XH
      - :math:`b_H \cdot \eta_{h,end} \cdot \frac{K_{h,O2}}{SO + K_{h,O2}} \cdot \frac{SNO}{SNO + K_{h,NO3}} \cdot XH`
    * - :math:`r_8`
      - Aerobic respiration of internal cell storage products
      - :math:`b_H \cdot \frac{SO}{SO + K_{h,O2}} \cdot \frac{XSTO}{XSTO + K_{h,STO}} \cdot XH`
    * - :math:`r_9`
      - Anoxic respiration of internal cell storage products
      - :math:`b_H \cdot \eta_{h,end} \cdot \frac{K_{h,O2}}{SO + K_{h,O2}} \cdot \frac{SNO}{SNO + K_{h,NO3}} \cdot \frac{XSTO}{XSTO + K_{h,STO}} \cdot XH`
    * - :math:`r_{10}`
      - Aerobic growth of XA
      - :math:`\mu_{AUT} \cdot \frac{SO}{SO + K_{N,O2}} \cdot \frac{SNH}{SNH + K_{N,NH4}} \cdot \frac{SALK}{SALK + K_{N,ALK}} \cdot XA`
    * - :math:`r_{11}`
      - Aerobic endogenous respiration of XA
      - :math:`b_{AUT} \cdot \frac{SO}{SO + K_{N,O2}} \cdot XA`
    * - :math:`r_{12}`
      - Anoxic endogenous respiration of XA
      - :math:`b_{AUT} \cdot \eta_{N,end} \cdot \frac{K_{N,O2}}{SO + K_{N,O2}} \cdot \frac{SNO}{SNO + K_{h,NO3}} \cdot XA`
    * - :math:`r_{13}`
      - Aeration
      - :math:`\frac{iO_2}{1000}`
    * - :math:`r_{14}`
      - Anoxic storage of SS\ :sub:`ad` (optional)
      - :math:`k_{STO} \cdot \eta_{h,NO3} \cdot \frac{K_{h,O2}}{SO + K_{h,O2}} \cdot \frac{SS_{ad}}{SS + SS_{ad} + K_{h,SS}} \cdot \frac{SNO}{SNO + K_{h,NO3}} \cdot XH`
    * - :math:`r_{15}`
      - Aerobic storage of SS\ :sub:`ad` (optional)
      - :math:`k_{STO} \cdot \frac{SO}{SO + K_{h,O2}} \cdot \frac{SS_{ad}}{SS + SS_{ad} + K_{h,SS}} \cdot XH`

Reactions :math:`r_{14}` and :math:`r_{15}` are only active when 15 components are defined (including SS\ :sub:`ad` and SI\ :sub:`ad`).


Stoichiometric Matrix
~~~~~~~~~~~~~~~~~~~~~

The stoichiometric matrix :math:`S` encodes the mass balance coefficients for each biochemical reaction:

.. list-table::
    :header-rows: 1
    :widths: 10 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7

    * - Component\\Reaction
      - :math:`r_1`
      - :math:`r_2`
      - :math:`r_3`
      - :math:`r_4`
      - :math:`r_5`
      - :math:`r_6`
      - :math:`r_7`
      - :math:`r_8`
      - :math:`r_9`
      - :math:`r_{10}`
      - :math:`r_{11}`
      - :math:`r_{12}`
      - :math:`r_{13}`
      - :math:`r_{14}^*`
      - :math:`r_{15}^*`
    * - :math:`SO`
      - 0
      - :math:`Y_{STO,aer} - 1`
      - :math:`Y_{STO,anox} - 1`
      - :math:`1 - 1/Y_{H,aer}`
      - :math:`1 - 1/Y_{H,anox}`
      - :math:` f_{XI} - 1`
      - :math:`f_{XI} - 1`
      - -1
      - -1
      - :math:`-(64/14) \cdot 1/Y_A + 1`
      - :math:`f_{XI} - 1`
      - :math:`f_{XI} - 1`
      - 1
      - 0
      - 0
    * - :math:`SS`
      - :math:`(1 - f_{SI})(1 - f_{SS,ad})`
      - -1
      - -1
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
    * - :math:`SS_{ad}^*`
      - :math:`(1 - f_{SI}) \cdot f_{SS,ad}`
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - -1
      - -1
    * - :math:`SNH`
      - :math:`i_{N,XS} - i_{N,SI} \cdot f_{SI} - (1 - f_{SI}) \cdot i_{N,SS}`
      - :math:`i_{N,SS}`
      - :math:`i_{N,SS}`
      - :math:`-i_{N,BM}`
      - :math:`-i_{N,BM}`
      - :math:`i_{N,BM} - f_{XI} \cdot i_{N,XI}`
      - :math:`i_{N,BM} -f_{XI} \cdot i_{N,XI}`
      - 0
      - 0
      - :math:`-1/Y_A - i_{N,BM}`
      - :math:`f_{XI} \cdot i_{N,XI} + i_{N,BM}`
      - :math:`f_{XI} \cdot i_{N,XI} + i_{N,BM}`
      - 0
      - 0
      - 0
    * - :math:`SNO`
      - 0
      - 0
      - :math:`(Y_{STO,anox} - 1) / (40/14)`
      - 0
      - :math:`(1 - 1/Y_{H,anox}) / (40/14)`
      - 0
      - :math:`(f_{XI} - 1) / (40/14)`
      - 0
      - :math:`-14/40`
      - :math:`1/Y_A`
      - 0
      - :math:`(f_{XI} - 1) / (40/14)`
      - 0
      - 0
      - 0
    * - :math:`SN2`
      - 0
      - 0
      - :math:`(1 - Y_{STO,anox}) / (40/14)`
      - 0
      - :math:`(1/Y_{H,anox} - 1) / (40/14)`
      - 0
      - :math:`(-f_{XI} - 1) / (40/14)`
      - 0
      - :math:`--14/40`
      - 0
      - 0
      - :math:`(1 + f_{XI} ) / (40/14)`
      - 0
      - 0
      - 0
    * - :math:`SALK`
      - :math:`(i_{N,XS} - i_{N,SI} \cdot f_{SI} - (1 - f_{SI}) \cdot i_{N,SS})/14`
      - :math:`i_{N,SS} / 14`
      - :math:`i_{N,SS} / 14 - (Y_{STO,anox} - 1) / 40`
      - :math:`-i_{N,BM} / 14`
      - :math:`-i_{N,BM} / 14 - (1 - 1/Y_{H,anox}) / 40`
      - :math:`(-f_{XI} \cdot i_{N,XI} + i_{N,BM}) / 14`
      - :math:`(-f_{XI} \cdot i_{N,XI} + i_{N,BM}) / 14  (-f_{XI} + 1) / 40`
      - 0
      - :math:`1/40`
      - :math:`-(1/Y_A) 1/7- i_{N,BM}/14`
      - :math:`f_{XI} \cdot i_{N,XI} + i_{N,BM} / 14`
      - :math:`(-f_{XI} \cdot i_{N,XI} + i_{N,BM})/14 - (f_{XI} - 1) / 40`
      - 0
      - 0
      - 0
    * - :math:`SI`
      - :math:`f_{SI}(1 - f_{SI,ad})`
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
    * - :math:`SI_{ad}^*`
      - :math:`f_{SI} \cdot f_{SI,ad}`
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
    * - :math:`XI`
      - 0
      - 0
      - 0
      - 0
      - 0
      - :math:`f_{XI}`
      - :math:`f_{XI}`
      - 0
      - 0
      - 0
      - :math:`f_{XI}`
      - :math:`f_{XI}`
      - 0
      - 0
      - 0
    * - :math:`XS`
      - -1
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
    * - :math:`XH`
      - 0
      - 0
      - 0
      - 1
      - 1
      - -1
      - -1
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
    * - :math:`XSTO`
      - 0
      - :math:`Y_{STO,aer}`
      - :math:`Y_{STO,anox}`
      - :math:`-1/Y_{H,aer}`
      - :math:`-1/Y_{H,anox}`
      - 0
      - 0
      - -1
      - -1
      - 0
      - 0
      - 0
      - 0
      - :math:`Y_{STO,anox}`
      - :math:`Y_{STO,aer}`
    * - :math:`XA`
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 0
      - 1
      - -1
      - -1
      - 0
      - 0
      - 0
    * - :math:`XMI`
      - 0
      - 0
      - 0
      - 0
      - 0
      - :math:`f_{XMI,BI}`
      - :math:`f_{XMI,BI}`
      - 0
      - 0
      - 0
      - :math:`f_{XMI,BI}`
      - :math:`f_{XMI,BI}`
      - 0
      - 0
      - 0

:sup:`*` Optional: only active when 15 components are defined (SS\ :sub:`ad` and SI\ :sub:`ad`).


.. figure:: images/ASM3h+.png
    :width: 100%
    :align: center

    Network representation of the ASM3h model.

Model Parameters
~~~~~~~~~~~~~~~~

Stoichiometric Parameters
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 60 20

   * - Symbol
     - Description
     - Unit
   * - :math:`i_{N,SI}`
     - Nitrogen content of soluble inert organics
     - :math:`g~N~(g~COD)^{-1}`
   * - :math:`i_{N,SS}`
     - Nitrogen content of readily biodegradable substrate
     - :math:`g~N~(g~COD)^{-1}`
   * - :math:`i_{N,XI}`
     - Nitrogen content of inert particulate organics
     - :math:`g~N~(g~COD)^{-1}`
   * - :math:`i_{N,XS}`
     - Nitrogen content of slowly biodegradable substrate
     - :math:`g~N~(g~COD)^{-1}`
   * - :math:`i_{N,BM}`
     - Nitrogen content of biomass
     - :math:`g~N~(g~COD)^{-1}`
   * - :math:`i_{VSS,XI}`
     - VSS to COD ratio for inert particulates
     - :math:`g~VSS~(g~COD)^{-1}`
   * - :math:`i_{VSS,XS}`
     - VSS to COD ratio for slowly biodegradable substrate
     - :math:`g~VSS~(g~COD)^{-1}`
   * - :math:`i_{VSS,STO}`
     - VSS to COD ratio for storage products
     - :math:`g~VSS~(g~COD)^{-1}`
   * - :math:`i_{VSS,BM}`
     - VSS to COD ratio for biomass
     - :math:`g~VSS~(g~COD)^{-1}`
   * - :math:`i_{TSS/VSS,BM}`
     - TSS to VSS ratio for biomass
     - :math:`g~TSS~(g~VSS)^{-1}`
   * - :math:`f_{SI}`
     - Fraction of SI from hydrolysis
     - dimensionless
   * - :math:`f_{XI}`
     - Fraction of XI from endogenous respiration
     - dimensionless
   * - :math:`f_{ISS,BM,prod}`
     - Fraction of ISS from biomass production
     - dimensionless
   * - :math:`Y_{H,aer}`
     - Aerobic yield of heterotrophs
     - :math:`g~COD~(g~COD)^{-1}`
   * - :math:`Y_{H,anox}`
     - Anoxic yield of heterotrophs
     - :math:`g~COD~(g~COD)^{-1}`
   * - :math:`Y_{STO,aer}`
     - Aerobic yield of storage
     - :math:`g~COD~(g~COD)^{-1}`
   * - :math:`Y_{STO,anox}`
     - Anoxic yield of storage
     - :math:`g~COD~(g~COD)^{-1}`
   * - :math:`Y_A`
     - Yield of autotrophs
     - :math:`g~COD~(g~N)^{-1}`

Kinetic Parameters
------------------

.. list-table::
   :header-rows: 1
   :widths: 20 60 20

   * - Symbol
     - Description
     - Unit
   * - :math:`k_H`
     - Hydrolysis rate constant (20 °C)
     - :math:`d^{-1}`
   * - :math:`K_{STO}`
     - Maximum storage rate (20 °C)
     - :math:`d^{-1}`
   * - :math:`\mu_H`
     - Heterotrophic max. growth rate (20 °C)
     - :math:`d^{-1}`
   * - :math:`b_H`
     - Rate constant for lysis and decay (20 °C)
     - :math:`d^{-1}`
   * - :math:`\mu_{AUT}`
     - Autotrophic max. growth rate (20 °C)
     - :math:`d^{-1}`
   * - :math:`b_{AUT}`
     - Rate constant for decay of autotrophs (20 °C)
     - :math:`d^{-1}`

Anoxic Reduction Factors
------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 60 20

   * - Symbol
     - Description
     - Unit
   * - :math:`\eta_{HNO_3}`
     - Anoxic reduction factor for heterotrophic growth
     - dimensionless
   * - :math:`\eta_{H,end}`
     - Anoxic reduction factor for endogenous respiration of XH
     - dimensionless
   * - :math:`\eta_{N,end}`
     - Anoxic reduction factor for endogenous respiration of XA
     - dimensionless

Saturation/Inhibition Coefficients
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 60 20

   * - Symbol
     - Description
     - Unit
   * - :math:`K_X`
     - Saturation coefficient for particulate COD
     - :math:`g~COD~(g~COD)^{-1}`
   * - :math:`K_{H,O_2}`
     - Saturation coefficient for oxygen (heterotrophs)
     - :math:`g~O_2~m^{-3}`
   * - :math:`K_{H,SS}`
     - Saturation coefficient for readily biodegradable substrates
     - :math:`g~COD~m^{-3}`
   * - :math:`K_{H,NO_3}`
     - Saturation coefficient for nitrate (heterotrophs)
     - :math:`g~N~m^{-3}`
   * - :math:`K_{H,NH_4}`
     - Saturation coefficient for ammonium (nutrient)
     - :math:`g~N~m^{-3}`
   * - :math:`K_{H,ALK}`
     - Saturation coefficient for alkalinity (heterotrophs)
     - :math:`mol~m^{-3}`
   * - :math:`K_{H,STO}`
     - Saturation coefficient for storage products
     - :math:`g~COD~(g~COD)^{-1}`
   * - :math:`K_{N,O_2}`
     - Saturation coefficient for oxygen (autotrophs)
     - :math:`g~O_2~m^{-3}`
   * - :math:`K_{N,NH_4}`
     - Saturation coefficient for ammonium (substrate)
     - :math:`g~N~m^{-3}`
   * - :math:`K_{N,ALK}`
     - Saturation coefficient for alkalinity (autotrophs)
     - :math:`mol~m^{-3}`

Temperature Dependence
----------------------

Kinetic rates are corrected for temperature using Arrhenius-type expressions:

.. math::

   k(T) = k_{20} \cdot \exp\left(-\theta \cdot (20 - T)\right)

with temperature coefficients:

- :math:`\theta = 0.04` for hydrolysis
- :math:`\theta = 0.06952` for heterotrophic processes (:math:`K_{STO}`, :math:`\mu_H`, :math:`b_H`)
- :math:`\theta = 0.105` for autotrophic processes (:math:`\mu_{AUT}`, :math:`b_{AUT}`)

Optional Parameters
-------------------

.. list-table::
   :header-rows: 1
   :widths: 20 60 20

   * - Symbol
     - Description
     - Unit
   * - :math:`f_{SI,ad}`
     - Fraction of adsorbable SI
     - dimensionless
   * - :math:`f_{SS,ad}`
     - Fraction of adsorbable SS
     - dimensionless
   * - :math:`iO_2`
     - Aeration oxygen input rate
     - :math:`g~O_2~m^{-3}~d^{-1}`
   * - :math:`V`
     - Reactor volume (for aeration)
     - :math:`m^3`

Configuration in CADET
~~~~~~~~~~~~~~~~~~~~~~

For configuration details, see :ref:`activated_sludge_model_config`.

Aeration Strategy
-----------------

The ASM3h model includes an aeration reaction (:math:`r_{13}`) that is activated when both ``ASM3_IO2`` and ``ASM3_V`` are set.
The aeration rate is computed as :math:`iO_2 / 1000`. To disable aeration, omit these parameters.
For more flexible aeration handling, model oxygen input explicitly via the Inlet unit operation (see :ref:`inlet_operation`).

Substrate Fractionation
-----------------------

The organic substrates can be fractionated into adsorbable and non-adsorbable components via ``ASM3_FSI_AD`` and ``ASM3_FSS_AD``.
This requires at least 15 components (including ``SI_ad`` and ``SS_ad``) and enables coupling with binding models for 
adsorptive interactions during biological treatment. If fewer than 15 components are defined, adsorbed components are disabled.

When fractionation is enabled, reactions :math:`r_{14}` and :math:`r_{15}` become active for aerobic and anoxic storage of adsorbed substrate SS\ :sub:`ad`.


Implementation Notes
--------------------

- Component indices can be customized via ``ASM3_COMP_IDX`` if components are ordered differently or offset.
- A minimum threshold of :math:`X_H = 0.1` is applied to prevent division by zero in reactions involving :math:`X_S/X_H` and :math:`X_{STO}/X_H` ratios.


References
~~~~~~~~~~

TODO
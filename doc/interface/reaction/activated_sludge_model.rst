.. _activated_sludge_model_config:

Activated Sludge Model (ASM3h)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX(/particle_type_YYY)/phase_reaction_ZZZ/ - TYPE = ACTIVATED_SLUDGE_MODEL3**

For information on model equations, refer to :ref:`activated_sludge_model`.

Notes
-----

- ``phase_reaction_ZZZ`` refers to one of the phase-specific reaction groups listed in :ref:`FFReaction`.
  This model acts on the liquid phase components, so only ``liquid_reaction_ZZZ`` is supported;
  configuring it in a solid phase is rejected.

**Component configuration**

``ASM3_COMP_IDX``

   Optional component indexes. Set this in case the relevant components start at a certain offset or are provided
   in a different order than listed below:

   ================  =========
   **Component**     **Index**
   ================  =========
   SO                0
   SS                1
   SNH               2
   SNO               3
   SN2               4
   SALK              5
   SI                6
   XI                7
   XS                8
   XH                9
   XSTO              10
   XA                11
   XMI               12
   SI_ad (optional)  13
   SS_ad (optional)  14
   ================  =========

   The two adsorbed components are only part of the list when 15 or more components are defined,
   in which case all 15 indices have to be given. Otherwise the field holds exactly 13 indices.

   ================  ===============================  ========================================================
   **Type:** int     **Range:** :math:`\mathbb{N}_0`   **Length:** 13 / 15
   ================  ===============================  ========================================================

**Numerical parameters**

``XH_MIN_THRESHOLD``

   (Optional) Lower-bound parameter :math:`c_{X_H,\mathrm{min}}` for regularizing the heterotrophic biomass concentration. Defaults to :math:`10^{-6}`.

   **Unit:** Same concentration unit as :math:`X_H`.

   ================  =======================  ===============
   **Type:** double  **Range:** finite, > 0   **Length:** 1
   ================  =======================  ===============

``XH_SMOOTH_EPS``

   (Optional) Smoothing parameter :math:`\varepsilon_X` for the lower bound on the heterotrophic biomass concentration. Defaults to :math:`10^{-4}`.

   **Unit:** Same concentration unit as :math:`X_H`.

   ================  =======================  ===============
   **Type:** double  **Range:** finite, > 0   **Length:** 1
   ================  =======================  ===============

   See :ref:`activated_sludge_model` for the lower-bound equation.

**Environment/process parameters**

``ASM3_T``

   Temperature :math:`T`.

   **Unit:** :math:`°C`

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\ge 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_V``

   (Optional) Reactor volume :math:`V`, used by the aeration reaction. Aeration is only active when both
   ``ASM3_V`` and ``ASM3_IO2`` are set; if either is missing, aeration is not considered in the model. The
   volume is not taken from the unit operation, so it is left to the user to ensure consistency between the
   volume set here and the unit operation volume.

   **Unit:** :math:`m^3`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_IO2``

   (Optional) Aeration oxygen input :math:`iO_2`. If not set, aeration is not considered in the model.
   
   **Unit:** :math:`g~O_2~m^{-3}~d^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\ge 0`       **Length:** 1
   ================  =============================  ========================================================

**Stoichiometric parameters**

The stoichiometric parameters below are read once at configuration time and built into the stoichiometric matrix.
They are not registered as model parameters, so they cannot be changed after configuration and no parameter
sensitivities can be computed for them. The kinetic parameters do not have this restriction.

``ASM3_INSI``

   Nitrogen content of soluble inert organics :math:`i_{N,SI}`.
   
   **Unit:** :math:`g~N~(g~COD)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\ge 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_INSS``

   Nitrogen content of readily biodegradable substrate :math:`i_{N,SS}`.
   
   **Unit:** :math:`g~N~(g~COD)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\ge 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_INXI``

   Nitrogen content of inert particulate organics :math:`i_{N,XI}`.
   
   **Unit:** :math:`g~N~(g~COD)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\ge 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_INXS``

   Nitrogen content of slowly biodegradable substrate :math:`i_{N,XS}`.
   
   **Unit:** :math:`g~N~(g~COD)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\ge 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_INBM``

   Nitrogen content of biomass :math:`i_{N,BM}`.
   
   **Unit:** :math:`g~N~(g~COD)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\ge 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_IVSS_BM``

   VSS to COD ratio for biomass :math:`i_{VSS,BM}`.
   
   **Unit:** :math:`g~VSS~(g~COD)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_ITSS_VSS_BM``

   TSS to VSS ratio for biomass :math:`i_{TSS/VSS,BM}`.
   
   **Unit:** :math:`g~TSS~(g~VSS)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_FISS_BM_PROD``

   Fraction of ISS from biomass production :math:`f_{ISS,BM,prod}`.
   
   **Unit:** dimensionless
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`[0, 1]`      **Length:** 1
   ================  =============================  ========================================================

``ASM3_FSI``

   Fraction of SI from hydrolysis :math:`f_{SI}`.
   
   **Unit:** dimensionless
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`[0, 1]`      **Length:** 1
   ================  =============================  ========================================================

``ASM3_FXI``

   Fraction of XI from endogenous respiration :math:`f_{XI}`.
   
   **Unit:** dimensionless
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`[0, 1]`      **Length:** 1
   ================  =============================  ========================================================

**Yield coefficients**

``ASM3_YH_AER``

   Aerobic yield of heterotrophs :math:`Y_{H,aer}`.
   
   **Unit:** :math:`g~COD~(g~COD)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_YH_ANOX``

   Anoxic yield of heterotrophs :math:`Y_{H,anox}`.
   
   **Unit:** :math:`g~COD~(g~COD)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_YSTO_AER``

   Aerobic yield of storage :math:`Y_{STO,aer}`.
   
   **Unit:** :math:`g~COD~(g~COD)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_YSTO_ANOX``

   Anoxic yield of storage :math:`Y_{STO,anox}`.
   
   **Unit:** :math:`g~COD~(g~COD)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_YA``

   Yield of autotrophs :math:`Y_A`.
   
   **Unit:** :math:`g~COD~(g~N)^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

**Adsorption fractionation (optional)**

``ASM3_FSI_AD``

   (Optional) Fraction of adsorbable SI :math:`f_{SI,ad}`. Only used if 15 components are defined.
   
   **Unit:** dimensionless
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`[0, 1]`      **Length:** 1
   ================  =============================  ========================================================

``ASM3_FSS_AD``

   (Optional) Fraction of adsorbable SS :math:`f_{SS,ad}`. Only used if 15 components are defined.
   
   **Unit:** dimensionless
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`[0, 1]`      **Length:** 1
   ================  =============================  ========================================================

**Maximum rates**

``ASM3_KH20``

   Hydrolysis rate constant :math:`k_H` at 20 °C.
   
   **Unit:** :math:`d^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KSTO20``

   Maximum storage rate :math:`K_{STO}` at 20 °C.

   **Unit:** :math:`d^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_MU_H20``

   Heterotrophic max. growth rate :math:`\mu_{H}` at 20 °C.

   **Unit:** :math:`d^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_BH20``

   Rate constant for lysis and decay :math:`b_H` at 20 °C.

   **Unit:** :math:`d^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_MU_AUT20``

   Autotrophic max. growth rate :math:`\mu_{AUT}` at 20 °C.

   **Unit:** :math:`d^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_BAUT20``

   Rate constant :math:`b_{AUT}` for decay of autotrophs at 20 °C.

   **Unit:** :math:`d^{-1}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

**Anoxic reduction factors**

``ASM3_ETA_HNO3``

   Anoxic reduction factor :math:`\eta_{HNO_3}` – heterotrophic growth.

   **Unit:** dimensionless
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_ETAH_END``

   Anoxic reduction factor :math:`\eta_{H_{end}}` – endogenous respiration XH.

   **Unit:** dimensionless
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_ETAN_END``

   Anoxic reduction factor :math:`\eta_{N_{end}}` – endogenous respiration XA.

   **Unit:** dimensionless
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

**Saturation/inhibition coefficients**

``ASM3_KX``

   Saturation/inhibition coefficient :math:`KX` for particulate COD.

   **Unit:** :math:`g~COD~m_{SP}^{-3}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHO2``

   Saturation/inhibition coefficient :math:`KHO_2` for oxygen, heterotrophic growth.

   **Unit:** :math:`g~O_2~m_{SP}^{-3}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHSS``

   Saturation/inhibition coefficient :math:`KHSS` for readily biodegradable substrates.

   **Unit:** :math:`g~COD~m_{SP}^{-3}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHNO3``

   Saturation/inhibition coefficient :math:`KHNO_3` for nitrate.

   **Unit:** :math:`g~N~m_{SP}^{-3}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHNH4``

   Saturation/inhibition coefficient :math:`KHNH_4` for ammonium (nutrient).

   **Unit:** :math:`g~N~m_{SP}^{-3}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHALK``

   Saturation coefficient :math:`KH_{ALK}` for alkalinity (HCO3-).

   **Unit:** :math:`mol~m_{SP}^{-3}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHSTO``

   Saturation coefficient :math:`KH_{STO}` for storage products.

   **Unit:** dimensionless
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KNO2``

   Saturation coefficient :math:`K_{NO_2}` for oxygen, autotrophic growth.

   **Unit:** :math:`g~O_2~m_{SP}^{-3}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KNNH4``

   Saturation coefficient :math:`K_{NNH_4}` for ammonium (substrate), autotrophic growth.

   **Unit:** :math:`g~N~m_{SP}^{-3}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KNALK``

   Saturation coefficient :math:`K_{NALK}` for alkalinity (HCO3-), autotrophic growth.

   **Unit:** :math:`mol~m_{SP}^{-3}`
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

**Example**

.. code-block:: python3

	# Setting up an ASM3 reaction in the liquid phase of unit 000, with example values

	model.root.input.model.unit_000.nreac_liquid = 1
	rxn = model.root.input.model.unit_000.liquid_reaction_000
	rxn.type = 'ACTIVATED_SLUDGE_MODEL3'

	# Stoichiometric parameters
	rxn.asm3_insi = 0.01
	rxn.asm3_inss = 0.03
	rxn.asm3_inxi = 0.04
	rxn.asm3_inxs = 0.03
	rxn.asm3_inbm = 0.07
	rxn.asm3_ivss_bm = 0.704225352
	rxn.asm3_itss_vss_bm = 1.086956522
	rxn.asm3_fiss_bm_prod = 1
	rxn.asm3_fsi = 0
	rxn.asm3_fxi = 0.2
	rxn.asm3_yh_aer = 0.8
	rxn.asm3_yh_anox = 0.65
	rxn.asm3_ysto_aer = 0.8375
	rxn.asm3_ysto_anox = 0.7
	rxn.asm3_ya = 0.24

	# Kinetic parameters
	rxn.asm3_kh20 = 9
	rxn.asm3_kx = 1
	rxn.asm3_ksto20 = 12
	rxn.asm3_mu_h20 = 3
	rxn.asm3_bh20 = 0.33
	rxn.asm3_mu_aut20 = 1.12
	rxn.asm3_baut20 = 0.18

	# Anoxic reduction factors
	rxn.asm3_eta_hno3 = 0.5
	rxn.asm3_etah_end = 0.5
	rxn.asm3_etan_end = 0.5

	# Saturation coefficients
	rxn.asm3_kho2 = 0.2
	rxn.asm3_khss = 10
	rxn.asm3_khno3 = 0.5
	rxn.asm3_khnh4 = 0.01
	rxn.asm3_khalk = 0.1
	rxn.asm3_khsto = 0.1
	rxn.asm3_kno2 = 0.5
	rxn.asm3_knnh4 = 0.7
	rxn.asm3_knalk = 0.5

	# Temperature, and the optional aeration
	rxn.asm3_t = 12
	rxn.asm3_v = 1000.0
	rxn.asm3_io2 = 0.0

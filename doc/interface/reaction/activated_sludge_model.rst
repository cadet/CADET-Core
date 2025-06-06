.. _activated_sludge_model_config:

Activated Sludge Model (ASM3h)
~~~~~~~~~~~~~~~~~~~~~~~~~

**Group /input/model/unit_XXX/reaction - REACTION_MODEL = ACTIVATED_SLUDGE_MODEL3**

For information on model equations, refer to :ref:`activated_sludge_model`.

Environment/process parameters
----

``ASM3_T``

   Temperature :math:`T`.
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_V``

   Reactor volume :math:`V`. Set this to the column volume since the model can't compute it itself.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_IO2``

   Aeration oxygen input :math:`iO_2`.
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\ge 0`       **Length:** 1
   ================  =============================  ========================================================

Maximum rates
----

``ASM3_KH20``

   Hydrolysis rate constant :math:`k_H` at 20 °C.
   
   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KSTO20``

   Maximum storage rate :math:`K_{STO}` at 20 °C.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\gt 0`       **Length:** 1
   ================  =============================  ========================================================

``ASM3_MU_H20``

   Heterotrophic max. growth rate :math:`\mu_{H}` at 20 °C.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_BH20``

   Rate constant for lysis and decay :math:`b_H` at 20 °C.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_MU_AUT20``

   Autotrophic max. growth rate :math:`\mu_{AUT}` at 20 °C.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_BAUT20``

   Rate constant :math:`b_{AUT}` for decay of autotrophs at 20 °C.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

Anoxic reduction factors
----

``ASM3_ETA_HNO3``

   Anoxic reduction factor :math:`\eta_{HNO_3}` – heterotrophic growth.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_ETAH_END``

   Anoxic reduction factor :math:`\eta_{H_{end}}` – endogenous respiration XH.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_ETAN_END``

   Anoxic reduction factor :math:`\eta_{N_{end}}` – endogenous respiration XA.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

Saturation/inhibition coefficients
----

``ASM3_KX``

   Saturation/inhibition coefficient :math:`KX` for particulate COD.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHO2``

   Saturation/inhibition coefficient :math:`KHO_2` for oxygen, heterotrophic growth.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHSS``

   Saturation/inhibition coefficient :math:`KHSS` for readily biodegradable substrates.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHNO3``

   Saturation/inhibition coefficient :math:`KHNO_3` for nitrate.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHNH4``

   Saturation/inhibition coefficient :math:`KHNH_4` for ammonium (nutrient).

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHALK``

   Saturation coefficient :math:`KH_{ALK}` for alkalinity (HCO3-).

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KHSTO``

   Saturation coefficient :math:`KH_{STO}` for storage products.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KNO2``

   Saturation coefficient :math:`K_{NO_2}` for oxygen, autotrophic growth.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KNNH4``

   Saturation coefficient :math:`K_{NNH_4}` for ammonium (substrate), autotrophic growth.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

``ASM3_KNALK``

   Saturation coefficient :math:`K_{NALK}` for alkalinity (HCO3-), autotrophic growth.

   ================  =============================  ========================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** 1
   ================  =============================  ========================================================

Example
-------

.. python::

	# Example of setting up an ASM3 reaction model in a unit operation with bulk reaction

	# Setup ASM3 reaction for unit 000 with example values
	model.root.input.model.unit_000.reaction_model = 'ACTIVATED_SLUDGE_MODEL3'
	model.root.input.model.unit_000.reaction_bulk.asm3_insi = 0.01
	model.root.input.model.unit_000.reaction_bulk.asm3_inss = 0.03
	model.root.input.model.unit_000.reaction_bulk.asm3_inxi = 0.04
	model.root.input.model.unit_000.reaction_bulk.asm3_inxs = 0.03
	model.root.input.model.unit_000.reaction_bulk.asm3_inbm = 0.07
	model.root.input.model.unit_000.reaction_bulk.asm3_ivss_xi = 0.751879699 # not used
	model.root.input.model.unit_000.reaction_bulk.asm3_ivss_xs = 0.555555556 # not used
	model.root.input.model.unit_000.reaction_bulk.asm3_ivss_sto = 0.6 # not used
	model.root.input.model.unit_000.reaction_bulk.asm3_ivss_bm = 0.704225352 
	model.root.input.model.unit_000.reaction_bulk.asm3_itss_vss_bm = 1.086956522


	model.root.input.model.unit_000.reaction_bulk.asm3_fiss_bm_prod = 1
	model.root.input.model.unit_000.reaction_bulk.asm3_fsi = 0
	model.root.input.model.unit_000.reaction_bulk.asm3_yh_aer = 0.8
	model.root.input.model.unit_000.reaction_bulk.asm3_yh_anox = 0.65

	model.root.input.model.unit_000.reaction_bulk.asm3_ysto_aer = 0.8375
	model.root.input.model.unit_000.reaction_bulk.asm3_ysto_anox = 0.7
	model.root.input.model.unit_000.reaction_bulk.asm3_fxi = 0.2
	model.root.input.model.unit_000.reaction_bulk.asm3_ya = 0.24
	model.root.input.model.unit_000.reaction_bulk.asm3_kh20 = 9
	model.root.input.model.unit_000.reaction_bulk.asm3_kx = 1
	model.root.input.model.unit_000.reaction_bulk.asm3_ksto20 = 12
	model.root.input.model.unit_000.reaction_bulk.asm3_mu_h20 = 3
	model.root.input.model.unit_000.reaction_bulk.asm3_bh20 = 0.33
	model.root.input.model.unit_000.reaction_bulk.asm3_eta_hno3 = 0.5
	model.root.input.model.unit_000.reaction_bulk.asm3_khO2 = 0.2
	model.root.input.model.unit_000.reaction_bulk.asm3_khss = 10
	model.root.input.model.unit_000.reaction_bulk.asm3_khno3 = 0.5
	model.root.input.model.unit_000.reaction_bulk.asm3_khnh4 = 0.01
	model.root.input.model.unit_000.reaction_bulk.asm3_khalk = 0.1
	model.root.input.model.unit_000.reaction_bulk.asm3_khsto = 0.1
	model.root.input.model.unit_000.reaction_bulk.asm3_mu_aut20 = 1.12
	model.root.input.model.unit_000.reaction_bulk.asm3_baut20 = 0.18
	model.root.input.model.unit_000.reaction_bulk.asm3_etah_end = 0.5
	model.root.input.model.unit_000.reaction_bulk.asm3_etan_end = 0.5
	model.root.input.model.unit_000.reaction_bulk.asm3_kno2 = 0.5
	model.root.input.model.unit_000.reaction_bulk.asm3_knnh4 = 0.7
	model.root.input.model.unit_000.reaction_bulk.asm3_knalk = 0.5
	model.root.input.model.unit_000.reaction_bulk.asm3_t = 12


	model.root.input.model.unit_000.reaction_bulk.asm3_v = 1000.0
	model.root.input.model.unit_000.reaction_bulk.asm3_io2 = 0.0

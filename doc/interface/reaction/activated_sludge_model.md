(activated-sludge-model-config)=

# Activated Sludge Model (ASM3h)

**Group /input/model/unit_XXX(/particle_type_YYY)/phase_reaction_ZZZ/ - TYPE = ACTIVATED_SLUDGE_MODEL3**

For information on model equations, refer to [Activated Sludge Model (ASM3h)](../../modelling/reaction/activated_sludge_model.md).

## Notes

- `phase_reaction_ZZZ` refers to one of the phase-specific reaction groups listed in [Reaction models](index.md).
  This model acts on the liquid phase components, so only `liquid_reaction_ZZZ` is supported;
  configuring it in a solid phase is rejected.

> **Warning**
>
> When using ASM3h in the liquid phase of a CSTR, the bulk liquid volume must remain constant throughout the simulation.
> To prescribe this volume, set the CSTR parameter `LIQUID_VOLUME`. In every section, the inlet flow rate must equal the sum of the outlet and filter flow rates.

**Component configuration**

`ASM3_COMP_IDX`

Optional component indexes. Set this in case the relevant components start at a certain offset or are provided
in a different order than listed below:

| **Component** | **Index** |
| --- | --- |
| SO | 0 |
| SS | 1 |
| SNH | 2 |
| SNO | 3 |
| SN2 | 4 |
| SALK | 5 |
| SI | 6 |
| XI | 7 |
| XS | 8 |
| XH | 9 |
| XSTO | 10 |
| XA | 11 |
| XMI | 12 |
| SI_ad (optional) | 13 |
| SS_ad (optional) | 14 |

The two adsorbed components are only part of the list when 15 or more components are defined,
in which case all 15 indices have to be given. Otherwise the field holds exactly 13 indices.

| **Type:** int | **Range:** $\mathbb{N}_0$ | **Length:** 13 / 15 |
| --- | --- | --- |

**Numerical parameters**

`XH_MIN_THRESHOLD`

(Optional) Lower-bound parameter $c_{X_H,\mathrm{min}}$ for regularizing the heterotrophic biomass concentration. Defaults to $10^{-6}$.

**Unit:** Same concentration unit as $X_H$.

| **Type:** double | **Range:** finite, > 0 | **Length:** 1 |
| --- | --- | --- |

`XH_SMOOTH_EPS`

(Optional) Smoothing parameter $\varepsilon_X$ for the lower bound on the heterotrophic biomass concentration. Defaults to $10^{-4}$.

**Unit:** Same concentration unit as $X_H$.

| **Type:** double | **Range:** finite, > 0 | **Length:** 1 |
| --- | --- | --- |

See [Activated Sludge Model (ASM3h)](../../modelling/reaction/activated_sludge_model.md) for the lower-bound equation.

**Environment/process parameters**

`ASM3_T`

Temperature $T$.

**Unit:** $°C$

| **Type:** double | **Range:** $\ge 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_V`

(Optional) Reactor volume $V$, used by the aeration reaction. Aeration is only active when both
`ASM3_V` and `ASM3_IO2` are set; if either is missing, aeration is not considered in the model. The
volume is not taken from the unit operation, so it is left to the user to ensure consistency between the
volume set here and the unit operation volume.

**Unit:** $m^3$

| **Type:** double | **Range:** $\gt 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_IO2`

(Optional) Aeration oxygen input $iO_2$. If not set, aeration is not considered in the model.

**Unit:** $g~O_2~m^{-3}~d^{-1}$

| **Type:** double | **Range:** $\ge 0$ | **Length:** 1 |
| --- | --- | --- |

**Stoichiometric parameters**

The stoichiometric parameters below are read once at configuration time and built into the stoichiometric matrix.
They are not registered as model parameters, so they cannot be changed after configuration and no parameter
sensitivities can be computed for them. The kinetic parameters do not have this restriction.

`ASM3_INSI`

Nitrogen content of soluble inert organics $i_{N,SI}$.

**Unit:** $g~N~(g~COD)^{-1}$

| **Type:** double | **Range:** $\ge 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_INSS`

Nitrogen content of readily biodegradable substrate $i_{N,SS}$.

**Unit:** $g~N~(g~COD)^{-1}$

| **Type:** double | **Range:** $\ge 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_INXI`

Nitrogen content of inert particulate organics $i_{N,XI}$.

**Unit:** $g~N~(g~COD)^{-1}$

| **Type:** double | **Range:** $\ge 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_INXS`

Nitrogen content of slowly biodegradable substrate $i_{N,XS}$.

**Unit:** $g~N~(g~COD)^{-1}$

| **Type:** double | **Range:** $\ge 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_INBM`

Nitrogen content of biomass $i_{N,BM}$.

**Unit:** $g~N~(g~COD)^{-1}$

| **Type:** double | **Range:** $\ge 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_IVSS_BM`

VSS to COD ratio for biomass $i_{VSS,BM}$.

**Unit:** $g~VSS~(g~COD)^{-1}$

| **Type:** double | **Range:** $\gt 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_ITSS_VSS_BM`

TSS to VSS ratio for biomass $i_{TSS/VSS,BM}$.

**Unit:** $g~TSS~(g~VSS)^{-1}$

| **Type:** double | **Range:** $\gt 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_FISS_BM_PROD`

Fraction of ISS from biomass production $f_{ISS,BM,prod}$.

**Unit:** dimensionless

| **Type:** double | **Range:** $[0, 1]$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_FSI`

Fraction of SI from hydrolysis $f_{SI}$.

**Unit:** dimensionless

| **Type:** double | **Range:** $[0, 1]$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_FXI`

Fraction of XI from endogenous respiration $f_{XI}$.

**Unit:** dimensionless

| **Type:** double | **Range:** $[0, 1]$ | **Length:** 1 |
| --- | --- | --- |

**Yield coefficients**

`ASM3_YH_AER`

Aerobic yield of heterotrophs $Y_{H,aer}$.

**Unit:** $g~COD~(g~COD)^{-1}$

| **Type:** double | **Range:** $\gt 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_YH_ANOX`

Anoxic yield of heterotrophs $Y_{H,anox}$.

**Unit:** $g~COD~(g~COD)^{-1}$

| **Type:** double | **Range:** $\gt 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_YSTO_AER`

Aerobic yield of storage $Y_{STO,aer}$.

**Unit:** $g~COD~(g~COD)^{-1}$

| **Type:** double | **Range:** $\gt 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_YSTO_ANOX`

Anoxic yield of storage $Y_{STO,anox}$.

**Unit:** $g~COD~(g~COD)^{-1}$

| **Type:** double | **Range:** $\gt 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_YA`

Yield of autotrophs $Y_A$.

**Unit:** $g~COD~(g~N)^{-1}$

| **Type:** double | **Range:** $\gt 0$ | **Length:** 1 |
| --- | --- | --- |

**Adsorption fractionation (optional)**

`ASM3_FSI_AD`

(Optional) Fraction of adsorbable SI $f_{SI,ad}$. Only used if 15 components are defined.

**Unit:** dimensionless

| **Type:** double | **Range:** $[0, 1]$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_FSS_AD`

(Optional) Fraction of adsorbable SS $f_{SS,ad}$. Only used if 15 components are defined.

**Unit:** dimensionless

| **Type:** double | **Range:** $[0, 1]$ | **Length:** 1 |
| --- | --- | --- |

**Maximum rates**

`ASM3_KH20`

Hydrolysis rate constant $k_H$ at 20 °C.

**Unit:** $d^{-1}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_KSTO20`

Maximum storage rate $K_{STO}$ at 20 °C.

**Unit:** $d^{-1}$

| **Type:** double | **Range:** $\gt 0$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_MU_H20`

Heterotrophic max. growth rate $\mu_{H}$ at 20 °C.

**Unit:** $d^{-1}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_BH20`

Rate constant for lysis and decay $b_H$ at 20 °C.

**Unit:** $d^{-1}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_MU_AUT20`

Autotrophic max. growth rate $\mu_{AUT}$ at 20 °C.

**Unit:** $d^{-1}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_BAUT20`

Rate constant $b_{AUT}$ for decay of autotrophs at 20 °C.

**Unit:** $d^{-1}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

**Anoxic reduction factors**

`ASM3_ETA_HNO3`

Anoxic reduction factor $\eta_{HNO_3}$ – heterotrophic growth.

**Unit:** dimensionless

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_ETAH_END`

Anoxic reduction factor $\eta_{H_{end}}$ – endogenous respiration XH.

**Unit:** dimensionless

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_ETAN_END`

Anoxic reduction factor $\eta_{N_{end}}$ – endogenous respiration XA.

**Unit:** dimensionless

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

**Saturation/inhibition coefficients**

`ASM3_KX`

Saturation/inhibition coefficient $KX$ for particulate COD.

**Unit:** $g~COD~m_{SP}^{-3}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_KHO2`

Saturation/inhibition coefficient $KHO_2$ for oxygen, heterotrophic growth.

**Unit:** $g~O_2~m_{SP}^{-3}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_KHSS`

Saturation/inhibition coefficient $KHSS$ for readily biodegradable substrates.

**Unit:** $g~COD~m_{SP}^{-3}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_KHNO3`

Saturation/inhibition coefficient $KHNO_3$ for nitrate.

**Unit:** $g~N~m_{SP}^{-3}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_KHNH4`

Saturation/inhibition coefficient $KHNH_4$ for ammonium (nutrient).

**Unit:** $g~N~m_{SP}^{-3}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_KHALK`

Saturation coefficient $KH_{ALK}$ for alkalinity (HCO3-).

**Unit:** $mol~m_{SP}^{-3}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_KHSTO`

Saturation coefficient $KH_{STO}$ for storage products.

**Unit:** dimensionless

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_KNO2`

Saturation coefficient $K_{NO_2}$ for oxygen, autotrophic growth.

**Unit:** $g~O_2~m_{SP}^{-3}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_KNNH4`

Saturation coefficient $K_{NNH_4}$ for ammonium (substrate), autotrophic growth.

**Unit:** $g~N~m_{SP}^{-3}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

`ASM3_KNALK`

Saturation coefficient $K_{NALK}$ for alkalinity (HCO3-), autotrophic growth.

**Unit:** $mol~m_{SP}^{-3}$

| **Type:** double | **Range:** $\mathbb{R}$ | **Length:** 1 |
| --- | --- | --- |

**Example**

```python3
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
```

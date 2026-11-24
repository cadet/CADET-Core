(colloidal-particle-adsorption-config)=

# Colloidal Particle Adsorption

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = COLLOIDAL_PARTICLE_ADSORPTION**

For information on model equations, refer to {ref}`colloidal-particle-adsorption-model`.

`IS_KINETIC`

: Selects kinetic or quasi-stationary adsorption mode: 1 = kinetic, 0 =
  quasi-stationary. If a single value is given, the mode is set for all
  bound states. Otherwise, the adsorption mode is set for each bound
  state separately.

```{list-table}

* - **Type:** int
  - **Range:** {0,1}
  - **Length:** 1/NTOTALBND
```

`CPA_TEMPERATURE`

: Absolute temperature $T$

**Unit:** $\mathrm{K}$

```{list-table}

* - **Type:** double
  - **Range:** $> 0$
  - **Length:** 1
```

`CPA_IONIC_STRENGTH`

: Ionic strength $I_m$ of the mobile phase. Required but ignored if
  `CPA_IONIC_VALENCE` is set (ionic strength is then computed
  from the pore-phase concentrations).

**Unit:** $\mathrm{mol \, m^{-3}}$

```{list-table}

* - **Type:** double
  - **Range:** $> 0$
  - **Length:** 1
```

`CPA_PERMITTIVITY`

: Relative permittivity $\varepsilon$ of the solvent

```{list-table}

* - **Type:** double
  - **Range:** $> 0$
  - **Length:** 1
```

`CPA_LIGAND_DENSITY`

: Ligand surface density $\Gamma_L$

**Unit:** $\mathrm{mol \, m^{-2}}$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** 1
```

`CPA_LIGAND_CHARGE_FULL`

: Charge of the fully protonated ligand $\zeta_L$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** 1
```

`CPA_LIGAND_PK`

: Dissociation constant $\mathrm{p}K_L$ of the ligand

```{list-table}

* - **Type:** double
  - **Range:** $> 0$
  - **Length:** 1
```

`CPA_PH_REF`

: Reference pH value $\mathrm{pH}_{\mathrm{ref}}$ for the
  protein charge polynomial

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** 1
```

`CPA_SPECIFIC_SURFACE_AREA`

: Specific adsorber surface area per skeleton volume $A_{s,i}$.
  Entries of non-binding components are not used and may be set to zero.

**Unit:** $\mathrm{m^{-1}}$

```{list-table}

* - **Type:** double
  - **Range:** $> 0$
  - **Length:** NCOMP
```

`CPA_RADIUS`

: Hydrodynamic radius $a_i$ of each protein component.
  Entries of non-binding components are not used and may be set to zero.

**Unit:** $\mathrm{m}$

```{list-table}

* - **Type:** double
  - **Range:** $> 0$
  - **Length:** NCOMP
```

`CPA_EFFECTIVE_CHARGE_COEF`

: Matrix of coefficients $z_{i,k}$ for the pH-dependent protein
  net charge. Matrix rows correspond to increasing polynomial powers
  $k=0,\ldots,P$, and columns correspond to components. The first row
  therefore contains the constant charges $z_{i,0}$ at
  $\mathrm{pH}=\mathrm{pH}_{\mathrm{ref}}$. The net charge is

  $$
  Z_i(\mathrm{pH}) = z_{i,0} + \sum_{k=1}^{P} z_{i,k}
  \left(\mathrm{pH}-\mathrm{pH}_{\mathrm{ref}}\right)^k.
  $$

  The matrix is supplied in polynomial-order-row-major ordering.
  Every row must contain one value for every component, including non-binding
  components.
  Coefficients that do not apply should be set to zero.

```{list-table}

* - **Type:** double
  - **Range:** $\mathbb{R}$
  - **Length:** NCOMP * (P + 1)
```

`CPA_LAT_CHARGE`

: Lateral charge $Z_{\mathrm{lat},i}$ used for computing
  the pairwise Yukawa coefficient $\beta_{ij}$ in the lateral
  protein–protein interaction

```{list-table}

* - **Type:** double
  - **Range:** $\mathbb{R}$
  - **Length:** NCOMP
```

`CPA_DELTA_REF`

: Reference interaction layer thickness $\delta_{i,\mathrm{ref}}$
  Entries of non-binding components are not used and may be set to zero.

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NCOMP
```

`CPA_DELTA_LIN`

: Linear coefficient $\delta_{i,\mathrm{lin}}$ for the
  logarithmic interaction layer thickness parameterisation

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NCOMP
```

`CPA_KKIN`

: Kinetic prefactor $k^*_{\mathrm{kin},i}$ used to calculate the
  kinetic rate constant $k_{\mathrm{kin},i}$.

**Unit:** $\mathrm{s^{-1}}$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NCOMP
```

`CPA_PROTON_IDX`

: 0-based index of the proton/pH component in the mobile phase
  (optional, defaults to 0). The corresponding component must be
  non-binding (`NBOUND = 0`).

```{list-table}

* - **Type:** int
  - **Range:** $\ge 0$
  - **Length:** 1
```

`CPA_IONIC_VALENCE`

: Integer valence charge $z_i$ for each component (optional).
  If provided, the ionic strength is computed from the pore-phase
  concentrations, and the Davies activity correction is applied to the proton
  component when computing pH. The vector must contain one entry per
  component (`NCOMP` values). Components that do not contribute to the ionic
  strength, should be assigned a charge of zero. Note that a charge of zero also
  disables the Davies activity correction for that component, so assigning
  $z = 0$ to the proton component yields $\gamma_{\mathrm{H}^+} = 1$.

```{list-table}

* - **Type:** int
  - **Range:** $\mathbb{Z}^+$
  - **Length:** NCOMP
```

`CPA_MAXITER`

: Maximum number of Newton iterations for solving the adsorber surface
  potential $\psi_{0,A}$ (optional, defaults to 100). If the solver
  does not converge, a warning is emitted and the last iterate is used.

```{list-table}

* - **Type:** int
  - **Range:** $\ge 1$
  - **Length:** 1
```

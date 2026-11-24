(multi-component-colloidal-config)=

# Multi Component Colloidal

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = MULTI_COMPONENT_COLLOIDAL**

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

`COL_PHI`

: Phase ratio $\Phi$

**Unit:** $m^{2} m_{s}^{-3}$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** 1
```

`COL_KAPPA_EXP`

: Screening term exponent factor $\kappa_{ef}$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** 1
```

`COL_KAPPA_FACT`

: Screening term factor $\kappa_{f}$

**Unit:** $m \cdot mM^{-1}$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** 1
```

`COL_KAPPA_CONST`

: Screening term constant $\kappa_{c}$

**Unit:** $m$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** 1
```

`COL_CORDNUM`

: Coordination number $n$

```{list-table}

* - **Type:** int
  - **Range:** $\ge 0$
  - **Length:** 1
```

`COL_LOGKEQ_PH_EXP`

: Protein-resin interaction $K_{e,i}$ equilibrium: Constant exponent $k_{e,i}$ for pH
  If pH is not considered, this value will be not be used but must still be specified, i.e. can be any number.

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_LOGKEQ_SALT_POWEXP`

: Protein-resin interaction $K_{e,i}$ equilibrium: Constant pre-factor $k_{a,i}$ for salt power

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_LOGKEQ_SALT_POWFACT`

: Protein-resin interaction $K_{e,i}$ equilibrium: Constant exponent $k_{b,i}$ for salt power

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_LOGKEQ_SALT_EXPFACT`

: Protein-resin interaction $K_{e,i}$ equilibrium: Constant pre-factor $k_{c,i}$ for e-function with salt power

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_LOGKEQ_SALT_EXPARGMULT`

: Protein-resin interaction $K_{e,i}$ equilibrium: Constant power factor $k_{d,i}$ for salt in e-function

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_BPP_PH_EXP`

: Protein-protein interaction $b_{pp,i}$: Constant power term $b_{e,i}$ for pH.
  If pH is not considered, this value will be not be used but must still be specified, i.e. can be any number.

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_BPP_SALT_POWFACT`

: Protein-protein interaction $b_{pp,i}$: Constant power pre-factor $b_{a,i}$ for salt

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_BPP_SALT_POWEXP`

: Protein-protein interaction $b_{pp,i}$: Constant power $b_{b,i}$ for salt

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_BPP_SALT_EXPFACT`

: Protein-protein interaction $b_{pp,i}$: Constant pre-factor $b_{c,i}$ e-function with salt power

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_BPP_SALT_EXPARGMULT`

: Protein-protein interaction $b_{pp,i}$: Constant power factor $b_{d,i}$ for salt in e-function

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_PROTEIN_RADIUS`

: Protein radius $r_i$

**Unit:** $m$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_KKIN`

: Adsorption rate constants $K_\text{kin}$ in state-major ordering

**Unit:** $s^{-1}$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NTOTALBND
```

`COL_LINEAR_THRESHOLD`

: Threshold concentration $c_\epsilon$ for switching to linear approximation

```{list-table}

* - **Type:** double
  - **Range:** $> 0$
  - **Length:** 1
```

`COL_USE_PH`

: Selects if pH is included in the model or not: 1 = yes, 0 = no.

```{list-table}

* - **Type:** int
  - **Range:** {0,1}
  - **Length:** 1
```

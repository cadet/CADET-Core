(multi-component-bi-langmuir-ldf-config)=

# Multi Component Bi-Langmuir LDF

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption – ADSORPTION_MODEL = MULTI_COMPONENT_BILANGMUIR_LDF**

For information on model equations, refer to {ref}`multi-component-bi-langmuir-ldf-model`.

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

`MCBLLDF_KEQ`

: Equilibrium loading constants in state-major ordering (see {ref}`ordering-multi-dimensional-data`)

**Unit:** $m_{MP}^3~mol^{-1}$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NCOMP $\cdot$ NSTATES
```

`MCBLLDF_KKIN`

: Linear driving force coefficients in state-major ordering

**Unit:** $s^{-1}$

```{list-table}

* - **Type:** double
  - **Range:** $\ge 0$
  - **Length:** NCOMP $\cdot$ NSTATES
```

`MCBLLDF_QMAX`

: Maximum adsorption capacities in state-major ordering

**Unit:** $mol~m_{SP}^{-3}$

```{list-table}

* - **Type:** double
  - **Range:** $\gt 0$
  - **Length:** NCOMP $\cdot$ NSTATES
```

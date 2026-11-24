(michaelis-menten-kinetics-config)=

# Michaelis Menten kinetics

**Group /input/model/unit_XXX(/particle_type_YYY)/phase_reaction_ZZZ/ - TYPE = MICHAELIS_MENTEN**

For information on model equations, refer to {ref}`michaelis-menten-kinetics-model`.

## Notes

- `phase_reaction_ZZZ` refers to one of the phase-specific reaction groups listed in {ref}`ffreaction`, e.g., `liquid_reaction_ZZZ`, `solid_reaction_ZZZ`.
- Some dimensions below depend on the hosting phase of this model instance:
  : - Bulk phase or particle liquid phase: `NVAR = NCOMP`
    - Particle solid phase: `NVAR = NTOTALBOUND` (total number of bound states across all components)

`MM_STOICHIOMETRY`

: Stoichiometric matrix $S$.
  This matrix defines the quantitative relationships between reactants and products for each reaction in the system.
  Each entry $S_{i,j}$ specifies the stoichiometric coefficient for component $i$ in reaction $j$.
  Negative values indicate consumption (substrate), while positive values indicate production (products).
  Input as reaction index major.

  **Unit:** None

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NREACT} \cdot \texttt{NVAR}$
  ```

`MM_VMAX`

: Maximum reaction rate $v_{\mathrm{max},j}$ at substrate saturation for reaction $j$.
  This parameter defines the upper limit of the reaction rate when the substrate concentration is sufficiently high such that the enzyme is saturated.

  **Unit:** $mol~m^{-3}~s^{-1}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NREACT}$
  ```

`MM_KMM`

: Michaelis constant $K_{\mathrm{M}_{i,j}}$ for reaction $j$ and substrate $i$.
  This constant represents the substrate concentration at which the reaction rate is half of its maximum value.
  Input as reaction index major.

  **Unit:** $mol~m^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NREACT} \cdot \texttt{NVAR}$
  ```

`MM_KI_C`

: (Optional) Inhibition constant for competitive inhibition $K^{c}_{I_{k}}$.
  The index $k$ corresponds to the inhibitors acting on substrate $c_{i,j}$ in reaction $j$, i.e. $k = (j,i,k)$, where $k$ is the index of the inhibitor.
  If $K^{c}_{I_{k}} > 0$, the component inhibits the reaction.
  Input as reaction index major.

  **Unit:** $mol~m^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NREACT} \cdot \texttt{NVAR} \cdot \texttt{NVAR}$
  ```

`MM_KI_UC`

: (Optional) Inhibition constant for uncompetitive inhibition $K^{uc}_{I_{k}}$.
  The index $k$ corresponds to the inhibitors acting on substrate $c_{i,j}$ in reaction $j$, i.e. $k = (j,i,k)$, where $k$ is the index of the inhibitor.
  Input as reaction index major.

  **Unit:** $mol~m^{-3}$

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NREACT} \cdot \texttt{NVAR} \cdot \texttt{NVAR}$
  ```

`MM_PRE_K`

: (Optional) Factor ($\neq 0$) before component for multiplying with Michaelis-Menten kinetics.
  The index $k$ corresponds to the component acting on reaction $j$. Note that components acting as substrates cannot be defined as prefactor components.
  Input as reaction index major.

  **Unit:** None

  ```{list-table}

  * - **Type:** double
    - **Range:** $\mathbb{R}$
    - **Length:** $\texttt{NREACT} \cdot \texttt{NVAR}$
  ```

`MM_EXP`

: (Optional) Hill exponent $n_{i,j}$ for substrate $i$ in reaction $j$.
  Only $\neq 1$ for substrate without competitive or uncompetitive inhibition.
  Use $1$ for substrates with competitive or uncompetitive inhibition.
  Input as reaction index major.

  **Unit:** None

  ```{list-table}

  * - **Type:** double
    - **Range:** $> 0$
    - **Length:** $\texttt{NREACT} \cdot \texttt{NVAR}$
  ```

## CADET Python Interface Example

This example shows the configuration of one Michaelis-Menten reaction system in bulk liquid phase.
The system has two components A and B, where A is the substrate and B is the product.
In addition to that the model includes:

- A Michaelis constant `KM_a`,
- competitive inhibition constant of `KI_b_a` for B inhibiting A,
- and a maximum rate of `vmax`.

```python3
# Configure the reaction system
input.model.unit_000.NREAC_LIQUID = 1
input.model.unit_000.liquid_reaction_000.TYPE = 'MICHAELIS_MENTEN'

# Km values 2D array [reaction][components]
input.model.unit_000.liquid_reaction_000.MM_KMM = [
    [KM_a, 0.0]  # A is substrate
]

# Competitive inhibition constants - 3D array [reaction][components][components]
input.model.unit_000.liquid_reaction_000.MM_KI_C = [
    [
        [0.0, KI_b_a],  # Inhibition constant for A (Product inhibition B inhibits A)
        [0.0, 0.0],     # Inhibition constant for B (not active)
    ]
]

# Uncompetitive inhibition constants - 3D array [reaction][components][components]
input.model.unit_000.liquid_reaction_000.MM_KI_UC = [
    [
        [0.0, 0.0],  # Inhibition constant for A (not active)
        [0.0, 0.0],  # Inhibition constant for B (not active)
    ]
]

# Vmax values 1D array [reaction]
input.model.unit_000.liquid_reaction_000.MM_VMAX = [vmax]

# Stoichiometry matrix 2D array [components][reaction]
input.model.unit_000.liquid_reaction_000.MM_STOICHIOMETRY = [
    [-1],
    [1]  # A -> B
]
```

(neural-network-config)=

# Neural Network Binding

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption - ADSORPTION_MODEL = NEURAL_NETWORK**

For information on model equations, refer to {ref}`neural-network`.

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

`NN_KKIN`

: Linear-driving-force coefficients in component-major ordering.
  Controls the rate at which the solid-phase loading approaches the
  neural network-predicted equilibrium.

**Unit:** $s^{-1}$

```{list-table}

* - **Type:** double
  - **Range:** $\geq 0$
  - **Length:** NTOTALBND
```

`NLAYERS`

: Number of hidden layers $N$ in the neural network architecture.

```{list-table}

* - **Type:** int
  - **Range:** $\geq 1$
  - **Length:** 1
```

`NNODES`

: Number of nodes per hidden layer.
  All hidden layers share the same width.

```{list-table}

* - **Type:** int
  - **Range:** $\geq 1$
  - **Length:** 1
```

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_XXX**

`NORM_FACTOR`

: Normalization factors applied element-wise to the pore-phase concentration
  vector before it is fed into the neural network, i.e. the network is
  evaluated on $x = c^p \odot \alpha_{\mathrm{norm}}$. Each factor is
  therefore the reciprocal of a characteristic concentration scale, for example
  $1 / \max\left(c^p\right)$, chosen such that the network inputs are of
  order one. The network must have been trained on exactly this scaling.

  A single value is applied to all input dimensions; otherwise one value per
  component has to be given. Each bound state has its own factors because each
  network is trained independently.

```{list-table}

* - **Type:** double
  - **Range:** $\neq 0$
  - **Length:** 1/NCOMP
```

`POROSITY_FACTOR`

: Optional scaling factor applied to the neural network prediction. This can
  be used to account for porosity differences or unit conversions between
  training data and simulation conditions. Defaults to 1.0 if not given.

```{list-table}

* - **Type:** double
  - **Range:** $> 0$
  - **Length:** 1
```

Neural network weights and biases are organized hierarchically by bound state and layer.
For each bound state, there are `NLAYERS + 1` layer groups:
`layer_0` through `layer_{NLAYERS-1}` are the hidden layers and
`layer_{NLAYERS}` is the output layer. All weight matrices must be stored in
column-major (Fortran) order.

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_XXX/layer_0**

First hidden layer parameters. Always present.

`KERNEL`

: Weight matrix $W_1$ connecting the input to the first hidden layer.
  Shape: (NNODES x NCOMP). Stored in column-major order.

```{list-table}

* - **Type:** double
  - **Range:** $\mathbb{R}$
  - **Length:** NNODES * NCOMP
```

`BIAS`

: Bias vector $b_1$ for the first hidden layer.

```{list-table}

* - **Type:** double
  - **Range:** $\mathbb{R}$
  - **Length:** NNODES
```

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_XXX/layer_l**
**(for l = 1, ..., NLAYERS-1)**

Intermediate hidden layer parameters. Present only when `NLAYERS >= 2`.
Layer `layer_l` connects hidden layer $l$ to hidden layer $l+1$.

`KERNEL`

: Weight matrix $W_{l+1}$ of shape (NNODES x NNODES). Stored in column-major order.

```{list-table}

* - **Type:** double
  - **Range:** $\mathbb{R}$
  - **Length:** NNODES\*NNODES
```

`BIAS`

: Bias vector $b_{l+1}$ for hidden layer $l+1$.

```{list-table}

* - **Type:** double
  - **Range:** $\mathbb{R}$
  - **Length:** NNODES
```

**Group /input/model/unit_XXX/particle_type_ZZZ/adsorption/bound_state_XXX/layer_NLAYERS**

Output layer parameters. Always present. Index equals `NLAYERS`.

`KERNEL`

: Weight matrix $W_{N+1}$ connecting the last hidden layer to the scalar output.
  Shape: (1 x NNODES). Stored in column-major order.

```{list-table}

* - **Type:** double
  - **Range:** $\mathbb{R}$
  - **Length:** NNODES
```

`BIAS`

: Bias scalar $b_{N+1}$ for the output layer.

```{list-table}

* - **Type:** double
  - **Range:** $\mathbb{R}$
  - **Length:** 1
```

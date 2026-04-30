(affinity-energy-distribution)=

# Affinity Energy Distribuition

Unlike the multi component Langmuir model and the multi component Bi-Langmuir model that presume a fixed number of discrete binding sites, the affinity energy distribution (AED) model considers a continuous distribution of binding site types on a heterogeneous surface on the energy space $\ln(K_{eq}) \in \mathbb{R}^{>0}$ (see e.g. {cite}`Guiochon2006`).
Assuming a single component system, the AED isotherm reads {cite}`Guiochon2006`:

$$
\begin{aligned}
    q(c) = \int_0^{\infty} F(\ln(K_{eq})) \theta(K_{eq}, c) \, \mathrm{d} \ln(K_{eq}),
\end{aligned}
$$

where $F\colon \mathbb{R}^{>0} \mapsto \mathbb{R}$ is the affinity energy distribution and $\theta\colon  \mathbb{R}^{>0} \times \mathbb{R}^{>0} \mapsto \mathbb{R}$ is the kernel function.
In general, the kernel can be chosen as any isotherm model that governs the adsorption and desorption behavior of each individual site, however, the most common and also well-tested kernel is the Langmuir kernel:

$$
\begin{aligned}
    \theta (K_{eq}, c) =  \frac{K_{eq}c}{1 + K_{eq}c},
\end{aligned}
$$

meaning that each binding site is governed by a Langmuir isotherm model.

When solving the AED model in a unit operation, the integral is disretized on an equidistant logarithmic grid between $\ln(K_{min})$ and $\ln(K_{max})$ to give:

> $$
> \begin{aligned}
>     q(c) = \sum_{i=1}^{N_k} F(\ln(K_{eq, i})) \theta (K_{eq,i}, c) \Delta \ln(K_{eq}),
> \end{aligned}
> $$

where $K_{min}$ is the mimimum binding equilibrium constant, $K_{max}$ is the maximum binding equilibrium constant, $N_k$ is the total number of of nodes (sites) and $\Delta \ln(K_{eq})$ is the distance between two nodes.
Hence

> $$
> \begin{aligned}
>     F(\ln(K_{eq, i})) \Delta \ln(K_{eq}) = q_{\text{max}, i},
> \end{aligned}
> $$

where $q_{max, i}$ is the binding capacity for the $i^{\text{th}}$ binding site. Inserting the Langmuir kernel, the above equation is equivilent to:

> $$
> \begin{aligned}
>     q(c) = \sum_{i=1}^{N_k} \frac{q_{\text{max}, i} K_{eq, i}c}{1 + K_{eq, i}c},
> \end{aligned}
> $$

indicating that the overall bound concentration $q$ is the summation of the bound concentration for each binding site. The kinetic isotherm for each individual binding site is given by:

$$
\begin{aligned}
    \frac{\mathrm{d} q_i}{\mathrm{d} t} &=  k_{a,i}\: c \: q_{\text{max}, i} \left( 1 - \frac{q_i}{ q_{\text{max}, i} } \right) - k_{d, i} q_{i} & i = 0, \dots, N_k - 1.
\end{aligned}
$$

The above discretized AED isotherm model is equivilent to the [](#multi-component-bi-langmuir-model) with a single component and $N_k$ binding sites. Therefore, the AED isotherm model is not registered in CADET as a standalone isotherm model and the [](#multi-component-bi-langmuir-model) is used instead.
Please refer to [](#multi-component-bi-langmuir-config) for information on the configuration of the CADET file format.

We emphasize that while the Multi-component Bi-langmuir isotherm implementation in CADET allows multiple components, the AED isotherm is originally derived and used for a single-component system.

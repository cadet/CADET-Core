(multi-component-langmuir-ldf-model)=

# Multi Component Langmuir LDF

This a linear driving force model variant of the [](#multi-component-langmuir-model) model.
It is based on the equilibrium concentration $q^*$ for a given liquid phase concentration $c$ (see also [](#ldf-model)).

$$
\begin{aligned}
    q_i^*=\frac{q_{m,i} k_{eq,i} c_i}{1 + \sum_{j=1}^{n_{comp}}{k_{eq,j} c_j}} && i = 0, \dots, N_{\text{comp}} - 1.
\end{aligned}
$$

For more information on model parameters required to define in CADET file format, see [](#multi-component-langmuir-ldf-config).

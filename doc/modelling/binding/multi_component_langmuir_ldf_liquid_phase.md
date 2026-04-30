(multi-component-langmuir-ldf-liquid-phase-model)=

# Multi Component Langmuir LDF Liquid Phase

This a linear driving force model variant of the [](#multi-component-langmuir-model) model.
It is based on the equilibrium concentration $c^*$ for a given solid phase concentration $q$ (see also [](#ldf-model)).

$$
\begin{aligned}
    c_i^*=\frac{q_{i}}{k_{eq,i} q_{m,i} \left(1 - \sum_{j=1}^{N_{\text{comp}}} \frac{q_j}{q_{m,j}}\right) } && i = 0, \dots, N_{\text{comp}} - 1.
\end{aligned}
$$

For more information on model parameters required to define in CADET file format, see [](#multi-component-langmuir-ldf-liquid-phase-config).

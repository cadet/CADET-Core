(multi-component-bi-langmuir-ldf-model)=

# Multi Component Bi-Langmuir LDF

This a linear driving force model variant of the [](#multi-component-bi-langmuir-model) model.
It is based on the equilibrium concentration $q^*$ for a given liquid phase concentration $c$ (see also [](#ldf-model)).

$$
\begin{aligned}
q_{i,j}^*=\frac{q_{m,i,j} k_{eq,i,j} c_i}{1 + \sum_{k=1}^{N_{comp}}{k_{eq,k,j} c_k}} \quad i = 0, \dots, N_{\text{comp}} - 1, \: j = 0, \dots, M - 1.
\end{aligned}
$$

For more information on model parameters required to define in CADET file format, see [](#multi-component-bi-langmuir-ldf-config).

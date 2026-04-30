(multi-component-ldf-freundlich-model)=

# Multi Component Linear Driving Force Freundlich

A multi-component extension to the classical Freundlich adsorption model.
A linear driving force approach is applied to obtain a kinetic form.

$$
\begin{aligned}
    \frac{\mathrm{d} q_i}{\mathrm{d} t} &= k_{\text{ldf},i}\: \left(q_i^* - q_i \right) & i = 0, \dots, N_{\text{comp}} - 1 \\
    q_i^* &= k_{F,i} c_{p,i} \left( \sum_j a_{ij} c_{p,j} + \tau \right)^{1 / n_i - 1}.
\end{aligned}
$$

Here, $\tau > 0$ is a small constant that ensures numerical stability.

In a rapid-equilibrium setting with a diagonal matrix (i.e., $a_{ii} = 1$ and $a_{ij} = 0$ for $j \neq i$), the traditional Freundlich isotherm is recovered.

For more information on model parameters required to define in CADET file format, see [](#multi-component-ldf-freundlich-config).

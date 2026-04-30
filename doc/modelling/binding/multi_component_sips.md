(multi-component-sips-model)=

# Multi Component Sips

The Sips binding model is a combination of the [Freundlich](#freundlich-ldf-model) and the [Langmuir adsorption model](#multi-component-langmuir-model).

$$
\begin{aligned}
    \frac{\mathrm{d} q_i}{\mathrm{d} t} = k_{a,i}\: \left( \frac{c_{p,i}}{ c_{\text{ref}} }\right)^{1 / n_i}\: q_{\text{max},i} \left( 1 - \sum_{j=0}^{N_{\text{comp}} - 1} \frac{q_j}{q_{\text{max},j}} \right) - k_{d,i} \left( \frac{q_i}{q_{\text{ref}}} \right) && i = 0, \dots, N_{\text{comp}} - 1.
\end{aligned}
$$

Here, $c_{\text{ref}}$ is a [reference concentration](#reference-concentrations), $n_i$ is the Freundlich exponent, $k_{a,i}, k_{d,i}$ are the adsorption and desorption rates, and $q_{\text{max},j}$ is the adsorption capacity.

As for the [Freundlich](#freundlich-ldf-model) isotherm, the first order Jacobian $\left(\frac{dq^*}{dc_p}\right)$ tends to infinity as $c_{p} \rightarrow 0$ for $n>1$.
Additionally, the isotherm is undefined for $c_{p} < 0$ if $\frac{1}{n_i}$ can be expressed as $\frac{p}{q}$ with $p,q \in \mathbb{N}$ where $q$ is an even number.
Negative concentrations can arise during simulations due to numerical fluctuations.
To address these issues an approximation of the isotherm is considered below a threshold concentration $c_p < \varepsilon$.
This approximation matches the isotherm in such a way that $q=0$ at $c_p=0$ and also matches the value and the first derivative of the istotherm at $c_p = \varepsilon$, where $\varepsilon$ is a very small number, for example $1e-10$.
The form of approximation and its derivative is given below for $c_p < \varepsilon$:

$$
\begin{aligned}
    c_{p,i,lin}  &= \left(\frac{\varepsilon}{c_{p,\text{ref}}}\right)^{\frac{1}{n_i} - 2} \frac{c_{p,i}}{{c_{p,\text{ref}}}^2} \left( \left(2-\frac{1}{n_i}\right)\varepsilon + c_{p,i}\left(\frac{1}{n}-1\right) \right)  \\
    \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} c_{p,i,lin} q_{\text{max},i} \left( 1 - \sum_j \frac{q_j}{q_{\text{max},j}} \right) - k_{d,i} \frac{q_i}{q_{i,\text{ref}}}
\end{aligned}
$$

For more information on model parameters required to define in CADET file format, see [](#multi-component-sips-config).
For more information on the model and its origin, please refer to {cite}`sips1948`.

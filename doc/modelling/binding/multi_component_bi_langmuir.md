(multi-component-bi-langmuir-model)=

# Multi Component Bi-Langmuir

The multi component Bi-Langmuir model {cite}`Guiochon2006` adds $M - 1$ additional types of binding sites $q_{i,j}$ ($0 \leq j \leq M - 1$) to the Langmuir model (see Section [](#multi-component-langmuir-model)) without allowing an exchange between the different sites $q_{i,j}$ and $q_{i,k}$ ($k \neq j$).
Therefore, there are no competitivity effects between the different types of binding sites and they have independent capacities.

$$
\begin{aligned}
    \frac{\mathrm{d} q_{i,j}}{\mathrm{d} t} &=  k_{a,i}^{(j)}\: c_{p,i}\: q_{\text{max},i}^{(j)} \left( 1 - \sum_{k=0}^{N_{\text{comp}} - 1} \frac{q_{k,j}}{q_{\text{max},k}^{(j)}}\right) - k_{d,i}^{(j)} q_{i,j} & i = 0, \dots, N_{\text{comp}} - 1, \: j = 0, \dots, M - 1.% (0 \leq i \leq N_{\text{comp}} - 1, \: 0 \leq j \leq M - 1).
\end{aligned}
$$

Note that all binding components must have exactly the same number of binding site types $M \geq 1$.
See the Section [](#multi-component-langmuir-model).

Originally, the Bi-Langmuir model is limited to two different binding site types.
Here, the model has been extended to arbitrary many binding site types.

For more information on model parameters required to define in CADET file format, see [](#multi-component-bi-langmuir-config).

(simplified-multi-state-steric-mass-action-model)=

# Simplified Multi-State Steric Mass Action

The simplified multi-state steric mass action is the same as the multi-state SMA model described above (see Section [](#multi-state-steric-mass-action-model)), but with additional assumptions:

- Molecules are only exchanged between two adjacent states, that is, no transfer from state $q_{i,1}$ to state $q_{i,3}$ is allowed.
- Characteristic charge $\nu_{i,j}$ and shielding factor $\sigma_{i,j}$ only depend on the index of the state $j$.

Thus, the exchange parameters $k^{(i)}_{j\ell}$, the characteristic charge $\nu_{i,j}$, and the shielding $\sigma_{i,j}$ can be parameterized with few degrees of freedom.
For all $i = 1,\dots,N_{\text{comp}} - 1$ and $j,\ell = 0,\dots,M_i - 1$ let

$$
\begin{aligned}
    k^{(i)}_{j\ell} &= \begin{cases}
    0, & \text{for } \left\lvert j-\ell\right\rvert \neq 1 \\
    K^{(i)}_{ws} + j K^{(i)}_{ws,\text{lin}} - K^{(i)}_{ws,\text{quad}} j(j - M_i+2), & \text{for } \ell = j+1 \\
    K^{(i)}_{sw} + \ell K^{(i)}_{sw,\text{lin}} - K^{(i)}_{sw,\text{quad}} \ell(\ell - M_i+2), & \text{for } \ell = j-1, \end{cases}\\
    \nu_{i,j} &= \nu_{\text{min},i} + \frac{j}{M_i-1} \left( \nu_{\text{max},i} - \nu_{\text{min},i} \right) - \nu_{\text{quad},i} j (j-M_i+1), \\
    \sigma_{i,j} &= \sigma_{\text{min},i} + \frac{j}{M_i-1} \left( \sigma_{\text{max},i} - \sigma_{\text{min},i} \right) - \sigma_{\text{quad},i} j (j-M_i+1).
\end{aligned}
$$

Note that the characteristic charge $\nu_{i,j}$ has to be monotonically non-decreasing in the second index $j$ and all other rates and the steric factor $\sigma_{i,j}$ have to be non-negative.

For more information on model parameters required to define in CADET file format, see [](#simplified-multi-state-steric-mass-action-config).

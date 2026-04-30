(michaelis-menten-kinetics-model)=

# Michaelis Menten kinetics

Implements Michaelis-Menten reaction kinetics of the form

$$
\begin{aligned}
    f_\text{react} = S \nu,
\end{aligned}
$$

where $S$ is the stoichiometric matrix and $\nu$ a flux vector with

$$
\begin{aligned}
    \nu_{j} = v_{\mathrm{max},j} \prod_{i = 1}^{N_{sub,j}} \nu_{i,j} = v_{\mathrm{max},j} \prod_{i = 1}^{N_{sub,j}} \frac{ c_{i,j}}{K_{\mathrm{M}_{i,j}} + c_{i,j}}
\end{aligned}
$$

where

- $\nu_{j}$ is the flux of reaction $j$ and $\nu_{i,j}$ is the flux of reaction $j$ with respect to substrate $i$,
- $N_{sub,j}$ is the number of substrates for reaction $j$,
- $v_{\mathrm{max},j}$ is the maximum reaction rate in reaction $j$,
- and $c_{i,j}$ is the concentration of substrate $i$ in reaction $j$.
- $K_{\mathrm{M}_{i,j}}$ is the Michaelis constant for substrate $i$ in reaction $j$.

The selection of which components act as substrates is controlled via the sign of the entries of the stoichiometric matrix.
For more information about the configuration can be found in [](#michaelis-menten-kinetics-config).

In addition, CADET supports three types of inhibition reactions.
In this case the flux $\nu_{i,j}$ can be modified as one of the following:

## Competitive Inhibition

In competitive inhibition, the inhibitor binds at the enzyme's active site. The modified flux expression is:

$$
\begin{aligned}
    \nu_{i,j} =  \frac{c_{i,j}}{K_{\mathrm{M}_{i,j}}\,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{c}_{I_{k}}}) + c_{i,j}},
\end{aligned}
$$

where
: - $c_{i,j}$ is the substrate component and $c_{k}$ is one inhibitor acting on substrate $c_{i,j}$,
  - $K^{c}_{I_{k}}$ is the inhibition constant with respect to inhibitor $c_{k}$ i.e if $K^{c}_{I_{k}} > 0$, component $c_{k}$ acts as an inhibitor to substrate $c_{i,j}$,
  - $\mathcal{I}^{c}_{i,j}$ is the index set of inhibitors for substrate $c_{i,j}$, i.e the indices $k$ where $K^{c}_{I_{k}} > 0$.

## Uncompetitive Inhibition

In an uncompetitive inhibition, the inhibitor binds to the enzyme-substrate complex, preventing the reaction from proceeding. The modified flux expression is:

$$
\begin{aligned}
    \nu_{i,j} = \frac{c_{i,j}}{K_{\mathrm{M}_{i,j}} + c_{i,j} \, (1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{uc}_{I_{k}}})},
\end{aligned}
$$

where
: - $c_{i,j}$ is the substrate component and $c_{k}$ is one inhibitor acting on substrate $c_{i,j}$.
  - $K^{uc}_{I_{k}}$ is the inhibition constant with respect to component $c_{k}$ in reaction $j$ i.e if $K^{uc}_{I_{k}} > 0$, component $c_{k}$ acts as an inhibitor to substrate $c_{i,j}$.
  - $\mathcal{I}^{uc}_{i,j}$ is the index set of inhibitors in reaction $j$, i.e the indices $k$ where $K^{uc}_{I_{k}} > 0$.

## Mixed Inhibition

In mixed inhibition, the inhibitor can bind to both the enzyme and the enzyme-substrate complex, preventing the reaction from proceeding. The modified flux expression is:

$$
\begin{aligned}
   \nu_{i,j} =  \frac{c_{i,j}}{ K_{\mathrm{M}_{i,j}} \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{c}_{I_{k}}}) + c_{i,j} \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{uc}_{I_{k}}})},
\end{aligned}
$$

where
: - $c_{i,j}$ is the substrate component and $c_{k}$ is one the inhibitor acting on substrate $c_{i,j}$.
  - $K^{c}_{I_{k}}$ and $K^{uc}_{I_{k}}$ are the inhibition constants with respect to component $c_{k}$ in reaction $j$ i.e if $K^{c}_{I_{k}} > 0$, component $c_{k}$ acts as an inhibitor to substrate $c_{i,j}$.
  - $\mathcal{I}^{c}_{i,j}$ is the index set of inhibitors in reaction $j$, i.e the indices $k$ where $K^{c}_{I_{k}} > 0$,
  - $\mathcal{I}^{uc}_{i,j}$ is the index set of inhibitors in reaction $j$, i.e the indices $k$ where $K^{uc}_{I_{k}} > 0$.

## Non-Competitive Inhibition

Non-competitive inhibition is a form of mixed inhibition where the inhibitor binds to both the enzyme and the enzyme-substrate complex with the **same affinity**.

$$
\begin{aligned}
   \nu_{i,j} =  \frac{c_{i,j}}{(K_{\mathrm{M}_{i,j}} + c_{i,j}) \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{n}_{I_{k}}})}
\end{aligned}
$$

where
: - $c_{i,j}$ is the substrate component and $c_{k}$ is one the inhibitor acting on substrate $c_{i,j}$.
  - $K^{n}_{I_{k}}$ is the inhibition constant with respect to component $c_{k}$ in reaction $j$ i.e if $K^{n}_{I_{k}} > 0$, component $c_{k}$ acts as an inhibitor to substrate $c_{i,j}$.
  - $\mathcal{I}_{i,j}$ is the index set of inhibitors in reaction $j$, i.e the indices $k$ where $K^{n}_{I_{k}} > 0$.

Note that the inhibition constant for the non-competitive inhibition is indirectly given if $K^{c}_{I_{k}} = K^{uc}_{I_{k}} = K^{n}_{I_{k}}$

For configuration information please refer to [](#michaelis-menten-kinetics-config).

## Monod Kinetics

The Michaelis-Menten model in CADET can also be used to represent Monod kinetics, as the mathematical formulation is identical.
The Monod equation is commonly used to describe microbial growth processes and corresponds to Michaelis-Menten kinetics with only one substrate.

The Monod equation for microbial growth has the form:

$$
\begin{aligned}
    \mu = \frac{\mu_{\mathrm{max}} \, c_S}{K_S + c_S}
\end{aligned}
$$

where:

- $\mu$ is the specific growth rate
- $\mu_{\mathrm{max}}$ is the maximum specific growth rate
- $c_S$ is the substrate concentration
- $K_S$ is the saturation constant (half-saturation constant)

By choosing a Michaelis-Menten kinetics configuration with one substrate and setting the the parameter accordingly,
i.e $\mu_{\mathrm{max}} = v_{\mathrm{max}}$, $K_S = K_{\mathrm{M}_{0,0}}$ and $c_S = c_{0,0}$,
the Monod equation can be expressed in the same form as the Michaelis-Menten kinetics.

## Literature

- Segel, I. H. (1993). Enzyme kinetics: Behavior and analysis of rapid equilibrium and steady-state enzyme systems. John Wiley & Sons.
- Monod, Jacques. 1949. “The Growth of Bacterial Cultures.” *Annual Review of Microbiology* 3 (1): 371–394. <https://doi.org/10.1146/annurev.mi.03.100149.002103>

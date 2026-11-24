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
More information about the configuration can be found in {ref}`michaelis-menten-kinetics-config`.

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

For configuration information please refer to {ref}`michaelis-menten-kinetics-config`.

## Prefactorial extension

In some cases a different component, like biomass, may have a prefactorial effect on the reaction rate, which can be represented by a prefactorial extension of the Michaelis-Menten kinetics:

$$
\begin{aligned}
    \nu_{i,j} = \prod_{q = 1}^{N_{q,j}} K_{q,j} c_{q,j} \cdot v_{\mathrm{max},j} \prod_{i = 1}^{N_{sub,j}} \frac{ c_{i,j}}{K_{\mathrm{M}_{i,j}} + c_{i,j}}
\end{aligned}
$$

where:

- $c_{q,j}$ is biomass concentration (or any other component that has a prefactorial effect on the reaction rate)
- $K_{q,j}$ is the prefactorial constant (if no prefactorial effect it is set as $0$)
- $N_{q,j}$ components that are **not** in $N_{sub,j}$ but have a prefactorial effect on the reaction rate

While those prefactorial components are not contributing to the reaction rate in the same way as substrates, they can still influence the overall reaction rate by acting as any form of inhibitor:

$$
\begin{aligned}
    \nu_{i,j} = \prod_{q = 1}^{N_{q,j}} K_{q,j} c_{q,j} \cdot v_{\mathrm{max},j} \prod_{i = 1}^{N_{sub,j}} \frac{c_{i,j}}{ K_{\mathrm{M}_{i,j}} \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{c}_{I_{k}}}) + c_{i,j} \,(1 + \sum_{k \in \mathcal{I}_{i,j}} \frac{c_{k}}{K^{uc}_{I_{k}}})},
\end{aligned}
$$

where:

- $N_{q,j}$ $\cap$ $N_{sub,j} = 0$, i.e. the prefactorial components which are not substrates
- $N_{q,j}$ $\cap$ $\mathcal{I}_{i,j} \neq 0$, i.e. the prefactorial components which **can** act as inhibitors

The Michaelis-Menten model with defined prefactor for biomass component in CADET can also be used to represent **Monod kinetics**, as the mathematical formulation is identical.

## Hill kinetics

For a substrate without competitive or uncompetitive inhibition, the optional Hill exponent $n_{i,j}$ changes its contribution to reaction $j$ to

$$
\nu_{i,j} = \frac{c_{i,j}^{n_{i,j}}}{K_{\mathrm{M}_{i,j}}^{n_{i,j}} + c_{i,j}^{n_{i,j}}}.
$$

The exponent must be positive. An exponent of one gives the usual Michaelis-Menten term. Exponents greater than one describe positive cooperativity; exponents between zero and one describe negative cooperativity. Different substrates in the same reaction may have different exponents. A substrate with inhibition must use an exponent of one.

## Literature

- Segel, I. H. (1993). Enzyme kinetics: Behavior and analysis of rapid equilibrium and steady-state enzyme systems. John Wiley & Sons.
- Monod, Jacques. 1949. “The Growth of Bacterial Cultures.” *Annual Review of Microbiology* 3 (1): 371–394. <https://doi.org/10.1146/annurev.mi.03.100149.002103>

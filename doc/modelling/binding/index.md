(binding-models)=

# Binding models

The following binding models are presented in dynamic binding mode.
By replacing all occurrences of $\mathrm{d}q / \mathrm{d}t$ with $0$, quasi-stationary (rapid-equlibrium) binding mode is achieved.
In quasi-stationary binding, it is assumed that ad- and desorption take place on a much faster time scale than the other transport processes such that bead liquid phase $c_{p,i}$ (or bulk liquid phase $c_i$ for certain unit operation models) are always in equilibrium with the solid phase $q_i$.

## Equilibrium constants

For the quasi-stationary binding mode, adsorption and desorption rate are no longer separate entities.
Instead, the quotient $k_{\text{eq}} = k_a / k_d$ of adsorption and desorption coefficient is the relevant parameter as shown for the linear binding model (see Section [](#linear-model)):

$$
\begin{aligned}
    \frac{\mathrm{d} q_i}{\mathrm{d} t} &= k_{a,i} c_{p,i} - k_{d,i} q_i \qquad \Rightarrow 0 = k_{a,i} c_{p,i} - k_{d,i} q_i \qquad \Leftrightarrow q_i = \frac{k_{a,i}}{k_{d,i}} c_{p,i} = k_{\text{eq},i} c_{p,i}.
\end{aligned}
$$

The equilibrium constant $k_{\text{eq},i}$ is used in CADET by setting $k_{d,i} = 1$ and $k_{a,i} = k_{\text{eq},i}$.

Note that adsorption rate $k_{a,i}$ and desorption rate $k_{d,i}$ are linearly correlated in both binding modes due to the form of the equilibrium constant $k_{\text{eq}}$:

$$
\begin{aligned}
    k_{a,i} = k_{\text{eq}} k_{d,i}.
\end{aligned}
$$

This correlation can potentially degrade performance of some optimization algorithms.
While in quasi-stationary binding mode this is prevented by using the technique above, a dynamic binding model has to be reparameterized in order to decouple parameters:

$$
\begin{aligned}
    \frac{\mathrm{d} q_i}{\mathrm{d} t} &= k_{a,i} c_{p,i} - k_{d,i} q_i = k_{d,i} \left[ k_{\text{eq},i} c_{p,i} - q_i \right] = k_{a,i} \left[ c_{p,i} - \frac{1}{k_{\text{eq},i}} q_i \right].
\end{aligned}
$$

This can be achieved by a (nonlinear) parameter transform

$$
\begin{aligned}
    F\left( k_{\text{eq},i}, k_{d,i} \right) = \begin{pmatrix} k_{\text{eq},i} k_{d,i} \\ k_{d,i} \end{pmatrix} \text{ with Jacobian } J_F\left( k_{\text{eq},i}, k_{d,i} \right) = \begin{pmatrix} k_{d,i} & k_{\text{eq},i} \\ 0 & 1 \end{pmatrix}.
\end{aligned}
$$

(ldf-model)=

## Linear Driving Force (LDF)

Some authors use the linear driving force (LDF) approximation instead of the native kinetic form of an isotherm.
All three approaches are equivalent in rapid equilibrium (`IS_KINETIC = 0`) but not equivalent when binding kinetics are considered (`IS_KINETIC = 1`).

1. In the native approach, $\frac{dq}{dt}$ is an explicit funtion of $c$ and $q$. For example $\frac{dq}{dt}=k_a c (q_m - q)-k_d q$ in the Langmuir model.

2\. A linear driving force approximation is based on the equilibrium concentration $q^*$ for given $c$.
For example $q^*= \frac{q_m k_{eq} c}{1 + k_{eq} c}$ in the Langmuir model.
Here, $\frac{dq}{dt}$ is proportional to the actual difference from equilibrium, i.e. $\frac{dq}{dt} = k_{kin}(q^*-q)$.
Note that, the sign of $\frac{dq}{dt}$ causes the resulting flux to act towards the equilibrium.
In CADET, this approach is denoted by `LDF`, for example in `MULTI_COMPONENT_LANGMUIR_LDF`.

3\. An alternative linear driving force approximation is based on the equilibrium concentration $c^*$ for given $q$.
For example $c^*=\frac{q}{k_{eq} \left(q_{m}-q\right)}$ in the Langmuir model.
Here, $\frac{dq}{dt}$ is proportional to the actual difference from equilibrium concentrations, i.e. $\frac{dq}{dt} = k_{kin}(c-c^*)$.
Note that, the sign of $\frac{dq}{dt}$ causes the resulting flux to act towards the equilibrium.
In CADET, this approach is denoted by `LDF_LIQUID_PHASE`, for example in `MULTI_COMPONENT_LANGMUIR_LDF_LIQUID_PHASE`.

In both LDF examples, the original rate constants $k_a$ and $k_d$ are replaced by the equilibrium contant $k_{eq}=\frac{k_a}{k_d}$.
The linear driving force approximations are based on a new kinetic constant, $k_{kin}$.

Note that some isotherms, such as the Freundlich model, don't have a native representation in the above sense.
LDF versions are availabe for some but not all binding models implemented in CADET.

(reference-concentrations)=

## Reference concentrations

Some binding models use reference concentrations $c_{\text{ref}}$ and $q_{\text{ref}}$ of the mobile phase modulator (e.g., salt) in the particle liquid and solid phase, respectively.
The reference values are mainly used for normalizing adsorption and desorption rates, but also for other parameters that appear with those concentrations.
They amount to a simple parameter transformation that is exemplified at one equation of the steric mass action binding model

$$
\begin{aligned}
    \frac{\mathrm{d} q_i}{\mathrm{d} t} = k_{a,i} c_{p,i}\bar{q}_0^{\nu_i} - k_{d,i} q_i c_{p,0}^{\nu_i},
\end{aligned}
$$

where $c_{p,0}$ denotes the mobile phase salt concentration and

$$
\begin{aligned}
    \bar{q}_0 = \Lambda - \sum_{j=1}^{N_{\text{comp}} - 1} \left( \nu_j + \sigma_j \right) q_j
\end{aligned}
$$

is the number of available binding sites which is related to the number of bound salt ions.
Using the parameter transformation

$$
\begin{aligned}
    k_{a,i} &= \tilde{k}_{a,i} q_{\text{ref}}^{-\nu_i}, \\
    k_{d,i} &= \tilde{k}_{d,i} c_{\text{ref}}^{-\nu_i},
\end{aligned}
$$

we obtain the modified model equation

$$
\begin{aligned}
    \frac{\mathrm{d} q_i}{\mathrm{d} t} = \tilde{k}_{a,i} c_{p,i} \left(\frac{\bar{q}_0}{q_{\text{ref}}}\right)^{\nu_i} - \tilde{k}_{d,i} q_i \left(\frac{c_{p,0}}{c_{\text{ref}}}\right)^{\nu_i}.
\end{aligned}
$$

This transformation serves as a (partial) nondimensionalization of the adsorption and desorption rates and, by properly choosing the reference concentrations $c_{\text{ref}}$ and $q_{\text{ref}}$, may improve the optimizer performance.

Recommended choices for $c_{\text{ref}}$ are the average or maximum inlet concentration of the mobile phase modifier $c_0$, and for $q_{\text{ref}}$ the ionic capacity $\Lambda$.
Note that setting the reference concentrations to $1.0$ each results in the original binding model.

(dependence-on-external-function-bind)=

## Dependence on external function

A binding model may depend on an external function or profile $T\colon \left[ 0, T_{\text{end}}\right] \times [0, L] \to \mathbb{R}$, where $L$ denotes the physical length of the unit operation, or $T\colon \left[0, T_{\text{end}}\right] \to \mathbb{R}$ if the unit operation model has no axial length.
By using an external profile, it is possible to account for effects that are not directly modeled in CADET (e.g., temperature).
The dependence of each parameter is modeled by a polynomial of third degree. For example, the adsorption rate $k_a$ is really given by

$$
\begin{aligned}
    k_a(T) &= k_{a,3} T^3 + k_{a,2} T^2 + k_{a,1} T + k_{a,0}.
\end{aligned}
$$

While $k_{a,0}$ is set by the original parameter `XXX_KA` of the file format (`XXX` being a placeholder for the binding model), the parameters $k_{a,3}$, $k_{a,2}$, and $k_{a,1}$ are given by `XXX_KA_TTT`, `XXX_KA_TT`, and `XXX_KA_T`, respectively.
The identifier of the externally dependent binding model is constructed from the original identifier by prepending `EXT_` (e.g., `MULTI_COMPONENT_LANGMUIR` is changed into `EXT_MULTI_COMPONENT_LANGMUIR`).
This pattern applies to all parameters and supporting binding models (see {numref}`MBFeatureMatrix`).
Note that the parameter units have to be adapted to the unit of the external profile by dividing with an appropriate power.

Each parameter of the externally dependent binding model can depend on a different external source.
The 0-based indices of the external source for each parameter is given in the dataset `EXTFUN`.
By assigning only one index to `EXTFUN`, all parameters use the same source.
The ordering of the parameters in `EXTFUN` is given by the ordering in the file format specification in Section [](#FFAdsorption).

(binding-model-feature)=

## Binding model feature matrix

A short comparison of the most prominent binding model features is given in {numref}`MBFeatureMatrix`.
The implemented binding models can be divided into two main classes: Single-state and multi-state binding.
While single-state models only have one bound state per component (or less), multi-state models provide multiple (possibly different) bound states for each component, which may correspond to different binding orientations or binding site types.
The models also differ in whether a mobile phase modifier (e.g., salt) is supported to modulate the binding behavior.

(mbfeaturematrix)=

| Binding model | Competitive | Mobile phase modifier | External function | Multi-state |
| --- | --- | --- | --- | --- |
| [Linear](#linear-model) | × | × | ✓ | × |
| [Multi-Component Langmuir](#multi-component-langmuir-model) | ✓ | × | ✓ | × |
| [Multi-Component Langmuir LDF](#multi-component-langmuir-ldf-model) | ✓ | × | ✓ | × |
| [Multi-Component Langmuir LDF Liquid Phase](#multi-component-langmuir-ldf-liquid-phase-model) | ✓ | × | ✓ | × |
| [Mobile Phase Modulator Langmuir](#mobile-phase-modulator-langmuir-model) | ✓ | ✓ | ✓ | × |
| [Extended Mobile Phase Modulator Langmuir](#extended-mobile-phase-modulator-langmuir-model) | ✓ | ✓ | ✓ | × |
| [Multi-Component Bi-Langmuir](#multi-component-bi-langmuir-model) | ✓ | × | ✓ | ✓ |
| [Affinity Energy Distribution](#affinity-energy-distribution) | ✓ | × | ✓ | ✓ |
| [Affinity Complex Titration](#affinity-complex-titration) | ✓ | × | ✓ | × |
| [Multi-Component Bi-Langmuir LDF](#multi-component-bi-langmuir-ldf-model) | ✓ | × | ✓ | ✓ |
| [Multi-Component Anti-Langmuir](#multi-component-anti-langmuir-model) | ✓ | × | ✓ | × |
| [Multi-Component Spreading](#multi-component-spreading-model) | ✓ | × | ✓ | ✓ |
| [Steric Mass Action](#steric-mass-action-model) | ✓ | ✓ | ✓ | × |
| [Multi-State Steric Mass Action](#multi-state-steric-mass-action-model) | ✓ | ✓ | ✓ | ✓ |
| [Simplified Multi-State Steric Mass Action](#simplified-multi-state-steric-mass-action-model) | ✓ | ✓ | × | ✓ |
| [Bi-Steric Mass Action](#bi-steric-mass-action-model) | ✓ | ✓ | ✓ | ✓ |
| [Multi-Component Colloidal](#multi-component-colloidal-model) | ✓ | ✓ | ✓ | × |
| [Generalized Ion Exchange](#generalized-ion-exchange-model) | ✓ | ✓ | ✓ | × |
| [Saska](#saska-model) | × | × | ✓ | × |
| [Self Association](#self-association-model) | ✓ | ✓ | ✓ | × |
| [Freundlich LDF](#freundlich-ldf-model) | × | × | ✓ | × |
| [HIC Water on Hydrophobic Surfaces](#hic-water-on-hydrophobic-surfaces-model) | ✓ | x | ✓ | x |
| [HIC Constant Water Activity](#hic-constant-water-activity-model) | ✓ | x | ✓ | x |
| [HIC Unified](#hic-unified-model) | ✓ | x | ✓ | x |
| [Multi-Component Sips](#multi-component-sips-model) | ✓ | × | ✓ | × |
| [Multi-Component LDF Freundlich](#multi-component-ldf-freundlich-model) | ✓ | × | ✓ | × |
| [Spline Interpolation](#spline-interpolation) | × | × | ✓ | ✓ |


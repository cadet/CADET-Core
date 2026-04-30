(lumped-rate-model-without-pores-model)=

# Lumped rate model without pores (LRM)

The lumped rate model without pores {cite}`Guiochon2006,Felinger2004` deviates from the lumped rate model with pores (see Section [](#lumped-rate-model-with-pores-model)) by neglecting pores completely.
The particle phase $c^p$ is removed and the porosity $\varepsilon_t$ is taken as total porosity

$$
\begin{aligned}
    \varepsilon_t = \varepsilon_c + \left( 1 - \varepsilon_c \right) \varepsilon_p.
\end{aligned}
$$ (TotalPorosity)

The phase ratio is denoted by $\beta_t = \varepsilon_t / (1 - \varepsilon_t)$ accordingly.
The model equations are given by

$$
\begin{aligned}
    \frac{\partial c^\ell_i}{\partial t} + \frac{1}{\beta_t} \frac{\partial}{\partial t} \sum_{m_i} c^s_{i,m_i} &= -u \frac{\partial c^\ell_i}{\partial z} + D_{\text{ax},i} \frac{\partial^2 c^\ell_i}{\partial z^2} + f_{\text{react},i}^\ell\left( c^\ell, c^s \right) + \frac{1}{\beta_t} f_{\text{react},i}^s\left( c^\ell, c^s \right),
\end{aligned}
$$ (ModelColumnLRM)

where $\beta_t = \varepsilon_t / (1 - \varepsilon_t)$ denotes the (total) phase ratio.
The equations are complemented by Danckwerts boundary conditions {cite}`Danckwerts1953`

$$
\begin{aligned}
    u c_{\text{in},i}(t) &= u c^\ell_i(t,0) - D_{\text{ax},i} \frac{\partial c^\ell_i}{\partial z}(t, 0) & \forall t > 0,\\
    \frac{\partial c^\ell_i}{\partial z}(t, L) &= 0 & \forall t > 0.
\end{aligned}
$$

Both quasi-stationary and dynamic binding models are supported:

$$
\begin{aligned}
    \text{quasi-stationary: }& & 0 &= f_{\text{ads}}\left( c^\ell, c^s\right), \\
    \text{dynamic: }& & \frac{\partial c^s}{\partial t} &= f_{\text{ads}}\left( c^\ell, c^s\right) + f_{\text{react}}^s\left( c^\ell, c^s \right).
\end{aligned}
$$

By default, the following initial conditions are applied for all $z \in [0,L]$:

$$
\begin{aligned}
    c^\ell_i(0, z) &= 0, & c^s_{i,m_i}(0,z) &= 0.
\end{aligned}
$$

Note that by setting $\varepsilon_t = 1$, removing all bound states by setting $N_{\text{bnd},i} = 0$ for all components $i$, and applying no binding model, a dispersive plug flow reactor (DPFR) is obtained.
For the specification of flow rate and direction, the same holds as for the general rate model (see Section [](#MUOPGRMflow)).

For information on model parameters see [](#axial-flow-column-1d-config) and [](#particle-model-config).

## Radial flow LRM

The radial flow LRM describes transport of solute molecules through the interstitial column volume by radial convective flow, band broadening caused by radial dispersion, and adsorption to the bead surfaces.

The main assumptions are:

- The shells of the column are homogenous in terms of interstitial volume, fluid flow, and distribution of components.
  Thus, only one spatial coordinate in radial direction $\rho$ is needed and axial transport is neglected in the column bulk volume.
- The bead radii $r_{p}$ are much smaller than the column radius $\mathrm{P}-\mathrm{P}_c$, with $\mathrm{P}$ and $\mathrm{P}_c$ being the inner and outer column radius respectively, and the column length $L$.
  Therefore, the beads can be seen as continuously distributed inside the column (i.e., at each point there is interstitial and bead volume).
- The fluids are incompressible, i.e. the velocity field $\mathrm{V} \colon \mathbb{R}^3 \to \mathbb{R}^3$ submits to $\operatorname{div}\left( \mathrm{V} \right) \equiv 0$.
  That is, the volumetric flow rate at the inner and outer column radius are the same.

Consider a hollow (double walled) column with inner column diameter $\mathrm{P}_c>0$ and outer diameter $\mathrm{P}>\mathrm{P}_c$, filled with spherical beads. The mass balance in the interstitial column volume is described by

$$
\begin{aligned}
    \frac{\partial c^\ell_i}{\partial t} + \frac{1}{\beta_t} \frac{\partial}{\partial t} \sum_{m_i} c^s_{i,m_i} &= -\frac{u}{\rho} \frac{\partial c^\ell_i}{\partial \rho} + D_{\text{rad},i} \frac{1}{\rho} \frac{\partial}{\partial \rho}  \left( \rho \frac{\partial c^\ell_i}{\partial \rho} \right) + f_{\text{react},i}^\ell\left( c^\ell, c^s \right) + \frac{1}{\beta_t} f_{\text{react},i}^s\left( c^\ell, c^s \right),
\end{aligned}
$$ (ModelRadialColumnLRM)

The equations are complemented by Danckwerts boundary conditions {cite}`Danckwerts1953`

$$
\begin{aligned}
    u c_{\text{in},i}(t) &= u c^\ell_i(t,0) - D_{\text{rad},i} \frac{\partial c^\ell_i}{\partial \rho}(t, 0) & \forall t > 0,\\
    \frac{\partial c^\ell_i}{\partial \rho}(t, \mathrm{P}) &= 0 & \forall t > 0.
\end{aligned}
$$

The complementing binding equations are described by the same equations as for the axial LRM.

For information on model parameters see [](#radial-flow-column-1d-config) and [](#particle-model-config).

(general-rate-model-model)=
# General rate model (GRM)

The general rate model is the most comprehensive model of mass transfer in column liquid chromatography, when only the axial coordinate in the column and the radial coordinate in the beads are considered {cite}`Kucera1965,Gu1995,Guiochon2006,Felinger2004`.

The main assumptions are:

- The cross sections of the column are homogeneous in terms of interstitial volume, fluid flow, and distribution of components.
  Thus, only one spatial coordinate in axial direction is needed and radial transport is neglected in the column bulk volume.
- The bead radii $r_{p}$ are much smaller than the column radius $\mathrm{P}$ and the column length $L$.
  Therefore, the beads can be seen as continuously distributed inside the column (i.e., at each point there is interstitial and bead volume).

(table-features)=
```{list-table} Variables and unknowns
:header-rows: 1

* - Variable
  - Domain
  - Description
* - $i$
  - $\left\{ 0, \dots, N_{\text{comp}} - 1 \right\}$
  - Component index
* - $j$
  - $\left\{ 0, \dots, N_{\text{partype}} - 1 \right\}$
  - Particle type index
* - $m_{j,i}$
  - $\left\{ 0, \dots, N_{\text{bnd},j,i} - 1 \right\}$
  - Bound state index of $i$th component in $j$th particle type
* - $m_j$
  - $\left\{ 0, \dots, \sum_{i=0}^{N_{\text{comp}}-1} N_{\text{bnd},j,i} - 1 \right\}$
  - Total bound state index in particle type $j$
* - $t$
  - $\left[0, T_{\text{end}}\right]$
  - Time coordinate
* - $z$
  - $\left[0, L\right]$
  - Axial coordinate
* - $r$
  - $\left[r_{c,j}, r_{p,j}\right]$
  - Generic bead radial coordinate
* - $c^\ell_{i}(t,z)$
  - $\left[0, T_{\text{end}}\right] \times [0, L]$
  - Interstitial concentration of the $i$th component
* - $c^p_{j,i}(t, z, r)$
  - $\left[0, T_{\text{end}}\right] \times [0, L] \times \left[r_{c,j}, r_{p,j}\right]$
  - Mobile phase concentration of the $i$th component in the $j$th particle type
* - $c^s_{j,i,m_{j,i}}(t, z, r)$
  - $\left[0, T_{\text{end}}\right] \times [0,L] \times \left[r_{c,j}, r_{p,j}\right]$
  - Solid phase concentration of the $i$th component's $m_{j,i}$th bound state in particles of type $j$
* - $j_{f,j,i}(t, z)$
  - $\left[0, T_{\text{end}}\right] \times [0, L]$
  - Flux of the $i$th component through stagnant film into the bead of type $j$
```

(modelgrmcolumn)=
:::{figure} column_bulk_model.png
Column bulk model
:::

The GRM describes transport of solute molecules through the interstitial column volume by convective flow, band broadening caused by axial dispersion, mass transfer resistance through a stagnant film around the beads, pore (and surface) diffusion in the porous beads {cite}`Ma1996,Schneider1968a,Miyabe2007`, and adsorption to the inner bead surfaces.

Consider a column of length $L>0$ filled with spherical beads of (possibly) multiple types with radius $r_{p,j} \ll L$ (see {numref}`modelgrmcolumn`), where $j$ is the particle type index. The mass balance in the interstitial column volume is described by

$$
\begin{aligned}
    \frac{\partial c^\ell_i}{\partial t} = -u \frac{\partial c^\ell_i}{\partial z} + D_{\text{ax},i} \frac{\partial^2 c^\ell_i}{\partial z^2} &- \frac{1}{\beta_c} \sum_j d_j \frac{3}{r_{p,j}} k_{f,j,i} \left[ c^\ell_i - c^p_{j,i}(\cdot, \cdot, r_{p,j}) \right] \\
    &+ f_{\text{react},i}^\ell\left(c^\ell\right).
\end{aligned}
$$ (ModelColumn)

Here, $c^\ell_i\colon \left[0, T_{\text{end}}\right] \times [0, L] \rightarrow \mathbb{R}^{\geq 0}$ denotes the concentration in the interstitial column volume, $c^p_{j,i}\colon \left[0, T_{\text{end}}\right] \times [0, L] \times [r_{c,j}, r_{p,j}] \rightarrow \mathbb{R}^{\geq 0}$ the liquid phase concentration in the beads, $k_{f,j,i}\geq 0$ the film diffusion coefficient, $D_{\text{ax},i}\geq 0$ the dispersion coefficient, $u\in\mathbb{R}$ the interstitial velocity, $d_j > 0$ the volume fraction of particle type $j$, and $\frac{1}{\beta_c} = (1 - \varepsilon_c) / \varepsilon_c$ the column phase ratio, where $\varepsilon_c\in (0,1]$ is the column porosity (ratio of interstitial volume to total column volume).
If reactions are considered, the term $f_{\text{react},i}^\ell\left(c^\ell\right)$ represents the net change of concentration $c_i$ due to reactions involving component $i$.

Danckwerts boundary conditions {cite}`Danckwerts1953` are applied to inlet and outlet of the column:

$$
\begin{aligned}
    u c_{\text{in},i}(t) &= u c^\ell_i(t,0) - D_{\text{ax},i} \frac{\partial c^\ell_i}{\partial z}(t, 0) & \forall t > 0,
\end{aligned}
$$ (BCOutlet)

$$
\begin{aligned}
    \frac{\partial c^\ell_i}{\partial z}(t, L) &= 0 & \forall t > 0.
\end{aligned}
$$ (BCInlet)

Note that the outlet boundary condition Eq. {eq}`BCOutlet` is also known as “do nothing” or natural outflow condition.

In the liquid phase of the porous beads (see {numref}`modelgrmcolumn`) the mass balance is given by

$$
\begin{aligned}
    \frac{\partial c^p_{j,i}}{\partial t} &+ \frac{1 - \varepsilon_{p,j}}{F_{\text{acc},j,i} \varepsilon_{p,j}} \frac{\partial}{\partial t} \sum_{m_{j,i}} c^s_{j,i,m_{j,i}} \\
    &= \underbrace{D_{p,j,i} \left[\frac{\partial^2}{\partial r^2} + \frac{2}{r} \frac{\partial}{\partial r} \right]c^p_{j,i}}_{\text{Pore diffusion}} \\
    &+ \underbrace{\frac{1 - \varepsilon_{p,j}}{F_{\text{acc},j,i} \varepsilon_{p,j}} D_{s,j,i} \left[\frac{\partial^2}{\partial r^2} + \frac{2}{r} \frac{\partial }{\partial r} \right] \sum_{m_{j,i}} c^s_{j,i,m_{j,i}} }_{\text{Surface diffusion}} \\
    &+ f_{\text{react},j,i}^p\left( c_j^p, c_j^s \right) + \frac{1 - \varepsilon_{p,j}}{F_{\text{acc},j,i} \varepsilon_{p,j}} f_{\text{react},j,i}^s\left( c_j^p, c_j^s \right),
\end{aligned}
$$ (ModelBead)

where $c^s_{j,i,m_{j,i}}\colon \left[0, T_{\text{end}}\right] \times [0,L] \times [r_{c,j}, r_{p,j}] \rightarrow \mathbb{R}^{\geq 0}$ denotes the solid phase concentration of the $i$th component’s $m_{j,i}$th bound state in the beads of $j$th type, $D_{p,j,i}>0$ the effective diffusion coefficient in the beads, $D_{s,j,i}\geq 0$ the surface diffusion coefficient, $F_{\text{acc},j,i}\geq 0 \in [0,1]$ the pore accessibility factor, and $\varepsilon_{p,j}\in (0,1]$ the particle porosity (ratio of pore volume to total bead volume).
The inner bead radius $r_{c,j} \in [0, r_{p,j})$ is assumed to be $0$ by default, but can be positive in order to account for core-shell particles that have an impermeable core.
Reaction terms in liquid and solid phase are collected in $f_{\text{react},j,i}^p( c_j^p, c_j^s)$ and $f_{\text{react},j,i}^s(c_j^p, c_j^s)$, respectively.

The GRM is used with both quasi-stationary (Eq. {eq}`REqBinding`) and dynamic (Eq. {eq}`DynBinding`) binding models.

$$
\begin{aligned}
    \text{quasi-stationary: } 0 &= f_{\text{ads},j}\left( c^p_j, c^s_j\right)
\end{aligned}
$$ (REqBinding)

$$
\begin{aligned}
    \text{dynamic: } \frac{\partial c^s_j}{\partial t} &= D_{s,j} \left[\frac{\partial^2}{\partial r^2} + \frac{2}{r} \frac{\partial }{\partial r} \right] c^s_{j} \\
    &+ f_{\text{ads},j}\left( c^p_j, c^s_j\right) + f_{\text{react},j}^s\left( c_j^p, c_j^s \right).
\end{aligned}
$$ (DynBinding)

Note that $c^p_j$ and $c^s_j$ denote the vector of all $c^p_{j,i}$ and $c^s_{j,i,m_{j,i}}$, respectively.

The boundary conditions of the bead model the film diffusion and are given for all ${t \in (0,\infty)}$ and $z \in [0,L]$ by

$$
\begin{aligned}
    k_{f,j,i}\left[ c^\ell_i - c^p_{j,i}(\cdot, \cdot, r_{p,j}) \right] &= F_{\text{acc},j,i} \varepsilon_{p,j} D_{p,j,i} \frac{\partial c^p_{j,i}}{\partial r}(\cdot, \cdot, r_{p,j}) \\
    &+ \left( 1 - \varepsilon_{p,j}\right) D_{s,j,i} \sum_{m_{j,i}} \frac{\partial c^s_{j,i,m_{j,i}}}{\partial r}(\cdot, \cdot, r_{p,j}),
\end{aligned}
$$ (BCBeadIn)

$$
\begin{aligned}
    \frac{\partial c^p_{j,i}}{\partial r}(\cdot, \cdot, r_{c,j}) &= 0.
\end{aligned}
$$ (BCBeadCenter)

By default, the following initial conditions are applied for all $z \in [0,L]$ and $r \in \left[r_{c,j}, r_{p,j}\right]$:

$$
\begin{aligned}
    c^\ell_i(0, z) &= 0, & c^p_{j,i}(0, z, r) &= 0, & c^s_{j,i,m_{j,i}}(0,z,r) &= 0.
\end{aligned}
$$ (InitialConditions)

(modelgrmbead)=
:::{figure} column_bead_model.png
Column bead model
:::

(modelgrmstates)=
:::{figure} multiple_bound_states.png
:width: 50%

Binding with multiple bound states
:::

This model follows the unified and composable modeling framework presented in {cite}`Breuer2026`.

See Table {ref}`axial-flow-column-1d-config` and {ref}`particle-model-config`.

(muopgrmmultiparticletypes)=
## Multiple particle types

A particle type has its own set of mass transfer parameters $\varepsilon_{p,j}$, $D_{p,j}$, $D_{s,j}$, etc. (see Eq. {eq}`ModelBead`) and its own binding model $f_{\mathrm{ads}}$ (including a possibly differing number of bound states).
This allows, for example, modeling of particle size distributions or potential applications with differently functionalized beads (e.g., immobilized enzymes).

The distribution of the particle types is governed by their volume fractions $d_j$ in Eq.
{eq}`ModelColumn`. The volume fractions have to sum to $1$:

$$
\begin{aligned}
    \sum_{j=0}^{N_{\text{partype}} - 1} d_j = 1.
\end{aligned}
$$

The particle type volume fractions can be spatially constant throughout the column, or depend on the position inside the column bulk volume.
In the latter case, the user can specify a set of volume fractions for each discretized finite volume cell.
This allows, for example, the placement of smaller particles near the frits.

(muopgrmparticlegeometry)=
## Particle Geometry

In the model above, spherical particles are considered.
Other supported particle forms are cylinders and slabs.
For cylinders, it is assumed that molecules can only enter through the lateral surface (i.e., the caps are sealed).
Slabs are assumed to have two large sides such that molecules enter through the two large faces (i.e., the remaining four small faces are sealed).

All particle forms support core-shell beads that have an impermeable core.
The particles are characterized by their (outer) "radius" $r_{p,j}$ and their (inner) core "radius" $r_{c,j} \in [0, r_{p,j})$.
See {numref}`modelgrmparticlegeometries`.

(modelgrmparticlegeometries)=
:::{figure} column_particle_geometries.png
Particle geometries
:::

For cylinders, the factor $3 / r_{p,j}$ in Eq. ({eq}`ModelColumn`) changes to $2 / r_{p,j}$ and the diffusion operator in Eq. ({eq}`ModelBead`) and Eq. ({eq}`DynBinding`) changes as

$$
\begin{aligned}
    \left[\frac{\partial^2}{\partial r^2} + \frac{2}{r} \frac{\partial }{\partial r} \right] \quad \rightarrow \quad \left[\frac{\partial^2}{\partial r^2} + \frac{1}{r} \frac{\partial }{\partial r} \right].
\end{aligned}
$$

For slabs, the factor $3 / r_{p,j}$ in (see Eq. ({eq}`ModelColumn`)) changes to $1 / r_{p,j}$ and the diffusion operator in Eq. ({eq}`ModelBead`) and Eq. ({eq}`DynBinding`) changes as

$$
\begin{aligned}
    \left[\frac{\partial^2}{\partial r^2} + \frac{2}{r} \frac{\partial }{\partial r} \right] \quad \rightarrow \quad \frac{\partial^2}{\partial r^2}.
\end{aligned}
$$

(muopgrmsizeexclusion)=
## Size exclusion chromatography

The general rate model can be used to simulate size exclusion chromatography (SEC) {cite}`Gu1995`.
The particle porosity $\varepsilon_{p,j}$ on the mobile phase side of the transport equations is replaced by a component-dependent accessible porosity

$$
\begin{aligned}
    \varepsilon_{p,j,i} = F_{\text{acc},j,i} \varepsilon_{p,j},
\end{aligned}
$$

where the pore accessibility factor $F_{\text{acc},j,i}$ ranges in $(0, 1]$.

Small molecules that can enter any pore have $F_{\text{acc},j,i} = 1$, whereas larger molecules that can enter some, but not small pores, have values $0 < F_{\text{acc},j,i} < 1$.
The other extreme is given by molecules so large that they cannot enter any pore and, consequently, $F_{\text{acc},j,i} = 0$.
Note that $F_{\text{acc},j,i} = 0$ is not allowed in a simulation, which can be circumvented by setting $k_{f,j,i} = 0$.

By default, $F_{\text{acc},j,i} = 1$ for all components $i$ and all particle types $j$, which disables size exclusion chromatography.

Note that the solid phase weighting $\left(1 - \varepsilon_{p,j}\right) / \left(F_{\text{acc},j,i} \varepsilon_{p,j}\right)$ in Eq. {eq}`ModelBead` retains the full solid volume fraction $1 - \varepsilon_{p,j}$, that is, the accessible porosity replaces $\varepsilon_{p,j}$ in the pore liquid transport and reaction terms only.
The pore accessibility factor can hence only be combined with a binding model, if the binding capacities remain well-defined.
Physically, a solute that reaches only part of the pore space will also reach only part of the immobilized ligands, so its saturation capacity (e.g., $q_{\text{max}}$ of Langmuir-type binding models) is reduced accordingly.
As long as the affected components do not compete for the same bound states, this reduction is simply expressed in the respective component's own capacity parameter.
If, however, components with different accessible porosities compete for one shared population of binding sites (e.g., multi component Langmuir binding, whose saturation term is shared by all components), no single set of capacity parameters can represent their different accessible ligand fractions consistently.
Binding models should therefore only be combined with the pore accessibility factor if this case is excluded.
In order to simulate pure SEC, binding is disabled by setting $N_{\text{bnd},i} = 0$ for all components $i$ and applying no binding model.

Note that multiple particle types can also be used to aid in modeling size exclusion effects, see Section {ref}`muopgrmmultiparticletypes`.

(muopgrmflow)=
## Specification of flow rate / velocity and direction

Since volumetric flow rates are specified for each network connection, the unit operation always infers its interstitial velocity via

$$
\begin{aligned}
    u = u_{\text{int}} = \frac{F_{\text{in}}}{A \varepsilon_c},
\end{aligned}
$$

where $F_{\text{in}}$ denotes the volumetric flow rate and $A$ the cross section area.
Note that without the bulk porosity $\varepsilon_c$, the superficial velocity would be obtained.

The direction of flow inside the unit operation is governed by the field `FORWARD_FLOW`, where the value $1$ means flow from the larger towards the smaller radius for the radial flow and frustum geometries.
Note that in case of reversed flow, the chromatogram is returned at the unit operation’s `INLET`, which may not be returned from simulation by default.

For information on model parameters see {ref}`axial-flow-column-1d-config` and {ref}`particle-model-config`.

(muopgrmcolumngeometry)=
# Column Geometry

In the model above, a cylindrical axial flow column is considered, see figure {numref}`geometrygrmaxialcolumn`.
Other geometries are being used, and CADET-Core additionally supports cylindrical radial flow columns, frustum columns, and columns whose cross section area varies smoothly along the flow path and is prescribed by the user.

(geometrygrmaxialcolumn)=
:::{figure} geometry_axial_flow_column.png
Axial flow cylindrical column geometry
:::

(muopgrmradialflow)=
## Radial flow GRM

The radial flow GRM describes transport of solute molecules through the interstitial column volume by radial convective flow, band broadening caused by radial dispersion, mass transfer resistance through a stagnant film around the beads, pore (and surface) diffusion in the porous beads {cite}`Ma1996,Schneider1968a,Miyabe2007`, and adsorption to the inner bead surfaces.
Figure {numref}`geometrygrmradialcolumn` shows the geometry of the radial flow column.

(geometrygrmradialcolumn)=
:::{figure} geometry_radial_flow_column.png
Radial flow cylindrical column geometry
:::

The main assumptions are:

- The cylindrical shells of the column are homogeneous in terms of interstitial volume, fluid flow, and distribution of components.
  Thus, only one spatial coordinate in radial direction $\rho$ is needed and axial transport is neglected in the column bulk volume.
- The bead radii $r_{p}$ are much smaller than the column radius $\mathrm{P}-\mathrm{P}_c$, with $\mathrm{P}$ and $\mathrm{P}_c$ being the outer and inner column radius, respectively, and the column length $L$.
  Therefore, the beads can be seen as continuously distributed inside the column (i.e., at each point there is interstitial and bead volume).
- The fluids are incompressible, i.e. the velocity field $\mathrm{V} \colon \mathbb{R}^3 \to \mathbb{R}^3$ satisfies $\operatorname{div}\left( \mathrm{V} \right) \equiv 0$.
  That is, the volumetric flow rate at the inner and outer column radius are the same.

Consider a hollow (double walled) column with inner column diameter $\mathrm{P}_c>0$ and outer diameter $\mathrm{P}>\mathrm{P}_c$, filled with spherical beads of (possibly) multiple types with radius $r_{p,j} \ll L$ (see {numref}`modelgrmcolumn`), where $j$ is the particle type index.
Alternatively, consider a wedge of this cylinder.
The mass balance in the interstitial column volume is described by

$$
\begin{aligned}
    \frac{\partial c^\ell_i}{\partial t} = -\frac{u}{\rho} \frac{\partial c^\ell_i}{\partial \rho} + D_{\text{rad},i} \frac{1}{\rho} \frac{\partial}{\partial \rho} \left(\rho \frac{\partial c^\ell_i}{\partial \rho} \right) &- \frac{1}{\beta_c} \sum_j d_j \frac{3}{r_{p,j}} k_{f,j,i} \left[ c^\ell_i - c^p_{j,i}(\cdot, \cdot, r_{p,j}) \right] \\
    &+ f_{\text{react},i}^\ell\left(c^\ell\right).
\end{aligned}
$$ (ModelRadialColumnGRM)

Here, $c^\ell_i\colon \left[0, T_{\text{end}}\right] \times [\mathrm{P}_c, \mathrm{P}] \rightarrow \mathbb{R}^{\geq 0}$ denotes the concentration in the interstitial column volume, $c^p_{j,i}\colon \left[0, T_{\text{end}}\right] \times [P_c, P] \times [r_{c,j}, r_{p,j}] \rightarrow \mathbb{R}^{\geq 0}$ the liquid phase concentration in the beads, $k_{f,j,i}\geq 0$ the film diffusion coefficient, $D_{\text{rad},i}\geq 0$ the dispersion coefficient, $u>0$ the interstitial velocity, $d_j>0$ the volume fraction of particle type $j$, and $\beta_c = \varepsilon_c / (1 - \varepsilon_c)$ the column phase ratio, where $\varepsilon_c\in(0,1)$ is the column porosity (ratio of interstitial volume to total column volume).
If reactions are considered, the term $f_{\text{react},i}^\ell\left(c^\ell\right)$ represents the net change of concentration $c_i$ due to reactions involving component $i$.

Danckwerts boundary conditions {cite}`Danckwerts1953` are applied to inlet and outlet of the column:

$$
\begin{aligned}
    u c_{\text{in},i}(t) &= u c^\ell_i(t,\mathrm{P}) + D_{\text{rad},i}\, \mathrm{P} \frac{\partial c^\ell_i}{\partial \rho}(t, \mathrm{P}) & \forall t > 0,
\end{aligned}
$$ (BCInletRadial)

$$
\begin{aligned}
    \frac{\partial c^\ell_i}{\partial \rho}(t, \mathrm{P}_c) &= 0 & \forall t > 0.
\end{aligned}
$$ (BCOutletRadial)

Under the default flow direction (`FORWARD_FLOW` $= 1$) the fluid enters at the outer radius $\mathrm{P}$ and leaves at the inner radius $\mathrm{P}_c$, so the inlet condition Eq. {eq}`BCInletRadial` is imposed at $\mathrm{P}$ and the outlet condition Eq. {eq}`BCOutletRadial` at $\mathrm{P}_c$; setting `FORWARD_FLOW` $= 0$ swaps the two boundaries.
Note that the radial coordinate $\rho$, and with it the ordering of the bulk solution, always runs from the inner to the outer radius, i.e. it is oriented against the default flow direction.
Note further that the outlet boundary condition Eq. {eq}`BCOutletRadial` is also known as “do nothing” or natural outflow condition.

The complementing mass transport and binding equations for the liquid and solid phases of the porous beads are described by the same equations as for the axial GRM.

For information on model parameters see {ref}`radial-flow-column-1d-config` and {ref}`particle-model-config`.

## Frustum device GRM

The frustum GRM describes transport of solute molecules through the interstitial column volume by convective flow from the larger to the smaller radius of the frustum device ("axial" direction), band broadening caused by dispersion, mass transfer resistance through a stagnant film around the beads, pore (and surface) diffusion in the porous beads {cite}`Ma1996,Schneider1968a,Miyabe2007`, and adsorption to the inner bead surfaces.
Figure {numref}`geometrygrmfrustumcolumn` shows the geometry of the frustum column.

(geometrygrmfrustumcolumn)=
:::{figure} geometry_frustum_column.png
Frustum column geometry
:::

The main assumptions are:

- The frustum shells of the column are homogeneous in terms of interstitial volume, fluid flow, and distribution of components.
  Thus, only one spatial coordinate in axial direction $x$ is needed.
- The bead radii $r_{p}$ are much smaller than the column length $H$ and radius $r(x) \forall x\in(0,H)$.
  Therefore, the beads can be seen as continuously distributed inside the column (i.e., at each point there is interstitial and bead volume).
- The fluids are incompressible, i.e. the velocity field $\mathrm{V} \colon \mathbb{R}^3 \to \mathbb{R}^3$ satisfies $\operatorname{div}\left( \mathrm{V} \right) \equiv 0$.
  That is, the volumetric flow rate at the inner and outer column radius are the same.
- A cross-sectional average $\overline{u}(x) := \int_A \vec{u} dA$ of the velocity is considered. This is a necessary assumption to derive the one-dimensionalized model.

We consider axial flow through a conical frustum, where the axial coordinate from bottom to top is denoted by $x \in (0, H)$.
The frustum is filled with spherical beads of (possibly) multiple types with radius $r_{p,j} \ll H$ (see {numref}`modelgrmcolumn`), where $j$ is the particle type index.

The flow-directional cross-section area is circular and varies with axial position $A(x) = \pi r(x)^2$, with $r(x) = r_0 + \frac{x}{H} (r_H - r_0)$ being the frustum radius at axial position $x$, and $r_H, r_0$ being the radii at the frustum top and bottom, respectively.

The *averaged bulk velocity* at position $x$ is then given by $\overline{u}(x) = \frac{Q}{A(x)} = \frac{Q}{\pi r(x)^2}$.

We define the *velocity coefficient* $u := \frac{Q}{\pi}$, which has units of $\mathrm{m}^3/\mathrm{s}$, and expresses the scaled volumetric flow rate.

The mass balance in the interstitial column volume is described by

$$
\begin{aligned}
    \frac{\partial c}{\partial t}
    = -\frac{u}{r(x)^2} \frac{\partial c}{\partial x}
    + \frac{1}{r(x)^2} \frac{\partial}{\partial x} \left( D(x) r(x)^2 \frac{\partial c}{\partial x} \right)
    &- \frac{1}{\beta_c} \sum_j d_j \frac{3}{r_{p,j}} k_{f,j,i} \left[ c^\ell_i - c^p_{j,i}(\cdot, \cdot, r_{p,j}) \right] \\
    &+ f_{\text{react},i}^\ell\left(c^\ell\right).
\end{aligned}
$$ (ModelFrustumColumnGRM)

Here, $c^\ell_i\colon \left[0, T_{\text{end}}\right] \times (0, H) \rightarrow \mathbb{R}^{\geq 0}$ denotes the concentration in the interstitial column volume, $c^p_{j,i}\colon \left[0, T_{\text{end}}\right] \times (0, H) \times (r_{c,j}, r_{p,j}) \rightarrow \mathbb{R}^{\geq 0}$ the liquid phase concentration in the beads, $k_{f,j,i}\geq 0$ the film diffusion coefficient, $D_{i}\geq 0$ the dispersion coefficient, $u>0$ the interstitial velocity coefficient, $d_j>0$ the volume fraction of particle type $j$, and $\beta_c = \varepsilon_c / (1 - \varepsilon_c)$ the column phase ratio, where $\varepsilon_c\in(0,1)$ is the column porosity (ratio of interstitial volume to total column volume).
If reactions are considered, the term $f_{\text{react},i}^\ell\left(c^\ell\right)$ represents the net change of concentration $c_i$ due to reactions involving component $i$.

Danckwerts boundary conditions {cite}`Danckwerts1953` are applied to inlet and outlet of the column:

$$
\begin{aligned}
    u c_{\text{in},i}(t) &= u c^\ell_i(t,0) - D_{i}\, r(0)^2 \frac{\partial c^\ell_i}{\partial x}(t, 0) & \forall t > 0,
\end{aligned}
$$ (BCInletFrustum)

$$
\begin{aligned}
    \frac{\partial c^\ell_i}{\partial x}(t, H) &= 0 & \forall t > 0.
\end{aligned}
$$ (BCOutletFrustum)

Under the default flow direction (`FORWARD_FLOW` $= 1$) the fluid enters at the large end $x = 0$ and leaves at the small end $x = H$; setting `FORWARD_FLOW` $= 0$ swaps the two boundaries.
Note that the outlet boundary condition Eq. {eq}`BCOutletFrustum` is also known as “do nothing” or natural outflow condition.

The complementing mass transport and binding equations for the liquid and solid phases of the porous beads are described by the same equations as for the axial GRM.

For information on model parameters see {ref}`frustum-flow-column-1d-config` and {ref}`particle-model-config`.

(muopgrmsmoothlyvarying)=
## Smoothly varying cross section GRM

The radial flow and frustum geometries both prescribe the cross section area $A(x)$ through a small number of geometric parameters, namely the end radii and, for the radial flow column, the cylinder height.
The smoothly varying cross section GRM instead lets the user prescribe $A(x)$ directly, which covers columns and devices whose shape is neither a cylinder shell nor a frustum.

We consider axial flow along the coordinate $x \in (0, H)$, with a cross section area $A \colon [0, H] \to \mathbb{R}^{>0}$ that is smooth and strictly positive, and a constant volumetric flow rate $Q$.
The *averaged bulk velocity* at position $x$ is then $\overline{u}(x) = \frac{Q}{A(x)}$, and the mass balance in the interstitial column volume reads

$$
\begin{aligned}
    \frac{\partial c}{\partial t}
    = -\frac{Q}{A(x)} \frac{\partial c}{\partial x}
    + \frac{1}{A(x)} \frac{\partial}{\partial x} \left( D(x) A(x) \frac{\partial c}{\partial x} \right)
    &- \frac{1}{\beta_c} \sum_j d_j \frac{3}{r_{p,j}} k_{f,j,i} \left[ c^\ell_i - c^p_{j,i}(\cdot, \cdot, r_{p,j}) \right] \\
    &+ f_{\text{react},i}^\ell\left(c^\ell\right).
\end{aligned}
$$ (ModelSmoothlyVaryingColumnGRM)

Danckwerts boundary conditions are applied as before, with $A(0)$ in place of $r(0)^2$ in Eq. {eq}`BCInletFrustum`.
Under the default flow direction (`FORWARD_FLOW` $= 1$) the fluid enters at $x = 0$ and leaves at $x = H$.

The area is supplied as the field `CROSS_SECTIONAL_AREA_AT_NODES`, i.e. by its values at the nodes of the DG discretization, and is represented within the scheme by the nodal Lagrange interpolant of those values on each element.
The area therefore lives in the same polynomial space as the solution, and the integrals weighted by $A$ (and by $A D$) that the discretization requires are evaluated by Gauss quadrature that is exact for the resulting integrands.
For a smooth area profile, the interpolation error is of the same order as the discretization error of the scheme itself, so this representation does not reduce the order of convergence.
Since the representation is tied to the DG grid, this geometry is only available for the DG spatial discretization, and the field has to be recomputed whenever the discretization is refined.

For information on model parameters see {ref}`smoothly-varying-column-1d-config` and {ref}`particle-model-config`.

(muopgrmvaryingporosity)=
## Axially varying bulk porosity

A packing whose density varies along the flow path, $\varepsilon_c = \varepsilon_c(x)$, is independent of the column geometry: a cylinder, a cylinder shell, a frustum and a user prescribed cross section may each be packed non-uniformly.
The structure of the mass balance is preserved in either case.
Convection and dispersion both act on the pore space rather than on the geometric cross section, so both fluxes carry the weight $\varepsilon_c A$ in place of $A$.
Introducing the *flow through area*

$$
w(x) = \varepsilon_c(x) A(x),
$$ (FlowThroughArea)

the quantity that is constant along the column is the delivered volumetric flow rate $\dot V = w(x) \overline{u}(x)$, and the mass balance in the interstitial column volume becomes

$$
\begin{aligned}
    \frac{\partial c}{\partial t}
    = -\frac{\dot V}{w(x)} \frac{\partial c}{\partial x}
    + \frac{1}{w(x)} \frac{\partial}{\partial x} \left( D(x) w(x) \frac{\partial c}{\partial x} \right)
    &- \frac{1 - \varepsilon_c(x)}{\varepsilon_c(x)} \sum_j d_j \frac{3}{r_{p,j}} k_{f,j,i} \left[ c^\ell_i - c^p_{j,i}(\cdot, \cdot, r_{p,j}) \right] \\
    &+ f_{\text{react},i}^\ell\left(c^\ell\right).
\end{aligned}
$$ (ModelColumnGRMVariablePorosity)

This is Eq. {eq}`ModelSmoothlyVaryingColumnGRM` with $A$ replaced by $w$ and $Q$ by $\dot V$, and it returns to it for a constant porosity.
A spatially varying bulk porosity is therefore not a separate model, but the same model evaluated with the flow through area in place of the geometric cross section area.

The bulk transport operator depends on the porosity only through the product $w$.
`COL_POROSITY` is prescribed at the DG nodes, in the same ordering as `CROSS_SECTIONAL_AREA_AT_NODES`, and CADET-Core multiplies it with the geometric cross section area of the chosen geometry, so that no other field has to be adjusted.
Since $w$ is then no longer the closed-form polynomial that the analytic mass matrices of the cylinder, the cylinder shell and the frustum are derived from, the element integrals are evaluated by the same Gauss quadrature as for the smoothly varying cross section, for every geometry.
Note that this represents $w$ by its nodal Lagrange interpolant, so a frustum needs $\texttt{POLYDEG} \geq 2$ to represent its own quadratic area without an additional interpolation error.
The rule is sized to integrate a weight of degree `FLOW_AREA_QUADRATURE_POLYDEG` exactly, which defaults to `POLYDEG` and is therefore exact for that interpolant.

What does not follow from the weight is the film transfer term.
Weighting the bulk balance with $w$ gives that term the weight $w (1 - \varepsilon_c) / \varepsilon_c = A (1 - \varepsilon_c)$, the solid cross section area, whereas the time derivative carries $w$.
The two are no longer proportional once the porosity, or through the velocity the film mass transfer coefficient, varies within an element, so that the prefactor cannot be taken out of the element integral.
Inverting the mass matrix of the time derivative then leaves the element wise operator

$$
\begin{aligned}
    T^{(e)} = \left( \mathcal{M}^{w} \right)^{-1} \tilde{\mathcal{M}}^{w (1 - \varepsilon_c) k_{f} / \varepsilon_c},
\end{aligned}
$$ (FilmTransferOperator)

which takes the place of the scalar prefactor, so that the film transfer flux at one discrete point contributes to the bulk residual at every point of its element.
For a porosity and a film mass transfer coefficient that are constant on an element, $T^{(e)}$ is that scalar times the identity; CADET-Core detects this case and keeps evaluating the term node locally.
The operator is assembled wherever $w$ is represented by its nodal interpolant, i.e. for a prescribed cross section and for an axially varying porosity on any geometry.
With a constant porosity on one of the analytic geometries, $w$ is instead integrated in closed form and the term stays node local, even if the film mass transfer coefficient varies within an element through its velocity dependence.
All nodal fields are represented by their Lagrange interpolant, so a smooth porosity profile does not reduce the order of convergence of the scheme, just as a smooth area profile does not.

Three restrictions apply.
The porosity has to be constant over each cross section, since a porosity that varies within the cross section contributes further correlation terms to the cross sectional averaging that the one dimensional model rests on.
It is prescribed at the DG nodes, like the cross section area itself, and is therefore only available for the DG discretization and has to be recomputed whenever the discretization is refined.
And it routes the unit to the variable coefficient operator, so `USE_COLLOCATION_DG` has to be $0$, which is the default for every geometry other than the axial flow cylinder.

For information on model parameters see {ref}`axial-flow-column-1d-config` and {ref}`particle-model-config`.

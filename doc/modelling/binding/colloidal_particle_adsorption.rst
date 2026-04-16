.. _colloidal_particle_adsorption_model:

Colloidal Particle Adsorption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The colloidal particle adsorption (CPA) model describes protein adsorption on ion exchange resins based on colloidal interaction theory :cite:`Briskot2021`.
The model captures three key contributions to adsorption: electrostatic protein-adsorber interactions, lateral protein-protein interactions on the surface, and steric exclusion effects via scaled-particle theory (hard-disc available surface function).

A designated proton component :math:`c_{p,\mathrm{pH}}` (by default the first component, index configurable via ``CPA_PROTON_IDX``) acts as a non-binding pH state, where :math:`\mathrm{pH} = \log_{10}(c_{p,\mathrm{pH}})`.
Optionally, a salt component (configurable via ``CPA_SALT_IDX``) can provide the ionic strength :math:`I_m` from the mobile phase.
Both the proton and salt components must be non-binding.

The kinetic formulation reads for each binding component :math:`i`:

.. math::

    \frac{\mathrm{d} q_{v,i}}{\mathrm{d} t} = k_{\mathrm{kin},i} \left( K_{v,i} \, c_{p,i} - q_{v,i} \right),

where :math:`q_{v,i}` is the volumetric solid phase concentration, :math:`c_{p,i}` is the pore liquid phase concentration, :math:`K_{v,i}` is the volumetric equilibrium constant, and :math:`k_{\mathrm{kin},i}` is the kinetic rate constant.

Multiple bound states per component are not supported.


Inverse Debye length
^^^^^^^^^^^^^^^^^^^^

The inverse Debye length :math:`\kappa` characterises the range of electrostatic interactions:

.. math::

    \kappa = e \sqrt{\frac{2 \, I_m \, N_A}{k_B \, T \, \varepsilon \, \varepsilon_0}},

where :math:`e` is the elementary charge, :math:`I_m` the ionic strength, :math:`N_A` Avogadro's number, :math:`k_B` the Boltzmann constant, :math:`T` the absolute temperature, :math:`\varepsilon` the relative permittivity, and :math:`\varepsilon_0` the vacuum permittivity.


Adsorber surface potential
^^^^^^^^^^^^^^^^^^^^^^^^^^

The adsorber surface potential :math:`\psi_{0,A}` is obtained by solving the electroneutrality condition

.. math::

    \sigma_{I,A}(\psi_{0,A}) = \sigma_D(\psi_{0,A})

via Newton's method, where the ionisable surface charge density is

.. math::

    \sigma_{I,A} = e \, N_A \, \Gamma_L \left( \zeta_L - \frac{1}{1 + 10^{\mathrm{p}K_L - \mathrm{pH}_0}} \right)

with the surface pH

.. math::

    \mathrm{pH}_0 = \mathrm{pH} + \frac{e \, \psi_{0,A}}{\ln(10) \, k_B \, T},

and the diffuse layer charge density is

.. math::

    \sigma_D = 2 \, \varepsilon \, \varepsilon_0 \, \kappa \, \frac{k_B T}{e} \, \sinh\!\left( \frac{e \, \psi_{0,A}}{2 \, k_B \, T} \right).

Here, :math:`\Gamma_L` is the ligand surface density, :math:`\zeta_L` the charge of the fully protonated ligand, and :math:`\mathrm{p}K_L` the dissociation constant of the ligand.


Protein net charge and surface potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The protein net charge :math:`Z_i` depends on pH via a quadratic relation:

.. math::

    Z_i(\mathrm{pH}) = Z_{i,\mathrm{ref}} + Z_{i,\mathrm{lin}} \left(\mathrm{pH} - \mathrm{pH}_{\mathrm{ref}}\right) + Z_{i,\mathrm{quad}} \left(\mathrm{pH} - \mathrm{pH}_{\mathrm{ref}}\right)^2.

The protein surface potential :math:`\psi_{0,i}` is obtained from the Debye–Hückel approximation for a sphere:

.. math::

    \psi_{0,i} = \frac{2 \, k_B T}{e} \operatorname{arcsinh}\!\left( \frac{Z_i \, e^2}{8 \pi \, a_i^2 \, \varepsilon \, \varepsilon_0 \, \kappa \, k_B T} \right),

where :math:`a_i` is the protein hydrodynamic radius.


Distance of closest approach
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The distance of closest approach :math:`\delta_{m,i}` between protein :math:`i` and the adsorber surface is computed analytically from the superposition of the two surface potentials:

.. math::

    \delta_{m,i} = -\frac{1}{\kappa} \ln\!\left( \frac{-2 \, \psi_{0,A} \, \psi_{0,i}}{\psi_{0,A}^2 + \psi_{0,i}^2} \right).


Interaction layer thickness
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The interaction layer thickness :math:`\delta_i` is parameterised in terms of the protein surface charge density :math:`\sigma_{I,i} = Z_i \, e / (4\pi a_i^2)`:

.. math::

    \log_{10}(\delta_i) = \delta_{i,\mathrm{ref}} + \delta_{i,\mathrm{lin}} \left( |\sigma_{I,i}| - |\sigma_{I,i}^{\mathrm{ref}}| \right),

where :math:`\sigma_{I,i}^{\mathrm{ref}} = Z_{i,\mathrm{ref}} \, e / (4\pi a_i^2)`.
The effective adsorption distance is then :math:`d_i^* = \delta_{m,i} + \delta_i / A_{s,i}`.


Protein–adsorber interaction energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The electrostatic interaction energy between protein :math:`i` and the adsorber at distance :math:`\delta_{m,i}` follows the Hogg–Healy–Fuerstenau model:

.. math::

    u_{A,i}(\delta_{m,i}) = \pi \, a_i \, \varepsilon \, \varepsilon_0 \left[ 2 \, \psi_{0,A} \, \psi_{0,i} \ln\!\left(\frac{1 + e^{-\kappa \delta_{m,i}}}{1 - e^{-\kappa \delta_{m,i}}}\right) - \left(\psi_{0,A}^2 + \psi_{0,i}^2\right) \ln\!\left(1 - e^{-2\kappa \delta_{m,i}}\right) \right].


Henry coefficient
^^^^^^^^^^^^^^^^^

The Henry adsorption coefficient :math:`K_{H,i}` is derived from the interaction potential:

.. math::

    K_{H,i} = \frac{k_B T}{u_{A,i}} \left( 1 - \exp\!\left(-\frac{u_{A,i}}{k_B T}\right) \right).


Available surface function (steric blocking)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The available surface function :math:`B_i(\Theta)` follows scaled-particle theory for hard discs on a surface:

.. math::

    B_i(\Theta) = (1 - \Theta) \exp\!\left( -\frac{\pi a_i^2 \, N_A \sum_j \tilde{q}_j + 2\pi a_i \, N_A \sum_j a_j \tilde{q}_j}{1 - \Theta} - \frac{\pi a_i^2 \left(N_A \sum_j a_j \tilde{q}_j\right)^2}{(1 - \Theta)^2} \right),

where :math:`\tilde{q}_j = q_{v,j} / A_{s,j}` denotes the surface concentration, and the total surface coverage is

.. math::

    \Theta = \pi \, N_A \sum_j a_j^2 \, \tilde{q}_j.


Lateral protein–protein interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lateral interactions between adsorbed proteins are modelled via a screened Coulomb (Yukawa) potential on a hexagonal lattice with spacing

.. math::

    D_{\mathrm{hex}} = \sqrt{\frac{2\sqrt{3}}{3 \, N_A \sum_j \tilde{q}_j}}.

The lateral interaction energy for component :math:`i` reads

.. math::

    u_{\mathrm{lat},i} = \frac{3\sqrt{3} \, D_{\mathrm{hex}} \, N_A \, e^{-\kappa D_{\mathrm{hex}}}}{1 - \exp\!\left(-\frac{3\sqrt{3}}{2\pi} \kappa D_{\mathrm{hex}}\right)} \sum_j \tilde{q}_j \, \beta_{ij},

with the pairwise Yukawa coefficient

.. math::

    \beta_{ij} = \frac{Z_{\mathrm{lat},i} \, Z_{\mathrm{lat},j} \, e^2}{4\pi \, \varepsilon \, \varepsilon_0} \cdot \frac{\exp\bigl(\kappa(a_i + a_j)\bigr)}{(1 + \kappa a_i)(1 + \kappa a_j)},

where :math:`Z_{\mathrm{lat},i}` is the lateral charge of component :math:`i`.


Volumetric equilibrium constant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Combining all contributions, the volumetric equilibrium constant is

.. math::

    K_{v,i} = A_{s,i} \left(d_i^* - \delta_{m,i}\right) K_{H,i} \, B_i(\Theta) \, \exp\!\left(-\frac{u_{\mathrm{lat},i}}{k_B T}\right).


Kinetic rate constant
^^^^^^^^^^^^^^^^^^^^^

The kinetic rate constant :math:`k_{\mathrm{kin},i}` is related to the pore diffusion coefficient :math:`D_i`:

.. math::

    k_{\mathrm{kin},i} = \frac{D_i}{2 \left(d_i^* - \delta_{m,i}\right)^2} \cdot \frac{\left(u_{A,i} / (k_B T)\right)^2}{\cosh\!\left(u_{A,i} / (k_B T)\right) - 1}.


Model assumptions and limitations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- One component must serve as a non-binding proton/pH state (index configurable, default 0).
- An optional non-binding salt component can provide the ionic strength from the mobile phase.
- Multiple bound states per component are not supported.
- Physical constants (:math:`e`, :math:`N_A`, :math:`k_B`, :math:`\varepsilon_0`) are hard-coded to CODATA 2018 values.

For more information on model parameters required to define in CADET file format, see :ref:`colloidal_particle_adsorption_config`.

.. _colloidal_particle_adsorption_model:

Colloidal Particle Adsorption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The colloidal particle adsorption (CPA) model describes protein adsorption based on colloidal interaction theory.
The model captures three key contributions to adsorption:

- electrostatic protein-adsorber interactions, 
- lateral protein-protein interactions on the surface, 
- steric exclusion effects via scaled-particle theory (hard-disc available surface function).

A defined proton component :math:`c^p_{\mathrm{H}+}` (by default the first component, index configurable via ``CPA_PROTON_IDX``) acts as a non-binding pH state.
The proton component must be non-binding.

The kinetic formulation reads for each binding component :math:`i`:

.. math::

    \frac{\mathrm{d} c^{s}_{i}}{\mathrm{d} t} = k_{\mathrm{kin},i} \left( K_{v,i} \, c^p_i - c^{s}_{i} \right),

where :math:`c^{s}_{i}` is the volumetric solid phase concentration, :math:`c^p_i` is the pore liquid phase concentration, :math:`K_{v,i}` is the volumetric equilibrium constant, and :math:`k_{\mathrm{kin},i}` is the kinetic rate constant :cite:`Briskot2021_2`.

In rapid-equilibrium mode, the corresponding bound-state equation is algebraic:

.. math::

    0 = c^{s}_{i} - K_{v,i} \, c^p_i.

The adsorption mode is selected by ``IS_KINETIC``.

Multiple bound states per component are not supported.


Ionic strength and activity coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ionic strength :math:`I_m` can be supplied in two ways:

- **Fixed parameter** (default): :math:`I_m` is read from ``CPA_IONIC_STRENGTH`` and the pH is computed directly from the proton component concentration.
  Since the standard definition of pH is based on concentration in mol/L, :math:`\mathrm{mol}/\mathrm{m}^3` is converted to mol/L via the factor :math:`10^{-3}`:

  .. math::

      \mathrm{pH} = -\log_{10}\!\left(c^p_{\mathrm{H}+} \cdot 10^{-3}\right).

  This corresponds to the setting of :cite:`Briskot2021_3`, where pH and :math:`I_m` are prescribed properties of the mobile phase.
  Note that ``CPA_IONIC_STRENGTH`` must always be given, even when it is ignored because ``CPA_IONIC_VALENCE`` is set.

- **Computed from concentrations**: if component charges :math:`z_i` are provided via ``CPA_IONIC_VALENCE``, the ionic strength is computed from the pore-phase concentrations at each time step, and the Davies activity correction is applied to the proton activity,

  .. math::

      I_m = \frac{1}{2} \sum_i z_i^2 \, c^p_i.

  Components that should not contribute to the ionic strength, should be assigned a charge of zero.
  This also switches off the activity correction for that component; in particular, assigning :math:`z = 0` to the proton component reduces :math:`\gamma_{\mathrm{H}^+}` to one.
  For the Davies model, this value is converted from :math:`\mathrm{mol}/\mathrm{m}^3` to :math:`\mathrm{mol\,/ L}` as :math:`I_M = 10^{-3} I_m`. The activity correction is then

  .. math::

      \log_{10}\gamma_i = -0.509 \, z_i^2 \left( \frac{\sqrt{I_M}}{1 + \sqrt{I_M}} - 0.3 \, I_M \right),

  so that

  .. math::

      \mathrm{pH} = -\log_{10}\!\left(\gamma_{\mathrm{H}^+} \, c^p_{\mathrm{H}+} \cdot 10^{-3}\right).

  **The constant 0.509 is the Debye–Hückel slope for water at 25 °C**.  The factor :math:`10^{-3}` converts the pore-phase concentration from :math:`\mathrm{mol}/\mathrm{m}^3` to mol/L, as required by the pH definition.

In both cases, a proton activity less than or equal to :math:`10^{-14}\mathrm{mol}/\mathrm{m}^3` is mapped to pH 14.

Inverse Debye length
^^^^^^^^^^^^^^^^^^^^

The inverse Debye length :math:`\kappa` characterises the range of electrostatic interactions:

.. math::

    \kappa = e \sqrt{\frac{2 \, I_m \, N_A}{k_B \, T \, \varepsilon \, \varepsilon_0}},

where :math:`e` is the elementary charge, :math:`I_m` the ionic strength, :math:`N_A` Avogadro's number, :math:`k_B` the Boltzmann constant, :math:`T` the absolute temperature, :math:`\varepsilon` the relative permittivity, and :math:`\varepsilon_0` the vacuum permittivity.


Adsorber surface potential
^^^^^^^^^^^^^^^^^^^^^^^^^^

The adsorber surface potential :math:`\psi_{0,A}` is obtained by solving the electroneutrality condition :cite:`Briskot2020,Briskot2021_3`

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

The protein surface potential :math:`\psi_{0,i}` follows from the same electroneutrality condition, applied to the protein surface charge density :math:`\sigma_{I,i} = Z_i \, e / (4 \pi a_i^2)`.
Since :math:`\sigma_{I,i}` does not depend on :math:`\psi_{0,i}`, the diffuse layer relation can be inverted analytically,

.. math::

    \psi_{0,i} = \frac{2 \, k_B T}{e} \operatorname{arcsinh}\!\left( \frac{Z_i \, e^2}{8 \pi \, a_i^2 \, \varepsilon \, \varepsilon_0 \, \kappa \, k_B T} \right),

where :math:`a_i` is the protein hydrodynamic radius.


Protein net charge and surface potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The protein net charge :math:`Z_i` depends on pH through a polynomial of
arbitrary degree :math:`P`:

.. math::

    Z_i(\mathrm{pH}) = z_{i,0} + \sum_{k=1}^{P} z_{i,k} \left(\mathrm{pH} - \mathrm{pH}_{\mathrm{ref}}\right)^k.

The coefficients are supplied by ``CPA_EFFECTIVE_CHARGE_COEF`` as a polynomial-order-row-major matrix. Its rows correspond to increasing powers, its columns to components, and :math:`z_{i,0}` is the charge at
:math:`\mathrm{pH}_{\mathrm{ref}}`. The number of matrix rows determines the polynomial degree.


Distance of closest approach
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The distance of closest approach :math:`\delta_{m,i}` between protein :math:`i` and the adsorber surface is defined as the positive distance at which the protein--adsorber interaction potential :math:`u_{A,i}(z)` attains its minimum:

.. math::

    \delta_{m,i} := \operatorname*{arg\,min}_{z>0} u_{A,i}(z).

Setting :math:`\partial u_{A,i}(z) / \partial z = 0` gives the analytical expression

.. math::

    \delta_{m,i} = -\frac{1}{\kappa} \ln\!\left( \frac{-2 \, \psi_{0,A} \, \psi_{0,i}}{\psi_{0,A}^2 + \psi_{0,i}^2} \right).

The analytical minimum requires opposite signs of :math:`\psi_{0,A}` and :math:`\psi_{0,i}` and a logarithm argument strictly between zero and one.


Interaction layer thickness
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The dimensionless interaction layer parameter :math:`\delta_i` corresponds to :math:`\Delta_i = A_{s,i} \left( d_i^* - \delta_{m,i} \right)` of :cite:`Briskot2021_3,Briskot2021_2`, i.e. the thickness of the interaction boundary layer scaled by the specific adsorber surface area.
It is parameterised in terms of the protein surface charge density :math:`\sigma_{I,i} = Z_i \, e / (4\pi a_i^2)`:

.. math::

    \log_{10}(\delta_i) = \log_{10}(\delta_{i,\mathrm{ref}}) + \delta_{i,\mathrm{lin}} \left( |\sigma_{I,i}| - |\sigma_{I,i}^{\mathrm{ref}}| \right),

where :math:`\delta_{i,\mathrm{ref}} > 0` and :math:`\sigma_{I,i}^{\mathrm{ref}} = z_{i,0} \, e / (4\pi a_i^2)`.
The outer edge of the interaction layer is then :math:`d_i^* = \delta_{m,i} + \delta_i / A_{s,i}`, so that :math:`A_{s,i} \left( d_i^* - \delta_{m,i} \right) = \delta_i`.


Protein–adsorber interaction energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The electrostatic interaction energy between protein :math:`i` and the adsorber at distance :math:`\delta_{m,i}` follows from the Derjaguin approximation applied to the constant-charge plate--plate interaction free energy of :cite:`Parsegian1972`, as described in :cite:`Briskot2020,Briskot2021_3`:

.. math::

    u_{A,i}(\delta_{m,i}) = \pi \, a_i \, \varepsilon \, \varepsilon_0 \left[ 2 \, \psi_{0,A} \, \psi_{0,i} \ln\!\left(\frac{1 + e^{-\kappa \delta_{m,i}}}{1 - e^{-\kappa \delta_{m,i}}}\right) - \left(\psi_{0,A}^2 + \psi_{0,i}^2\right) \ln\!\left(1 - e^{-2\kappa \delta_{m,i}}\right) \right].

Note the sign of the second term: it is the constant-charge counterpart of the constant-potential (Hogg--Healy--Fuerstenau) expression, and it is what makes :math:`u_{A,i}` attain a finite minimum at :math:`\delta_{m,i} > 0` for oppositely charged surfaces.


Henry coefficient
^^^^^^^^^^^^^^^^^

The Henry adsorption coefficient :math:`K_{H,i}` is obtained by linearising :math:`u_{A,i}(z)` on :math:`z \in [\delta_{m,i}, d_i^*]` and averaging the Boltzmann factor over the interaction layer :cite:`Briskot2020`:

.. math::

    K_{H,i} = \frac{k_B T}{u_{A,i}} \left( 1 - \exp\!\left(-\frac{u_{A,i}}{k_B T}\right) \right).

In the limit :math:`u_{A,i} \to 0`, :math:`K_{H,i} \to 1`.


Available surface function (steric blocking)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The available surface function :math:`B_i(\Theta)` follows the two-dimensional scaled-particle theory for hard discs of :cite:`Talbot1994`, as used in :cite:`Briskot2021_3`:

.. math::

    B_i(\Theta) = (1 - \Theta) \exp\!\left( -\frac{\pi a_i^2 \, N_A \sum_j \tilde{c}^{s}_{j} + 2\pi a_i \, N_A \sum_j a_j \tilde{c}^{s}_{j}}{1 - \Theta} - \frac{\pi^2 a_i^2 \left(N_A \sum_j a_j \tilde{c}^{s}_{j}\right)^2}{(1 - \Theta)^2} \right),

where :math:`\tilde{c}^{s}_{j} = c^{s}_{j} / A_{s,j}` denotes the surface concentration, and the total surface coverage is

.. math::

    \Theta = \pi \, N_A \sum_j a_j^2 \, \tilde{c}^{s}_{j}.

Note that Eq. (28) of :cite:`Briskot2021_3` carries :math:`\pi a_i^2` instead of :math:`\pi^2 a_i^2` in the second term of the exponent.
This is a typographical error in the publication: only with :math:`\pi^2 a_i^2` does the single-component limit reduce to the classical hard-disc result
:math:`B(\Theta) = (1-\Theta) \exp\!\left( -3\Theta/(1-\Theta) - \Theta^2/(1-\Theta)^2 \right)` of :cite:`Talbot1994`.
CADET-Core implements the corrected form given above.


Lateral protein–protein interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adsorbed proteins are assumed to be arranged on a hexagonal lattice whose centre-to-centre distance :math:`D_{\mathrm{hex}}` follows from the maximum packing density :cite:`Briskot2021_3`,

.. math::

    D_{\mathrm{hex}} = \sqrt{\frac{2\sqrt{3}}{3 \, N_A \sum_j \tilde{c}^{s}_{j}}}.

Summing the screened Coulomb (Yukawa) pair potential over the lattice gives the lateral interaction energy for component :math:`i`,

.. math::

    u_{\mathrm{lat},i} = \frac{3\sqrt{3} \, D_{\mathrm{hex}} \, N_A \, e^{-\kappa D_{\mathrm{hex}}}}{1 - \exp\!\left(-\frac{3\sqrt{3}}{2\pi} \kappa D_{\mathrm{hex}}\right)} \sum_j \tilde{c}^{s}_{j} \, \beta_{ij},

with

.. math::

    \beta_{ij} = \frac{Z_{\mathrm{lat},i} \, Z_{\mathrm{lat},j} \, e^2}{4\pi \, \varepsilon \, \varepsilon_0} \cdot \frac{\exp\bigl(\kappa(a_i + a_j)\bigr)}{(1 + \kappa a_i)(1 + \kappa a_j)},

where :math:`Z_{\mathrm{lat},i}` is the lateral charge of component :math:`i`.


Volumetric equilibrium constant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Combining all contributions, the volumetric equilibrium constant is :cite:`Briskot2021_2`

.. math::

    K_{v,i} = A_{s,i} \left(d_i^* - \delta_{m,i}\right) K_{H,i} \, B_i(\Theta) \, \exp\!\left(-\frac{u_{\mathrm{lat},i}}{k_B T}\right)
            = \delta_i \, K_{H,i} \, B_i(\Theta) \, \exp\!\left(-\frac{u_{\mathrm{lat},i}}{k_B T}\right).


Kinetic rate constant
^^^^^^^^^^^^^^^^^^^^^

The kinetic rate constant :math:`k_{\mathrm{kin},i}` is calculated from the kinetic prefactor :math:`k^*_{\mathrm{kin},i}` supplied by ``CPA_KKIN``:

.. math::

    k_{\mathrm{kin},i} = \frac{k^*_{\mathrm{kin},i}}{2} \cdot \frac{\left(u_{A,i} / (k_B T)\right)^2}{\cosh\!\left(u_{A,i} / (k_B T)\right) - 1}.

This expression is taken from :cite:`Briskot2021_2`; it approaches :math:`k^*_{\mathrm{kin},i}` in the limit :math:`u_{A,i} \to 0`.
The kinetic rate and ``CPA_KKIN`` are not evaluated for bound states in rapid-equilibrium mode.


Model assumptions and limitations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- One component must serve as a non-binding proton/pH state (index configurable, default 0).
- If ``CPA_IONIC_VALENCE`` is provided, the ionic strength is computed from the pore-phase concentrations and the Davies activity correction is applied to the proton activity. Otherwise, ``CPA_IONIC_STRENGTH`` is used as a fixed parameter.
- The degree of the protein charge polynomial is inferred from the number of rows in ``CPA_EFFECTIVE_CHARGE_COEF`` and is the same for all components.
- Kinetic and rapid-equilibrium adsorption can be selected globally or per bound state through ``IS_KINETIC``.
- Multiple bound states per component are not supported.
- Physical constants (:math:`e`, :math:`N_A`, :math:`k_B`, :math:`\varepsilon_0`) are hard-coded to CODATA 2018 values.
- The Davies activity correction uses the Debye--Hückel slope for water at 25 °C. Setting ``CPA_TEMPERATURE`` to a different value while ``CPA_IONIC_VALENCE`` is given triggers a warning.
- Treating the pH as a transported state derived from a proton component is a CADET-Core extension; in :cite:`Briskot2020,Briskot2021_3,Briskot2021_2` the pH is a prescribed property of the mobile phase.

For more information on model parameters required to configure in CADET-Core, see :ref:`colloidal_particle_adsorption_config`.


Literature
^^^^^^^^^^
- :cite:`Briskot2020` — charge regulation, protein--adsorber interaction energy, and the Henry coefficient (linear range).
- :cite:`Briskot2021_3` — extension to the nonlinear range: available surface function, lateral interactions, and the pH dependence of :math:`Z_i` and :math:`\delta_i`.
- :cite:`Briskot2021_2` — kinetic formulation, :math:`k_{\mathrm{kin},i}` and :math:`K_{v,i}`.
- :cite:`Briskot2021_1` — application to protein mixtures.
- :cite:`Talbot1994` — scaled-particle theory for the hard-disc available surface function.
- :cite:`Parsegian1972` — constant-charge double-layer interaction free energy.

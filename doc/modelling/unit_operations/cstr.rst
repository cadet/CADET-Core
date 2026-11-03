.. _cstr_model:

Continuous stirred tank reactor model (CSTR)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: ../equations.rst

The continuous stirred tank reactor model is a basic building block in unit operation networks and often used to model holdup volume.
It describes a well mixed vessel that may hold a distribution of particle types, which makes it equally suited to the batch uptake experiment, in which a known volume of solution is contacted with a known amount of adsorbent and the decay of the liquid phase concentration is recorded.

Assuming that the liquid in the vessel is perfectly mixed, the bulk mass balance is written on the amount of substance held by the bulk liquid rather than on its concentration,

.. math::
    :label: ModelCSTRBulk

    \begin{aligned}
        V^\ell \frac{\mathrm{d} c^\ell_i}{\mathrm{d} t} + \frac{\mathrm{d} V^\ell}{\mathrm{d} t} c^\ell_i &=
        F_{\text{in}} c^\ell_{\text{in},i} - F_{\text{out}} c^\ell_i
        + \Phi_i
        + V^\ell f_{\text{react},i}^\ell\left( c^\ell \right), \\
        \frac{\mathrm{d} V^\ell}{\mathrm{d} t} &= F_{\text{in}} - F_{\text{out}} - F_{\text{filter}},
    \end{aligned}

where :math:`\Phi_i` is what the particles take up, see below, and

- :math:`c^\ell_i` is the bulk liquid phase concentration of component :math:`i`,
- :math:`V^\ell` is the volume of the bulk liquid, which is a state of the model,
- :math:`F_{\text{in}}`, :math:`F_{\text{out}}` and :math:`F_{\text{filter}}` are the volumetric flow rates into the vessel, out of it, and of pure liquid drawn without components.

The additional parameter :math:`F_{\text{filter}}` can be used to model a filtering unit.
Note that it is the user's duty to make sure that the liquid volume does not fall below 0.
If it does, the simulation may fail to run or may produce unreasonable (e.g., unphysical) results.

A vessel of constant volume is obtained by prescribing :math:`V^\ell` instead of its initial value, which requires the flow rates to cancel, :math:`F_{\text{in}} = F_{\text{out}} + F_{\text{filter}}` in every section.
A closed vessel, that is, a batch uptake experiment, is obtained for :math:`F_{\text{in}} = F_{\text{out}} = F_{\text{filter}} = 0`.
Removing all bound states by setting :math:`N_{\text{bnd},j,i} = 0` for all components :math:`i` and particle types :math:`j`, and applying no binding model, results in a simple tank.

Particles in rapid equilibrium with the bulk liquid
"""""""""""""""""""""""""""""""""""""""""""""""""""

A particle type without pores and without a film diffusion resistance is in rapid equilibrium with the bulk liquid.
It has no liquid phase of its own, so its bound states take up directly what leaves the bulk liquid,

.. math::

    \begin{aligned}
        \Phi_i &= - V^{s} \sum_j d_j \sum_{m_{j,i}}^{N_{\text{bnd},j,i}} \frac{\mathrm{d}}{\mathrm{d}t} c^s_{j,i,m_{j,i}}
        + V^{s}\sum_j d_j f_{\text{react},j,i}^s\left( c^\ell, c_j^s \right),
    \end{aligned}

where:

- |component_sij|
- |volume_solid|
- |volume_fac_pat|

The second term corrects the first: the uptake term assumes that whatever the bound states gain comes out of the bulk liquid, which a reaction between bound states does not.
Since such a particle type has no liquid phase of its own, the liquid phase part of a cross-phase reaction of that type acts on the bulk liquid and enters :math:`f_{\text{react},i}^\ell` above, weighted with :math:`V^\ell`.

The binding behaviour of :math:`c_i^{\ell}` and :math:`c^s_j` is described by

.. math::

    \begin{aligned}
        \text{quasi-stationary: }& & 0 &= f_{\text{ads},j}\left( c_i^\ell, c^s_j\right), \\
        \text{dynamic: }& & \frac{\partial c^s_j}{\partial t} &= f_{\text{ads},j}\left( c_i^\ell, c^s_j\right) + f_{\text{react},j}^s\left( c_i^\ell, c_j^s \right),
    \end{aligned}

where :math:`f_{\text{ads},j}` is the adsorption rate in the solid phase of particle type :math:`j`.

Particles with a film diffusion resistance
""""""""""""""""""""""""""""""""""""""""""

A particle type that has a film diffusion resistance carries a liquid phase of its own, and the transport into it is driven by the concentration difference at the particle surface,

.. math::

    \begin{aligned}
        \Phi_i &= - V^p \sum_{j} d_j a_{p,j} k_{f,j,i} \left[ c^\ell_i - c^p_{j,i}\left(\cdot, r_{p,j}\right) \right],
    \end{aligned}

where:

- :math:`c^p_{j,i}` is the particle liquid phase concentration of component :math:`i` in particle type :math:`j`, evaluated at the particle surface :math:`r_{p,j}`,
- :math:`V^p` is the constant total volume of the particles, **including the liquid held in their pores**,
- :math:`a_{p,j}` is the surface to volume ratio of particle type :math:`j`, which is :math:`3 / r_{p,j}` for a sphere without a solid core,
- :math:`k_{f,j,i}` is the film diffusion coefficient.

**The two volumes are disjoint**: :math:`V^\ell` is the liquid outside the particles only, and the liquid in the particle pores counts towards :math:`V^p`.
The bulk porosity follows as :math:`\varepsilon_b = V^\ell / \left( V^\ell + V^p \right)`, which changes over time whenever the liquid volume does.

Writing the balance on the amount of substance rather than on the concentration is what keeps a changing volume compatible with the particle models: the film transfer term carries the constant :math:`V^p`, whereas dividing by :math:`V^\ell` would turn its prefactor :math:`\left( 1 - \varepsilon_b \right) / \varepsilon_b = V^p / V^\ell` into a quantity that depends on a state.

The particle equations are the same as in the general rate model (see Section :ref:`general_rate_model_model`) and in the lumped rate model with pores (see Section :ref:`lumped_rate_model_with_pores_model`): a particle type is spatially resolved if it has pore or surface diffusion and lumped otherwise.
Both quasi-stationary and dynamic binding models are supported:

.. math::

    \begin{aligned}
        \text{quasi-stationary: }& & 0 &= f_{\text{ads},j}\left( c^p_j, c^s_j\right), \\
        \text{dynamic: }& & \frac{\partial c^s_j}{\partial t} &= f_{\text{ads},j}\left( c^p_j, c^s_j\right) + f_{\text{react},j}^s\left( c_j^p, c_j^s \right).
    \end{aligned}

Particle types
""""""""""""""

:ref:`MUOPGRMMultiParticleTypes` types are supported, but particles in rapid equilibrium with the bulk liquid must not be mixed with particle types that have pores, as the two give the tank a different state layout.

By default, the following initial conditions are applied:

.. math::

    \begin{aligned}
        c^\ell_i(0) &= 0, & c^p_{j,i}(0) &= 0, & c^s_{j,i,m_{j,i}}(0) &= 0.
    \end{aligned}

For information on model parameters see :ref:`cstr_config` and :ref:`particle_model_config`.

.. _finite_bath_model:

Finite bath
~~~~~~~~~~~

The finite bath describes a well mixed vessel of constant volume that holds a distribution of particle types.
It combines the homogeneous bulk model of the :ref:`cstr_model` with the particle models of the :ref:`general_rate_model_model` and the :ref:`lumped_rate_model_with_pores_model`, that is, the transport into the particles is limited by film diffusion.
Its primary use is the batch uptake experiment, in which a known volume of solution is contacted with a known amount of adsorbent and the decay of the liquid phase concentration is recorded.
Since the vessel has an inlet and an outlet, the model can equally be used as a stirred tank with film diffusion limited particles inside a unit operation network.

Assuming that the liquid in the vessel is perfectly mixed, the bulk mass balance per unit bulk liquid volume :math:`V^\ell` reads

.. math::
    :label: ModelFiniteBathBulk

    \begin{aligned}
        \frac{\mathrm{d} c^\ell_i}{\mathrm{d} t} &= \frac{1}{V^\ell} \left( F_{\text{in}} c^\ell_{\text{in},i} - F_{\text{out}} c^\ell_i \right)
        - \frac{1 - \varepsilon_b}{\varepsilon_b} \sum_{j} d_j a_{p,j} k_{f,j,i} \left[ c^\ell_i - c^p_{j,i}\left(\cdot, r_{p,j}\right) \right]
        + f_{\text{react},i}^\ell\left( c^\ell \right),
    \end{aligned}

where:

- :math:`c^\ell_i` is the bulk liquid phase concentration of component :math:`i`,
- :math:`c^p_{j,i}` is the particle liquid phase concentration of component :math:`i` in particle type :math:`j`, evaluated at the particle surface :math:`r_{p,j}`,
- :math:`V^\ell` is the constant volume of the bulk liquid,
- :math:`\varepsilon_b = V^\ell / \left( V^\ell + V^p \right)` is the bulk porosity, that is, the ratio of bulk liquid volume to total volume, with :math:`V^p` the total volume of the particles,
- :math:`d_j` is the volume fraction of particle type :math:`j`,
- :math:`a_{p,j}` is the surface to volume ratio of particle type :math:`j`, which is :math:`3 / r_{p,j}` for a sphere without a solid core,
- :math:`k_{f,j,i}` is the film diffusion coefficient,
- :math:`F_{\text{in}}` and :math:`F_{\text{out}}` are the volumetric flow rates into and out of the vessel.

Because the volume is constant, the flow rates have to cancel, :math:`F_{\text{in}} = F_{\text{out}}`.
A closed vessel, that is, a batch uptake experiment, is obtained for :math:`F_{\text{in}} = F_{\text{out}} = 0`.

The particle equations are the same as in the general rate model (see Section :ref:`general_rate_model_model`) and in the lumped rate model with pores (see Section :ref:`lumped_rate_model_with_pores_model`): a particle type is spatially resolved if it has pore or surface diffusion and lumped otherwise.
Note that this is the only difference between the finite bath and the :ref:`cstr_model`, whose particles are in rapid equilibrium with the bulk liquid, that is, without pores and without a film diffusion resistance.
The unit operation is selected accordingly, see :ref:`finite_bath_config`.

Both quasi-stationary and dynamic binding models are supported:

.. math::

    \begin{aligned}
        \text{quasi-stationary: }& & 0 &= f_{\text{ads},j}\left( c^p_j, c^s_j\right), \\
        \text{dynamic: }& & \frac{\partial c^s_j}{\partial t} &= f_{\text{ads},j}\left( c^p_j, c^s_j\right) + f_{\text{react},j}^s\left( c_j^p, c_j^s \right).
    \end{aligned}

By default, the following initial conditions are applied:

.. math::

    \begin{aligned}
        c^\ell_i(0) &= 0, & c^p_{j,i}(0) &= 0, & c^s_{j,i,m_{j,i}}(0) &= 0.
    \end{aligned}

:ref:`MUOPGRMMultiParticleTypes` types are supported.

For information on model parameters see :ref:`finite_bath_config` and :ref:`particle_model_config`.

.. _kumar_langmuir_model:

Kumar-Langmuir
~~~~~~~~~~~~~~

This extension of the Langmuir isotherm (see Section :ref:`multi_component_langmuir_model`) developed in :cite:`Kumar2015` was used to model charge variants of monoclonal antibodies in ion-exchange chromatography.
A non-binding salt component :math:`c_{p,0}` is added to modulate the ad- and desorption process.

.. math::

    \begin{aligned}
        \frac{\mathrm{d} q_i}{\mathrm{d} t} &= k_{a,i} \exp\left( \frac{k_{\text{act},i}}{T} \right) c_{p,i} q_{\text{max},i} \left( 1 - \sum_{j=1}^{N_{\text{comp}} - 1} \frac{q_j}{q_{\text{max},j}} \right) - k_{d,i} \left( c_{p,0} \right)^{\nu_i} q_i && i = 1, \dots, N_{\text{comp}} - 1
    \end{aligned}

In this model, the true adsorption rate :math:`k_{a,i,\text{true}}` is governed by the Arrhenius law in order to take temperature into account

.. math::

    \begin{aligned}
        k_{a,i,\text{true}} = k_{a,i} \exp\left( \frac{k_{\text{act},i}}{T} \right).
    \end{aligned}

Here, :math:`k_{a,i}` is the frequency or pre-exponential factor, :math:`k_{\text{act},i} = E / R` is the activation temperature (:math:`E` denotes the activation energy and :math:`R` the Boltzmann gas constant), and :math:`T` is the temperature.
The characteristic charge :math:`\nu` of the protein is taken into account by the power law.  


For more information on model parameters required to define in CADET file format, see :ref:`kumar_langmuir_config`.

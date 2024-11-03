.. _rate_constant_distribution_theory:

Rate Constant Distribution Theory
---------------------------------

The rate constant distribution (RCD) theory describes a heterogeneous system by assuming a continuum of different adsorption sites, in contrast to the single site presumed in the Thomas model. 
The traditional RCD model assumes that each binding site is governed by a Langmuir isotherm model, while in theory any single-component isotherm models can be used. 

The RCD model reads:

.. math::

    \begin{aligned}
        \frac{\mathrm{d} q}{\mathrm{d} t} = \int_{0}^{\infty} \int_{0}^{\infty} k_a c(t) (\mathbf{q_{\text{max}}}(k_a, k_d) - \mathbf{q}(k_a, k_d, c, t) ) \mathrm{d} k_a \mathrm{d} k_d  - \int_{0}^{\infty} \int_{0}^{\infty} k_d \mathbf{q}(k_a, k_d, c, t) \mathrm{d} k_a \mathrm{d} k_d.
    \end{aligned}

where :math:`q` is the total solid phase concentration, :math:`\mathbf{q_{\text{max}}}(k_a, k_d)` and :math:`\mathbf{q}(k_a, k_d, c, t)` are the binding capacity distribution and the solid phase concentration distribution in the :math:`k_a` - :math:`k_d` domain, respectively.

To solve the above integral-differntial equation, the :math:`k_a` and :math:`k_d` domains can be separately discretized on an equidistant logarithmic grid: 
:math:`\ln k_a^{\text{min}} = \ln k_a^{1} < \ln k_a^{2} < ...< \ln k_a^{N_{ka}-1} < \ln k_a^{N_{ka}} = \ln k_a^{\text{max}}`
and :math:`\ln k_d^{\text{min}} = \ln k_d^{1} < \ln k_d^{2} < ...< \ln k_d^{N_{kd}-1} < \ln k_d^{N_{kd}} = \ln k_d^{\text{max}}`, 
where :math:`k_a^{\text{max}}`, :math:`k_a^{\text{min}}` and :math:`k_d^{\text{max}}`, :math:`k_d^{\text{min}}` are the maximum and minimum adsorption and desorption rate constants considered. :math:`N_{ka}` and :math:`N_{kd}` are the number of nodes considered for the :math:`k_a` and :math:`k_d` domains, respectively.

The discretized RCD model reads:

.. math::

    \begin{aligned}
        \frac{\mathrm{d} q}{\mathrm{d} t} = \sum_{i=1}^{N_{ka}} \sum_{j=1}^{N_{kd}} k_{a}^i c (q_{\text{max}}^{i, j} - q^{i,j}) - \sum_{i=1}^{N_{ka}} \sum_{j=1}^{N_{kd}} k_{d}^{j} q^{i,j}, 
    \end{aligned}

where :math:`q^{i,j}` and :math:`q_{\text{max}}^{i, j}` are the solid phase concentration and binding capacity for binding site :math:`(i,j)`. 

By definition, the total solid phase concentration :math:`q` is a collective consequence of the amount absorbed by each site, rendering

.. math::

    \begin{aligned}
        q = \sum_{i=1}^{N_{ka}} \sum_{j=1}^{N_{kd}} q^{i, j}. 
    \end{aligned}

Therefore, we obtain a governing equation for each discretized binding site:

.. math::

    \begin{aligned}
        \frac{\mathrm{d} q^{i,j}}{\mathrm{d} t} = k_{a}^i c (q_{max}^{i, j} - q^{i,j}) - k_{d}^j q^{i,j}.
    \end{aligned}

This equation is essentially a single-component Langmuir isotherm model. Hence, to configure the entire RCD model, one can utilize the :ref:`multi_component_bi_langmuir_model` implementation and the RCD model is not implemented as a standalone model.

The associated mass balance equation is given by:

.. math::

    \begin{aligned}
        \frac{\mathrm{d} c}{\mathrm{d} t} =  -\frac{1}{\beta} \sum_{i=1}^{N_{ka}} \sum_{j=1}^{N_{kd}} \frac{\mathrm{d} q^{i,j}}{\mathrm{d} t}.
    \end{aligned}

Similar to the :ref:`thomas_model`, the RCD model can also be configured in different unit operation models like :ref:`cstr_model` and :ref:`lumped_rate_model_without_pores_model`.

A major problem with the RCD model lies in its ill-posedness: it is usually difficult to uniquely determine many model parameters in the RCD model. Users should use caution when using the RCD model. 

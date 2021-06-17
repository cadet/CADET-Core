.. _reference_concentrations:

Reference concentrations
~~~~~~~~~~~~~~~~~~~~~~~~

Some binding models use reference concentrations :math:`c_{\text{ref}}` and :math:`q_{\text{ref}}` of the mobile phase modulator (e.g., salt) in the particle liquid and solid phase, respectively.
The reference values are mainly used for normalizing adsorption and desorption rates, but also for other parameters that appear with those concentrations.
They amount to a simple parameter transformation that is exemplified at one equation of the steric mass action binding model

.. math::

    \begin{aligned}
        \frac{\mathrm{d} q_i}{\mathrm{d} t} = k_{a,i} c_{p,i}\bar{q}_0^{\nu_i} - k_{d,i} q_i c_{p,0}^{\nu_i},
    \end{aligned}

where :math:`c_{p,0}` denotes the mobile phase salt concentration and

.. math::

    \begin{aligned}
        \bar{q}_0 = \Lambda - \sum_{j=1}^{N_{\text{comp}} - 1} \left( \nu_j + \sigma_j \right) q_j
    \end{aligned}

is the number of available binding sites which is related to the number of bound salt ions.
Using the parameter transformation

.. math::

    \begin{aligned}
        k_{a,i} &= \tilde{k}_{a,i} q_{\text{ref}}^{-\nu_i}, \\
        k_{d,i} &= \tilde{k}_{d,i} c_{\text{ref}}^{-\nu_i},
    \end{aligned}

we obtain the modified model equation

.. math::

    \begin{aligned}
        \frac{\mathrm{d} q_i}{\mathrm{d} t} = \tilde{k}_{a,i} c_{p,i} \left(\frac{\bar{q}_0}{q_{\text{ref}}}\right)^{\nu_i} - \tilde{k}_{d,i} q_i \left(\frac{c_{p,0}}{c_{\text{ref}}}\right)^{\nu_i}.
    \end{aligned}

This transformation serves as a (partial) nondimensionalization of the adsorption and desorption rates and, by properly choosing the reference concentrations :math:`c_{\text{ref}}` and :math:`q_{\text{ref}}`, may improve the optimizer performance.

Recommended choices for :math:`c_{\text{ref}}` are the average or maximum inlet concentration of the mobile phase modifier :math:`c_0`, and for :math:`q_{\text{ref}}` the ionic capacity :math:`\Lambda`.
Note that setting the reference concentrations to :math:`1.0` each results in the original binding model.


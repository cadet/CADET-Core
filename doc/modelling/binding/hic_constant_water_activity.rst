.. _hic_constant_water_activity_model:

HIC Constant Water Activity
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This model implements the HIC isotherm assuming a constant water activity as described by Jäpel and Buyel :cite:`Jaepel2022`.

.. math::
    \begin{align}
        \beta &= \beta_0 e^{c_{p,0}\beta_1} \\
        \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} c_{p,i} \left( 1 - \sum_j \frac{q_j}{q_{max,j}} \right)^{\nu_i} - k_{d,i} q_i 0.1^{\nu_i \beta}
    \end{align}

- Component :math:`c_0` is assumed to be salt without a bound state.
- Multiple bound states are not supported.
- Components without bound state (i.e., salt and non-binding components) are supported.

For more information on model parameters required to define in CADET file format, see :ref:`hic_constant_water_activity_config`.


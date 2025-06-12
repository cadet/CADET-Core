.. _hic_unified_model:

HIC Unified
~~~~~~~~~~~
This model implements a unified hydrophobic interaction chromatography (HIC) isotherm as described by Jäpel et al. :cite:`Jaepel2025`.
It combines the isotherms proposed by Mollerup :cite:`Mollerup2006`, Deitcher :cite:`Deitcher2010` and the SWA isotherm by Jäpel et al. into one isotherm.
In addition to the first component :math:`c_{p,0}`, which represents salt, the second component :math:`c_{p,1}` represents another non-binding modifier (e.g., pH).

.. math::
    \begin{align}
        k_{a,i} &= k_{a,i,0} \exp\left( k_{a,i,lin} ({c_{p,1}}-{c}_{1,ref})\right)\\
        \nu_i &= \nu_{i,0}  + \nu_{i,lin} ({c_{p,1}}-{c}_{1,ref})\\
        \beta &= \beta_0 \exp\left(\beta_1 c_{p,0}\right) \\
        \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} c_{p,i}\exp\left({k_{p,i} c_{p,i}+k_{s,i} c_{p,0}}\right)\left( 1 - \sum_j \frac{q_j}{q_{max,j}} \right)^{\nu_i}
        -k_{d,i} q_{p,i}(1+\epsilon q_{p,i}) \exp\left(({\rho c_{p,0}})^{\nu_i \beta} \right)
    \end{align}

Model Assumptions and Limitations:

- Component :math:`c_0` is assumed to be salt with no bound state.
- Component :math:`c_1` is assumed to be a modifying component with no bound state. It is measured relative to a reference concentration of :math:`c_{1,ref}`, which defaults to 0.0.
- Multiple bound states are not supported.
- Components without bound state (i.e., salt and additional non-binding components) are supported.

For more information on model parameters required to define in CADET file format, see :ref:`hic_unified_config`.


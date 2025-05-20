.. _affinity_complex_titration:

Affinity Complex Titration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The affinity complex titration (ACT) isotherm is a modified Langmuir isotherm where both the binding capacity and equilibrium constant are dependent on pH via a Hill-type relationship :cite:`Zhang2025ACT`. 
pH, rather than the proton concentration, is treated as a mobile phase modulator with a concentration. Multiple bound state is not supported. 
The current implementation requires the first component to be pH and non-binding. Although the original derivation and the equation shown below is based on pH, the mobile phase modulator can also be any type of salt. 
The ACT isotherm reads:

.. math::

    \begin{aligned}
        \frac{\mathrm{d}q_i}{\mathrm{d}t} = k_{a,i} q_{\text{max},i} \left( f_{A, i}-\sum_{j=1}^{N_{\text{comp}}} \frac{q_j}{q_{\text{max},j}} \right) f_{G,i} c_i - k_{d,i}q_i, 
    \end{aligned}

where :math:`f_{A, i}` is the modification factor for the binding capacity :math:`q_{\text{max}, i}`, and :math:`f_{G,i}` is the modification factor for the equilibrium constant :math:`K_{eq, i} = k_{a,i} / k_{d,i}`.
The modification factors are defined by:

.. math::

    \begin{aligned}
        f_{A, i} =\frac{1}{1+10^{\eta_{A, i} (\mathrm{p}Ka_{A, i}-\mathrm{pH})}} , \quad f_{G, i} =\frac{1}{1+10^{\eta_{G, i} (\mathrm{p}Ka_{G, i}-\mathrm{pH})}}, 
    \end{aligned}

where :math:`\eta_{A, i}` and :math:`\eta_{G, i}` denote the Hill-type coefficients that control the slope of the :math:`q_{max, i}` and :math:`K_{eq, i}` responses as a function of the pH, respectively,
while :math:`\mathrm{p}Ka_{A, i}` and :math:`\mathrm{p}Ka_{G, i}` denote the center of their responses. respectively. For more details and interpretations on these parameters, please refer to :cite:`Zhang2025ACT`. 

For more information on model parameters required to define in CADET file format, see :ref:`affinity_complex_titration_config`.
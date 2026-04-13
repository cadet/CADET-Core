.. _affinity_complex_titration:

Affinity Complex Titration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The affinity complex titration (ACT) isotherm is a modified Langmuir isotherm where both the binding capacity and equilibrium constant are dependent on a mobile phase modulator via a Hill-type relationship in :cite:`Zhang2025ACT`. 
CADET-Core implements a general form of ACT that can use either a raw ion concentration or a negative log ion concentration. 
If the ion of interest is proton and users choose to use negative log proton concentration, then it becomes pH which is the default case originally derived in the paper.
The current implementation requires the first component to be the non-binding mobile phase modulator and does not support multiple bound states. 

The ACT isotherm is given by:

.. math::

    \begin{aligned}
        \frac{\mathrm{d}q_i}{\mathrm{d}t} = k_{a,i} q_{\text{max},i} \left( f_{A, i}-\sum_{j=1}^{N_{\text{comp}}} \frac{q_j}{q_{\text{max},j}} \right) f_{G,i} c_i - k_{d,i}q_i, 
    \end{aligned}

where :math:`f_{A, i}` is the modification factor for the binding capacity :math:`q_{\text{max}, i}`, and :math:`f_{G,i}` is the modification factor for the equilibrium constant :math:`K_{eq, i} = k_{a,i} / k_{d,i}`.

The following two approaches result in the same mathematical model; they only differ in how the modulator is provided.

If pH is used, the modification factors are defined by:

.. math::

    \begin{aligned}
        f_{A, i} =\frac{1}{1+10^{\eta_{A, i} (\mathrm{p}Ka_{A, i}-\mathrm{pH})}} , \quad f_{G, i} =\frac{1}{1+10^{\eta_{G, i} (\mathrm{p}Ka_{G, i}-\mathrm{pH})}}, 
    \end{aligned}

where :math:`\eta_{A, i}` and :math:`\eta_{G, i}` denote the Hill-type coefficients that control the slope of the :math:`q_{max, i}` and :math:`K_{eq, i}` responses as a function of the pH, respectively,
while :math:`\mathrm{p}Ka_{A, i}` and :math:`\mathrm{p}Ka_{G, i}` denote the center of their responses. respectively. For more details and interpretations on these parameters, please refer to :cite:`Zhang2025ACT`.

If ion concentration :math:`c_\mathrm{ion}` is used, the modification factors become:

.. math::

    \begin{aligned}
        f_{A, i} =\frac{1}{1+\left(\frac{c_\mathrm{ion}}{c_{\mathrm{mid},A,i}}\right)^{\eta_{A, i}}} , \quad f_{G, i} =\frac{1}{1+\left(\frac{c_\mathrm{ion}}{c_{\mathrm{mid},G,i}}\right)^{\eta_{G, i}}}.
    \end{aligned}

In this form, :math:`c_{\mathrm{mid},A,i}` and :math:`c_{\mathrm{mid},G,i}` are the midpoint ion concentrations corresponding to ion concentraion changes. 

For more information on model parameters required to define in CADET file format, see :ref:`affinity_complex_titration_config`.
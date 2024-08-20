.. _michaelis_menten_kinetics_model:

Michaelis Menten kinetics
-------------------------

Implements liquid phase Michaelis-Menten reaction kinetics of the form

.. math::

    \begin{aligned}
        f_\text{react} = S \mathbf{\nu},
    \end{aligned}

where :math:`S` is the stoichiometric matrix and the fluxes are given by

.. math::

    \begin{aligned}
        \nu_j = \frac{\mu_{\mathrm{max},j} \, c_S}{k_{\mathrm{MM},j} + c_S},
    \end{aligned}

where

- :math:`j` is the reaction index,
- :math:`c_S` is the substrate component,
- :math:`\mu_{\mathrm{max},j}`, is the limiting rate approached by the system at saturation,
- :math:`k_{\mathrm{MM},j}` is the Michaelis constant, which is defined as the concentration of substrate at which the reaction rate is half ov :math:`\mu_{\mathrm{max},j}`.


In addition, the reaction might be inhibited by other components.
In this case, the flux has the form

.. math::

    \begin{aligned}
        \nu_j = \frac{\mu_{\mathrm{max},j} \, c_S}{k_{\mathrm{MM},j} + c_S} \prod_i \frac{k_{\mathrm{I},j,i}}{k_{\mathrm{I},j,i} + c_{\mathrm{I},i}}.
    \end{aligned}

Here, :math:`k_{\mathrm{I},j,i}` is the inhibition constant w.r.t component :math:`i` and reaction :math:`j`.
If :math:`k_{\mathrm{I},j,i} <= 0`, component :math:`i` does not inhibit the reaction.

.. _fragmentation_model:

Fragmentation Models
~~~~~~~~~~~~~~~~~~~~

For detailed information on the crystallization models implemented in CADET, including aggregation, please refer to :cite:`Zhang2025`.

The fragmentation/breakage model considered here can be combined with :ref:`pbm_model` and / or :ref:`aggregation_model`.
Further, it can be applied in any of the unit operations, specifically in a tank or DPFR.

The fragmentation/breakage crystallization model describes the evolution of the particle number density :math:`n` driven by particle fragmentation.
Here, we consider multiple fragmentation, i.e. the general breakage of a particle into a particle size distribution, based on particle size :math:`x`, which is called internal coordinate.

Size-based fragmentation is governed by the integro-differential equation

.. math::
    :label: FragmentationSizeBased

    \begin{aligned}
        \frac{\partial n(x)}{\partial t} &= 
        - S(x) n(x) + \int_x^{x_{\mathrm{max}}} S(\lambda) b(x | \lambda) n(\lambda) d\lambda.
    \end{aligned}


Here, :math:`x_{\mathrm{end}}` is the maximal considered particle size, :math:`b(x | \lambda)` is the probability density function for the generation of a particle of size :math`:`x` from breakage of a particle of size :math:`\lambda`,
and :math:`S` is the selection function which determines the rate of fragmentation.


.. math::
    :label: selectionFunction

    \begin{aligned}
        S(x) := S_0 x^{3 \alpha},
    \end{aligned}

where :math:`\alpha > 0` reckons the breakage rate as a function of particle volume.


The propability breakage function is defined as

.. math::
    :label: propabilityDensityFunc

    \begin{aligned}
        b(x | \lambda) := 3 x^2 \frac{\gamma}{\lambda^3} \left( \frac{x^3}{\lambda^3} \right)^{\gamma - 2},
    \end{aligned}

where :math:`\gamma > 1` determines the average number of daughter particles into which a mother particle breaks.
Further, :math:`b` satisfies

.. math::
    :label: propabilityDensityFuncConstraints

    \begin{aligned}
        \int_{x_\mathrm{c}}^\lambda x^3 b(x | \lambda) dx &= \lambda^3,
        \\
        N(\lambda) = \int_{x_\mathrm{c}}^\lambda b(x | \lambda) dx &\geq 1,
    \end{aligned}

where :math:`N(\lambda)` is the total number of daughter particles that a mother particle of size :math:`\lambda` generates on average.

For information on model parameters and how to specify the model interface, see :ref:`pbm_config`.


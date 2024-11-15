.. _thomas_model:

Thomas Model
------------

The so-called Thomas model ignores all mass transfer resistances from the resin bead and treats the beads as a homogeneous binding site which follows a single-component isotherm model. 

The Thomas model reads:

.. math::

    \begin{aligned}
        \frac{\mathrm{d} q}{\mathrm{d} t} = f_{\text{ads}}, \quad  \frac{\mathrm{d} c}{\mathrm{d} t} = -\frac{1}{\beta} \frac{\mathrm{d} q}{\mathrm{d} t}, 
    \end{aligned}

where :math:`\beta = V_{\text{L}} / V_{\text{S}}` is the liquid-solid phase ratio, :math:`f_{\text{ads}}` can be any single-component isotherm model, see :ref:`binding_models`. 
The first equation describes adsorption and the second equation accounts for the mass balance of the system. 
The classic Thomas model uses an Langmuir isotherm model which follows the mass action law and, therefore, is considered as a reaction model. 

.. math::

    \begin{aligned}
        f_{\text{ads}} = k_a c q_{\text{max}} (1 - \frac{q}{q_{\text{max}}}) - k_d q . 
    \end{aligned}

The above equations can be solved in CADET using :ref:`cstr_model` and configuring :ref:`multi_component_langmuir_model`.  

The above concept of treating the bead as a binding site can also be applied in other unit operation models. 
Note that these unit operation models usually vary in the complexity of describing the mass transfer resistances in the bead.
Since all mass transfer resistances within the beads are ignored in the Thomas model, it is reasonable to choose :ref:`lumped_rate_model_without_pores_model` and assign an isotherm model. 

We note that the Thomas model was successfully applied in simple systems. For more complex systems, more advanced :ref:`reaction_models` should be considered. 
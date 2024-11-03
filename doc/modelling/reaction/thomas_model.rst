.. _thomas_model:

Thomas Model
------------

The so-called Thomas model ignores all mass transfer resistances from the resin bead and treats the beads as a homogeneous binding site which follows a single-component isotherm model. 

The Thomas model reads:

.. math::

    \begin{aligned}
        \frac{\mathrm{d} q}{\mathrm{d} t} = f_{\text{ads}}, \quad  \frac{\mathrm{d} c}{\mathrm{d} t} = -\frac{1}{\beta} \frac{\mathrm{d} q}{\mathrm{d} t}, 
    \end{aligned}

where :math:`\beta = V_{\text{Liquid}} / V_{\text{Solid phase}}` is the phase ratio, :math:`f_{\text{ads}}` can be any single-component isotherm model defined in :ref:`binding_models`. 
The first equation considers the isotherm while the second equation accounts for the mass balance for the bulk phase. 
The classic Thomas model uses an Langmuir isotherm model which follows the mass action law and, therefore, is considered as a reaction model. 

.. math::

    \begin{aligned}
        f_{\text{ads}} = k_a c q_{\text{max}} (1 - \frac{q}{q_{\text{max}}}) - k_d q . 
    \end{aligned}

The above equations can be solved in CADET using :ref:`cstr_model` and configuring :ref:`multi_component_langmuir_model`.  

The same concept of treating the bead as a binding site can also be applied in other unit operation models. 
Note that these unit operation models usually vary in the complexity of describing the mass transfer resistances in the bead.
Since all mass transfer resistances are ignored in the Thomas model, it is more reasonable to choose :ref:`lumped_rate_model_without_pores_model` and assign an isotherm model. 

The Thomas model has some success in simple systems, but is limited for more advanced systems, in which cases one should consider using more advanced modern models in CADET. 
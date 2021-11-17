.. _multi_component_langmuir_model_ldf:

Multi Component Langmuir LDF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Langmuir :ref:`ldf_model` binding model also includes a saturation term and takes into account the capacity of the resin :cite:`Langmuir1916,Guiochon2006`. 
All components compete for the same binding sites.

.. math::

    \begin{aligned}
        q_i^*=\frac{q_{m,i} k_{eq,i} c_i}{1 + \sum_{j=1}^{n_{comp}}{k_{eq,j} c_j}} && i = 0, \dots, N_{\text{comp}} - 1.
    \end{aligned}


For more information on model parameters required to define in CADET file format, see :ref:`multi_component_langmuir_ldf_config`.

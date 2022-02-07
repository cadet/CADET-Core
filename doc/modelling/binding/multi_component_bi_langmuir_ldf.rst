.. _multi_component_bi_langmuir_model_ldf:

Multi Component Bi-Langmuir LDF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This model is a :ref:`ldf_model` approximation of  multi component Bi-Langmuir model :cite:`Guiochon2006`, that adds :math:`M - 1` *additional* types of binding sites :math:`q_{i,j}` (:math:`0 \leq j \leq M - 1`) to the LDF based Langmuir model (see Section :ref:`multi_component_langmuir_model_ldf`). The implementation follows the same principle as followed in :ref:`multi_component_bi_langmuir_model`.
Adsorbed phase concnetration averaged over the entire bulk volume is given as:

.. math::
    \begin{aligned}
        q_{i,j}^*=\frac{q_{m,i,j} k_{eq,i,j} c_i}{1 + \sum_{k=1}^{N_{comp}}{k_{eq,k,j} c_k}} && i = 0, \dots, N_{\text{comp}} - 1, \: j = 0, \dots, M - 1.% 	           (0 \leq i \leq N_{\text{comp}} - 1, \: 0 \leq j \leq M - 1). 
    \end{aligned}


Note that all binding components must have exactly the same number of binding site types :math:`M \geq 1`.
See the Section :ref:`multi_component_langmuir_model_ldf`.

Originally, the Bi-Langmuir model is limited to two different binding site types.
Here, the model has been extended to arbitrary many binding site types.


For more information on model parameters required to define in CADET file format, see :ref:`multi_component_bi_langmuir_ldf_config`.

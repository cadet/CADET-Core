.. _freundlich_ldf_model:

Freundlich Linear Driving Force (LDF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A freundlich binding model, which is often employed for low concentrations or in analytic settings :cite:`Singh2016` :cite:`Benedikt2019`.

.. math::

    \begin{aligned}
        \frac{\mathrm{d} q_i}{\mathrm{d} t} = kkin_{,i} c_{p,i}^{\frac{1}{n}} && i = 0, \dots, N_{\text{comp}} - 1.
    \end{aligned}

One of the limitation of this isotherm is the undefined first order Jacobian :math:`\left(\frac{dq}{dc_p}\right)` at :math:`c_{p}=0`. To address this issue an approximation of isotherm is considered near the origin. This approximation matches the isotherm in such a way that  $q=0$ at $c_p=0$ and also matches the first derivative of the istotherm at :math:`c_p = \epsilon`, where :math:`\epsilon` is a very small number close to zero, for example :math:`1e-14`. The form of approximation and its derivative is given below for :math:`c_p < \varepsilon`:

.. math::

	\begin{aligned} 
		q = \alpha_0+\alpha_1 c+\alpha_2 c_p^2  
	\end{aligned}
	
	\begin{aligned} 
		\frac{dq}{dc_p} = \alpha_1+ 2 \alpha_2 c_p 
	\end{aligned}

where :math:`\alpha_0` can be determined from the condition :math:`q=0 ~at~ c_p=0`, while :math:`\alpha_1` and :math:`\alpha_2` are determined from :math:`\alpha_1 \varepsilon+\alpha_2 \varepsilon^2 = k_f \varepsilon^{1/n}` and :math:`\alpha_1 + 2 \alpha_2 \varepsilon = \frac{1}{n}k_f \varepsilon^{\frac{1-n}{n}}`.

.. math::

	\begin{aligned}
		\alpha_0 = 0.0
	\end{aligned}	
.. math::
	\begin{aligned}
		\alpha_1 = \frac{2 n-1}{n}k_f \varepsilon^{\frac{1-n}{n}}
	\end{aligned}
.. math::
	\begin{aligned}
		\alpha_2 = \frac{1-n}{n}k_f \varepsilon^{\frac{1-2 n}{n}}
	\end{aligned}

This approximation can be used for any pore phase cocentration $c_p < \epsilon$, including all negative ones.
For more information on model parameters required to define in CADET file format, see :ref:`freundlich_ldf_config`.

.. bibliography::
	:style: unsrt



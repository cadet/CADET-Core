.. _multi_channel_transport_model_config:

Multichannel Transport model (MCT model) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``UNIT_TYPE``

   Specifies the type of unit operation model
   
   ================  ===================================================  =============
   **Type:** string  **Range:** :math:`\texttt{MULTI_CHANNEL_TRANSPORT}`  **Length:** 1
   ================  ===================================================  =============
   
``NCOMP``

   Number of chemical components in the chromatographic medium
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============
   
``REACTION_MODEL_BULK``

   Specifies the type of reaction model of the bulk volume. The model is configured in the subgroup :math:`\texttt{reaction_bulk}`.
   
   ================  ========================================  =============
   **Type:** string  **Range:** See Section :ref:`FFReaction`  **Length:** 1
   ================  ========================================  =============
   
``INIT_C``

   Initial concentrations for each component in all channels of the bulk mobile phase (length :math:`\texttt{NCOMP}`), or for each component in each channel (length :math:`\texttt{NCOMP} \cdot \texttt{NCHANNEL}`, ordering channel-major)

   **Unit:** :math:`\mathrm{mol}\,\mathrm{m}^{-3}`
   
   ================  =========================  =========================================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** :math:`\texttt{NCOMP} / \texttt{NCOMP} \cdot \texttt{NCHANNEL}`
   ================  =========================  =========================================================================
   
``INIT_STATE``

   Full state vector for initialization (optional, :math:`\texttt{INIT_C}` will be ignored; if length is :math:`2\texttt{NDOF}`, then the second half is used for time derivatives)

   **Unit:** :math:`various`
   
   ================  =============================  ==================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** :math:`\texttt{NDOF} / 2\texttt{NDOF}`
   ================  =============================  ==================================================
   
``COL_DISPERSION``

   Axial dispersion coefficient.  In case of a spatially inhomogeneous setting, the :math:`\texttt{SENS_PARTYPE}` field is used for indexing the channel when specifying parameter sensitivities.

   **Unit:** :math:`\mathrm{m}^{2}\,\mathrm{s}^{-1}`
   
   ================  =========================  =========================================================
   **Type:** double  **Range:** :math:`\geq 0`  **Length:** see :math:`\texttt{COL_DISPERSION_MULTIPLEX}`
   ================  =========================  =========================================================
   
``COL_DISPERSION_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{COL_DISPERSION}`. Determines whether :math:`\texttt{COL_DISPERSION}` is treated as component-, channel-, and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{COL_DISPERSION}`.  Valid modes are: 

  1. Component-independent, channel-independent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is 1 
  2. Component-independent, channel-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCHANNEL}` 
  3. Component-dependent, channel-independent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCOMP}` 
  4. Component-dependent, channel-dependent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCOMP} \cdot \texttt{NCHANNEL}`; ordering is channel-major 
  5. Component-independent, channel-independent, section-dependent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NSEC}` 
  6. Component-independent, channel-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCHANNEL} \cdot \texttt{NSEC}`; ordering is section-major 
  7. Component-dependent, channel-independent, section-independent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCOMP} \cdot \texttt{NSEC}`; ordering is section-major 
  8. Component-dependent, channel-dependent, section-dependent; length of :math:`\texttt{COL_DISPERSION}` is :math:`\texttt{NCOMP} \cdot \texttt{NCHANNEL} \cdot \texttt{NSEC}`; ordering is section-channel-major 
   
   =============  ===================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, 7 \}`  **Length:** 1
   =============  ===================================  =============
   
``COL_LENGTH``

   Column length

   **Unit:** :math:`\mathrm{m}`
   
   ================  ======================  =============
   **Type:** double  **Range:** :math:`> 0`  **Length:** 1
   ================  ======================  =============
   
``VELOCITY``

   Velocity

   **Unit:** :math:`\mathrm{m}\,\mathrm{s}^{-1}`

   Indicates flow direction in each channel (forward if value is positive, backward if value is negative), see Section :ref:`MUOPGRMflow2D`).  In case of a spatially inhomogeneous setting, the :math:`\texttt{SENS_PARTYPE}` field is used for indexing the channel when specifying parameter sensitivities.
   
   ================  =============================  ===================================================
   **Type:** double  **Range:** :math:`\mathbb{R}`  **Length:** see :math:`\texttt{VELOCITY_MULTIPLEX}`
   ================  =============================  ===================================================
   
``VELOCITY_MULTIPLEX``

   Multiplexing mode of :math:`\texttt{VELOCITY}`. Determines whether :math:`\texttt{VELOCITY}` is treated as channel- and/or section-independent.  This field is optional. When left out, multiplexing behavior is inferred from the length of :math:`\texttt{VELOCITY}`.  Valid modes are: 

  1. Channel-independent, section-independent; length of :math:`\texttt{VELOCITY}` is 1 
  2. Channel-dependent, section-independent; length of :math:`\texttt{VELOCITY}` is :math:`\texttt{NCHANNEL}` 
  3. Section-dependent; length of :math:`\texttt{VELOCITY}` is :math:`\texttt{NSEC}` 
  4. Channel-dependent, section-dependent; length of :math:`\texttt{VELOCITY}` is :math:`\texttt{NCHANNEL} \cdot \texttt{NSEC}`; ordering is section-major 
   
   =============  ===================================  =============
   **Type:** int  **Range:** :math:`\{0, \dots, 3 \}`  **Length:** 1
   =============  ===================================  =============

``EXCHANGE_MATRIX``

   Exchange matrix

   **Unit:** :math:`\mathrm{s}^{-1}`

   Ordered list containing all exchange rates :math:`e^k_{ij}` for component :math:`k` from compartment :math:`i` to :math:`j` based on the exchange matrix :math:`E^k`. The vector ordering is source channel - destination channel - component (i.e. i-j-k) major.

   .. math::
    
    E^k=\begin{bmatrix} 
    0 & e^k_{12} & \dots & e^k_{1N} \\
    e^k_{21} & \ddots & & \vdots\\
    \vdots & & \ddots & e^k_{(N-1)N}\\
    e^k_{N1} & \dots & e^k_{N(N-1)} & 0 
    \end{bmatrix}    

   ================  ========================  ===============================================
   **Type:** double  **Range:** :math:`[0,1]`  **Length:** :math:`\texttt{NCHANNEL}*\texttt{NCHANNEL}*\texttt{NCOMP}`
   ================  ========================  ===============================================

``NCOL``

   Number of axial column discretization cells
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============
   
``NCHANNEL``

   Number of channels :math:`ij` 
   
   =============  =========================  =============
   **Type:** int  **Range:** :math:`\geq 1`  **Length:** 1
   =============  =========================  =============
      
   
``CHANNEL_CROSS_SECTION_AREAS``

   Cross section areas

   **Unit:** :math:`\mathrm{m}^{2}`

   Defines the cross section area of each channel
   
   ================  ====================== ======================================
   **Type:** double  **Range:** :math:`> 0`  **Length:** :math:`\texttt{NCHANNEL}`
   ================  ====================== ======================================
   
``USE_ANALYTIC_JACOBIAN``

   Determines whether analytically computed Jacobian matrix (faster) is used (value is 1) instead of Jacobians generated by algorithmic differentiation (slower, value is 0)
   
   =============  ===========================  =============
   **Type:** int  **Range:** :math:`\{0, 1\}`  **Length:** 1
   =============  ===========================  =============
   
``LINEAR_SOLVER_BULK``

   Linear solver used for the sparse column bulk block. This field is optional, the best available method is selected (i.e., sparse direct solver if possible).  Valid values are: 

  - :math:`\texttt{DENSE}` Converts the sparse matrix into a banded matrix and uses regular LAPACK. Slow and memory intensive, but always available. 
  - :math:`\texttt{UMFPACK}` Uses the UMFPACK sparse direct solver (LU decomposition) from SuiteSparse. Fast, but has to be enabled when compiling and requires UMFPACK library. 
  - :math:`\texttt{SUPERLU}` Uses the SuperLU sparse direct solver (LU decomposition). Fast, but has to be enabled when compiling and requires SuperLU library. 
   
   ================  =======================================================================  =============
   **Type:** string  **Range:** :math:`\{\texttt{DENSE},\texttt{UMFPACK},\texttt{SUPERLU}\}`  **Length:** 1
   ================  =======================================================================  =============

``RECONSTRUCTION``

   Type of reconstruction method for fluxes
   
   ================  ================================  =============
   **Type:** string  **Range:** :math:`\texttt{WENO}`  **Length:** 1
   ================  ================================  =============
   


.. _langmuir_exchange_model:

Langmuir Exchange Model
~~~~~~~~~~~~~~~~~~~~~~~

The Langmuir Exchange model extends the exchange framework of the Multi-Channel Transport (MCT) model by incorporating capacity-limited, competitive exchange kinetics between channels. 
This model is particularly relevant for systems where the exchange process is governed by a finite capacity in the destination phase.

Mathematical Formulation
^^^^^^^^^^^^^^^^^^^^^^^^

The Langmuir Exchange model modifies the standard MCT exchange term by introducing a Langmuir-type capacity limitation. For a component :math:`i` in channel :math:`l`, the exchange contribution to the mass balance equation becomes:

.. math::
   :label: langmuir_exchange_term

   \frac{\partial c_{i,l}}{\partial t}\Big|_{\text{exchange}} = \sum_{k=1, k \neq l}^{N_k} J_{k \rightarrow l}^i - J_{l \rightarrow k}^i

where the exchange flux :math:`J_{l \rightarrow k}^i` from channel :math:`l` to channel :math:`k` for component :math:`i` is given by:

.. math::
   :label: langmuir_exchange_flux

   J_{l \rightarrow k}^i = k_{lk}^i \cdot c_{i,l} \cdot q_{\max,k}^i \cdot \left(1 - \sum_{j=1}^{N_c} \frac{c_{j,k}}{q_{\max,k}^j}\right) \cdot \frac{A_l}{A_k}

Here:

- :math:`k_{lk}^i` is the exchange rate constant for component :math:`i` from channel :math:`l` to :math:`k` [s⁻¹]
- :math:`c_{i,l}` is the concentration of component :math:`i` in source channel :math:`l` [mol m⁻³]
- :math:`q_{\max,k}^i` is the maximum capacity for component :math:`i` in destination channel :math:`k` [mol m⁻³]
- :math:`\sum_{j=1}^{N_c} \frac{c_{j,k}}{q_{\max,k}^j}` is the total normalized occupancy in destination channel :math:`k` [-]
- :math:`A_l, A_k` are the cross section areas of channels :math:`l` and :math:`k` [m²]

Limiting Cases
^^^^^^^^^^^^^^

The Langmuir Exchange model reduces to simpler cases under specific conditions:

1. **Linear Exchange** (:math:`q_{\max,k}^i \rightarrow \infty`):
   
   .. math::
      J_{l \rightarrow k}^i \rightarrow k_{lk}^i \cdot c_{i,l} \cdot \frac{A_l}{A_k}

2. **Single Component** (:math:`N_c = 1`):
   
   .. math::
      J_{l \rightarrow k}^i = k_{lk}^i \cdot c_{i,l} \cdot q_{\max,k}^i \cdot \left(1 - \frac{c_{i,k}}{q_{\max,k}^i}\right) \cdot \frac{A_l}{A_k}


Physical Interpretation
^^^^^^^^^^^^^^^^^^^^^^^

The Langmuir Exchange model can represent various physical phenomena:

**Liquid-Liquid Chromatography**: The model describes partitioning between two immiscible liquid phases where the stationary phase has finite capacity for solute molecules.

**Biological Transport**: Saturable transport processes in biological systems, such as protein-mediated transport across membranes.

**Phase Separation**: Mass transfer between phases with limited solubility or capacity.

Implementation Notes
^^^^^^^^^^^^^^^^^^^^

The model is implemented within the MCT framework and inherits all transport properties (convection, dispersion) from the base MCT model. The exchange terms are added to the residual formulation and the corresponding Jacobian contributions are computed analytically for efficient solution.

The implementation ensures:

- Mass conservation across all channels
- Numerical stability through appropriate bounds checking
- Efficient computation through pre-calculated factors
- Compatibility with the existing MCT discretization schemes

Example Application
^^^^^^^^^^^^^^^^^^^

Consider a two-channel system representing aqueous and organic phases in liquid-liquid extraction:

- Channel 1: Aqueous phase with limited binding sites
- Channel 2: Organic phase with high capacity
- Component A: Hydrophilic, prefers aqueous phase
- Component B: Hydrophobic, prefers organic phase

The exchange matrix might be configured as:

.. math::
   K = \begin{pmatrix}
   0 & k_{12}^A & k_{12}^B \\
   k_{21}^A & 0 & k_{21}^B \\
   \end{pmatrix}

with :math:`k_{12}^A < k_{21}^A` (A prefers aqueous) and :math:`k_{12}^B > k_{21}^B` (B prefers organic), and appropriate capacity values reflecting the physical properties of each phase.

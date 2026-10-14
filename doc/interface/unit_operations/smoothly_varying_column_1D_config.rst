.. _smoothly_varying_column_1D_config:

Smoothly Varying Cross Section Column 1D
========================================

Group /input/model/unit_XXX - UNIT_TYPE - COLUMN_MODEL_1D
------------------------------------------------------------

This geometry describes an axial flow column whose cross section area varies smoothly along
the flow path and is prescribed by the user, see :ref:`MUOPGRMsmoothlyVarying`.
It is selected by ``GEOMETRY`` :math:`= \texttt{SMOOTHLY_VARYING}` and is only available for the
DG spatial discretization, since the area is given at the DG nodes.

Only the geometry fields differ from the axial flow cylinder; all remaining fields of the unit
(``NCOMP``, ``COL_POROSITY``, ``TOTAL_POROSITY``, ``NPARTYPE``, ``PAR_TYPE_VOLFRAC``,
``COL_DISPERSION``, ``INIT_C``, ``INIT_STATE``, the reaction and particle subgroups, and the
discretization group) are identical to :ref:`axial_flow_column_1D_config`.

``UNIT_TYPE``

   Specifies the type of unit operation model

   ================  ===========================================  =============
   **Type:** string  **Range:** :math:`\texttt{COLUMN_MODEL_1D}`  **Length:** 1
   ================  ===========================================  =============

``GEOMETRY``

   Column geometry

   ================  ============================================  =============
   **Type:** string  **Range:** :math:`\texttt{SMOOTHLY_VARYING}`  **Length:** 1
   ================  ============================================  =============

``CROSS_SECTIONAL_AREA_AT_NODES``

   Cross section area at every DG node, in the ordering of the bulk solution, i.e. element by
   element from the column inlet at :math:`x = 0` to :math:`x = \texttt{BED_LENGTH}`, with the
   nodes of each element in ascending order.
   Since neighbouring elements share their interface node, that node appears twice, once as the
   last entry of the left and once as the first entry of the right element, and the two entries
   must be equal.
   The area has to be strictly positive.

   **Unit:** :math:`\mathrm{m}^2`

   ================  ======================  ===============================================================
   **Type:** double  **Range:** :math:`> 0`  **Length:** :math:`(\texttt{POLYDEG} + 1) \cdot \texttt{NELEM}`
   ================  ======================  ===============================================================

   The nodes are the Legendre-Gauss-Lobatto nodes of each element, which are also written out by
   ``WRITE_COORDINATES``, so a convenient way to obtain them is to run the same discretization
   once with any other geometry and evaluate the intended area profile at the reported
   coordinates.

``BED_LENGTH``

   Column length. NOT optional for this geometry.

   **Unit:** :math:`\mathrm{m}`

   ================  ======================  =============
   **Type:** double  **Range:** :math:`> 0`  **Length:** 1
   ================  ======================  =============

``FORWARD_FLOW``

   Flow direction of the column. If set to 1, the flow runs from :math:`x = 0` towards
   :math:`x = \texttt{BED_LENGTH}`, i.e. the inlet is the end whose area is the first entry of
   ``CROSS_SECTIONAL_AREA_AT_NODES``.

   ==============  ==========================  ==================
   **Type:** bool  **Range:** :math:`\{0,1\}`  **Length:** NSEC
   ==============  ==========================  ==================

Group /input/model/unit_XXX/discretization - UNIT_TYPE - COLUMN_MODEL_1D
-------------------------------------------------------------------------

The discretization group is the one of :ref:`axial_flow_column_1D_config`, with two restrictions:

- ``SPATIAL_METHOD`` must be :math:`\texttt{DG}`. There is no Finite Volume variant of this
  geometry, because the area is prescribed at the DG nodes.
- ``USE_COLLOCATION_DG`` must be :math:`0`, which is the default for every geometry other than
  the axial flow cylinder.

Note that ``CROSS_SECTIONAL_AREA_AT_NODES`` depends on ``POLYDEG`` and ``NELEM``, so it has to be
recomputed whenever the discretization is refined.

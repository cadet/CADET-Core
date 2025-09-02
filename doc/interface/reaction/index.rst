.. _FFReaction:

Reaction models
===============

Reactions can take place on different phases in the model. Depending on the unit operation, reactions can be defined for the bulk phase, solid phase, particle phase, or as cross-phase reactions between liquid and solid phases.
In general the following phases are supported with reactions

- Bulk - main fluid phase in GRM and LRMP units
- Solid - solid phase of particles
- Pore - liquid within particles
- Liquid - main fluid phase in LRM and CSTR

The following table provides an overview of reaction phase support across different unit operations:


.. list-table:: Unit Operation Reaction Phase Support
    :header-rows: 1
    :widths: 40 20 20 20

    * - Unit Operation
      - Bulk/Liquid
      - Pore/Particle Liquid
      - Solid/Particle Solid
    * - GeneralRateModel (1D/2D/DG)
      - bulk
      - pore (particle liquid)
      - solid (particle solid)
    * - LumpedRateModelWithPores
      - bulk
      - liquid (particle liquid)
      - solid (particle solid)
    * - LumpedRateModelWithoutPores
      - liquid (particle liquid / bulk)
      - --
      - solid (particle solid)
    * - CSTR
      - liquid (particle liquid / bulk)
      - --
      - solid (particle solid)
    * - MCTM
      - liquid (bulk)
      - --
      - --


**Single-Phase Reaction Models** (for bulk/liquid, solid, and particle phases):

.. toctree::
    :maxdepth: 2

    mass_action_law
    michaelis_menten_kinetics


Cross Phase Reactions
---------------------

Cross-phase reactions enable chemical processes that occur at the interface between different phases, 
such as dissolution, precipitation, or catalytic reactions on solid surfaces. These reactions 
are different from single-phase reactions as they involve mass transfer between phases.
They are defined between a solid phase and a liquid phase i.e bulk/particle liquid and solid/particle solid phase.

**Cross-Phase Reaction Models**:

.. toctree::
    :maxdepth: 2

    mass_action_law_cross_phase

.. note::
   Cross-phase reaction models can **only** be used in the cross-phase formulation
   For reactions within a single phase, use the standard reaction models listed above.


Particle Reactions
------------------

If the reaction takes place in a phase associated with a particle, the reaction group is defined under the corresponding particle phase i.e. ``particle_type_XXX``.
Each particle can have its own reaction model configuration.

.. note::
   - For reactions within a single phase, use the standard reaction models listed above.
   - For reactions between a liquid and solid phase, use cross-phase reaction models.


Multiple Reaction Models per Phase
----------------------------------

Multiple reaction models can be defined for each phase of a unit operation.
The number of reaction models is specified by the parameter ``NREAC_PHASE`` (e.g., ``NREAC_BULK`` for the bulk phase).
If your reaction model is defined in a particle phase, the input group is defined under the corresponding particle type (e.g., ``particle_type_000/NREAC_SOLID``).


Step by Step Guide to Define Reactions
--------------------------------------

1. Choose the appropriate reaction model based on the kinetics of your system (e.g., Mass Action Law, Michaelis-Menten).
2. Check if your reaction occurs within a single phase or between phases (cross-phase).
3. Choose a phase for your reaction (bulk/liquid, solid, particle liquid, particle solid). **Note** that for particle reactions, the reaction model should be defined under the specific particle type.
4. Choose how many reaction types you want to define for the selected phase and set up the corresponding groups in the configuration file (e.g., ``NREAC_BULK``, ``NREAC_SOLID``, ``NREAC_CROSS_PHASE``).
  **Note** here that the input group depends on whether the reaction is defined in a particle phase or not. 
5. Define the reaction parameters according to the chosen model's specifications. Refer to the specific reaction model documentation for details on required parameters.


Externally dependent reaction models
------------------------------------

Some reaction models have a variant that can use external sources as specified :ref:`/input/model/external/<FFModelExternalSourceLinInterp>` (also see Section :ref:`dependence-on-external-function_react`).
For the sake of brevity, only the standard variant of those reaction models is specified below.
In order to obtain the format for the externally dependent variant, first replace the reaction model name ``XXX`` by ``EXT_XXX``.
Each parameter :math:`p` (except for stoichiometric and exponent matrices) depends on a (possibly distinct) external source in a polynomial way:

.. math::

    \begin{aligned}
        p(T) &= p_{\texttt{TTT}} T^3 + p_{\texttt{TT}} T^2 + p_{\texttt{T}} T + p.
    \end{aligned}

Thus, a parameter ``XXX_YYY`` of the standard reaction model variant is replaced by the four parameters ``EXT_XXX_YYY``, ``EXT_XXX_YYY_T``, ``EXT_XXX_YYY_TT``, and ``EXT_XXX_YYY_TTT``.
Since each parameter can depend on a different external source, the dataset ``EXTFUN`` (not listed in the standard variants below) should contain a vector of 0-based integer indices of the external source of each parameter.
The ordering of the parameters in ``EXTFUN`` is given by the ordering in the standard variant.
However, if only one index is passed in ``EXTFUN``, this external source is used for all parameters.

Note that parameter sensitivities with respect to column radius, column length, particle core radius, and particle radius may be wrong when using externally dependent reaction models.
This is caused by not taking into account the derivative of the external profile with respect to column position.

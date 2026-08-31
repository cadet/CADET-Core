// SPDX-License-Identifier: AGPL-3.0-or-later
// =================================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Affero General Public
//  License v3.0 (or, at your option, any later version).
// =================================================================================
 
 #include "model/particle/ParticleModel.hpp"
#include "model/ColumnModel1D.hpp"
#ifdef ENABLE_2D_MODELS
	#include "model/ColumnModel2D.hpp"
	#include "model/GeneralRateModel2D.hpp"
#endif
#include "model/GeneralRateModel.hpp"
#include "model/LumpedRateModelWithoutPores.hpp"
#include "model/LumpedRateModelWithoutPoresDG.hpp"
#include "model/LumpedRateModelWithPores.hpp"
#include "CompileTimeConfig.hpp"
#include "LoggingUtils.hpp"
#include "Logging.hpp"


namespace cadet
{
namespace model
{

#ifdef ENABLE_2D_MODELS
	IUnitOperation* select2DColumnUnitOperation(UnitOpIdx uoId, IParameterProvider& paramProvider)
	{
		std::string uoType = paramProvider.getString("UNIT_TYPE");

		paramProvider.pushScope("discretization");
		const std::string discName = paramProvider.getString("SPATIAL_METHOD");
		paramProvider.popScope(); // discretization

		if (discName == "DG")
			return new ColumnModel2D(uoId);
		else if (discName == "FV")
		{
			if (uoType == "GENERAL_RATE_MODEL_2D")
				return new GeneralRateModel2D(uoId);

			paramProvider.pushScope("particle_type_000");
			const bool filmDiffusion = paramProvider.getBool("HAS_FILM_DIFFUSION");
			const bool poreDiffusion = paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false;
			const bool surfaceDiffusion = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;
			paramProvider.popScope();

			const std::string particleType = ParticleModel(filmDiffusion, poreDiffusion, surfaceDiffusion).getParticleTransportType();

			if (particleType == "GENERAL_RATE_PARTICLE")
				return new GeneralRateModel2D(uoId);
			else
				LOG(Error) << "This particle Type (check HAS_FILM_DIFFUSION, HAS_PORE_DIFFUSION, HAS_SURFACE_DIFFUSION) is not implemented for FV discretization of the bulk phase, but was specified as such for unit " << uoId;
		}
		else
			LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;

		return nullptr;
	}
#endif

	/*
	 * @brief Selects the appropriate column unit operation based on the given parameters
	 * @detail We have different unit operation types based on: numerical discretization, particle model, column geometry and spatial resolution of the column.
	 *		   We have a total of 3 different column geometries (axial, radial, frustum) and 2 different bulk discretization methods (FV, DG) and 3 different particle models (EPM, LRMP, GRM).
	 *		   Some combinations yield the same unit operation due to various generalizations.
	 */
	IUnitOperation* selectColumnUnitOperation(UnitOpIdx uoId, IParameterProvider& paramProvider)
	{
		std::string uoType = paramProvider.getString("UNIT_TYPE");

#ifdef ENABLE_2D_MODELS
		if (uoType.find("_2D") != std::string::npos)
			return select2DColumnUnitOperation(uoId, paramProvider);
#endif

		IUnitOperation* model = nullptr;

		std::string colGeometry = paramProvider.getString("GEOMETRY");

		if (!(colGeometry == "AXIAL_FLOW_CYLINDER" || colGeometry == "RADIAL_FLOW_CYLINDER_SHELL" || colGeometry == "AXIAL_FLOW_FRUSTUM"))
		{
			throw InvalidParameterException("Unsupported column geometry " + colGeometry + " was specified for unit " + std::to_string(uoId));
		}

		paramProvider.pushScope("discretization");
		
		const std::string discName = paramProvider.getString("SPATIAL_METHOD");

		// ARROW_HEAD_OPTIMIZATION defaults to true for FV, preserving the existing block-structured
		// (arrow-head) Jacobian solver in the dedicated FV unit operation classes.
		// Set it to false to route FV through ColumnModel1D (global sparse solver).
		const bool arrowHeadOpt = paramProvider.exists("FV_ARROW_HEAD_OPTIMIZATION")
			? paramProvider.getBool("FV_ARROW_HEAD_OPTIMIZATION")
			: discName == "FV";

		if (arrowHeadOpt && discName != "FV")
		{
			throw InvalidParameterException("FV_ARROW_HEAD_OPTIMIZATION is only available for FV discretization but " + discName + " was specifiedfor unit " + std::to_string(uoId));
		}

		bool axialCollocationDG = (colGeometry == "AXIAL_FLOW_CYLINDER");
		if (discName == "DG")
		{
			// collocation DG default for axial flow, not implemented for other geometries
			if (paramProvider.exists("USE_COLLOCATION_DG"))
				axialCollocationDG = paramProvider.exists("USE_COLLOCATION_DG") ? paramProvider.getBool("USE_COLLOCATION_DG") : colGeometry == "AXIAL_FLOW_CYLINDER";

			if (colGeometry != "AXIAL_FLOW_CYLINDER" && axialCollocationDG)
			{
				throw InvalidParameterException("USE_COLLOCATION_DG is only available for AXIAL_FLOW_CYLINDER geometry but " + colGeometry + " was specified for unit " + std::to_string(uoId));
			}
		}

		paramProvider.popScope(); // discretization

		if (uoType == "COLUMN_MODEL_1D")
		{
			if (paramProvider.getInt("NPARTYPE") > 0)
			{
				paramProvider.pushScope("particle_type_000");

				bool filmDiffusion = true;

				if (paramProvider.exists("HAS_FILM_DIFFUSION"))
					filmDiffusion = paramProvider.getBool("HAS_FILM_DIFFUSION");

				const bool poreDiffusion = paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false;
				const bool surfaceDiffusion = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;

				const std::string particleType = ParticleModel(filmDiffusion, poreDiffusion, surfaceDiffusion).getParticleTransportType();

				if (particleType == "EQUILIBRIUM_PARTICLE")
				{
					if (discName == "DG")
						model = axialCollocationDG ? createAxialLRMDG(uoId) : createVariableCrossSectionLRMDG(uoId);
					else if (discName == "FV")
					{
						if (colGeometry == "AXIAL_FLOW_CYLINDER")
							model = createAxialFVLRM(uoId);
						else if (colGeometry == "RADIAL_FLOW_CYLINDER_SHELL")
							model = createRadialFVLRM(uoId);
						else if (colGeometry == "AXIAL_FLOW_FRUSTUM")
							model = createFrustumFVLRM(uoId);
					}
					else
						LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;
				}
				else if (discName == "FV" && arrowHeadOpt)
				{
					if (particleType == "HOMOGENEOUS_PARTICLE")
					{
						if (colGeometry == "AXIAL_FLOW_CYLINDER")
							model = createAxialFVLRMP(uoId);
						else if (colGeometry == "RADIAL_FLOW_CYLINDER_SHELL")
							model = createRadialFVLRMP(uoId);
						else if (colGeometry == "AXIAL_FLOW_FRUSTUM")
							model = createFrustumFVLRMP(uoId);
					}
					else if (particleType == "GENERAL_RATE_PARTICLE")
					{
						if (!arrowHeadOpt)
							LOG(Info) << "FV_ARROW_HEAD_OPTIMIZATION is set to false for a general rate model unit operation, probably resulting in a less efficient computation, we recommend using the arrow head optimization for the combination of FV and GRM";
						if (colGeometry == "AXIAL_FLOW_CYLINDER")
							model = createAxialFVGRM(uoId);
						else if (colGeometry == "RADIAL_FLOW_CYLINDER_SHELL")
							model = createRadialFVGRM(uoId);
						else if (colGeometry == "AXIAL_FLOW_FRUSTUM")
							model = createFrustumFVGRM(uoId);
					}
				}
				else if (discName == "FV" && !arrowHeadOpt)
				{
					if (colGeometry == "AXIAL_FLOW_CYLINDER")
						model = createAxialCol1DFV(uoId);
					else if (colGeometry == "RADIAL_FLOW_CYLINDER_SHELL")
						model = createRadialCol1DFV(uoId);
					else if (colGeometry == "AXIAL_FLOW_FRUSTUM")
						model = createFrustumCol1DFV(uoId);
				}
				else if (discName == "DG")
					model = axialCollocationDG ? createAxialCol1DDG(uoId) : createVariableCrossSectionCol1DDG(uoId);

				paramProvider.popScope();
			}
			else if (discName == "DG")
			{
				model = axialCollocationDG ? createAxialCol1DDG(uoId) : createVariableCrossSectionCol1DDG(uoId);
			}
			else if (discName == "FV")
			{
				if (arrowHeadOpt)
				{
					if (colGeometry == "AXIAL_FLOW_CYLINDER")
						model = createAxialFVLRM(uoId);
					else if (colGeometry == "RADIAL_FLOW_CYLINDER_SHELL")
						model = createRadialFVLRM(uoId);
					else if (colGeometry == "AXIAL_FLOW_FRUSTUM")
						model = createFrustumFVLRM(uoId);
				}
				else
				{
					if (colGeometry == "AXIAL_FLOW_CYLINDER")
						model = createAxialCol1DFV(uoId);
					else if (colGeometry == "RADIAL_FLOW_CYLINDER_SHELL")
						model = createRadialCol1DFV(uoId);
					else if (colGeometry == "AXIAL_FLOW_FRUSTUM")
						model = createFrustumCol1DFV(uoId);
				}
			}
			else
			{
				throw InvalidParameterException("Unknown bulk discretization type " + discName + " for unit " + std::to_string(uoId));
			}
		}
		else // feature: we support legacy unit operation names
		{
			if (paramProvider.getInt("NPARTYPE") < 1)
				throw InvalidParameterException("NPARTYPE must be at least 1 for unit operation " + uoType);

			if (discName == "DG")
			{
				if (uoType == "LUMPED_RATE_MODEL_WITHOUT_PORES")
					model = axialCollocationDG ? createAxialLRMDG(uoId) : createVariableCrossSectionLRMDG(uoId);
				else if (uoType == "LUMPED_RATE_MODEL_WITH_PORES" || uoType == "GENERAL_RATE_MODEL")
					model = axialCollocationDG ? createAxialCol1DDG(uoId) : createVariableCrossSectionCol1DDG(uoId);
			}
			else if (discName == "FV")
			{
				if(colGeometry == "AXIAL_FLOW_CYLINDER")
				{
					if (uoType == "LUMPED_RATE_MODEL_WITHOUT_PORES")
						model = arrowHeadOpt ? createAxialFVLRM(uoId) : createAxialCol1DFV(uoId);
					else if (uoType == "LUMPED_RATE_MODEL_WITH_PORES")
						model = arrowHeadOpt ? createAxialFVLRMP(uoId) : createAxialCol1DFV(uoId);
					else if (uoType == "GENERAL_RATE_MODEL")
						model = arrowHeadOpt ? createAxialFVGRM(uoId) : createAxialCol1DFV(uoId);
				}
				else if(colGeometry == "RADIAL_FLOW_CYLINDER_SHELL")
				{
					if (uoType == "LUMPED_RATE_MODEL_WITHOUT_PORES")
						model = arrowHeadOpt ? createRadialFVLRM(uoId) : createRadialCol1DFV(uoId);
					else if (uoType == "LUMPED_RATE_MODEL_WITH_PORES")
						model = arrowHeadOpt ? createRadialFVLRMP(uoId) : createRadialCol1DFV(uoId);
					else if (uoType == "GENERAL_RATE_MODEL")
						model = arrowHeadOpt ? createRadialFVGRM(uoId) : createRadialCol1DFV(uoId);
				}
				else if(colGeometry == "AXIAL_FLOW_FRUSTUM")
				{
					if (uoType == "LUMPED_RATE_MODEL_WITHOUT_PORES")
						model = arrowHeadOpt ? createFrustumFVLRM(uoId) : createFrustumCol1DFV(uoId);
					else if (uoType == "LUMPED_RATE_MODEL_WITH_PORES")
						model = arrowHeadOpt ? createFrustumFVLRMP(uoId) : createFrustumCol1DFV(uoId);
					else if (uoType == "GENERAL_RATE_MODEL")
						model = arrowHeadOpt ? createFrustumFVGRM(uoId) : createFrustumCol1DFV(uoId);
				}
			}
			else
				LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;
		}

		return model;
	}

	void registerColumnModel(std::unordered_map<std::string, std::function<IUnitOperation* (UnitOpIdx, IParameterProvider&)>>& models)
	{
		models[ColumnModel2D::identifier()] = selectColumnUnitOperation;
		models[GeneralRateModel2D::identifier()] = selectColumnUnitOperation;

		models["COLUMN_MODEL_1D"] = selectColumnUnitOperation;
		models["GENERAL_RATE_MODEL"] = selectColumnUnitOperation;
		models["GRM"] = selectColumnUnitOperation;
		models["LUMPED_RATE_MODEL_WITH_PORES"] = selectColumnUnitOperation;
		models["LRMP"] = selectColumnUnitOperation;
		models["LUMPED_RATE_MODEL_WITHOUT_PORES"] = selectColumnUnitOperation;
		models["LRM"] = selectColumnUnitOperation;
		models["DISPERSIVE_PLUG_FLOW_REACTOR"] = selectColumnUnitOperation;
		models["DPFR"] = selectColumnUnitOperation;

		typedef ColumnModel1D<parts::VariableCrossSectionConvectionDispersionOperatorBaseDG> VarCrossSectionCol1DDG;
		typedef ColumnModel1D<parts::AxialConvectionDispersionOperatorBaseCollocationDG> AxialCol1DCollocationDG;
		typedef ColumnModel1D<parts::AxialConvectionDispersionOperatorBaseFV> AxialCol1DCollocationFV;
		typedef ColumnModel1D<parts::RadialConvectionDispersionOperatorBaseFV> RadialCol1DCollocationFV;
		typedef ColumnModel1D<parts::FrustumConvectionDispersionOperatorBaseFV> FrustumCol1DCollocationFV;
		
		models[VarCrossSectionCol1DDG::identifier()] = selectColumnUnitOperation;
		models[AxialCol1DCollocationDG::identifier()] = selectColumnUnitOperation;
		models[AxialCol1DCollocationFV::identifier()] = selectColumnUnitOperation;
		models[RadialCol1DCollocationFV::identifier()] = selectColumnUnitOperation;
		models[FrustumCol1DCollocationFV::identifier()] = selectColumnUnitOperation;

		typedef LumpedRateModelWithoutPoresDG<parts::VariableCrossSectionConvectionDispersionOperatorBaseDG> VarCrossSectionLRMDG;
		typedef LumpedRateModelWithoutPoresDG<parts::AxialConvectionDispersionOperatorBaseCollocationDG> AxialLRMCOLLOCATIONDG;

		models[AxialLRMCOLLOCATIONDG::identifier()] = selectColumnUnitOperation;
		models[VarCrossSectionLRMDG::identifier()] = selectColumnUnitOperation;

		typedef GeneralRateModel<parts::AxialConvectionDispersionOperatorFV> AxialGRMFV;
		typedef GeneralRateModel<parts::RadialConvectionDispersionOperatorFV> RadialGRMFV;
		typedef GeneralRateModel<parts::FrustumConvectionDispersionOperatorFV> FrustumGRMFV;

		models[AxialGRMFV::identifier()] = selectColumnUnitOperation;
		models[RadialGRMFV::identifier()] = selectColumnUnitOperation;
		models[FrustumGRMFV::identifier()] = selectColumnUnitOperation;

		typedef LumpedRateModelWithPores<parts::AxialConvectionDispersionOperatorFV> AxialLRMPFV;
		typedef LumpedRateModelWithPores<parts::RadialConvectionDispersionOperatorFV> RadialLRMPFV;

		models[AxialLRMPFV::identifier()] = selectColumnUnitOperation;
		models[RadialLRMPFV::identifier()] = selectColumnUnitOperation;

		typedef LumpedRateModelWithoutPores<parts::AxialConvectionDispersionOperatorBaseFV> AxialLRMFV;
		typedef LumpedRateModelWithoutPores<parts::RadialConvectionDispersionOperatorBaseFV> RadialLRMFV;

		models[AxialLRMFV::identifier()] = selectColumnUnitOperation;
		models[RadialLRMFV::identifier()] = selectColumnUnitOperation;
	}

}  // namespace model

}  // namespace cadet
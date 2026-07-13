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

	IUnitOperation* selectAxialFlowColumnUnitOperation(UnitOpIdx uoId, IParameterProvider& paramProvider)
	{
		IUnitOperation* model = nullptr;

		std::string uoType = paramProvider.getString("UNIT_TYPE");

		paramProvider.pushScope("discretization");
		
		const std::string discName = paramProvider.getString("SPATIAL_METHOD");

		// ARROW_HEAD_OPTIMIZATION defaults to true for FV, preserving the existing block-structured
		// (arrow-head) Jacobian solver in the dedicated FV unit operation classes.
		// Set it to false to route FV through ColumnModel1D (unified path, global sparse solver).
		const bool arrowHeadOpt = paramProvider.exists("FV_ARROW_HEAD_OPTIMIZATION")
			? paramProvider.getBool("FV_ARROW_HEAD_OPTIMIZATION")
			: discName == "FV";

		if (discName == "FV")
		{
			if (!arrowHeadOpt)
				LOG(Info) << "FV_ARROW_HEAD_OPTIMIZATION is set to false, possibly resulting in a less efficient computation";
		}
		else if (arrowHeadOpt)
		{
			throw InvalidParameterException("FV_ARROW_HEAD_OPTIMIZATION is only available for FV discretization but " + discName + " was specified for " + uoType);
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
						model = createAxialLRMDG(uoId) ;
					else if (discName == "FV")
						model = arrowHeadOpt ? createAxialFVLRM(uoId) : createAxialCol1DFV(uoId);
					else
						LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;
				}
				else if (discName == "FV" && arrowHeadOpt)
				{
					if (particleType == "HOMOGENEOUS_PARTICLE")
						model = createAxialFVLRMP(uoId);
					else if (particleType == "GENERAL_RATE_PARTICLE")
						model = createAxialFVGRM(uoId);
				}
				else if (discName == "FV" && !arrowHeadOpt)
					model = createAxialCol1DFV(uoId);
				else if (discName == "DG")
					model = createAxialCol1DDG(uoId);

				paramProvider.popScope();
			}
			else if (discName == "DG")
			{
				model = createAxialCol1DDG(uoId);
			}
			else if (discName == "FV")
			{
				if(arrowHeadOpt)
					model = createAxialFVLRM(uoId);
				else
				model = createAxialCol1DFV(uoId);
			}
			else
			{
				throw InvalidParameterException("Unknown bulk discretization type " + discName + " for unit " + std::to_string(uoId));
			}
		}
#ifdef ENABLE_2D_MODELS
		else if (uoType.find("_2D") != std::string::npos)
		{
			if (discName == "DG")
				model = new ColumnModel2D(uoId);
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
					model = new GeneralRateModel2D(uoId);
				else
					LOG(Error) << "This particle Type (check HAS_FILM_DIFFUSION, HAS_PORE_DIFFUSION, HAS_SURFACE_DIFFUSION) is not implemented for FV discretization of the bulk phase, but was specified as such for unit " << uoId;
			}
			else
				LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;
		}
#endif
		else
		{
			if (paramProvider.getInt("NPARTYPE") < 1)
				throw InvalidParameterException("NPARTYPE must be at least 1 for unit operation " + uoType);

			if (discName == "DG")
			{
				if (uoType == "LUMPED_RATE_MODEL_WITHOUT_PORES")
					model = createAxialLRMDG(uoId);
				else if (uoType == "LUMPED_RATE_MODEL_WITH_PORES" || uoType == "GENERAL_RATE_MODEL")
					model = createAxialCol1DDG(uoId);
			}
			else if (discName == "FV")
			{
				if (uoType == "LUMPED_RATE_MODEL_WITHOUT_PORES")
					model = arrowHeadOpt ? createAxialFVLRM(uoId) : createAxialCol1DFV(uoId);
				else if (uoType == "LUMPED_RATE_MODEL_WITH_PORES")
					model = arrowHeadOpt ? createAxialFVLRMP(uoId) : createAxialCol1DFV(uoId);
				else if (uoType == "GENERAL_RATE_MODEL")
					model = arrowHeadOpt ? createAxialFVGRM(uoId) : createAxialCol1DFV(uoId);
			}
			else
				LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;
		}

		return model;
	}

	IUnitOperation* selectRadialFlowColumnUnitOperation(UnitOpIdx uoId, IParameterProvider& paramProvider)
	{
		IUnitOperation* model = nullptr;

		std::string uoType = paramProvider.getString("UNIT_TYPE");

		paramProvider.pushScope("discretization");

		const std::string discName = paramProvider.getString("SPATIAL_METHOD");

		// ARROW_HEAD_OPTIMIZATION defaults to true, preserving the existing block-structured
		// (arrow-head) Jacobian solver in the dedicated FV unit operation classes.
		// Set it to false to route FV through ColumnModel1D (unified path, global sparse solver).
		const bool arrowHeadOpt = paramProvider.exists("FV_ARROW_HEAD_OPTIMIZATION")
			? paramProvider.getBool("FV_ARROW_HEAD_OPTIMIZATION")
			: discName == "FV";

		if (discName == "FV")
		{
			if (!arrowHeadOpt)
				LOG(Info) << "FV_ARROW_HEAD_OPTIMIZATION is set to false, possibly resulting in a less efficient computation";
		}
		else if (arrowHeadOpt)
		{
			throw InvalidParameterException("FV_ARROW_HEAD_OPTIMIZATION is only available for FV discretization but " + discName + " was specified for " + uoType);
		}

		paramProvider.popScope(); // discretization

		if (uoType == "RADIAL_COLUMN_MODEL_1D")
		{
			if (paramProvider.getInt("NPARTYPE") > 0)
			{
				paramProvider.pushScope("particle_type_000");

				const bool filmDiffusion = paramProvider.getBool("HAS_FILM_DIFFUSION");
				const bool poreDiffusion = paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false;
				const bool surfaceDiffusion = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;

				const std::string particleType = ParticleModel(filmDiffusion, poreDiffusion, surfaceDiffusion).getParticleTransportType();

				if (discName == "DG")
				{
					if (particleType == "EQUILIBRIUM_PARTICLE")
						model = createRadialLRMDG(uoId);
					else
						model = createRadialCol1DDG(uoId);
				}
				else if (discName == "FV")
				{
					if (particleType == "EQUILIBRIUM_PARTICLE")
						model = arrowHeadOpt ? createRadialFVLRM(uoId) : createRadialCol1DFV(uoId);
					else if (particleType == "HOMOGENEOUS_PARTICLE")
						model = arrowHeadOpt ? createRadialFVLRMP(uoId) : createRadialCol1DFV(uoId);
					else if (particleType == "GENERAL_RATE_PARTICLE")
						model = arrowHeadOpt ? createRadialFVGRM(uoId) : createRadialCol1DFV(uoId);
					else
						LOG(Error) << "Unknown particle type " << particleType << " for unit " << uoId;
				}
				else
					LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;

				paramProvider.popScope();
			}
			else
			{
				if (discName == "DG")
					model = createRadialCol1DDG(uoId);
				else
					model = arrowHeadOpt ? createRadialFVLRM(uoId) : createRadialCol1DFV(uoId);
			}
		}
		else
		{
			if (discName == "DG")
			{
				if (uoType == "RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES")
					model = createRadialLRMDG(uoId);
				else if (uoType == "RADIAL_LUMPED_RATE_MODEL_WITH_PORES" || uoType == "RADIAL_GENERAL_RATE_MODEL")
					model = createRadialCol1DDG(uoId);
				else
					LOG(Error) << "Radial DG only supports LRM, LRMP, and GRM currently for unit " << uoId;
			}
			else if (discName == "FV")
			{
				if (uoType == "RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES")
					model = createRadialFVLRM(uoId);
				else if (uoType == "RADIAL_LUMPED_RATE_MODEL_WITH_PORES")
					model = createRadialFVLRMP(uoId);
				else if (uoType == "RADIAL_GENERAL_RATE_MODEL")
					model = createRadialFVGRM(uoId);
			}
			else
				LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;
		}

		return model;
	}

	IUnitOperation* selectFrustumFlowColumnUnitOperation(UnitOpIdx uoId, IParameterProvider& paramProvider)
	{
		IUnitOperation* model = nullptr;

		std::string uoType = paramProvider.getString("UNIT_TYPE");

		paramProvider.pushScope("discretization");

		const std::string discName = paramProvider.getString("SPATIAL_METHOD");

		// ARROW_HEAD_OPTIMIZATION defaults to true, preserving the existing block-structured
		// (arrow-head) Jacobian solver in the dedicated FV unit operation classes.
		// Set it to false to route FV through ColumnModel1D (unified path, global sparse solver).
		const bool arrowHeadOpt = paramProvider.exists("FV_ARROW_HEAD_OPTIMIZATION")
			? paramProvider.getBool("FV_ARROW_HEAD_OPTIMIZATION")
			: discName == "FV";

		if (discName != "FV")
			throw InvalidParameterException("Only FV discretization supported for Frustum geometry but " + discName + " was asked for unit " + std::to_string(uoId));
		else if (!arrowHeadOpt)
			throw InvalidParameterException("FV_ARROW_HEAD_OPTIMIZATION was set to false for unit " + std::to_string(uoId) + " but only arrow head implementation is available");

		paramProvider.popScope(); // discretization

		if (uoType == "FRUSTUM_COLUMN_MODEL_1D")
		{
			if (paramProvider.exists("particle_type_000"))
			{
				if (discName == "DG")
					LOG(Error) << "Frustum flow not implemented for DG discretization yet, was called for unit " << uoId;
				else if (discName != "FV")
					LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;

				paramProvider.pushScope("particle_type_000");

				const bool filmDiffusion = paramProvider.getBool("HAS_FILM_DIFFUSION");
				const bool poreDiffusion = paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false;
				const bool surfaceDiffusion = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;

				const std::string particleType = ParticleModel(filmDiffusion, poreDiffusion, surfaceDiffusion).getParticleTransportType();

				if (particleType == "EQUILIBRIUM_PARTICLE")
					model = createFrustumFVLRM(uoId);
				if (particleType == "HOMOGENEOUS_PARTICLE")
					model = createFrustumFVLRMP(uoId);
				else if (particleType == "GENERAL_RATE_PARTICLE")
					model = createFrustumFVGRM(uoId);

				paramProvider.popScope();
			}
			else
				model = createFrustumFVLRM(uoId); // LRMP used for npartype = 0
		}
		else
		{
			if (discName == "DG")
				LOG(Error) << "Frustum flow not implemented for DG discretization yet, was called for unit " << uoId;
			else if (discName == "FV")
			{
				if (uoType == "FRUSTUM_LUMPED_RATE_MODEL_WITHOUT_PORES")
					model = createFrustumFVLRM(uoId);
				else if (uoType == "FRUSTUM_LUMPED_RATE_MODEL_WITH_PORES")
					model = createFrustumFVLRMP(uoId);
				else if (uoType == "FRUSTUM_GENERAL_RATE_MODEL")
					model = createFrustumFVGRM(uoId);
			}
			else
				LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;

		}

		return model;
	}

	void registerColumnModel(std::unordered_map<std::string, std::function<IUnitOperation* (UnitOpIdx, IParameterProvider&)>>& models)
	{
		models[ColumnModel1D<>::identifier()] = selectAxialFlowColumnUnitOperation;
		models["COLUMN_MODEL_1D"] = selectAxialFlowColumnUnitOperation;
		models["RADIAL_COLUMN_MODEL_1D"] = selectRadialFlowColumnUnitOperation;
		models["FRUSTUM_COLUMN_MODEL_1D"] = selectFrustumFlowColumnUnitOperation;

		models[ColumnModel2D::identifier()] = selectAxialFlowColumnUnitOperation;
		models["COLUMN_MODEL_2D"] = selectAxialFlowColumnUnitOperation;

		typedef LumpedRateModelWithoutPoresDG<parts::RadialConvectionDispersionOperatorBaseDG> RadialLRMDG;

		models[LumpedRateModelWithoutPoresDG<>::identifier()] = selectAxialFlowColumnUnitOperation;
		models["LRM_DG"] = selectAxialFlowColumnUnitOperation;

		models[RadialLRMDG::identifier()] = selectRadialFlowColumnUnitOperation;
		models["RLRM_DG"] = selectRadialFlowColumnUnitOperation;

		typedef GeneralRateModel<parts::AxialConvectionDispersionOperator> AxialGRM;
		typedef GeneralRateModel<parts::RadialConvectionDispersionOperator> RadialGRM;
		typedef GeneralRateModel<parts::FrustumConvectionDispersionOperator> FrustumGRM;

		models[AxialGRM::identifier()] = selectAxialFlowColumnUnitOperation;
		models["GRM"] = selectAxialFlowColumnUnitOperation;

		models[RadialGRM::identifier()] = selectRadialFlowColumnUnitOperation;
		models["RGRM"] = selectRadialFlowColumnUnitOperation;

		models[FrustumGRM::identifier()] = selectFrustumFlowColumnUnitOperation;
		models["FGRM"] = selectFrustumFlowColumnUnitOperation;

		models[GeneralRateModel2D::identifier()] = selectAxialFlowColumnUnitOperation;
		models["GRM2D"] = selectAxialFlowColumnUnitOperation;

		typedef LumpedRateModelWithPores<parts::AxialConvectionDispersionOperator> AxialLRMP;
		typedef LumpedRateModelWithPores<parts::RadialConvectionDispersionOperator> RadialLRMP;

		models[AxialLRMP::identifier()] = selectAxialFlowColumnUnitOperation;
		models["LRMP"] = selectAxialFlowColumnUnitOperation;

		models[RadialLRMP::identifier()] = selectRadialFlowColumnUnitOperation;
		models["RLRMP"] = selectRadialFlowColumnUnitOperation;

		typedef LumpedRateModelWithoutPores<parts::AxialConvectionDispersionOperatorBaseFV> AxialLRM;
		typedef LumpedRateModelWithoutPores<parts::RadialConvectionDispersionOperatorBaseFV> RadialLRM;

		models[AxialLRM::identifier()] = selectAxialFlowColumnUnitOperation;
		models["LRM"] = selectAxialFlowColumnUnitOperation;
		models["DPFR"] = selectAxialFlowColumnUnitOperation;

		models[RadialLRM::identifier()] = selectRadialFlowColumnUnitOperation;
		models["RLRM"] = selectRadialFlowColumnUnitOperation;
	}

}  // namespace model

}  // namespace cadet
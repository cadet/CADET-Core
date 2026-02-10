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

		if (uoType == "COLUMN_MODEL_1D")
		{
			if (paramProvider.exists("particle_type_000"))
			{
				paramProvider.pushScope("discretization");
				const std::string discName = paramProvider.getString("SPATIAL_METHOD");
				paramProvider.popScope();

				paramProvider.pushScope("particle_type_000");

				bool filmDiffusion = true;

				if (paramProvider.exists("HAS_FILM_DIFFUSION"))
					filmDiffusion = paramProvider.getBool("HAS_FILM_DIFFUSION");

				const bool poreDiffusion = paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false;
				const bool surfaceDiffusion = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;

				const std::string particleType = ParticleModel(filmDiffusion, poreDiffusion, surfaceDiffusion).getParticleTransportType();

				if(particleType == "EQUILIBRIUM_PARTICLE")
				{
					if (discName == "DG")
						model = new LumpedRateModelWithoutPoresDG(uoId);
					else if (discName == "FV")
						model = createAxialFVLRM(uoId);
					else
						LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;
				}
				else if (discName == "FV")
				{
					if(particleType == "HOMOGENEOUS_PARTICLE")
						model = createAxialFVLRMP(uoId);
					else if (particleType == "GENERAL_RATE_PARTICLE")
						model = createAxialFVGRM(uoId);
				}
				else
					model = new ColumnModel1D(uoId);

				paramProvider.popScope();
			}
			else
				model = new ColumnModel1D(uoId);
		}
#ifdef ENABLE_2D_MODELS
		else if (uoType.find("_2D") != std::string::npos)
		{
			paramProvider.pushScope("discretization");
			const std::string discName = paramProvider.getString("SPATIAL_METHOD");
			paramProvider.popScope();

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
			paramProvider.pushScope("discretization");
			const std::string discName = paramProvider.getString("SPATIAL_METHOD");
			paramProvider.popScope();

			if (discName == "DG")
			{
				if (uoType == "LUMPED_RATE_MODEL_WITHOUT_PORES")
					model = new LumpedRateModelWithoutPoresDG(uoId);
				else if (uoType == "LUMPED_RATE_MODEL_WITH_PORES" || uoType == "GENERAL_RATE_MODEL")
					model = new ColumnModel1D(uoId);
			}
			else if (discName == "FV")
			{
				if (uoType == "LUMPED_RATE_MODEL_WITHOUT_PORES")
					model = createAxialFVLRM(uoId);
				else if (uoType == "LUMPED_RATE_MODEL_WITH_PORES")
					model = createAxialFVLRMP(uoId);
				else if (uoType == "GENERAL_RATE_MODEL")
					model = createAxialFVGRM(uoId);
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

		if (uoType == "RADIAL_COLUMN_MODEL_1D")
		{
			if (paramProvider.exists("particle_type_000"))
			{
				paramProvider.pushScope("discretization");
				const std::string discName = paramProvider.getString("SPATIAL_METHOD");
				paramProvider.popScope();

				if (discName == "DG")
					LOG(Error) << "Radial flow not implemented for DG discretization yet, was called for unit " << uoId;
				else if (discName != "FV")
					LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;

				paramProvider.pushScope("particle_type_000");

				const bool filmDiffusion = paramProvider.getBool("HAS_FILM_DIFFUSION");
				const bool poreDiffusion = paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false;
				const bool surfaceDiffusion = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;

				const std::string particleType = ParticleModel(filmDiffusion, poreDiffusion, surfaceDiffusion).getParticleTransportType();

				if (particleType == "EQUILIBRIUM_PARTICLE")
					model = createRadialFVLRM(uoId);
				if (particleType == "HOMOGENEOUS_PARTICLE")
					model = createRadialFVLRMP(uoId);
				else if (particleType == "GENERAL_RATE_PARTICLE")
					model = createRadialFVGRM(uoId);

				paramProvider.popScope();
			}
			else
				model = createRadialFVLRMP(uoId); // LRMP used for npartype = 0
		}
		else
		{
			paramProvider.pushScope("discretization");
			const std::string discName = paramProvider.getString("SPATIAL_METHOD");

			if (discName == "DG")
				LOG(Error) << "Radial flow not implemented for DG discretization yet, was called for unit " << uoId;
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

			paramProvider.popScope();
		}

		return model;
	}

	IUnitOperation* selectFrustumFlowColumnUnitOperation(UnitOpIdx uoId, IParameterProvider& paramProvider)
	{
		IUnitOperation* model = nullptr;

		std::string uoType = paramProvider.getString("UNIT_TYPE");

		if (uoType == "FRUSTUM_COLUMN_MODEL_1D")
		{
			if (paramProvider.exists("particle_type_000"))
			{
				paramProvider.pushScope("discretization");
				const std::string discName = paramProvider.getString("SPATIAL_METHOD");
				paramProvider.popScope();

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
					LOG(Error) << "Frustum flow LRM not implemented, please use GRM, was called for unit " << uoId;
				if (particleType == "HOMOGENEOUS_PARTICLE")
					LOG(Error) << "Frustum flow LRMP not implemented, please use GRM, was called for unit " << uoId;
				else if (particleType == "GENERAL_RATE_PARTICLE")
					model = createFrustumFVGRM(uoId);

				paramProvider.popScope();
			}
			else
				model = createRadialFVLRMP(uoId); // LRMP used for npartype = 0
		}
		else
		{
			paramProvider.pushScope("discretization");
			const std::string discName = paramProvider.getString("SPATIAL_METHOD");

			if (discName == "DG")
				LOG(Error) << "Frustum flow not implemented for DG discretization yet, was called for unit " << uoId;
			else if (discName == "FV")
			{
				if (uoType == "FRUSTUM_LUMPED_RATE_MODEL_WITHOUT_PORES")
					LOG(Error) << "Frustum flow LRM not implemented, please use GRM, was called for unit " << uoId;
				else if (uoType == "FRUSTUM_LUMPED_RATE_MODEL_WITH_PORES")
					LOG(Error) << "Frustum flow LRMP not implemented, please use GRM, was called for unit " << uoId;
				else if (uoType == "FRUSTUM_GENERAL_RATE_MODEL")
					model = createFrustumFVGRM(uoId);
			}
			else
				LOG(Error) << "Unknown bulk discretization type " << discName << " for unit " << uoId;

			paramProvider.popScope();
		}

		return model;
	}

	void registerColumnModel(std::unordered_map<std::string, std::function<IUnitOperation* (UnitOpIdx, IParameterProvider&)>>& models)
	{
		models[ColumnModel1D::identifier()] = selectAxialFlowColumnUnitOperation;
		models["COLUMN_MODEL_1D"] = selectAxialFlowColumnUnitOperation;
		models["RADIAL_COLUMN_MODEL_1D"] = selectRadialFlowColumnUnitOperation;
		models["FRUSTUM_COLUMN_MODEL_1D"] = selectFrustumFlowColumnUnitOperation;

		models[ColumnModel2D::identifier()] = selectAxialFlowColumnUnitOperation;
		models["COLUMN_MODEL_2D"] = selectAxialFlowColumnUnitOperation;

		models[LumpedRateModelWithoutPoresDG::identifier()] = selectAxialFlowColumnUnitOperation;
		models["LRM_DG"] = selectAxialFlowColumnUnitOperation;

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

		typedef LumpedRateModelWithoutPores<parts::AxialConvectionDispersionOperatorBase> AxialLRM;
		typedef LumpedRateModelWithoutPores<parts::RadialConvectionDispersionOperatorBase> RadialLRM;

		models[AxialLRM::identifier()] = selectAxialFlowColumnUnitOperation;
		models["LRM"] = selectAxialFlowColumnUnitOperation;
		models["DPFR"] = selectAxialFlowColumnUnitOperation;

		models[RadialLRM::identifier()] = selectRadialFlowColumnUnitOperation;
		models["RLRM"] = selectRadialFlowColumnUnitOperation;
	}

}  // namespace model

}  // namespace cadet
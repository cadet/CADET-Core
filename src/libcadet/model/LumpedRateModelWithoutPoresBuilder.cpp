#include "model/LumpedRateModelWithoutPores.hpp"
#include "CompileTimeConfig.hpp"
#ifdef ENABLE_DG
	#include "model/LumpedRateModelWithoutPoresDG.hpp"
#endif
#include "LoggingUtils.hpp"
#include "Logging.hpp"


namespace cadet
{
namespace model
{

	IUnitOperation* selectAxialFlowDiscretizationLRM(UnitOpIdx uoId, IParameterProvider& paramProvider)
	{
		IUnitOperation* model = nullptr;

		paramProvider.pushScope("discretization");

		if (paramProvider.exists("SPATIAL_METHOD")) {

			const std::string discName = paramProvider.getString("SPATIAL_METHOD");

#ifdef ENABLE_DG
			if(discName == "DG")
				model = new LumpedRateModelWithoutPoresDG(uoId);
			else if (discName == "FV")
#else
			if (discName == "FV")
#endif
				model = createAxialFVLRM(uoId);
			else
			{
				LOG(Error) << "Unknown discretization type " << discName << " for unit " << uoId;
			}
		}
		else {
			model = createAxialFVLRM(uoId);
		}

		paramProvider.popScope();

		return model;
	}

	IUnitOperation* selectRadialFlowDiscretizationLRM(UnitOpIdx uoId, IParameterProvider& paramProvider)
	{
		IUnitOperation* model = nullptr;

		paramProvider.pushScope("discretization");

		if (paramProvider.exists("SPATIAL_METHOD")) {

			const std::string discName = paramProvider.getString("SPATIAL_METHOD");

			if (discName == "DG")
			{
				LOG(Error) << "Radial flow not implemented for DG discretization yet, was called for unit " << uoId;
			}
			else if (discName == "FV")
				model = createRadialFVLRM(uoId);
			else
			{
				LOG(Error) << "Unknown discretization type " << discName << " for unit " << uoId;
			}
		}
		else {
			model = createRadialFVLRM(uoId);
		}

		paramProvider.popScope();

		return model;
	}

void registerLumpedRateModelWithoutPores(std::unordered_map<std::string, std::function<IUnitOperation* (UnitOpIdx, IParameterProvider&)>>& models)
{
	typedef LumpedRateModelWithoutPores<parts::AxialConvectionDispersionOperatorBase> AxialLRM;
	typedef LumpedRateModelWithoutPores<parts::RadialConvectionDispersionOperatorBase> RadialLRM;

#ifdef ENABLE_DG
	models[LumpedRateModelWithoutPoresDG::identifier()] = selectAxialFlowDiscretizationLRM;
	models["LRM_DG"] = selectAxialFlowDiscretizationLRM;
#endif

	models[AxialLRM::identifier()] = selectAxialFlowDiscretizationLRM;
	models["LRM"] = selectAxialFlowDiscretizationLRM;
	models["DPFR"] = selectAxialFlowDiscretizationLRM;

	models[RadialLRM::identifier()] = selectRadialFlowDiscretizationLRM;
	models["RLRM"] = selectRadialFlowDiscretizationLRM;
}

}  // namespace model

}  // namespace cadet
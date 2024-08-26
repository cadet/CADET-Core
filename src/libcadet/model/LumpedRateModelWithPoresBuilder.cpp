#include "model/LumpedRateModelWithPores.hpp"
#include "CompileTimeConfig.hpp"
#ifdef ENABLE_DG
#include "model/LumpedRateModelWithPoresDG.hpp"
#endif
#include "LoggingUtils.hpp"
#include "Logging.hpp"

namespace cadet
{
namespace model
{

IUnitOperation* selectAxialFlowDiscretizationLRMP(UnitOpIdx uoId, IParameterProvider& paramProvider)
{
	IUnitOperation* model = nullptr;

	paramProvider.pushScope("discretization");

	if (paramProvider.exists("SPATIAL_METHOD"))
	{

		const std::string discName = paramProvider.getString("SPATIAL_METHOD");

#ifdef ENABLE_DG
		if (discName == "DG")
			model = new LumpedRateModelWithPoresDG(uoId);
		else if (discName == "FV")
#else
		if (discName == "FV")
#endif
			model = createAxialFVLRMP(uoId);
		else
		{
			LOG(Error) << "Unknown discretization type " << discName << " for unit " << uoId;
		}
	}
	else
	{
		model = createAxialFVLRMP(uoId);
	}

	paramProvider.popScope();

	return model;
}

IUnitOperation* selectRadialFlowDiscretizationLRMP(UnitOpIdx uoId, IParameterProvider& paramProvider)
{
	IUnitOperation* model = nullptr;

	paramProvider.pushScope("discretization");

	if (paramProvider.exists("SPATIAL_METHOD"))
	{

		const std::string discName = paramProvider.getString("SPATIAL_METHOD");

		if (discName == "DG")
		{
			LOG(Error) << "Radial flow not implemented for DG discretization yet, was called for unit " << uoId;
		}
		else if (discName == "FV")
			model = createRadialFVLRMP(uoId);
		else
		{
			LOG(Error) << "Unknown discretization type " << discName << " for unit " << uoId;
		}
	}
	else
	{
		model = createRadialFVLRMP(uoId);
	}

	paramProvider.popScope();

	return model;
}

void registerLumpedRateModelWithPores(
	std::unordered_map<std::string, std::function<IUnitOperation*(UnitOpIdx, IParameterProvider&)>>& models)
{
	typedef LumpedRateModelWithPores<parts::AxialConvectionDispersionOperator> AxialLRMP;
	typedef LumpedRateModelWithPores<parts::RadialConvectionDispersionOperator> RadialLRMP;

#ifdef ENABLE_DG
	models[LumpedRateModelWithPoresDG::identifier()] = selectAxialFlowDiscretizationLRMP;
	models["LRMP_DG"] = selectAxialFlowDiscretizationLRMP;
#endif
	models[AxialLRMP::identifier()] = selectAxialFlowDiscretizationLRMP;
	models["LRMP"] = selectAxialFlowDiscretizationLRMP;

	models[RadialLRMP::identifier()] = selectRadialFlowDiscretizationLRMP;
	models["RLRMP"] = selectRadialFlowDiscretizationLRMP;
}

} // namespace model

} // namespace cadet

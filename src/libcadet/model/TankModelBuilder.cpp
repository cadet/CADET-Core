#include "model/StirredTankModel.hpp"
#include "model/StirredTankModelpH.hpp"
#include "CompileTimeConfig.hpp"
#include "LoggingUtils.hpp"
#include "Logging.hpp"


namespace cadet
{
namespace model
{
	IUnitOperation* selectTankUnitOperation(UnitOpIdx uoId, IParameterProvider& paramProvider)
	{
		IUnitOperation* model = nullptr;

		std::string uoType = paramProvider.getString("UNIT_TYPE");

		if (uoType != "CSTR")
			throw InvalidParameterException("Unknown tank type " + uoType + " for unit " + std::to_string(uoId));

		if (paramProvider.exists("PH_CONTROL") && paramProvider.getBool("PH_CONTROL"))
			model = new CSTRpH(uoId);
		else
			model = new CSTRModel(uoId);

		return model;
	}

	void registerTankModel(std::unordered_map<std::string, std::function<IUnitOperation* (UnitOpIdx, IParameterProvider&)>>& models)
	{
		models[CSTRpH::identifier()] = selectTankUnitOperation;
		models[CSTRModel::identifier()] = selectTankUnitOperation;
	}

}  // namespace model

}  // namespace cadet
#include "model/GeneralRateModel.hpp"
#include "CompileTimeConfig.hpp"
#include "model/ColumnModel1D.hpp"
#include "LoggingUtils.hpp"
#include "Logging.hpp"


namespace cadet
{
	namespace model
	{

		IUnitOperation* selectAxialFlowDiscretizationGRM(UnitOpIdx uoId, IParameterProvider& paramProvider)
		{
			IUnitOperation* model = nullptr;

			paramProvider.pushScope("discretization");

			if (paramProvider.exists("SPATIAL_METHOD")) {

				const std::string discName = paramProvider.getString("SPATIAL_METHOD");
				if (discName == "DG")
					model = new ColumnModel1D(uoId);
				else if (discName == "FV")
					model = createAxialFVGRM(uoId);
				else
				{
					LOG(Error) << "Unknown discretization type " << discName << " for unit " << uoId;
				}
			}
			else {
				model = createAxialFVGRM(uoId);
			}

			paramProvider.popScope();

			return model;
		}

		IUnitOperation* selectRadialFlowDiscretizationGRM(UnitOpIdx uoId, IParameterProvider& paramProvider)
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
					model = createRadialFVGRM(uoId);
				else
				{
					LOG(Error) << "Unknown discretization type " << discName << " for unit " << uoId;
				}
			}
			else {
				model = createRadialFVGRM(uoId);
			}

			paramProvider.popScope();

			return model;
		}

		void registerGeneralRateModel(std::unordered_map<std::string, std::function<IUnitOperation* (UnitOpIdx, IParameterProvider&)>>& models)
		{
			typedef GeneralRateModel<parts::AxialConvectionDispersionOperator> AxialGRM;
			typedef GeneralRateModel<parts::RadialConvectionDispersionOperator> RadialGRM;

			models[ColumnModel1D::identifier()] = selectAxialFlowDiscretizationGRM;
			models["COLUMN_MODEL_1D"] = selectAxialFlowDiscretizationGRM;

			models[AxialGRM::identifier()] = selectAxialFlowDiscretizationGRM;
			models["GRM"] = selectAxialFlowDiscretizationGRM;

			models[RadialGRM::identifier()] = selectRadialFlowDiscretizationGRM;
			models["RGRM"] = selectRadialFlowDiscretizationGRM;
		}

	}  // namespace model

}  // namespace cadet
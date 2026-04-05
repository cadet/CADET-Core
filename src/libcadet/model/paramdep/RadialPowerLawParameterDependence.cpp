// =============================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/paramdep/ParameterDependenceBase.hpp"
#include "SimulationTypes.hpp"
#include "cadet/ParameterId.hpp"

#include <cmath>

namespace cadet
{
namespace model
{

/**
 * @brief Defines a radial power law parameter parameter dependence
 * @details The parameter dependence is defined as
 * @f[\begin{align}
	p_{new}(p, r) = p \cdot \alpha \cdot (1 + \beta \cdot \hat{r}^k),
   \end{align} @f]
   where @f$ \hat{r} \in [0, 1] @f$ is the normalized radial position,
   @f$ \alpha @f$ is a base factor, @f$ \beta @f$ is the radial factor,
   and @f$ k @f$ is the exponent.
 */
class RadialPowerLawParameterParameterDependence : public ParameterParameterDependenceBase
{
public:

	RadialPowerLawParameterParameterDependence() { }
	virtual ~RadialPowerLawParameterParameterDependence() CADET_NOEXCEPT { }

	static const char* identifier() { return "RADIAL_POWER_LAW"; }
	virtual const char* name() const CADET_NOEXCEPT { return RadialPowerLawParameterParameterDependence::identifier(); }

	CADET_PARAMETERPARAMETERDEPENDENCE_BOILERPLATE

protected:
	active _base;
	active _factor;
	active _exponent;

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx, BoundStateIdx bndIdx, const std::string& name)
	{
		const std::string baseName = name + "_BASE";
		const std::string factorName = name + "_FACTOR";
		const std::string expName = name + "_EXPONENT";

		if (paramProvider.exists(baseName))
			_base = paramProvider.getDouble(baseName);
		else
			_base = 1.0;

		if (paramProvider.exists(factorName))
			_factor = paramProvider.getDouble(factorName);
		else
			_factor = 1.0;

		if (paramProvider.exists(expName))
			_exponent = paramProvider.getDouble(expName);
		else
			_exponent = 1.0;

		_parameters[makeParamId(hashStringRuntime(baseName), unitOpIdx, CompIndep, parTypeIdx, bndIdx, ReactionIndep, SectionIndep)] = &_base;
		_parameters[makeParamId(hashStringRuntime(factorName), unitOpIdx, CompIndep, parTypeIdx, bndIdx, ReactionIndep, SectionIndep)] = &_factor;
		_parameters[makeParamId(hashStringRuntime(expName), unitOpIdx, CompIndep, parTypeIdx, bndIdx, ReactionIndep, SectionIndep)] = &_exponent;

		return true;
	}

	template <typename ParamType>
	ParamType getValueImpl(const IModel& model, const ColumnPosition& colPos, int comp, int parType, int bnd) const
	{
		return 0.0;
	}

	template <typename ParamType>
	ParamType getValueImpl(const IModel& model, const ColumnPosition& colPos, int comp, int parType, int bnd, ParamType val) const
	{
		using std::pow;

		const ParamType relPos = static_cast<ParamType>(colPos.axial);

		return val * static_cast<ParamType>(_base) * (1.0 + static_cast<ParamType>(_factor) * pow(relPos, static_cast<ParamType>(_exponent)));
	}

};


namespace paramdep
{
	void registerRadialPowerLawParamDependence(std::unordered_map<std::string, std::function<model::IParameterParameterDependence*()>>& paramDeps)
	{
		paramDeps[RadialPowerLawParameterParameterDependence::identifier()] = []() { return new RadialPowerLawParameterParameterDependence(); };
	}
}  // namespace paramdep

}  // namespace model

}  // namespace cadet

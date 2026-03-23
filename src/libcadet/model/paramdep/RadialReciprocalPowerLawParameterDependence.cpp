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
 * @brief Defines a radial reciprocal power law parameter parameter dependence
 * @details The parameter dependence is defined as
 * @f[\begin{align}
	p_{new}(p, r) = p \cdot \alpha \cdot \left(\frac{\rho_{in}}{\rho_{in} + L \cdot \hat{r}}\right)^k,
   \end{align} @f]
   where @f$ \hat{r} \in [0, 1] @f$ is the normalized radial position,
   @f$ \rho_{in} @f$ is the inner radius, @f$ L = \rho_{out} - \rho_{in} @f$
   is the column length in the radial direction, @f$ \alpha @f$ is a base
   factor, and @f$ k @f$ is the exponent.

   This models physically motivated radial dependencies where parameters
   scale with the interstitial velocity u(rho) = Q / (2*pi*rho*L_col):
   - D_rad proportional to u:      exponent = 1   =>  D(rho) = D(rin) * (rin/rho)
   - k_f proportional to u^(1/3):  exponent = 1/3 =>  kf(rho) = kf(rin) * (rin/rho)^(1/3)

   Configuration:
   - <PARAM>_BASE:     base factor alpha (default 1.0)
   - <PARAM>_EXPONENT: power law exponent k (default 1.0)
   - <PARAM>_RINNER:   inner radius rho_in (required)
   - <PARAM>_LENGTH:   radial column length L (required)
 */
class RadialReciprocalPowerLawParameterParameterDependence : public ParameterParameterDependenceBase
{
public:

	RadialReciprocalPowerLawParameterParameterDependence() { }
	virtual ~RadialReciprocalPowerLawParameterParameterDependence() CADET_NOEXCEPT { }

	static const char* identifier() { return "RADIAL_RECIPROCAL_POWER_LAW"; }
	virtual const char* name() const CADET_NOEXCEPT { return RadialReciprocalPowerLawParameterParameterDependence::identifier(); }

	CADET_PARAMETERPARAMETERDEPENDENCE_BOILERPLATE

protected:
	active _base;
	active _exponent;
	active _rInner;
	active _length;

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx, BoundStateIdx bndIdx, const std::string& name)
	{
		const std::string baseName = name + "_BASE";
		const std::string expName = name + "_EXPONENT";
		const std::string rinName = name + "_RINNER";
		const std::string lenName = name + "_LENGTH";

		if (paramProvider.exists(baseName))
			_base = paramProvider.getDouble(baseName);
		else
			_base = 1.0;

		if (paramProvider.exists(expName))
			_exponent = paramProvider.getDouble(expName);
		else
			_exponent = 1.0;

		if (paramProvider.exists(rinName))
			_rInner = paramProvider.getDouble(rinName);
		else
			_rInner = 0.0;

		if (paramProvider.exists(lenName))
			_length = paramProvider.getDouble(lenName);
		else
			_length = 1.0;

		_parameters[makeParamId(hashStringRuntime(baseName), unitOpIdx, CompIndep, parTypeIdx, bndIdx, ReactionIndep, SectionIndep)] = &_base;
		_parameters[makeParamId(hashStringRuntime(expName), unitOpIdx, CompIndep, parTypeIdx, bndIdx, ReactionIndep, SectionIndep)] = &_exponent;
		_parameters[makeParamId(hashStringRuntime(rinName), unitOpIdx, CompIndep, parTypeIdx, bndIdx, ReactionIndep, SectionIndep)] = &_rInner;
		_parameters[makeParamId(hashStringRuntime(lenName), unitOpIdx, CompIndep, parTypeIdx, bndIdx, ReactionIndep, SectionIndep)] = &_length;

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

		// colPos.axial contains the normalized radial position [0, 1]
		const ParamType relPos = static_cast<ParamType>(colPos.axial);

		// rho = rInner + length * relPos
		// p_new = val * base * (rInner / rho)^exponent
		//       = val * base * (rInner / (rInner + length * relPos))^exponent
		const ParamType rho = static_cast<ParamType>(_rInner) + static_cast<ParamType>(_length) * relPos;
		const ParamType ratio = static_cast<ParamType>(_rInner) / rho;

		return val * static_cast<ParamType>(_base) * pow(ratio, static_cast<ParamType>(_exponent));
	}

};


namespace paramdep
{
	void registerRadialReciprocalPowerLawParamDependence(std::unordered_map<std::string, std::function<model::IParameterParameterDependence*()>>& paramDeps)
	{
		paramDeps[RadialReciprocalPowerLawParameterParameterDependence::identifier()] = []() { return new RadialReciprocalPowerLawParameterParameterDependence(); };
	}
}  // namespace paramdep

}  // namespace model

}  // namespace cadet
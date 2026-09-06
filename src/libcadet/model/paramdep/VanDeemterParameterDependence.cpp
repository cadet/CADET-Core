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
 * @brief Defines a van Deemter parameter parameter dependence
 * @details Reproduces a classical van Deemter plate-height correlation,
 *   H(v) = A + B / |v| + C * |v|,
 *   applied here to a dispersion coefficient (e.g. COL_DISPERSION) so that
 *   the *effective* dependence factor returned is
 *   p_new(p, v) = p * H(v) * |v| / 2 = p * (A * |v| + B + C * v^2) / 2,
 *   i.e., with the associated base parameter p left at 1 (the default,
 *   dimensionless placeholder), this directly evaluates the standard
 *   equilibrium-dispersive-model relation D_ax(v) = H(v) * v / 2 for a
 *   local, velocity-dependent (van Deemter) plate height H(v) -- allowing
 *   a genuine A + B/v + C*v flow-rate dependence (as opposed to the single
 *   monomial v^k achievable with POWER_LAW) to be used with
 *   COL_DISPERSION_DEP for velocity-varying bulk transport models
 *   (e.g. GEOMETRY=AXIAL_FLOW_FRUSTUM/RADIAL_FLOW_CYLINDER_SHELL).
 *   A, B, and C are (meta-)parameters.
 */
class VanDeemterParameterParameterDependence : public ParameterParameterDependenceBase
{
public:

	VanDeemterParameterParameterDependence() { }
	virtual ~VanDeemterParameterParameterDependence() CADET_NOEXCEPT { }

	static const char* identifier() { return "VAN_DEEMTER"; }
	virtual const char* name() const CADET_NOEXCEPT { return VanDeemterParameterParameterDependence::identifier(); }

	CADET_PARAMETERPARAMETERDEPENDENCE_BOILERPLATE

protected:
	active _a;
	active _b;
	active _c;

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx, BoundStateIdx bndIdx, const std::string& name)
	{
		const std::string aName = name + "_A";
		const std::string bName = name + "_B";
		const std::string cName = name + "_C";

		_a = paramProvider.getDouble(aName);
		_b = paramProvider.getDouble(bName);
		_c = paramProvider.getDouble(cName);

		_parameters[makeParamId(hashStringRuntime(aName), unitOpIdx, CompIndep, parTypeIdx, bndIdx, ReactionIndep, SectionIndep)] = &_a;
		_parameters[makeParamId(hashStringRuntime(bName), unitOpIdx, CompIndep, parTypeIdx, bndIdx, ReactionIndep, SectionIndep)] = &_b;
		_parameters[makeParamId(hashStringRuntime(cName), unitOpIdx, CompIndep, parTypeIdx, bndIdx, ReactionIndep, SectionIndep)] = &_c;

		return true;
	}

	template <typename ParamType>
	ParamType getValueImpl(const ColumnPosition& colPos, int comp, int parType, int bnd) const
	{
		return 0.0;
	}

	template <typename ParamType>
	ParamType getValueImpl(const ColumnPosition& colPos, int comp, int parType, int bnd, ParamType val) const
	{
		using std::abs;

		const ParamType v = abs(val);
		// H(v) * v / 2 = (A * v + B + C * v^2) / 2
		return (static_cast<ParamType>(_a) * v + static_cast<ParamType>(_b) + static_cast<ParamType>(_c) * v * v) * ParamType(0.5);
	}

};


namespace paramdep
{
	void registerVanDeemterParamDependence(std::unordered_map<std::string, std::function<model::IParameterParameterDependence*()>>& paramDeps)
	{
		paramDeps[VanDeemterParameterParameterDependence::identifier()] = []() { return new VanDeemterParameterParameterDependence(); };
	}
}  // namespace paramdep

}  // namespace model

}  // namespace cadet

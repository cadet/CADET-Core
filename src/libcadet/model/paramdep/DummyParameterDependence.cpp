// =============================================================================
//  CADET
//  
//  Copyright © 2008-2022: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/paramdep/ParameterDependenceBase.hpp"
#include "ParamReaderHelper.hpp"
#include "cadet/Exceptions.hpp"
#include "SimulationTypes.hpp"
#include "cadet/ParameterId.hpp"

#include <sstream>
#include <vector>
#include <iomanip>

namespace cadet
{

namespace model
{

/**
 * @brief Defines a dummy parameter dependence
 */
class DummyParameterStateDependence : public ParameterStateDependenceBase
{
public:

	DummyParameterStateDependence() { }
	virtual ~DummyParameterStateDependence() CADET_NOEXCEPT { }

	static const char* identifier() { return "DUMMY"; }
	virtual const char* name() const CADET_NOEXCEPT { return DummyParameterStateDependence::identifier(); }

	virtual int jacobianElementsPerRowLiquid() const CADET_NOEXCEPT { return 0; }
	virtual int jacobianElementsPerRowCombined() const CADET_NOEXCEPT { return 0; }

	virtual void analyticJacobianLiquidAdd(const ColumnPosition& colPos, double param, double const* y, int comp, double factor, int offset, int row, linalg::DoubleSparseMatrix& jac) const { }
	virtual void analyticJacobianCombinedAddLiquid(const ColumnPosition& colPos, double param, double const* yLiquid, double const* ySolid, int comp, double factor, int offset, int row, linalg::DoubleSparseMatrix& jac) const { }
	virtual void analyticJacobianCombinedAddSolid(const ColumnPosition& colPos, double param, double const* yLiquid, double const* ySolid, int bnd, double factor, int offset, int row, linalg::DoubleSparseMatrix& jac) const { }

	CADET_PARAMETERSTATEDEPENDENCE_BOILERPLATE

protected:

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx, const std::string& name)
	{
		return true;
	}

	template <typename StateType, typename ParamType>
	typename DoubleActivePromoter<StateType, ParamType>::type liquidParameterImpl(const ColumnPosition& colPos, const ParamType& param, StateType const* y, int comp) const
	{
		return 0.0;
	}

	template <typename RowIterator>
	void analyticJacobianLiquidAddImpl(const ColumnPosition& colPos, double param, double const* y, int comp, double factor, int offset, RowIterator jac) const { }

	template <typename StateType, typename ParamType>
	typename DoubleActivePromoter<StateType, ParamType>::type combinedParameterLiquidImpl(const ColumnPosition& colPos, const ParamType& param, StateType const* yLiquid, StateType const* ySolid, int comp) const
	{
		return 0.0;
	}

	template <typename RowIterator>
	void analyticJacobianCombinedAddLiquidImpl(const ColumnPosition& colPos, double param, double const* yLiquid, double const* ySolid, int comp, double factor, int offset, RowIterator jac) const { }

	template <typename StateType, typename ParamType>
	typename DoubleActivePromoter<StateType, ParamType>::type combinedParameterSolidImpl(const ColumnPosition& colPos, const ParamType& param, StateType const* yLiquid, StateType const* ySolid, int bnd) const
	{
		return 0.0;
	}

	template <typename RowIterator>
	void analyticJacobianCombinedAddSolidImpl(const ColumnPosition& colPos, double param, double const* yLiquid, double const* ySolid, int bnd, double factor, int offset, RowIterator jac) const { }
};


/**
 * @brief Defines a dummy parameter dependence
 */
class DummyParameterParameterDependence : public ParameterParameterDependenceBase
{
public:

	DummyParameterParameterDependence() { }
	virtual ~DummyParameterParameterDependence() CADET_NOEXCEPT { }

	static const char* identifier() { return "DUMMY"; }
	virtual const char* name() const CADET_NOEXCEPT { return DummyParameterParameterDependence::identifier(); }

	CADET_PARAMETERPARAMETERDEPENDENCE_BOILERPLATE

protected:

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx, BoundStateIdx bndIdx, const std::string& name)
	{
		return true;
	}

	template <typename ParamType>
	ParamType getValueImpl(UnitOpIdx unitOpIdx, const std::unordered_map<ParameterId, active*>& params, const ColumnPosition& colPos, int comp, int parType, int bnd) const
	{
		return 0.0;
	}

	template <typename ParamType>
	ParamType getValueImpl(UnitOpIdx unitOpIdx, const std::unordered_map<ParameterId, active*>& params, const ColumnPosition& colPos, int comp, int parType, int bnd, const ParamType& val) const
	{
		return 0.0;
	}

};


namespace paramdep
{
	void registerDummyParamDependence(std::unordered_map<std::string, std::function<model::IParameterStateDependence*()>>& paramDeps)
	{
		paramDeps[DummyParameterStateDependence::identifier()] = []() { return new DummyParameterStateDependence(); };
		paramDeps["NONE"] = []() { return new DummyParameterStateDependence(); };
	}

	void registerDummyParamDependence(std::unordered_map<std::string, std::function<model::IParameterParameterDependence*()>>& paramDeps)
	{
		paramDeps[DummyParameterParameterDependence::identifier()] = []() { return new DummyParameterParameterDependence(); };
		paramDeps["NONE"] = []() { return new DummyParameterParameterDependence(); };
	}
}  // namespace paramdep

}  // namespace model

}  // namespace cadet

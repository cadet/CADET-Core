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

/**
 * @file 
 * Defines a ReactionModel base classes.
 */

#ifndef LIBCADET_REACTIONMODELBASE_HPP_
#define LIBCADET_REACTIONMODELBASE_HPP_

#include "model/ReactionModel.hpp"
#include "ParamIdUtil.hpp"

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace cadet
{

namespace model
{

/**
 * @brief Defines a DynamicReactionModel base class that can be used to implement other dynamic reaction models
 * @details This base class can be used as a starting point for new dynamic reaction models.
 *          Some common parameter handling is provided using a hash map (std::unordered_map).
 */
class DynamicReactionModelBase : public IDynamicReactionModel
{
public:

	DynamicReactionModelBase();

	virtual ~DynamicReactionModelBase() CADET_NOEXCEPT;

	virtual bool requiresConfiguration() const CADET_NOEXCEPT { return true; }
	virtual bool usesParamProviderInDiscretizationConfig() const CADET_NOEXCEPT { return true; }
	virtual bool configure(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx);
	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset);

	virtual std::unordered_map<ParameterId, double> getAllParameterValues() const;
	virtual bool hasParameter(const ParameterId& pId) const;

	virtual bool setParameter(const ParameterId& pId, int value);
	virtual bool setParameter(const ParameterId& pId, double value);
	virtual bool setParameter(const ParameterId& pId, bool value);

	virtual active* getParameter(const ParameterId& pId);

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { }

	virtual bool isInEquilibrium(unsigned int reaction) const CADET_NOEXCEPT
	{
		return (reaction < _reactionQuasistationarity.size()) ? _reactionQuasistationarity[reaction] : false;
	}

	/**
	 * @brief Returns whether at least one reaction of this model is in rapid equilibrium
	 */
	bool hasQuasiStationaryReactions() const CADET_NOEXCEPT
	{
		return std::any_of(_reactionQuasistationarity.begin(), _reactionQuasistationarity.end(), [](bool qs) { return qs; });
	}

protected:
	int _nComp; //!< Number of components
	unsigned int const* _nBoundStates; //!< Array with number of bound states for each component
	unsigned int const* _boundOffset; //!< Array with offsets to the first bound state of each component
	int _nTotalBoundStates;

	std::vector<bool> _reactionQuasistationarity; //!< Determines whether each reaction is in rapid equilibrium (@c true) or kinetic (@c false)

	std::unordered_map<ParameterId, active*> _parameters; //!< Map used to translate ParameterIds to actual variables

	/**
	 * @brief Configures the reaction model
	 * @details This function implements the (re-)configuration of a reaction model. It is called when
	 *          the reaction model is configured or reconfigured. On call the _parameters map will always
	 *          be empty.
	 * @param [in] paramProvider Parameter provider
	 * @param [in] unitOpIdx Unit operation index
	 * @param [in] parTypeIdx Particle type index
	 * @return @c true if the configuration was successful, otherwise @c false
	 */
	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx) = 0;

	/**
	 * @brief Reads the reaction dynamics (kinetic or rapid equilibrium) of every reaction
	 * @details Reads the optional field @c IS_KINETIC, which is either a scalar applying to all
	 *          reactions or an array with one entry per reaction. Reactions default to kinetic.
	 *          Called after configureImpl(), since only then is numReactions() known.
	 * @param [in] paramProvider Parameter provider
	 * @throws InvalidParameterException if @c IS_KINETIC is given for a model without reactions or
	 *         has too few elements
	 */
	void readQuasiStationarity(IParameterProvider& paramProvider);
};


/**
 * @brief Inserts implementations of all residual() and analyticJacobian() method variants
 * @details An IDynamicReactionModel implementation has to provide residualFluxAdd(), residualCombinedAdd(),
 *          analyticJacobianAdd(), and analyticJacobianCombinedAdd() methods for different variants of state
 *          and parameter type. This macro saves some time by providing those implementations. It assumes that the 
 *          implementation provides templatized residualFluxImpl(), residualCombinedImpl(), jacobianFluxImpl(),
 *          and jacobianCombinedImpl() function that realize all required variants.
 *          
 *          The implementation is inserted inline in the class declaration.
 */
#define CADET_DYNAMICREACTIONMODEL_BOILERPLATE                                                                                                          \
	virtual int residualFluxAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, active const* y,                                         \
		active* res, const active& factor, LinearBufferAllocator workSpace) const                                                                       \
	{                                                                                                                                                   \
		return residualFluxImpl<active, active, double, active>(t, secIdx, colPos, nStates,y, res, factor, workSpace);                                        \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual int residualFluxAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, active const* y,                                         \
		active* res, double factor, LinearBufferAllocator workSpace) const                                                                              \
	{                                                                                                                                                   \
		return residualFluxImpl<active, active, double, double>(t, secIdx, colPos,nStates, y, res, factor, workSpace);                                        \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual int residualFluxAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,                                         \
		active* res, double factor, LinearBufferAllocator workSpace) const                                                                              \
	{                                                                                                                                                   \
		return residualFluxImpl<double, active, active, double>(t, secIdx, colPos, nStates,y, res, factor, workSpace);                                        \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual int residualFluxAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,                                         \
		double* res, double factor, LinearBufferAllocator workSpace) const                                                                              \
	{                                                                                                                                                   \
		return residualFluxImpl<double, double, double, double>(t, secIdx, colPos, nStates, y, res, factor, workSpace);                                        \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual int residualCombinedAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, active const* yLiquid,                                 \
		active const* ySolid, active* resLiquid, active* resSolid, double factor, LinearBufferAllocator workSpace) const                                \
	{                                                                                                                                                   \
		return residualCombinedImpl<active, active, double>(t, secIdx, colPos, yLiquid, ySolid, resLiquid, resSolid, factor, workSpace);                \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual int residualCombinedAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid,                                 \
		double const* ySolid, active* resLiquid, active* resSolid, double factor, LinearBufferAllocator workSpace) const                                \
	{                                                                                                                                                   \
		return residualCombinedImpl<double, active, active>(t, secIdx, colPos, yLiquid, ySolid, resLiquid, resSolid, factor, workSpace);                \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual int residualCombinedAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid,                                 \
		double const* ySolid, double* resLiquid, double* resSolid, double factor, LinearBufferAllocator workSpace) const                                \
	{                                                                                                                                                   \
		return residualCombinedImpl<double, double, double>(t, secIdx, colPos, yLiquid, ySolid, resLiquid, resSolid, factor, workSpace);                \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual void analyticJacobianAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,                                \
		double factor, linalg::BandMatrix::RowIterator jac, LinearBufferAllocator workSpace) const                                                      \
	{                                                                                                                                                   \
		jacobianFluxImpl(t, secIdx, colPos, nStates, y, factor, jac, workSpace);                                                                               \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual void analyticJacobianAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,                                \
		double factor, linalg::DenseBandedRowIterator jac, LinearBufferAllocator workSpace) const                                                       \
	{                                                                                                                                                   \
		jacobianFluxImpl(t, secIdx, colPos, nStates, y, factor, jac, workSpace);                                                                               \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual void analyticJacobianAdd(double t, unsigned int secIdx, const ColumnPosition& colPos,const unsigned int nStates, double const* y,                                \
		double factor, linalg::BandedSparseRowIterator jac, LinearBufferAllocator workSpace) const                                                      \
	{                                                                                                                                                   \
		jacobianFluxImpl(t, secIdx, colPos, nStates, y, factor, jac, workSpace);                                                                               \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual void analyticJacobianAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,                                \
		double factor, linalg::BandedEigenSparseRowIterator jac, LinearBufferAllocator workSpace) const                                                 \
	{                                                                                                                                                   \
		jacobianFluxImpl(t, secIdx, colPos, nStates, y, factor, jac, workSpace);                                                                               \
	}                                                                                                                                                   \
		                                                                                                                                                \
	virtual void analyticJacobianCombinedAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid,  \
		double factor, linalg::BandedEigenSparseRowIterator jacLiquid, linalg::DenseBandedRowIterator jacSolid, LinearBufferAllocator workSpace) const \
	{                                                                                                                                                   \
		jacobianCombinedImpl(t, secIdx, colPos, yLiquid, ySolid, factor, jacLiquid, jacSolid, workSpace);                                               \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual void analyticJacobianCombinedAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid,  \
		double factor, linalg::BandMatrix::RowIterator jacLiquid, linalg::BandMatrix::RowIterator jacSolid, LinearBufferAllocator workSpace) const      \
	{                                                                                                                                                   \
		jacobianCombinedImpl(t, secIdx, colPos, yLiquid, ySolid, factor, jacLiquid, jacSolid, workSpace);                                               \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual void analyticJacobianCombinedAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid,  \
		double factor, linalg::DenseBandedRowIterator jacLiquid, linalg::DenseBandedRowIterator jacSolid, LinearBufferAllocator workSpace) const        \
	{                                                                                                                                                   \
		jacobianCombinedImpl(t, secIdx, colPos, yLiquid, ySolid, factor, jacLiquid, jacSolid, workSpace);                                               \
	}                                                                                                                                                   \
	                                                                                                                                                    \
	virtual void analyticJacobianCombinedAdd(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid,  \
		double factor, linalg::BandMatrix::RowIterator jacLiquid, linalg::DenseBandedRowIterator jacSolid, LinearBufferAllocator workSpace) const       \
	{                                                                                                                                                   \
		jacobianCombinedImpl(t, secIdx, colPos, yLiquid, ySolid, factor, jacLiquid, jacSolid, workSpace);                                               \
	}																																					\
	virtual int residualEquilibriumFlux(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, active const* y,		\
		active* res, unsigned int& eqIdx, LinearBufferAllocator workSpace) const 																		\
	{																																					\
		return residualEquilibriumImpl<active, active, double>(t, secIdx, colPos, nStates, y, res, eqIdx , workSpace);									\
	}																																					\
	virtual int residualEquilibriumFlux(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,		\
		active* res, unsigned int& eqIdx, LinearBufferAllocator workSpace) const																		\
	{																																					\
		return residualEquilibriumImpl<double, active, active>(t, secIdx, colPos, nStates, y, res, eqIdx , workSpace);									\
	}																																					\
	virtual int residualEquilibriumFlux(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,		\
		double* res, unsigned int& eqIdx, LinearBufferAllocator workSpace) const																		\
	{																																					\
		return residualEquilibriumImpl<double, double, double>(t, secIdx, colPos, nStates, y, res, eqIdx , workSpace);\
	}\
	\
	virtual void analyticEquilibriumJacobian(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,\
		unsigned int& eqIdx, unsigned int eqRowOffset,linalg::BandMatrix::RowIterator jac, LinearBufferAllocator workSpace) const\
	{\
		jacobianEquilibriumImpl(t,secIdx,colPos, nStates, y, eqIdx, eqRowOffset, jac,workSpace);\
	}\
	virtual void analyticEquilibriumJacobian(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,\
		unsigned int& eqIdx,unsigned int eqRowOffset, linalg::DenseBandedRowIterator jac, LinearBufferAllocator workSpace) const\
	{\
		jacobianEquilibriumImpl(t,secIdx,colPos, nStates, y, eqIdx, eqRowOffset, jac,workSpace);\
	}\
	virtual void analyticEquilibriumJacobian(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,\
		unsigned int& eqIdx,unsigned int eqRowOffset, linalg::BandedSparseRowIterator jac, LinearBufferAllocator workSpace) const\
	{\
		jacobianEquilibriumImpl(t,secIdx,colPos, nStates, y, eqIdx, eqRowOffset, jac,workSpace);\
	}\
	virtual void analyticEquilibriumJacobian(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y,\
		unsigned int& eqIdx, unsigned int eqRowOffset,linalg::BandedEigenSparseRowIterator jac, LinearBufferAllocator workSpace) const \
	{\
		jacobianEquilibriumImpl(t,secIdx,colPos, nStates, y, eqIdx,eqRowOffset, jac,workSpace);\
	}

} // namespace model
} // namespace cadet

#endif  // LIBCADET_REACTIONMODELBASE_HPP_

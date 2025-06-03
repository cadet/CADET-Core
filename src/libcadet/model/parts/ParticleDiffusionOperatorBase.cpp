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

#include "model/parts/ParticleDiffusionOperatorBase.hpp"

namespace cadet
{

namespace model
{

namespace parts
{
	ParticleDiffusionOperatorBase::ParticleDiffusionOperatorBase() : _boundOffset(nullptr), _nBound(nullptr)
	{
	}

	bool ParticleDiffusionOperatorBase::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
	{
		_strideBulkComp = strideBulkComp;

		_parTypeIdx = parTypeIdx;
		_nComp = nComp;

		_filmDiffusion.resize(_nComp); // filled in notifyDiscontinuousSectionTransition
		_poreAccessFactor.resize(_nComp); // filled in notifyDiscontinuousSectionTransition
		_invBetaP.resize(_nComp); // filled in notifyDiscontinuousSectionTransition

		std::vector<int> nBound;
		const bool newNBoundInterface = paramProvider.exists("NBOUND");

		paramProvider.pushScope("discretization");

		if (!newNBoundInterface && paramProvider.exists("NBOUND")) // done here and in this order for backwards compatibility
			nBound = paramProvider.getIntArray("NBOUND");
		else
		{
			paramProvider.popScope();
			nBound = paramProvider.getIntArray("NBOUND");
			paramProvider.pushScope("discretization");
		}
		if (nBound.size() < _nComp)
			throw InvalidParameterException("Field NBOUND contains too few elements (NCOMP = " + std::to_string(_nComp) + " required)");

		std::vector<int> stridesParTypeBound(nParType + 1);
		std::vector<int> nBoundBeforeType(nParType);
		if (!_nBound)
			_nBound = new unsigned int[_nComp];

		if (nBound.size() < _nComp * nParType)
		{
			std::copy_n(nBound.begin(), _nComp, _nBound);

			stridesParTypeBound[0] = std::accumulate(nBound.begin(), nBound.begin() + _nComp, 0);
			nBoundBeforeType[0] = 0;

			for (int type = 1; type < nParType; type++)
			{
				stridesParTypeBound[type] = stridesParTypeBound[0];
				nBoundBeforeType[type] += nBoundBeforeType[type - 1] + _nBound[type - 1];
			}
		}
		else
		{
			std::copy_n(nBound.begin() + _parTypeIdx * _nComp, _nComp, _nBound);

			stridesParTypeBound[0] = std::accumulate(nBound.begin(), nBound.begin() + _nComp, 0);
			nBoundBeforeType[0] = 0;

			for (int type = 1; type < nParType; type++)
			{
				stridesParTypeBound[type] = std::accumulate(nBound.begin() + type * _nComp, nBound.begin() + (type + 1) * _nComp, 0);
				nBoundBeforeType[type] += nBoundBeforeType[type - 1] + nBound[type - 1];
			}
		}

		// Precompute offsets and total number of bound states (DOFs in solid phase)
		if (!_boundOffset)
			_boundOffset = new unsigned int[_nComp];

		_boundOffset[0] = 0.0;
		_strideBound = std::accumulate(_nBound, _nBound + _nComp, 0u);

		for (unsigned int i = 1; i < _nComp; ++i)
		{
			_boundOffset[i] = _boundOffset[i - 1] + _nBound[i - 1];
		}

		return nBoundBeforeType[_parTypeIdx];
	}


} // namespace parts
} // namespace model
} // namespace cadet
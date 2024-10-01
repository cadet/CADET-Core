// =============================================================================
//  CADET
//  
//  Copyright © 2008-2024: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================



#include "model/ExternalFunctionSupport.hpp"
#include "ParamIdUtil.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "Memory.hpp"
#include "model/PhaseTransitionModel.hpp"
#include "ParamReaderHelper.hpp"
#include "model/parts/MultiChannelConvectionDispersionOperator.hpp"

#include <vector>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <iterator>
#include <tuple>

namespace cadet
{

namespace model
{
/**
 * @brief Defines the linear exchange model
 * @details Implements the linear exchange model for the MCT model
 * The exchange is given by a matrix of exchange coefficients eij for each component i and j and the chross Ai section from channel i.
 * The exchange from channel i to all other channel j is given by \f$ \frac{\mathrm{d}ci}{\mathrm{d}t} = \sum_j eij cj Aj/Ai - eji ci \f$.
 */
class LangumirLDFExchangeBase : public IPhaseTransitionModel
{
public:

	LangumirLDFExchangeBase() : _nComp(0), _nChannel(0), _nCol(0) { }
	virtual ~LangumirLDFExchangeBase() CADET_NOEXCEPT { }

	static const char* identifier() { return "LANGMUIR_LDF_EX"; }
	virtual const char* name() const CADET_NOEXCEPT { return "LANGMUIR_LDF_EX"; }
	virtual bool requiresConfiguration() const CADET_NOEXCEPT { return true; }
	virtual bool usesParamProviderInDiscretizationConfig() const CADET_NOEXCEPT { return true; }

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nChannel, unsigned int nCol)
	{
		_nComp = nComp;
		_nChannel = nChannel;
		_nCol = nCol;

		return true;
	}

	virtual bool configure(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx)
	{
		_parameters.clear();
		readParameterMatrix(_exchangeMatrix, paramProvider, "EXCHANGE_MATRIX", _nChannel * _nChannel * _nComp, 1); // include parameterPeaderHelp in exchange modul
		_crossSections = paramProvider.getDoubleArray("CHANNEL_CROSS_SECTION_AREAS");
		readParameterMatrix(_capacityMatrix, paramProvider, "CAPACITY_MATRIX", _nComp * _nChannel, 1); // include parameterPeaderHelp in exchange modul

		return true;
	}

	virtual void fillExchangeInitialParameters(ParameterId* params, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx) const CADET_NOEXCEPT
	{
		unsigned int ctr = 0;
		for (int c = 0; c < _nComp; ++c)
		{
			for (unsigned int bp = 0; bp < _nChannel; ++bp, ++ctr)
				params[ctr] = makeParamId(hashString("INIT_C"), unitOpIdx, c, parTypeIdx, bp, ReactionIndep, SectionIndep);
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------// 
	virtual std::unordered_map<ParameterId, double> getAllParameterValues() const
	{
		std::unordered_map<ParameterId, double> data;
		std::transform(_parameters.begin(), _parameters.end(), std::inserter(data, data.end()),
			[](const std::pair<const ParameterId, active*>& p) { return std::make_pair(p.first, static_cast<double>(*p.second)); });

		return data;
	}

	virtual bool hasParameter(const ParameterId& pId) const
	{
		return _parameters.find(pId) != _parameters.end();
	}

	virtual bool setParameter(const ParameterId& pId, int value)
	{
		return false;
	}

	virtual bool setParameter(const ParameterId& pId, double value)
	{
		auto paramHandle = _parameters.find(pId);
		if (paramHandle != _parameters.end())
		{
			paramHandle->second->setValue(value);
			return true;
		}

		return false;
	}

	virtual bool setParameter(const ParameterId& pId, bool value)
	{
		return false;
	}

	virtual active* getParameter(const ParameterId& pId) {
		auto paramHandle = _parameters.find(pId);
		if (paramHandle != _parameters.end())
		{
			return paramHandle->second;
		}

		return nullptr;
	}

	virtual void configure(IExternalFunction** extFuns, unsigned int size) {}

	virtual int residual(active const* y, active* res, WithParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const
	{
		if (wantJac)
			return residualImpl<double, active, active, true>(reinterpret_cast<const double*>(y), res, jacBegin);
		else
			return residualImpl<active, active, active, false>(y, res, jacBegin);
	}

	virtual int residual(active const* y, active* res, WithoutParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const
	{	
		if (wantJac)
			return residualImpl<double, active, double, false>(reinterpret_cast<const double*>(y), res, jacBegin);
		else
			return residualImpl<active, active, double, false>(y, res, jacBegin);
	}
	
	virtual int residual(double const* y, active* res, WithParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const
	{
		if(wantJac)
			return residualImpl<double, active, active, true>(y, res, jacBegin);
		else
			return residualImpl<double, active, active, false>(y, res, jacBegin);
	}

	virtual int residual(double const* y, double* res, WithoutParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const
	{
		if (wantJac)
			return residualImpl<double, double, double, true>(y, res, jacBegin);
		else
			return residualImpl<double, double, double, false>(y, res, jacBegin);
	}


protected:
	int _nComp; //!< Number of components
	unsigned int _nChannel; //!< Total number of bound states
	unsigned int _nCol; //!< Number of columns

	std::vector<active> _exchangeMatrix; //!< Matrix of exchange coeffs for the linear inter-channel transport
	std::vector<double> _crossSections; //!< Cross sections of the channels
	std::vector<active> _capacityMatrix; //!< Capacity of the channels -> double  ncomp x nchannel
	parts::MultiChannelConvectionDispersionOperator _conDis; //!< Convection dispersion operator

	std::unordered_map<ParameterId, active*> _parameters; //!< Map used to translate ParameterIds to actual variables

	virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }


	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualImpl(StateType const* y, ResidualType* res, linalg::BandedSparseRowIterator jacBegin) const
	{	

		const unsigned int offsetC = _nChannel * _nComp;
		for (unsigned int col = 0; col < _nCol; ++col) // Collum
		{
			const unsigned int offsetColBlock = col * _nChannel * _nComp;
			ResidualType* const resColBlock = res + offsetC + offsetColBlock;
			StateType const* const yColBlock = y + offsetC + offsetColBlock;

			for (unsigned int rad_orig = 0; rad_orig < _nChannel; ++rad_orig)  // Channel orig
			{
				const unsigned int offsetToRadOrigBlock = rad_orig * _nComp;
				const unsigned int offsetColRadOrigBlock = offsetColBlock + offsetToRadOrigBlock;
				ResidualType* const resColRadOrigBlock = resColBlock + offsetToRadOrigBlock;
				StateType const* const yColRadOrigBlock = yColBlock + offsetToRadOrigBlock;

				for (unsigned int rad_dest = 0; rad_dest < _nChannel; ++rad_dest) // Channel dest
				{
					if (rad_orig == rad_dest)
						continue;

					const unsigned int offsetToRadDestBlock = rad_dest * _nComp;
					const unsigned int offsetColRadDestBlock = offsetColBlock + offsetToRadDestBlock;
					ResidualType* const resColRadDestBlock = resColBlock + offsetToRadDestBlock;
					StateType const* const yColRadDestBlock = yColBlock + offsetToRadDestBlock;

					for (unsigned int comp = 0; comp < _nComp; ++comp) // Component
					{
						const unsigned int offsetCur_orig = offsetColRadOrigBlock + comp; // component in orig channel
						const unsigned int offsetCur_dest = offsetColRadDestBlock + comp; // component in dest channel
						StateType const* const yCur_orig = yColRadOrigBlock + comp; // state of component in orig channel
						StateType const* const yCur_dest = yColRadDestBlock + comp; // state of component in dest channel
						ResidualType* const resCur_orig = resColRadOrigBlock + comp; // residual of component in orig channel
						ResidualType* const resCur_dest = resColRadDestBlock + comp; // residual of component in dest channel

						const ParamType exchange_orig_comp = static_cast<ParamType>(_exchangeMatrix[rad_orig * _nChannel * _nComp + rad_dest * _nComp + comp]);
						const Paramtype exchange_dest_comp = static_cast<ParamType>(_exchangeMatrix[rad_dest * _nChannel * _nComp + rad_orig * _nComp + comp]);

						double keq_comp_orig = exchange_orig_comp / exchange_dest_comp;

						const ParamType capacity_orig_comp = static_cast<ParamType>(_capacityMatrix[rad_orig * _nComp + comp]);
						const ParamType capacity_dest_dest = static_cast<ParamType>(_capacityMatrix[rad_dest * _nComp + comp]);

						if (cadet_likely(exchange_dest_comp > 0.0 ))
						{
							double compSum = 0.0;
							for (unsigned int component = 0; component < _nComp; component++) {

								StateType const* const compSum = yColRadDestBlock + component;
								double compSum_dest = static_cast<double>(compSum[0]);
								

								const ParamType exchange_orig_component = static_cast<ParamType>(_exchangeMatrix[rad_orig * _nChannel * _nComp + rad_dest * _nComp + component]);
								const Paramtype exchange_dest_component = static_cast<ParamType>(_exchangeMatrix[rad_dest * _nChannel * _nComp + rad_orig * _nComp + component]);

								double keq_comp  = exchange_orig_component / exchange_dest_component;
								
								if (cMax_destComp < 1e-16)
								{
									continue;
								}

								compSum += compSum_dest * keq_comp;

								//double chanSum_orig = 0.0;
								//double chanSum_dest = 0.0;
								//for (unsigned int channel = 0; channel < _nChannel; channel++) {
								// 	
								//	chanSum_orig += static_cast<double>(yColBlock[component + channel * _nChannel]) / cMax_origComp;
								//
								// 	chanSum_dest += static_cast<double>(yColBlock[component + channel * _nChannel]) / cMax_destComp;
								//}

							 }

							double qstar_nom = capacity_orig_comp * keq_comp_orig * yCur_orig[0];
							double qstar_denom = 1 + compSum
							*resCur_orig += (qstar_nom / qstar_denom) - yCur_orig[0];

							if (wantJac) {

								linalg::BandedSparseRowIterator jacorig;
								jacorig = jacBegin + offsetCur_orig;
								jacorig[0] += keq_comp_orig / qstar_denom;
							}


						}


					}

				}

			}
		}
		return 0;
	}
};

typedef LangumirLDFExchangeBase  LangumirExchange;

namespace exchange
{
	void registerLangumirExModel(std::unordered_map<std::string, std::function<model::IPhaseTransitionModel*()>>& exchange)
	{
		exchange[LangumirExchange::identifier()] = []() { return new  LangumirExchange(); };
	}
}  // namespace exchange

} // namespace model

}// namespace cadet

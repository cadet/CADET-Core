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

#include "model/reaction/ReactionModelBase.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "linalg/ActiveDenseMatrix.hpp"
#include "Memory.hpp"

#include <functional>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <iterator>
#include <vector>
#include <bitset>

namespace cadet
{
namespace model
{
namespace detail
{
	class CrystallizationMode {
	public:
		bool hasPBM()    const { return state_[0]; };
		bool hasAggregation()  const { return state_[1]; };
		bool hasFragmentation()  const { return state_[2]; };
		unsigned int getMode() const { return state_.to_ulong(); };

		void SetMode(unsigned int mode)
		{
			if (mode > 7 || mode < 1)
				throw InvalidParameterException("Crystallization mode specified as " + std::to_string(mode) + ", which is not defined.");

			state_ = std::bitset<3>(mode);
		}
	private:
		std::bitset<3> state_;
	};

	/**
	 * Numerical reconstruction type for solving the PBM
	*/
	enum class PBMReconstruction : int
	{
		/**
		 * Upwind scheme
		*/
		Upwind,

		/**
		 * High-resolution Koren-type reconstruction
		*/
		HRKoren,

		/**
		 * WENO23 reconstruction
		*/
		WENO23,

		/**
		 * WENO35 reconstruction
		*/
		WENO35
	};

	template <typename ParamType>
	struct ReconstructionParams
	{
		ParamType v_g;
		const ParamType k_g_times_s_g;
		const ParamType B_0;

		// WENO
		ParamType IS_0 = 0.0;
		ParamType IS_1 = 0.0;
		ParamType IS_2 = 0.0;
		ParamType alpha_0 = 0.0;
		ParamType alpha_1 = 0.0;
		ParamType alpha_2 = 0.0;
		ParamType W_0 = 0.0;
		ParamType W_1 = 0.0;
		ParamType W_2 = 0.0;
		ParamType q_0 = 0.0;
		ParamType q_1 = 0.0;
		ParamType q_2 = 0.0;

		// HRKoren
		ParamType r_x_i = 0.0;
		ParamType phi = 0.0;
		ParamType F_i = 0.0;

		ReconstructionParams(ParamType v_g, ParamType k_g_times_s_g, ParamType B_0,
			ParamType IS_0 = 0.0, ParamType IS_1 = 0.0, ParamType IS_2 = 0.0,
			ParamType alpha_0 = 0.0, ParamType alpha_1 = 0.0, ParamType alpha_2 = 0.0,
			ParamType W_0 = 0.0, ParamType W_1 = 0.0, ParamType W_2 = 0.0,
			ParamType q_0 = 0.0, ParamType q_1 = 0.0, ParamType q_2 = 0.0,
			ParamType r_x_i = 0.0, ParamType phi = 0.0, ParamType F_i = 0.0)
			: v_g(v_g), k_g_times_s_g(k_g_times_s_g), B_0(B_0),
			IS_0(IS_0), IS_1(IS_1), IS_2(IS_2),
			alpha_0(alpha_0), alpha_1(alpha_1), alpha_2(alpha_2),
			W_0(W_0), W_1(W_1), W_2(W_2),
			q_0(q_0), q_1(q_1), q_2(q_2),
			r_x_i(r_x_i), phi(phi), F_i(F_i) {}
	};

	struct AggCoefficients
	{
		std::pair<std::vector<unsigned short int>, std::vector<unsigned short int>> index;  // a pair of vectors to store combination of j and k
		std::vector<unsigned int> indexTracker;           // the idx can be a large number  
		std::vector<short int> iIndex;

		active sumVolumeCubeRoot = 0.0;
		unsigned int count = 0;

		// the sizes of j_index and k_index are not known a priori
		AggCoefficients(const std::vector<active>& binCenters, const std::vector<active>& bins, const std::vector<active>& binSizes)
			: index(), indexTracker(binCenters.size() + 1), iIndex(binCenters.size()* binCenters.size())
		{
			// make the first element of indexTracker 0
			indexTracker[0] = 0;

			// loop over daughter particle i
			for (unsigned int i = 0; i < binCenters.size(); ++i)
			{
				// source term, find lij in set A
				for (int j = 0; j < i; ++j)
				{
					for (int k = j; k < i + 1; ++k)
					{
						sumVolumeCubeRoot = pow(binCenters[k] * binCenters[k] * binCenters[k] + binCenters[j] * binCenters[j] * binCenters[j], 1.0 / 3.0);
						if ((sumVolumeCubeRoot > bins[i]) && (sumVolumeCubeRoot <= bins[i + 1]))
						{
							index.first.emplace_back(static_cast<short int>(j));
							index.second.emplace_back(static_cast<short int>(k));
							++count;
						}
					}
				}
				indexTracker[i + 1] = count;

				// sink term, loop over another mother particle j
				for (unsigned int j = 0; j < binCenters.size(); ++j)
				{
					sumVolumeCubeRoot = pow(binCenters[i] * binCenters[i] * binCenters[i] + binCenters[j] * binCenters[j] * binCenters[j], 1.0 / 3.0);
					int k = (i > j) ? i + 1 : j + 1;

					iIndex[i * binCenters.size() + j] = -1;        // -1 indicates there is no match
					for (k; k < binCenters.size(); ++k)
					{
						if ((sumVolumeCubeRoot > bins[k]) && (sumVolumeCubeRoot <= bins[k + 1]))
						{
							iIndex[i * binCenters.size() + j] = k;  // if there is a match, overwrite it
							break;
						}
					}
				}
			}
		}
	};

	struct FragCoefficients
	{
		// define fragmentation-related local parameters
		std::vector<active> upsilonSource;
		std::vector<active> upsilonSink;

		// constructor
		FragCoefficients(const std::vector<active>& binCenters, const std::vector<active>& bins, const active& fragKernelGamma)
			: upsilonSource(binCenters.size()), upsilonSink(binCenters.size())
		{
			const active N_j = fragKernelGamma / (fragKernelGamma - 1.0);

			active UpsilonBirthSum = 0.0;
			active UpsilonDeathSum = 0.0;
			active bIntegralBirth = 0.0;

			// calculate upsilon and store it
			for (int i = 0; i < binCenters.size(); ++i)
			{
				// reset the sum for each i
				UpsilonBirthSum = 0.0;
				UpsilonDeathSum = 0.0;
				const active x_i_3 = binCenters[i] * binCenters[i] * binCenters[i];
				for (int j = 0; j < i + 1; ++j)
				{
					const active x_j_3 = binCenters[j] * binCenters[j] * binCenters[j];
					if (cadet_likely(i != j))
					{
						bIntegralBirth = N_j * (pow(bins[j + 1], 3.0 * fragKernelGamma - 3.0) - pow(bins[j], 3.0 * fragKernelGamma - 3.0)) / pow(binCenters[i], 3.0 * fragKernelGamma - 3.0);
						UpsilonBirthSum += bIntegralBirth * (x_i_3 - x_j_3);
						UpsilonDeathSum += bIntegralBirth * x_j_3;
					}
					else
					{
						bIntegralBirth = N_j * (pow(binCenters[j], 3.0 * fragKernelGamma - 3.0) - pow(bins[j], 3.0 * fragKernelGamma - 3.0)) / pow(binCenters[i], 3.0 * fragKernelGamma - 3.0);
						UpsilonDeathSum += bIntegralBirth * x_j_3;
					}
				}
				if (cadet_likely(i > 0))
				{
					upsilonSource[i] = (N_j - 1.0) * x_i_3 / UpsilonBirthSum;
				}
				else
				{
					upsilonSource[i] = (N_j - 1.0) * x_i_3;
				}
				upsilonSink[i] = UpsilonDeathSum * upsilonSource[i] / x_i_3;
			}
		}
	};

} // namespace detail

/**
 * @brief Defines crystallization models as a reaction module. Details on the models and methods can be found in @cite Zhang2024 and @cite Zhang2025
*/
class CrystallizationReaction : public IDynamicReactionModel
{
public:

	CrystallizationReaction() : _nComp(0), _nBins(0), _bins(0), _binCenters(0), _binSizes(0), _agg(nullptr), _frag(nullptr), _reconstruction(nullptr), _jacParams(nullptr) { }
	virtual ~CrystallizationReaction() CADET_NOEXCEPT
	{
		if (_agg)
		{
			delete _agg;
			_agg = nullptr;
		}
		if (_frag)
		{
			delete _frag;
			_frag = nullptr;
		}
		if (_reconstruction)
		{
			delete _reconstruction;
			_reconstruction = nullptr;
		}
		if (_jacParams)
		{
			delete _jacParams;
			_jacParams = nullptr;
		}
	}

	static const char* identifier() { return "CRYSTALLIZATION"; }
	virtual const char* name() const CADET_NOEXCEPT { return identifier(); }

	virtual bool requiresConfiguration() const CADET_NOEXCEPT { return true; }
	virtual bool usesParamProviderInDiscretizationConfig() const CADET_NOEXCEPT { return false; }

	virtual bool configure(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		if (paramProvider.exists("CRY_MODE"))
			_mode.SetMode(paramProvider.getInt("CRY_MODE"));
		else
			_mode.SetMode(1); // For backwards compatibility: PBM as default

		if (_mode.hasPBM())
		{
			if (!paramProvider.exists("CRY_NUCLEI_MASS_DENSITY"))
						throw InvalidParameterException("Crystallization mode specified as " + std::to_string(_mode.getMode()) + ", i.e. crystallization including population mass balance, but field CRY_NUCLEI_MASS_DENSITY was not specified");

			// Comp 0 is substrate, last comp is equilibrium
			_nBins = _nComp - 2;

			if (_nBins < 1)
				throw InvalidParameterException("Expected at least 3 components (got " + std::to_string(_nComp) + ")");
		}
		else
		{
			_nBins = _nComp;

			if (_mode.hasAggregation())
			{
				if (!paramProvider.exists("CRY_AGGREGATION_RATE_CONSTANT"))
					throw InvalidParameterException("Crystallization mode specified as " + std::to_string(_mode.getMode()) + ", i.e. crystallization including aggregation, but field CRY_AGGREGATION_RATE_CONSTANT was not specified");

				if (!paramProvider.exists("CRY_AGGREGATION_INDEX"))
					throw InvalidParameterException("Crystallization mode specified as " + std::to_string(_mode.getMode()) + ", i.e. crystallization including aggregation, but field CRY_AGGREGATION_INDEX was not specified");
			}
			if (_mode.hasFragmentation() && !paramProvider.exists("CRY_FRAGMENTATION_RATE_CONSTANT"))
				throw InvalidParameterException("Crystallization mode specified as " + std::to_string(_mode.getMode()) + ", i.e. crystallization including fragmentation, but field CRY_FRAGMENTATION_RATE_CONSTANT was not specified");
		}

		readScalarParameterOrArray(_bins, paramProvider, "CRY_BINS", 1);

		if (_bins.size() != _nBins + 1)
			throw InvalidParameterException("Expected CRY_BINS to have " + std::to_string(_nBins + 1) + " elements (got " + std::to_string(_bins.size()) + ")");

		registerParam1DArray(_parameters, _bins, [=](bool multi, unsigned int idx) { return makeParamId(hashString("CRY_BINS"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, idx, SectionIndep); });

		_binCenters.resize(_nBins);
		_binSizes.resize(_nBins);
		_binCenterDists.resize(_nBins - 1);
		updateBinCoords();

		if (_mode.hasPBM())
		{
			_nucleiMassDensity = paramProvider.getDouble("CRY_NUCLEI_MASS_DENSITY");
			_parameters[makeParamId(hashString("CRY_NUCLEI_MASS_DENSITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_nucleiMassDensity;

			_volShapeFactor = paramProvider.getDouble("CRY_VOL_SHAPE_FACTOR");
			_parameters[makeParamId(hashString("CRY_VOL_SHAPE_FACTOR"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_volShapeFactor;

			_primaryNucleationRate = paramProvider.getDouble("CRY_PRIMARY_NUCLEATION_RATE");
			_parameters[makeParamId(hashString("CRY_PRIMARY_NUCLEATION_RATE"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_primaryNucleationRate;

			_secondaryNucleationRate = paramProvider.getDouble("CRY_SECONDARY_NUCLEATION_RATE");
			_parameters[makeParamId(hashString("CRY_SECONDARY_NUCLEATION_RATE"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_secondaryNucleationRate;

			_growthRateConstant = paramProvider.getDouble("CRY_GROWTH_RATE_CONSTANT");
			_parameters[makeParamId(hashString("CRY_GROWTH_RATE_CONSTANT"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_growthRateConstant;

			_growthConstant = paramProvider.getDouble("CRY_GROWTH_CONSTANT");
			_parameters[makeParamId(hashString("CRY_GROWTH_CONSTANT"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_growthConstant;

			_growthDispersionRate = paramProvider.getDouble("CRY_GROWTH_DISPERSION_RATE");
			_parameters[makeParamId(hashString("CRY_GROWTH_DISPERSION_RATE"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_growthDispersionRate;

			const int growthSchemeOrder = paramProvider.getInt("CRY_GROWTH_SCHEME_ORDER");

			if (!(growthSchemeOrder == 1 || growthSchemeOrder == 2 || growthSchemeOrder == 3 || growthSchemeOrder == 4))
				throw InvalidParameterException("CRY_GROWTH_SCHEME_ORDER needs to be an int between [1, 4]");

			_growthSchemeOrder = static_cast<detail::PBMReconstruction>(growthSchemeOrder - 1);

			_a = paramProvider.getDouble("CRY_A");
			_parameters[makeParamId(hashString("CRY_A"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_a;

			_b = paramProvider.getDouble("CRY_B");
			_parameters[makeParamId(hashString("CRY_B"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_b;

			_g = paramProvider.getDouble("CRY_G");
			_parameters[makeParamId(hashString("CRY_G"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_g;

			_p = paramProvider.getDouble("CRY_P");
			_parameters[makeParamId(hashString("CRY_P"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_p;

			_k = paramProvider.getDouble("CRY_K");
			_parameters[makeParamId(hashString("CRY_K"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_k;

			_u = paramProvider.getDouble("CRY_U");
			_parameters[makeParamId(hashString("CRY_U"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_u;

			// clear scheme coefficients
			if (_reconstruction)
			{
				delete _reconstruction;
				_reconstruction = nullptr;
				delete _jacParams;
				_jacParams = nullptr;
			}

			if (_growthSchemeOrder == detail::PBMReconstruction::HRKoren)
			{
				_reconstruction = new HRKorenParams<active>(_binSizes);
				_jacParams = new JacobianParamsHRKoren;
			}
			else if (_growthSchemeOrder == detail::PBMReconstruction::WENO23)
			{
				_reconstruction = new WENO23Params<active>(_binSizes);
				_jacParams = new JacobianParamsWENO23;
			}
			else if (_growthSchemeOrder == detail::PBMReconstruction::WENO35)
			{
				_reconstruction = new WENO35Params<active>(_binSizes);
				_jacParams = new JacobianParamsWENO35;
			}
			else
			{
				_reconstruction = new ReconstructionParams<active>();
				_jacParams = new JacobianParamsBase;
			}
		}
		if (_mode.hasAggregation())
		{
			_aggregationRateConstant = paramProvider.getDouble("CRY_AGGREGATION_RATE_CONSTANT");
			_parameters[makeParamId(hashString("CRY_AGGREGATION_RATE_CONSTANT"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_aggregationRateConstant;

			_aggregationIndex = paramProvider.getInt("CRY_AGGREGATION_INDEX");
			if (_aggregationIndex < 0 || _aggregationIndex > 4)
				throw InvalidParameterException("CRY_AGGREGATION_INDEX needs to be an integer in [0, 4]");

			if (_agg) // clear scheme coefficients
			{
				delete _agg;
				_agg = nullptr;
			}

			if (_aggregationRateConstant != 0.0)
				_agg = new detail::AggCoefficients(_binCenters, _bins, _binSizes);
		}
		if (_mode.hasFragmentation())
		{
			_fragRateConstant = paramProvider.getDouble("CRY_FRAGMENTATION_RATE_CONSTANT");
			_parameters[makeParamId(hashString("CRY_FRAGMENTATION_RATE_CONSTANT"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_fragRateConstant;

			_fragKernelGamma = paramProvider.getDouble("CRY_FRAGMENTATION_KERNEL_GAMMA");
			if (_fragKernelGamma <= 1.0)
				throw InvalidParameterException("CRY_FRAGMENTATION_KERNEL_GAMMA needs to be <= 1.0");
			_parameters[makeParamId(hashString("CRY_FRAGMENTATION_KERNEL_GAMMA"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_fragKernelGamma;

			_fragSelectionFunctionAlpha = paramProvider.getDouble("CRY_FRAGMENTATION_SELECTION_FUNCTION_ALPHA");
			_parameters[makeParamId(hashString("CRY_FRAGMENTATION_SELECTION_FUNCTION_ALPHA"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_fragSelectionFunctionAlpha;

			if (_frag) // clear scheme coefficients
			{
				delete _frag;
				_frag = nullptr;
			}

			if (_fragRateConstant != 0.0)
				_frag = new detail::FragCoefficients(_binCenters, _bins, _fragKernelGamma);
		}

		return true;
	}

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		_nComp = nComp;

		return true;
	}

	std::unordered_map<ParameterId, double> getAllParameterValues() const
	{
		std::unordered_map<ParameterId, double> data;
		std::transform(_parameters.begin(), _parameters.end(), std::inserter(data, data.end()),
			[](const std::pair<const ParameterId, active*>& p) { return std::make_pair(p.first, static_cast<double>(*p.second)); });
		return data;
	}

	bool hasParameter(const ParameterId& pId) const
	{
		return _parameters.find(pId) != _parameters.end();
	}

	bool setParameter(const ParameterId& pId, int value) { return false; }
	bool setParameter(const ParameterId& pId, bool value) { return false; }

	bool setParameter(const ParameterId& pId, double value)
	{
		auto paramHandle = _parameters.find(pId);
		if (paramHandle != _parameters.end())
		{
			paramHandle->second->setValue(value);

			// TODO: This does not handle parameter sensitivities wrt. to bin size
			if (pId.name == hashString("CRY_BINS"))
				updateBinCoords();

			return true;
		}

		return false;
	}

	active* getParameter(const ParameterId& pId)
	{
		auto paramHandle = _parameters.find(pId);
		if (paramHandle != _parameters.end())
		{
			return paramHandle->second;
		}

		return nullptr;
	}

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return false; }
	virtual bool requiresWorkspace() const CADET_NOEXCEPT { return false; }
	virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
	{
		return 0;
	}

	virtual unsigned int numReactionsLiquid() const CADET_NOEXCEPT { return 1; }
	virtual unsigned int numReactionsCombined() const CADET_NOEXCEPT { return 1; }

	CADET_DYNAMICREACTIONMODEL_BOILERPLATE

protected:

	template<typename ParamType>
	struct ReconstructionParams
	{
		ParamType v_g = 0.0;
		ParamType k_g_times_s_g = 0.0;
		ParamType B_0 = 0.0;

		ReconstructionParams(ParamType v_g = 0.0, ParamType k_g_times_s_g = 0.0, ParamType B_0 = 0.0) : v_g(v_g), k_g_times_s_g(k_g_times_s_g), B_0(B_0) {};

		void updateParams(ParamType k_g_times_s_g_, ParamType B_0_)
		{
			k_g_times_s_g = k_g_times_s_g_;
			B_0 = B_0_;
		}

		// Create a double view of this struct
		template <typename U>
		ReconstructionParams<U>& as() {
			//static_assert(std::is_same_v<U, double>, "Only double views are supported");
			//return reinterpret_cast<ReconstructionParams<U>&>(*this);
			return *ReconstructionParams<U>(
				static_cast<U>(v_g),
				static_cast<U>(k_g_times_s_g),
				static_cast<U>(B_0)
			);
		}

		// Overload for const objects
		template <typename U>
		const ReconstructionParams<U>& as() const {
			//static_assert(std::is_same_v<U, double>, "Only double views are supported");
			return reinterpret_cast<const ReconstructionParams<U>&>(*this);
		}
	};

	template<typename ParamType>
	struct HRKorenParams : public ReconstructionParams<ParamType>
	{
		ParamType r_x_i = 0.0;
		ParamType phi = 0.0;
		ParamType F_i = 0.0;

		std::vector<ParamType> A_coeff;
		std::vector<ParamType> R_coeff;

		HRKorenParams(const std::vector<active>& binSizes)
			: r_x_i(r_x_i), phi(phi), F_i(F_i), A_coeff(binSizes.size() - 1), R_coeff(binSizes.size() - 1)
		{
			// calculate the coefficients
			for (int i = 1; i + 1 < binSizes.size(); ++i)
			{
				const ParamType delta_xi_xip1 = static_cast<ParamType>(binSizes[i] + binSizes[i + 1]);
				A_coeff[i] = delta_xi_xip1 / (static_cast<ParamType>(binSizes[i] + binSizes[i - 1]));
				R_coeff[i] = delta_xi_xip1 / static_cast<ParamType>(binSizes[i]);
			}
		}

		// Create a double view of this struct
		template <typename U>
		HRKorenParams<U>& as() {
			return reinterpret_cast<HRKorenParams<U>&>(*this);
		}

		// Overload for const objects
		template <typename U>
		const HRKorenParams<U>& as() const {
			return reinterpret_cast<const HRKorenParams<U>&>(*this);
		}
	};

	/**
	 * Flux reconstruction scheme coefficients.For details see @cite Zhang2024
	*/
	template<typename ParamType>
	struct WENO23Params : public ReconstructionParams<ParamType>
	{
		std::vector<ParamType> q_1_right_coeff;
		std::vector<ParamType> q_0_right_coeff;
		std::vector<ParamType> C_right_coeff;
		std::vector<ParamType> IS_0_coeff;
		std::vector<ParamType> IS_1_coeff;

		ParamType IS_0 = 0.0;
		ParamType IS_1 = 0.0;
		ParamType alpha_0 = 0.0;
		ParamType alpha_1 = 0.0;
		ParamType W_0 = 0.0;
		ParamType W_1 = 0.0;
		ParamType q_0 = 0.0;
		ParamType q_1 = 0.0;

		WENO23Params(const std::vector<active>& binSizes)
			: q_1_right_coeff(binSizes.size() - 1), q_0_right_coeff(binSizes.size() - 1), C_right_coeff(binSizes.size() - 1), IS_0_coeff(binSizes.size() - 1), IS_1_coeff(binSizes.size() - 1)
		{
			// calculate the coefficients and store them, the first element is empty
			for (int i = 1; i + 1 < binSizes.size(); ++i)
			{
				const ParamType delta_sum = static_cast<ParamType>(binSizes[i] + binSizes[i - 1] + binSizes[i + 1]);
				const ParamType delta_left_sum = static_cast<ParamType>(binSizes[i] + binSizes[i - 1]);
				const ParamType delta_right_sum = static_cast<ParamType>(binSizes[i] + binSizes[i + 1]);
				q_0_right_coeff[i] = static_cast<ParamType>(binSizes[i + 1]) / delta_right_sum;
				q_1_right_coeff[i] = 1.0 + static_cast<ParamType>(binSizes[i]) / delta_left_sum;
				C_right_coeff[i] = delta_left_sum / delta_sum;
				IS_0_coeff[i] = (2.0 * static_cast<ParamType>(binSizes[i]) / delta_right_sum) * (2.0 * static_cast<ParamType>(binSizes[i]) / delta_right_sum);
				IS_1_coeff[i] = (2.0 * static_cast<ParamType>(binSizes[i]) / delta_left_sum) * (2.0 * static_cast<ParamType>(binSizes[i]) / delta_left_sum);
			}
		}

		// Create a double view of this struct
		template <typename U>
		WENO23Params<U>& as() {
			return reinterpret_cast<WENO23Params<U>&>(*this);
		}

		// Overload for const objects
		template <typename U>
		const WENO23Params<U>& as() const {
			return reinterpret_cast<const WENO23Params<U>&>(*this);
		}
	};

	template<typename ParamType>
	struct WENO35Params : public ReconstructionParams<ParamType>
	{
		std::vector<ParamType> q_2_coeff_1;
		std::vector<ParamType> q_2_coeff_2;
		std::vector<ParamType> q_1_coeff_1;
		std::vector<ParamType> q_1_coeff_2;
		std::vector<ParamType> q_0_coeff_1;
		std::vector<ParamType> q_0_coeff_2;
		std::vector<ParamType> C_0;
		std::vector<ParamType> C_1;
		std::vector<ParamType> C_2;
		std::vector<ParamType> IS_0_coeff_1;
		std::vector<ParamType> IS_0_coeff_2;
		std::vector<ParamType> IS_0_coeff_3;
		std::vector<ParamType> IS_1_coeff_1;
		std::vector<ParamType> IS_1_coeff_2;
		std::vector<ParamType> IS_1_coeff_3;
		std::vector<ParamType> IS_2_coeff_1;
		std::vector<ParamType> IS_2_coeff_2;
		std::vector<ParamType> IS_2_coeff_3;

		ParamType IS_0 = 0.0;
		ParamType IS_1 = 0.0;
		ParamType IS_2 = 0.0;
		ParamType alpha_0 = 0.0;
		ParamType alpha_1 = 0.0;
		ParamType alpha_2 = 0.0;
		ParamType W_0 = 0.0;
		ParamType W_1 = 0.0;
		ParamType W_2 = 0.0;
		ParamType q_0 = 0.0;
		ParamType q_1 = 0.0;
		ParamType q_2 = 0.0;

		// left side coefficients
		ParamType IS_0_coeff_weno3 = 0.0;
		ParamType IS_1_coeff_weno3 = 0.0;
		ParamType C_coeff1_weno3 = 0.0;
		ParamType C_coeff2_weno3 = 0.0;
		ParamType q0_coeff1_weno3 = 0.0;
		ParamType q0_coeff2_weno3 = 0.0;
		ParamType q1_coeff1_weno3 = 0.0;
		ParamType q1_coeff2_weno3 = 0.0;
		// right side coefficients
		ParamType IS_0_coeff_weno3_r = 0.0;
		ParamType IS_1_coeff_weno3_r = 0.0;
		ParamType C_coeff1_weno3_r = 0.0;
		ParamType C_coeff2_weno3_r = 0.0;
		ParamType q0_coeff1_weno3_r = 0.0;
		ParamType q0_coeff2_weno3_r = 0.0;
		ParamType q1_coeff1_weno3_r = 0.0;
		ParamType q1_coeff2_weno3_r = 0.0;

		WENO35Params(const std::vector<ParamType>& binSizes)
			: q_2_coeff_1(binSizes.size() - 2), q_2_coeff_2(binSizes.size() - 2), q_1_coeff_1(binSizes.size() - 2),
			q_1_coeff_2(binSizes.size() - 2), q_0_coeff_1(binSizes.size() - 2), q_0_coeff_2(binSizes.size() - 2),
			C_0(binSizes.size() - 2), C_1(binSizes.size() - 2), C_2(binSizes.size() - 2),
			IS_0_coeff_1(binSizes.size() - 2), IS_0_coeff_2(binSizes.size() - 2), IS_0_coeff_3(binSizes.size() - 2),
			IS_1_coeff_1(binSizes.size() - 2), IS_1_coeff_2(binSizes.size() - 2), IS_1_coeff_3(binSizes.size() - 2),
			IS_2_coeff_1(binSizes.size() - 2), IS_2_coeff_2(binSizes.size() - 2), IS_2_coeff_3(binSizes.size() - 2)
		{
			// calculate the coefficients and store them, the first and second elements are empty
			for (int i = 2; i + 2 < binSizes.size(); ++i)
			{
				const ParamType delta_sum = binSizes[i - 2] + binSizes[i - 1] + binSizes[i] + binSizes[i + 1] + binSizes[i + 2];
				const ParamType delta_sum_I0 = binSizes[i - 2] + binSizes[i - 1] + binSizes[i];
				const ParamType delta_sum_I1 = binSizes[i - 1] + binSizes[i] + binSizes[i + 1];
				const ParamType delta_sum_I2 = binSizes[i] + binSizes[i + 1] + binSizes[i + 2];
				q_0_coeff_1[i] = binSizes[i + 1] * (delta_sum_I2 - binSizes[i]) / (delta_sum_I2 - binSizes[i + 2]) / delta_sum_I2;
				q_0_coeff_2[i] = binSizes[i + 1] * binSizes[i] / delta_sum_I2 / (delta_sum_I2 - binSizes[i]);
				q_1_coeff_1[i] = binSizes[i] * (delta_sum_I1 - binSizes[i + 1]) / delta_sum_I1 / (delta_sum_I1 - binSizes[i - 1]);
				q_1_coeff_2[i] = binSizes[i] * binSizes[i + 1] / delta_sum_I1 / (delta_sum_I1 - binSizes[i + 1]);
				q_2_coeff_1[i] = binSizes[i] * (delta_sum_I0 - binSizes[i - 2]) / delta_sum_I0 / (delta_sum_I0 - binSizes[i]);
				q_2_coeff_2[i] = 1.0 + binSizes[i] / (binSizes[i - 1] + binSizes[i]) + binSizes[i] / delta_sum_I0;
				C_0[i] = delta_sum_I0 * (delta_sum_I0 - binSizes[i - 2]) / delta_sum / (delta_sum - binSizes[i - 2]);
				C_1[i] = delta_sum_I0 / delta_sum * (delta_sum_I2 - binSizes[i]) / (binSizes[i - 1] + delta_sum_I2) * (1.0 + (delta_sum - binSizes[i - 2]) / (delta_sum - binSizes[i + 2]));
				C_2[i] = binSizes[i + 1] * (binSizes[i + 1] + binSizes[i + 2]) / delta_sum / (delta_sum - binSizes[i + 2]);
				const ParamType IS_0_pre = 4.0 * (binSizes[i] / delta_sum_I2) * (binSizes[i] / delta_sum_I2);
				IS_0_coeff_1[i] = IS_0_pre * (10.0 * binSizes[i] * binSizes[i] + binSizes[i + 1] * (binSizes[i] + binSizes[i + 1])) / (binSizes[i + 1] + binSizes[i + 2]) / (binSizes[i + 1] + binSizes[i + 2]);
				IS_0_coeff_2[i] = IS_0_pre * (20.0 * binSizes[i] * binSizes[i] + 2.0 * binSizes[i + 1] * (binSizes[i] + binSizes[i + 1]) + (2.0 * binSizes[i + 1] + binSizes[i]) * delta_sum_I2) / (binSizes[i + 1] + binSizes[i + 2]) / (binSizes[i + 1] + binSizes[i]);
				IS_0_coeff_3[i] = IS_0_pre * (10.0 * binSizes[i] * binSizes[i] + (2.0 * delta_sum_I2 - binSizes[i + 2]) * (delta_sum_I2 + binSizes[i + 1])) / (binSizes[i] + binSizes[i + 1]) / (binSizes[i] + binSizes[i + 1]);
				const ParamType IS_1_pre = 4.0 * (binSizes[i] / delta_sum_I1) * (binSizes[i] / delta_sum_I1);
				IS_1_coeff_1[i] = IS_1_pre * (10.0 * binSizes[i] * binSizes[i] + binSizes[i + 1] * (binSizes[i] + binSizes[i + 1])) / (binSizes[i - 1] + binSizes[i]) / (binSizes[i - 1] + binSizes[i]);
				IS_1_coeff_2[i] = IS_1_pre * (20.0 * binSizes[i] * binSizes[i] - binSizes[i + 1] * binSizes[i - 1] - (binSizes[i] + binSizes[i + 1]) * (binSizes[i] + binSizes[i - 1])) / (binSizes[i] + binSizes[i + 1]) / (binSizes[i] + binSizes[i - 1]);
				IS_1_coeff_3[i] = IS_1_pre * (10.0 * binSizes[i] * binSizes[i] + binSizes[i - 1] * (binSizes[i - 1] + binSizes[i])) / (binSizes[i] + binSizes[i + 1]) / (binSizes[i] + binSizes[i + 1]);
				const ParamType IS_2_pre = 4.0 * (binSizes[i] / delta_sum_I0) * (binSizes[i] / delta_sum_I0);
				IS_2_coeff_1[i] = IS_2_pre * (10.0 * binSizes[i] * binSizes[i] + binSizes[i - 1] * (binSizes[i - 1] + binSizes[i])) / (binSizes[i - 2] + binSizes[i - 1]) / (binSizes[i - 2] + binSizes[i - 1]);
				IS_2_coeff_2[i] = IS_2_pre * (20.0 * binSizes[i] * binSizes[i] + 2.0 * binSizes[i - 1] * (binSizes[i - 1] + binSizes[i]) + delta_sum_I0 * (2.0 * binSizes[i - 1] + binSizes[i])) / (binSizes[i - 1] + binSizes[i]) / (binSizes[i - 2] + binSizes[i - 1]);
				IS_2_coeff_3[i] = IS_2_pre * (10.0 * binSizes[i] * binSizes[i] + (2.0 * delta_sum_I0 - binSizes[i - 2]) * (delta_sum_I0 + binSizes[i - 1])) / (binSizes[i - 1] + binSizes[i]) / (binSizes[i - 1] + binSizes[i]);
			}

			const int nBins = binSizes.size();
			IS_0_coeff_weno3 = 4.0 * binSizes[1] * binSizes[1] / (binSizes[1] + binSizes[2]) / (binSizes[1] + binSizes[2]);
			IS_1_coeff_weno3 = 4.0 * binSizes[1] * binSizes[1] / (binSizes[0] + binSizes[1]) / (binSizes[0] + binSizes[1]);
			C_coeff1_weno3 = (binSizes[0] + binSizes[1]) / (binSizes[0] + binSizes[1] + binSizes[2]);
			C_coeff2_weno3 = 1.0 - (C_coeff1_weno3);
			q0_coeff1_weno3 = binSizes[2] / (binSizes[1] + binSizes[2]);
			q0_coeff2_weno3 = 1.0 - (q0_coeff1_weno3);
			q1_coeff1_weno3 = 1.0 + binSizes[1] / (binSizes[0] + binSizes[1]);
			q1_coeff2_weno3 = 1.0 - (q1_coeff1_weno3);
			IS_0_coeff_weno3_r = 4.0 * binSizes[nBins - 2] * binSizes[nBins - 2] / (binSizes[nBins - 2] + binSizes[nBins - 1]) / (binSizes[nBins - 2] + binSizes[nBins - 1]);
			IS_1_coeff_weno3_r = 4.0 * binSizes[nBins - 2] * binSizes[nBins - 2] / (binSizes[nBins - 2] + binSizes[nBins - 3]) / (binSizes[nBins - 2] + binSizes[nBins - 3]);
			C_coeff1_weno3_r = (binSizes[nBins - 2] + binSizes[nBins - 1]) / (binSizes[nBins - 3] + binSizes[nBins - 2] + binSizes[nBins - 1]);
			C_coeff2_weno3_r = 1.0 - (C_coeff1_weno3_r);
			q0_coeff1_weno3_r = binSizes[nBins - 1] / (binSizes[nBins - 2] + binSizes[nBins - 1]);
			q0_coeff2_weno3_r = 1.0 - (q0_coeff1_weno3_r);
			q1_coeff1_weno3_r = 1.0 + binSizes[nBins - 2] / (binSizes[nBins - 2] + binSizes[nBins - 3]);
			q1_coeff2_weno3_r = 1.0 - (q1_coeff1_weno3_r);
		}

		// Create a double view of this struct
		template <typename U>
		WENO35Params<U>& as() {
			return reinterpret_cast<WENO35Params<U>&>(*this);
		}

		// Overload for const objects
		template <typename U>
		const WENO35Params<U>& as() const {
			return reinterpret_cast<const WENO35Params<U>&>(*this);
		}
	};

	struct JacobianParamsBase
	{
		// Constant parameters used in upwind scheme. Also included in parameters for HRKoren and WENO23 and WENO35 schemes
		double vG_factor;
		double dBp_dc;
		double dBp_dceq;
		double dBs_dc;
		double dBs_dceq;
		double dvG_dc_factor;
		double dvG_dceq_factor;
		double dBs_dni_factor;

		void updateParams(double vG_factor_, double dBp_dc_, double dBp_dceq_, double dBs_dc_, double dBs_dceq_, double dvG_dc_factor_, double dvG_dceq_factor_, double dBs_dni_factor_)
		{
			vG_factor = vG_factor_;
			dBp_dc = dBp_dc_;
			dBp_dceq = dBp_dceq_;
			dBs_dc = dBs_dc_;
			dBs_dceq = dBs_dceq_;
			dvG_dc_factor = dvG_dc_factor_;
			dvG_dceq_factor = dvG_dceq_factor_;
			dBs_dni_factor = dBs_dni_factor_;
		}
	};

	struct JacobianParamsHRKoren : public JacobianParamsBase
	{
		// HR Koren specific parameters
		const double epsilon = 1e-10;
		double r_x_i = 0.0;
		double phi = 0.0;
		double F_i = 0.0;
		double dvG_dc_right = 0.0;
		double dvG_dceq_right = 0.0;
		double vg_right = 0.0;

		double r_cluster = 0.0;
		double r_square_cluster = 0.0;
		double A_cluster = 0.0;
		double Ar_cluster = 0.0;
		double ni_difference = 0.0;

		double dFi_wrt_nim1 = 0.0;
		double dFi_wrt_ni = 0.0;
		double dFi_wrt_nip1 = 0.0;
	};

	struct JacobianParamsWENO23 : public JacobianParamsBase
	{
		// WENO23 related coefficients
		double dvG_dc_right = 0.0;
		double dvG_dceq_right = 0.0;
		double vg_right = 0.0;
		double IS_0_right = 0.0;
		double IS_1_right = 0.0;
		double alpha_0_right = 0.0;
		double alpha_1_right = 0.0;
		double W_0_right = 0.0;
		double W_1_right = 0.0;
		double q_0_right = 0.0;
		double q_1_right = 0.0;

		// Jacobian related coefficients
		double four_w1_B1_right = 0.0;
		double four_w0_B0_right = 0.0;
		double four_w0_w1_B1_right = 0.0;
		double four_w0_w0_B0_right = 0.0;
		double four_w1_w1_B1_right = 0.0;
		double four_w0_w1_B0_right = 0.0;

		double dw0_right_dni_m1 = 0.0;
		double dw0_right_dni = 0.0;
		double dw0_right_dni_p1 = 0.0;

		double dw1_right_dni_m1 = 0.0;
		double dw1_right_dni = 0.0;
		double dw1_right_dni_p1 = 0.0;

		double dq0_right_dni = 0.0;
		double dq0_right_dni_p1 = 0.0;

		double dq1_right_dni_m1 = 0.0;
		double dq1_right_dni = 0.0;
	};

	struct JacobianParamsWENO35 : public JacobianParamsBase
	{
		void updateParams(double vG_factor_, double dBp_dc_, double dBp_dceq_, double dBs_dc_, double dBs_dceq_, double dvG_dc_factor_, double dvG_dceq_factor_, double dBs_dni_factor_, const double* yCrystal, const int nComp, const std::vector<active>& binSizes, WENO35Params<active>& WENO5Params)
		{
			JacobianParamsBase::updateParams(vG_factor_, dBp_dc_, dBp_dceq_, dBs_dc_, dBs_dceq_, dvG_dc_factor_, dvG_dceq_factor_, dBs_dni_factor_);

			// WENO3 related coefficients, for bin N_x-2
			IS_0_weno3_right = static_cast<double>(WENO5Params.IS_0_coeff_weno3_r) * (yCrystal[nComp - 3] - yCrystal[nComp - 4]) * (yCrystal[nComp - 3] - yCrystal[nComp - 4]);
			IS_1_weno3_right = static_cast<double>(WENO5Params.IS_1_coeff_weno3_r) * (yCrystal[nComp - 4] - yCrystal[nComp - 5]) * (yCrystal[nComp - 4] - yCrystal[nComp - 5]);
			alpha_0_weno3_right = static_cast<double>(WENO5Params.C_coeff1_weno3_r) / (static_cast<double>(binSizes[nComp - 4]) + IS_0_weno3_right) / (static_cast<double>(binSizes[nComp - 4]) + IS_0_weno3_right);
			alpha_1_weno3_right = static_cast<double>(WENO5Params.C_coeff2_weno3_r) / (static_cast<double>(binSizes[nComp - 4]) + IS_1_weno3_right) / (static_cast<double>(binSizes[nComp - 4]) + IS_1_weno3_right);
			W_0_weno3_right = alpha_0_weno3_right / (alpha_0_weno3_right + alpha_1_weno3_right);
			W_1_weno3_right = 1.0 - W_0_weno3_right;
			q_0_weno3_right = static_cast<double>(WENO5Params.q0_coeff1_weno3_r) * yCrystal[nComp - 4] + static_cast<double>(WENO5Params.q0_coeff2_weno3_r) * yCrystal[nComp - 3];
			q_1_weno3_right = static_cast<double>(WENO5Params.q1_coeff1_weno3_r) * yCrystal[nComp - 4] + static_cast<double>(WENO5Params.q1_coeff2_weno3_r) * yCrystal[nComp - 5];

			// WENO3 jacobian related coefficients, for bin N_x-2
			four_w0_B0_weno3_right = 4.0 * W_0_weno3_right * static_cast<double>(WENO5Params.IS_0_coeff_weno3_r) * (yCrystal[nComp - 3] - yCrystal[nComp - 4]) / (static_cast<double>(binSizes[nComp - 4]) + IS_0_weno3_right);
			four_w1_B1_weno3_right = 4.0 * W_1_weno3_right * static_cast<double>(WENO5Params.IS_1_coeff_weno3_r) * (yCrystal[nComp - 4] - yCrystal[nComp - 5]) / (static_cast<double>(binSizes[nComp - 4]) + IS_1_weno3_right);
			four_w0_w0_B0_weno3_right = four_w0_B0_weno3_right * W_0_weno3_right;
			four_w0_w1_B1_weno3_right = four_w1_B1_weno3_right * W_0_weno3_right;
			four_w0_w1_B0_weno3_right = four_w0_B0_weno3_right * W_1_weno3_right;
			four_w1_w1_B1_weno3_right = four_w1_B1_weno3_right * W_1_weno3_right;

			dw0_right_dni_m1_weno3_right = -four_w0_w1_B1_weno3_right;
			dw0_right_dni_weno3_right = four_w0_B0_weno3_right - four_w0_w0_B0_weno3_right + four_w0_w1_B1_weno3_right;
			dw0_right_dni_p1_weno3_right = four_w0_w0_B0_weno3_right - four_w0_B0_weno3_right;

			dw1_right_dni_m1_weno3_right = four_w1_B1_weno3_right - four_w1_w1_B1_weno3_right;
			dw1_right_dni_weno3_right = -four_w1_B1_weno3_right - four_w0_w1_B0_weno3_right + four_w1_w1_B1_weno3_right;
			dw1_right_dni_p1_weno3_right = four_w0_w1_B0_weno3_right;

			dq0_right_dni_weno3_right = static_cast<double>(WENO5Params.q0_coeff1_weno3_r);
			dq0_right_dni_p1_weno3_right = static_cast<double>(WENO5Params.q0_coeff2_weno3_r);
			dq1_right_dni_m1_weno3_right = static_cast<double>(WENO5Params.q1_coeff2_weno3_r);
			dq1_right_dni_weno3_right = static_cast<double>(WENO5Params.q1_coeff1_weno3_r);

			// WENO3 related coefficients, for bin 2
			IS_0_weno3_left = static_cast<double>(WENO5Params.IS_0_coeff_weno3) * (yCrystal[2] - yCrystal[1]) * (yCrystal[2] - yCrystal[1]);
			const double IS_1_weno3_left = static_cast<double>(WENO5Params.IS_1_coeff_weno3) * (yCrystal[1] - yCrystal[0]) * (yCrystal[1] - yCrystal[0]);
			alpha_0_weno3_left = static_cast<double>(WENO5Params.C_coeff1_weno3) / (static_cast<double>(binSizes[1]) + IS_0_weno3_left) / (static_cast<double>(binSizes[1]) + IS_0_weno3_left);
			alpha_1_weno3_left = static_cast<double>(WENO5Params.C_coeff2_weno3) / (static_cast<double>(binSizes[1]) + IS_1_weno3_left) / (static_cast<double>(binSizes[1]) + IS_1_weno3_left);
			W_0_weno3_left = alpha_0_weno3_left / (alpha_0_weno3_left + alpha_1_weno3_left);
			W_1_weno3_left = 1.0 - W_0_weno3_left;
			q_0_weno3_left = static_cast<double>(WENO5Params.q0_coeff1_weno3) * yCrystal[1] + static_cast<double>(WENO5Params.q0_coeff2_weno3) * yCrystal[2];
			q_1_weno3_left = static_cast<double>(WENO5Params.q1_coeff1_weno3) * yCrystal[1] + static_cast<double>(WENO5Params.q1_coeff2_weno3) * yCrystal[0];

			// WENO3 Jacobian related coefficients, for bin 2
			four_w0_B0_weno3_left = 4.0 * W_0_weno3_left * static_cast<double>(WENO5Params.IS_0_coeff_weno3) * (yCrystal[2] - yCrystal[1]) / (static_cast<double>(binSizes[1]) + IS_0_weno3_left);
			four_w1_B1_weno3_left = 4.0 * W_1_weno3_left * static_cast<double>(WENO5Params.IS_1_coeff_weno3) * (yCrystal[1] - yCrystal[0]) / (static_cast<double>(binSizes[1]) + IS_1_weno3_left);
			four_w0_w0_B0_weno3_left = four_w0_B0_weno3_left * W_0_weno3_left;
			four_w0_w1_B1_weno3_left = four_w1_B1_weno3_left * W_0_weno3_left;
			four_w0_w1_B0_weno3_left = four_w0_B0_weno3_left * W_1_weno3_left;
			four_w1_w1_B1_weno3_left = four_w1_B1_weno3_left * W_1_weno3_left;

			dw0_right_dni_m1_weno3_left = -four_w0_w1_B1_weno3_left;
			dw0_right_dni_weno3_left = four_w0_B0_weno3_left - four_w0_w0_B0_weno3_left + four_w0_w1_B1_weno3_left;
			dw0_right_dni_p1_weno3_left = four_w0_w0_B0_weno3_left - four_w0_B0_weno3_left;

			dw1_right_dni_m1_weno3_left = four_w1_B1_weno3_left - four_w1_w1_B1_weno3_left;
			dw1_right_dni_weno3_left = -four_w1_B1_weno3_left - four_w0_w1_B0_weno3_left + four_w1_w1_B1_weno3_left;
			dw1_right_dni_p1_weno3_left = four_w0_w1_B0_weno3_left;

			dq0_right_dni_weno3_left = static_cast<double>(WENO5Params.q0_coeff1_weno3);
			dq0_right_dni_p1_weno3_left = static_cast<double>(WENO5Params.q0_coeff2_weno3);
			dq1_right_dni_m1_weno3_left = static_cast<double>(WENO5Params.q1_coeff2_weno3);
			dq1_right_dni_weno3_left = static_cast<double>(WENO5Params.q1_coeff1_weno3);
		}

		// WENO35 boundary WENO23 related coefficients, for bin N_x-2
		double IS_0_weno3_right = 0.0;
		double IS_1_weno3_right = 0.0;
		double alpha_0_weno3_right = 0.0;
		double alpha_1_weno3_right = 0.0;
		double W_0_weno3_right = 0.0;
		double W_1_weno3_right = 0.0;
		double q_0_weno3_right = 0.0;
		double q_1_weno3_right = 0.0;

		// WENO35 boundary WENO23 jacobian related coefficients, for bin N_x-2
		double four_w0_B0_weno3_right = 0.0;
		double four_w1_B1_weno3_right = 0.0;
		double four_w0_w0_B0_weno3_right = 0.0;
		double four_w0_w1_B1_weno3_right = 0.0;
		double four_w0_w1_B0_weno3_right = 0.0;
		double four_w1_w1_B1_weno3_right = 0.0;

		double dw0_right_dni_m1_weno3_right = 0.0;
		double dw0_right_dni_weno3_right = 0.0;
		double dw0_right_dni_p1_weno3_right = 0.0;

		double dw1_right_dni_m1_weno3_right = 0.0;
		double dw1_right_dni_weno3_right = 0.0;
		double dw1_right_dni_p1_weno3_right = 0.0;

		double dq0_right_dni_weno3_right = 0.0;
		double dq0_right_dni_p1_weno3_right = 0.0;
		double dq1_right_dni_m1_weno3_right = 0.0;
		double dq1_right_dni_weno3_right = 0.0;

		// WENO35 boundary WENO23 related coefficients, for bin 2
		double IS_0_weno3_left = 0.0;
		double IS_1_weno3_left = 0.0;
		double alpha_0_weno3_left = 0.0;
		double alpha_1_weno3_left = 0.0;
		double W_0_weno3_left = 0.0;
		double W_1_weno3_left = 0.0;
		double q_0_weno3_left = 0.0;
		double q_1_weno3_left = 0.0;

		// WENO35 boundary WENO23 Jacobian related coefficients, for bin 2
		double four_w0_B0_weno3_left;
		double four_w1_B1_weno3_left;
		double four_w0_w0_B0_weno3_left;
		double four_w0_w1_B1_weno3_left;
		double four_w0_w1_B0_weno3_left;
		double four_w1_w1_B1_weno3_left;

		double dw0_right_dni_m1_weno3_left;
		double dw0_right_dni_weno3_left;
		double dw0_right_dni_p1_weno3_left;

		double dw1_right_dni_m1_weno3_left;
		double dw1_right_dni_weno3_left;
		double dw1_right_dni_p1_weno3_left;

		double dq0_right_dni_weno3_left;
		double dq0_right_dni_p1_weno3_left;
		double dq1_right_dni_m1_weno3_left;
		double dq1_right_dni_weno3_left;

		// WENO35 related coefficients
		double IS_0 = 0.0;
		double IS_1 = 0.0;
		double IS_2 = 0.0;
		double alpha_0 = 0.0;
		double alpha_1 = 0.0;
		double alpha_2 = 0.0;
		double W_0 = 0.0;
		double W_1 = 0.0;
		double W_2 = 0.0;
		double q_0 = 0.0;
		double q_1 = 0.0;
		double q_2 = 0.0;

		double dvG_dc_right = 0.0;
		double dvG_dceq_right = 0.0;
		double vg_right = 0.0;

		// WENO35 Jacobian related coefficients
		double dIS0_wrt_dni_p2 = 0.0;
		double dIS0_wrt_dni_p1 = 0.0;
		double dIS0_wrt_dni = 0.0;
		double dIS1_wrt_dni_p1 = 0.0;
		double dIS1_wrt_dni = 0.0;
		double dIS1_wrt_dni_m1 = 0.0;
		double dIS2_wrt_dni = 0.0;
		double dIS2_wrt_dni_m1 = 0.0;
		double dIS2_wrt_dni_m2 = 0.0;
		double dq0_wrt_dni_p2 = 0.0;
		double dq0_wrt_dni_p1 = 0.0;
		double dq0_wrt_dni = 0.0;
		double dq1_wrt_dni_p1 = 0.0;
		double dq1_wrt_dni = 0.0;
		double dq1_wrt_dni_m1 = 0.0;
		double dq2_wrt_dni = 0.0;
		double dq2_wrt_dni_m1 = 0.0;
		double dq2_wrt_dni_m2 = 0.0;
		double w0_w2_is2 = 0.0;
		double w0_w1_is1 = 0.0;
		double w0_w0_is0 = 0.0;
		double w0_is0 = 0.0;
		double dw0_wrt_dni_m2 = 0.0;
		double dw0_wrt_dni_m1 = 0.0;
		double dw0_wrt_dni = 0.0;
		double dw0_wrt_dni_p1 = 0.0;
		double dw0_wrt_dni_p2 = 0.0;
		double w1_w2_is2 = 0.0;
		double w1_w1_is1 = 0.0;
		double w1_w0_is0 = 0.0;
		double w1_is1 = 0.0;
		double dw1_wrt_dni_m2 = 0.0;
		double dw1_wrt_dni_m1 = 0.0;
		double dw1_wrt_dni = 0.0;
		double dw1_wrt_dni_p1 = 0.0;
		double dw1_wrt_dni_p2 = 0.0;
		double w2_is2 = 0.0;
		double w2_w2_is2 = 0.0;
		double w2_w1_is1 = 0.0;
		double w2_w0_is0 = 0.0;
		double dw2_wrt_dni_m2 = 0.0;
		double dw2_wrt_dni_m1 = 0.0;
		double dw2_wrt_dni = 0.0;
		double dw2_wrt_dni_p1 = 0.0;
		double dw2_wrt_dni_p2 = 0.0;
	};

	std::unordered_map<ParameterId, active*> _parameters; //!< Map used to translate ParameterIds to actual variables
	detail::CrystallizationMode _mode; //!< Crystallization mode, i.e. specification of considered effects
	int _nComp; //!< Number of components
	int _nBins; //!< Number of crystal size bins

	std::vector<active> _bins;
	std::vector<active> _binCenters;
	std::vector<active> _binSizes;

	// PBM parameters
	std::vector<active> _binCenterDists;
	active _nucleiMassDensity; //!< rho
	active _volShapeFactor; //!< k_v
	active _primaryNucleationRate; //!< k_p
	active _secondaryNucleationRate; //!< k_b
	active _growthRateConstant; //!< k_g
	active _growthConstant; //!< gamma
	active _growthDispersionRate; //!< D_g
	detail::PBMReconstruction _growthSchemeOrder; // reconstruction type for PBM aka order of the growth scheme
	active _a; //!< System constant
	active _b; //!< System constant
	active _g; //!< System constant
	active _p; //!< System constant
	active _k; //!< System constant
	active _u; //!< System constant

	JacobianParamsBase* _jacParams; //!< Parameters used in Jacobian computation, dependent on the discretization scheme
	ReconstructionParams<active>* _reconstruction;

	// Fragmentation parameters
	active _fragRateConstant; // constant fragmentation rate constant
	active _fragKernelGamma; // gamma in the fragmentation kernel
	active _fragSelectionFunctionAlpha; // alpha in the selection function
	detail::FragCoefficients* _frag;

	// Aggregation parameters
	int _aggregationIndex; // determines which kernel to use
	active _aggregationRateConstant; // aggregation rate constant
	detail::AggCoefficients* _agg;

	void updateBinCoords() CADET_NOEXCEPT
	{
		for (int i = 0; i < _nBins; ++i)
		{
			if (cadet_likely(i + 1 < _nBins))
				_binCenterDists[i] = 0.5 * (_bins[i + 2] - _bins[i]);

			_binCenters[i] = 0.5 * (_bins[i] + _bins[i + 1]);
			_binSizes[i] = _bins[i + 1] - _bins[i];
		}
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	void upwindKernel(StateType const* yCrystal, ResidualType* resCrystal, const FactorType& factor, ReconstructionParams<active>& upwindParams, const int i) const
	{
		// Flux through left face
		if (cadet_likely((i > 0) && (i + 1 < _nBins)))
		{
			// flux through the left face, note that v_g is set from the last call
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * (static_cast<ResidualType>(upwindParams.v_g) * yCrystal[i - 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// flux through the right face
			upwindParams.v_g = static_cast<ResidualType>(upwindParams.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * (static_cast<ResidualType>(upwindParams.v_g) * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		else if (i == 0)
		{
			// Left boundary condition
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(upwindParams.B_0);
			// upwind
			upwindParams.v_g = static_cast<ResidualType>(upwindParams.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(upwindParams.v_g) * yCrystal[i] - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		else
		{
			// first order approximation
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(upwindParams.v_g) * yCrystal[i - 1] - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// no flux
		}
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	void HRKorenKernel(StateType const* yCrystal, ResidualType* resCrystal, const FactorType& factor, HRKorenParams<active>& HRKoren, const int i) const
	{
		if (cadet_likely((i > 1) && (i + 1 < _nBins)))
		{
			// Flux through left face, modified van Leer flux limiter
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(HRKoren.v_g) * static_cast<ResidualType>(HRKoren.F_i) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// Flux through left face, modified van Leer flux limiter, update r_x_i, R_i, and growth rate
			HRKoren.r_x_i = static_cast<ParamType>(HRKoren.A_coeff[i]) * (yCrystal[i] - yCrystal[i - 1] + 1e-10) / (yCrystal[i + 1] - yCrystal[i] + 1e-10);
			if (cadet_likely(static_cast<ResidualType>(HRKoren.r_x_i) > 0))
			{
				HRKoren.phi = static_cast<ResidualType>(HRKoren.r_x_i) / (static_cast<ParamType>(HRKoren.R_coeff[i]) - 1.0 + static_cast<ResidualType>(HRKoren.r_x_i));
			}
			else
			{
				HRKoren.phi = 0.0;
			}
			HRKoren.F_i = yCrystal[i] + static_cast<ResidualType>(HRKoren.phi) * (yCrystal[i + 1] - yCrystal[i]);
			HRKoren.v_g = static_cast<ResidualType>(HRKoren.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(HRKoren.v_g) * static_cast<ResidualType>(HRKoren.F_i) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		else if (i == 1)
		{
			// upwind to F_{1/2} and the diffusion
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(HRKoren.v_g) * yCrystal[i - 1] - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// HR scheme applied to F_{1+1/2}
			HRKoren.r_x_i = static_cast<ParamType>(HRKoren.A_coeff[i]) * (yCrystal[i] - yCrystal[i - 1] + 1e-10) / (yCrystal[i + 1] - yCrystal[i] + 1e-10);
			if (cadet_likely(static_cast<ResidualType>(HRKoren.r_x_i) > 0))
			{
				HRKoren.phi = static_cast<ResidualType>(HRKoren.r_x_i) / (static_cast<ParamType>(HRKoren.R_coeff[i]) - 1.0 + static_cast<ResidualType>(HRKoren.r_x_i));
			}
			else
			{
				HRKoren.phi = 0.0;
			}
			HRKoren.F_i = yCrystal[i] + static_cast<ResidualType>(HRKoren.phi) * (yCrystal[i + 1] - yCrystal[i]);
			HRKoren.v_g = static_cast<ResidualType>(HRKoren.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(HRKoren.v_g) * static_cast<ResidualType>(HRKoren.F_i) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		else if (i + 1 == _nBins)
		{
			// HR scheme applied to the influx of the last bin 
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(HRKoren.v_g) * static_cast<ResidualType>(HRKoren.F_i) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// no flux
		}
		else
		{
			// left boundary condition
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(HRKoren.B_0);
			// upwind to F_1
			HRKoren.v_g = static_cast<ResidualType>(HRKoren.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(HRKoren.v_g) * yCrystal[i] - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	void WENO23Kernel(StateType const* yCrystal, ResidualType* resCrystal, const FactorType& factor, WENO23Params<active>& WENO23Params, const int i) const
	{
		if (cadet_likely((i > 1) && (i + 1 < _nBins)))
		{
			// flux through left face. W_0, q_0, W_1 and q_1 are coming from the right face of bin (i-1).
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO23Params.v_g) * (static_cast<ResidualType>(WENO23Params.W_0) * static_cast<ResidualType>(WENO23Params.q_0) + static_cast<ResidualType>(WENO23Params.W_1) * static_cast<ResidualType>(WENO23Params.q_1)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// flux through right face, update IS, alpha, W, q and growth rate
			WENO23Params.IS_0 = static_cast<ParamType>(WENO23Params.IS_0_coeff[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]);
			WENO23Params.IS_1 = static_cast<ParamType>(WENO23Params.IS_1_coeff[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
			WENO23Params.alpha_0 = static_cast<ParamType>(WENO23Params.C_right_coeff[i]) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO23Params.IS_0)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO23Params.IS_0));
			WENO23Params.alpha_1 = (1.0 - static_cast<ParamType>(WENO23Params.C_right_coeff[i])) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO23Params.IS_1)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO23Params.IS_1));
			WENO23Params.W_0 = static_cast<ResidualType>(WENO23Params.alpha_0) / (static_cast<ResidualType>(WENO23Params.alpha_0) + static_cast<ResidualType>(WENO23Params.alpha_1));
			WENO23Params.W_1 = 1.0 - static_cast<ResidualType>(WENO23Params.W_0);
			WENO23Params.q_0 = static_cast<ParamType>(WENO23Params.q_0_right_coeff[i]) * yCrystal[i] + (1.0 - static_cast<ParamType>(WENO23Params.q_0_right_coeff[i])) * yCrystal[i + 1];
			WENO23Params.q_1 = static_cast<ParamType>(WENO23Params.q_1_right_coeff[i]) * yCrystal[i] + (1.0 - static_cast<ParamType>(WENO23Params.q_1_right_coeff[i])) * yCrystal[i - 1];
			WENO23Params.v_g = static_cast<ResidualType>(WENO23Params.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO23Params.v_g) * (static_cast<ResidualType>(WENO23Params.W_0) * static_cast<ResidualType>(WENO23Params.q_0) + static_cast<ResidualType>(WENO23Params.W_1) * static_cast<ResidualType>(WENO23Params.q_1)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		// boundary condition
		else if (i + 1 == _nBins)
		{
			// weno3 applied to the influx of the last bin
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO23Params.v_g) * (static_cast<ResidualType>(WENO23Params.W_0) * static_cast<ResidualType>(WENO23Params.q_0) + static_cast<ResidualType>(WENO23Params.W_1) * static_cast<ResidualType>(WENO23Params.q_1)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// no flux, regularity boundary condition
		}
		else if (i == 1)
		{
			// upwind applied to F_{1/2}
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO23Params.v_g) * yCrystal[i - 1] - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// weno23 applied to F_{1+1/2}
			WENO23Params.IS_0 = static_cast<ParamType>(WENO23Params.IS_0_coeff[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]);
			WENO23Params.IS_1 = static_cast<ParamType>(WENO23Params.IS_1_coeff[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
			WENO23Params.alpha_0 = static_cast<ParamType>(WENO23Params.C_right_coeff[i]) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO23Params.IS_0)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO23Params.IS_0));
			WENO23Params.alpha_1 = (1.0 - static_cast<ParamType>(WENO23Params.C_right_coeff[i])) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO23Params.IS_1)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO23Params.IS_1));
			WENO23Params.W_0 = static_cast<ResidualType>(WENO23Params.alpha_0) / (static_cast<ResidualType>(WENO23Params.alpha_0) + static_cast<ResidualType>(WENO23Params.alpha_1));
			WENO23Params.W_1 = 1.0 - static_cast<ResidualType>(WENO23Params.W_0);
			WENO23Params.q_0 = static_cast<ParamType>(WENO23Params.q_0_right_coeff[i]) * yCrystal[i] + (1.0 - static_cast<ParamType>(WENO23Params.q_0_right_coeff[i])) * yCrystal[i + 1];
			WENO23Params.q_1 = static_cast<ParamType>(WENO23Params.q_1_right_coeff[i]) * yCrystal[i] + (1.0 - static_cast<ParamType>(WENO23Params.q_1_right_coeff[i])) * yCrystal[i - 1];
			WENO23Params.v_g = static_cast<ResidualType>(WENO23Params.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO23Params.v_g) * (static_cast<ResidualType>(WENO23Params.W_0) * static_cast<ResidualType>(WENO23Params.q_0) + static_cast<ResidualType>(WENO23Params.W_1) * static_cast<ResidualType>(WENO23Params.q_1)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		else
		{
			// nucleation boundary condition
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO23Params.B_0);
			// upwind
			WENO23Params.v_g = static_cast<ResidualType>(WENO23Params.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO23Params.v_g) * yCrystal[i] - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	void WENO35Kernel(StateType const* yCrystal, ResidualType* resCrystal, const FactorType& factor, WENO35Params<active>& WENO35Params, const int i) const
	{
		if (cadet_likely((i > 2) && (i + 2 < _nBins)))
		{
			// flux through the left face 
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.v_g) * (static_cast<ResidualType>(WENO35Params.W_0) * static_cast<ResidualType>(WENO35Params.q_0) + static_cast<ResidualType>(WENO35Params.W_1) * static_cast<ResidualType>(WENO35Params.q_1) + static_cast<ResidualType>(WENO35Params.W_2) * static_cast<ResidualType>(WENO35Params.q_2)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// flux through the right face
			WENO35Params.IS_0 = static_cast<ParamType>(WENO35Params.IS_0_coeff_1[i]) * (yCrystal[i + 2] - yCrystal[i + 1]) * (yCrystal[i + 2] - yCrystal[i + 1]) + static_cast<ParamType>(WENO35Params.IS_0_coeff_2[i]) * (yCrystal[i + 2] - yCrystal[i + 1]) * (yCrystal[i] - yCrystal[i + 1]) + static_cast<ParamType>(WENO35Params.IS_0_coeff_3[i]) * (yCrystal[i] - yCrystal[i + 1]) * (yCrystal[i] - yCrystal[i + 1]);
			WENO35Params.IS_1 = static_cast<ParamType>(WENO35Params.IS_1_coeff_1[i]) * (yCrystal[i - 1] - yCrystal[i]) * (yCrystal[i - 1] - yCrystal[i]) + static_cast<ParamType>(WENO35Params.IS_1_coeff_2[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i - 1] - yCrystal[i]) + static_cast<ParamType>(WENO35Params.IS_1_coeff_3[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]);
			WENO35Params.IS_2 = static_cast<ParamType>(WENO35Params.IS_2_coeff_1[i]) * (yCrystal[i - 2] - yCrystal[i - 1]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(WENO35Params.IS_2_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(WENO35Params.IS_2_coeff_3[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
			WENO35Params.alpha_0 = static_cast<ParamType>(WENO35Params.C_0[i]) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_0)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_0));
			WENO35Params.alpha_1 = static_cast<ParamType>(WENO35Params.C_1[i]) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_1)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_1));
			WENO35Params.alpha_2 = static_cast<ParamType>(WENO35Params.C_2[i]) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_2)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_2));
			WENO35Params.W_0 = static_cast<ResidualType>(WENO35Params.alpha_0) / (static_cast<ResidualType>(WENO35Params.alpha_0) + static_cast<ResidualType>(WENO35Params.alpha_1) + static_cast<ResidualType>(WENO35Params.alpha_2));
			WENO35Params.W_1 = static_cast<ResidualType>(WENO35Params.alpha_1) / (static_cast<ResidualType>(WENO35Params.alpha_0) + static_cast<ResidualType>(WENO35Params.alpha_1) + static_cast<ResidualType>(WENO35Params.alpha_2));
			WENO35Params.W_2 = 1.0 - static_cast<ResidualType>(WENO35Params.W_0) - static_cast<ResidualType>(WENO35Params.W_1);
			WENO35Params.q_0 = yCrystal[i + 1] + static_cast<ParamType>(WENO35Params.q_0_coeff_1[i]) * (yCrystal[i] - yCrystal[i + 1]) + static_cast<ParamType>(WENO35Params.q_0_coeff_2[i]) * (yCrystal[i + 1] - yCrystal[i + 2]);
			WENO35Params.q_1 = yCrystal[i] + static_cast<ParamType>(WENO35Params.q_1_coeff_1[i]) * (yCrystal[i + 1] - yCrystal[i]) + static_cast<ParamType>(WENO35Params.q_1_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]);
			WENO35Params.q_2 = yCrystal[i - 1] + static_cast<ParamType>(WENO35Params.q_2_coeff_1[i]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(WENO35Params.q_2_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]);
			WENO35Params.v_g = static_cast<ResidualType>(WENO35Params.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.v_g) * (static_cast<ResidualType>(WENO35Params.W_0) * static_cast<ResidualType>(WENO35Params.q_0) + static_cast<ResidualType>(WENO35Params.W_1) * static_cast<ResidualType>(WENO35Params.q_1) + static_cast<ResidualType>(WENO35Params.W_2) * static_cast<ResidualType>(WENO35Params.q_2)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		else if (i == 2)
		{
			// weno23 applied to F_{2-1/2}
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.v_g) * (static_cast<ResidualType>(WENO35Params.W_0) * static_cast<ResidualType>(WENO35Params.q_0) + static_cast<ResidualType>(WENO35Params.W_1) * static_cast<ResidualType>(WENO35Params.q_1)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// weno35 applied to F_{2+1/2}
			WENO35Params.IS_0 = static_cast<ParamType>(WENO35Params.IS_0_coeff_1[i]) * (yCrystal[i + 2] - yCrystal[i + 1]) * (yCrystal[i + 2] - yCrystal[i + 1]) + static_cast<ParamType>(WENO35Params.IS_0_coeff_2[i]) * (yCrystal[i + 2] - yCrystal[i + 1]) * (yCrystal[i] - yCrystal[i + 1]) + static_cast<ParamType>(WENO35Params.IS_0_coeff_3[i]) * (yCrystal[i] - yCrystal[i + 1]) * (yCrystal[i] - yCrystal[i + 1]);
			WENO35Params.IS_1 = static_cast<ParamType>(WENO35Params.IS_1_coeff_1[i]) * (yCrystal[i - 1] - yCrystal[i]) * (yCrystal[i - 1] - yCrystal[i]) + static_cast<ParamType>(WENO35Params.IS_1_coeff_2[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i - 1] - yCrystal[i]) + static_cast<ParamType>(WENO35Params.IS_1_coeff_3[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]);
			WENO35Params.IS_2 = static_cast<ParamType>(WENO35Params.IS_2_coeff_1[i]) * (yCrystal[i - 2] - yCrystal[i - 1]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(WENO35Params.IS_2_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(WENO35Params.IS_2_coeff_3[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
			WENO35Params.alpha_0 = static_cast<ParamType>(WENO35Params.C_0[i]) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_0)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_0));
			WENO35Params.alpha_1 = static_cast<ParamType>(WENO35Params.C_1[i]) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_1)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_1));
			WENO35Params.alpha_2 = static_cast<ParamType>(WENO35Params.C_2[i]) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_2)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_2));
			WENO35Params.W_0 = static_cast<ResidualType>(WENO35Params.alpha_0) / (static_cast<ResidualType>(WENO35Params.alpha_0) + static_cast<ResidualType>(WENO35Params.alpha_1) + static_cast<ResidualType>(WENO35Params.alpha_2));
			WENO35Params.W_1 = static_cast<ResidualType>(WENO35Params.alpha_1) / (static_cast<ResidualType>(WENO35Params.alpha_0) + static_cast<ResidualType>(WENO35Params.alpha_1) + static_cast<ResidualType>(WENO35Params.alpha_2));
			WENO35Params.W_2 = 1.0 - static_cast<ResidualType>(WENO35Params.W_0) - static_cast<ResidualType>(WENO35Params.W_1);
			WENO35Params.q_0 = yCrystal[i + 1] + static_cast<ParamType>(WENO35Params.q_0_coeff_1[i]) * (yCrystal[i] - yCrystal[i + 1]) + static_cast<ParamType>(WENO35Params.q_0_coeff_2[i]) * (yCrystal[i + 1] - yCrystal[i + 2]);
			WENO35Params.q_1 = yCrystal[i] + static_cast<ParamType>(WENO35Params.q_1_coeff_1[i]) * (yCrystal[i + 1] - yCrystal[i]) + static_cast<ParamType>(WENO35Params.q_1_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]);
			WENO35Params.q_2 = yCrystal[i - 1] + static_cast<ParamType>(WENO35Params.q_2_coeff_1[i]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(WENO35Params.q_2_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]);
			WENO35Params.v_g = static_cast<ResidualType>(WENO35Params.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.v_g) * (static_cast<ResidualType>(WENO35Params.W_0) * static_cast<ResidualType>(WENO35Params.q_0) + static_cast<ResidualType>(WENO35Params.W_1) * static_cast<ResidualType>(WENO35Params.q_1) + static_cast<ResidualType>(WENO35Params.W_2) * static_cast<ResidualType>(WENO35Params.q_2)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		else if (i == 1)
		{
			// upwind applied to F_{1-1/2}
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.v_g) * yCrystal[i - 1] - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// weno23 applied to F_{1+1/2}
			WENO35Params.IS_0 = static_cast<ParamType>(WENO35Params.IS_0_coeff_weno3) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]);
			WENO35Params.IS_1 = static_cast<ParamType>(WENO35Params.IS_1_coeff_weno3) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
			WENO35Params.alpha_0 = static_cast<ParamType>(WENO35Params.C_coeff1_weno3) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_0)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_0));
			WENO35Params.alpha_1 = static_cast<ParamType>(WENO35Params.C_coeff2_weno3) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_1)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_1));
			WENO35Params.W_0 = static_cast<ResidualType>(WENO35Params.alpha_0) / (static_cast<ResidualType>(WENO35Params.alpha_0) + static_cast<ResidualType>(WENO35Params.alpha_1));
			WENO35Params.W_1 = 1.0 - static_cast<ResidualType>(WENO35Params.W_0);
			WENO35Params.q_0 = static_cast<ParamType>(WENO35Params.q0_coeff1_weno3) * yCrystal[i] + static_cast<ParamType>(WENO35Params.q0_coeff2_weno3) * yCrystal[i + 1];
			WENO35Params.q_1 = static_cast<ParamType>(WENO35Params.q1_coeff1_weno3) * yCrystal[i] + static_cast<ParamType>(WENO35Params.q1_coeff2_weno3) * yCrystal[i - 1];
			WENO35Params.v_g = static_cast<ResidualType>(WENO35Params.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.v_g) * (static_cast<ResidualType>(WENO35Params.W_0) * static_cast<ResidualType>(WENO35Params.q_0) + static_cast<ResidualType>(WENO35Params.W_1) * static_cast<ResidualType>(WENO35Params.q_1)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		else if (i == 0)
		{
			// boundary condition
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.B_0);
			// upwind to F_{1/2}
			WENO35Params.v_g = static_cast<ResidualType>(WENO35Params.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.v_g) * yCrystal[i] - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		else if (i + 2 == _nBins)
		{
			// weno35 applied to F_{nx-2-1/2}
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.v_g) * (static_cast<ResidualType>(WENO35Params.W_0) * static_cast<ResidualType>(WENO35Params.q_0) + static_cast<ResidualType>(WENO35Params.W_1) * static_cast<ResidualType>(WENO35Params.q_1) + static_cast<ResidualType>(WENO35Params.W_2) * static_cast<ResidualType>(WENO35Params.q_2)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// weno23 applied to F_{nx-2+1/2}
			WENO35Params.IS_0 = static_cast<ParamType>(WENO35Params.IS_0_coeff_weno3_r) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]);
			WENO35Params.IS_1 = static_cast<ParamType>(WENO35Params.IS_1_coeff_weno3_r) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
			WENO35Params.alpha_0 = static_cast<ParamType>(WENO35Params.C_coeff1_weno3_r) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_0)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_0));
			WENO35Params.alpha_1 = static_cast<ParamType>(WENO35Params.C_coeff2_weno3_r) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_1)) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ResidualType>(WENO35Params.IS_1));
			WENO35Params.W_0 = static_cast<ResidualType>(WENO35Params.alpha_0) / (static_cast<ResidualType>(WENO35Params.alpha_0) + static_cast<ResidualType>(WENO35Params.alpha_1));
			WENO35Params.W_1 = 1.0 - static_cast<ResidualType>(WENO35Params.W_0);
			WENO35Params.q_0 = static_cast<ParamType>(WENO35Params.q0_coeff1_weno3_r) * yCrystal[i] + static_cast<ParamType>(WENO35Params.q0_coeff2_weno3_r) * yCrystal[i + 1];
			WENO35Params.q_1 = static_cast<ParamType>(WENO35Params.q1_coeff1_weno3_r) * yCrystal[i] + static_cast<ParamType>(WENO35Params.q1_coeff2_weno3_r) * yCrystal[i - 1];
			WENO35Params.v_g = static_cast<ResidualType>(WENO35Params.k_g_times_s_g) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
			resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.v_g) * (static_cast<ResidualType>(WENO35Params.W_0) * static_cast<ResidualType>(WENO35Params.q_0) + static_cast<ResidualType>(WENO35Params.W_1) * static_cast<ResidualType>(WENO35Params.q_1)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
		}
		else
		{
			// weno23 to F_{nx-1-1/2}
			resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * static_cast<ResidualType>(WENO35Params.v_g) * (static_cast<ResidualType>(WENO35Params.W_0) * static_cast<ResidualType>(WENO35Params.q_0) + static_cast<ResidualType>(WENO35Params.W_1) * static_cast<ResidualType>(WENO35Params.q_1)) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
			// no flux BC
		}
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	void PBMKernel(StateType const* yCrystal, ResidualType* resCrystal, const FactorType& factor, const int i) const
	{
		switch (_growthSchemeOrder)
		{
		case detail::PBMReconstruction::Upwind:
		{
			upwindKernel<StateType, ResidualType, ParamType, FactorType>(yCrystal, resCrystal, factor, *_reconstruction, i);
		}
		break;
		case detail::PBMReconstruction::HRKoren:
		{
			HRKorenKernel<StateType, ResidualType, ParamType, FactorType>(yCrystal, resCrystal, factor, *static_cast<HRKorenParams<active>*>(_reconstruction), i);
		}
		break;
		case detail::PBMReconstruction::WENO23:
		{
			WENO23Kernel<StateType, ResidualType, ParamType, FactorType>(yCrystal, resCrystal, factor, *static_cast<WENO23Params<active>*>(_reconstruction), i);
		}
		break;
		case detail::PBMReconstruction::WENO35:
		{
			WENO35Kernel<StateType, ResidualType, ParamType, FactorType>(yCrystal, resCrystal, factor, *static_cast<WENO35Params<active>*>(_reconstruction), i);
		}
		break;
		}
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	int residualFluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, StateType const* y, ResidualType* res, const FactorType& factor, LinearBufferAllocator workSpace) const
	{
		ResidualType B_0 = 0.0;
		ResidualType k_g_times_s_g = 0.0;

		if (_mode.hasPBM())
		{
			// if we solve the mass balance, then we have the solubility entry (last state entry) and the solute entry (first position), which is why we advance the pointer.
			StateType const* const yCrystal = y + 1;
			ResidualType* const resCrystal = res + 1;

			const StateType sParam = (cadet_likely(y[0] / y[_nComp - 1] - 1.0 > 0)) ? y[0] / y[_nComp - 1] - 1.0 : StateType(0.0); // s = (c_0 - c_eq) / c_eq = c_0 / c_eq - 1, rewrite it to zero if s drops below 0
			const ParamType massDensityShapeFactor = static_cast<ParamType>(_nucleiMassDensity) * static_cast<ParamType>(_volShapeFactor);
			k_g_times_s_g = static_cast<ParamType>(_growthRateConstant) * pow(sParam, static_cast<ParamType>(_g));

			ResidualType MParam = 0.0;
			ResidualType substrateConversion = 0.0;

			for (int i = 0; i < _nBins; ++i)
			{
				MParam += yCrystal[i] * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binSizes[i]);
				substrateConversion += yCrystal[i] * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_binCenters[i]), static_cast<ParamType>(_p))) * static_cast<ParamType>(_binSizes[i]);
			}
			MParam *= massDensityShapeFactor;
			substrateConversion *= 3.0 * k_g_times_s_g;

			// B_0 = primary + secondary nucleation rate
			B_0 = static_cast<ParamType>(_primaryNucleationRate) * pow(sParam, static_cast<ParamType>(_u)) + static_cast<ParamType>(_secondaryNucleationRate) * pow(sParam, static_cast<ParamType>(_b)) * pow(MParam, static_cast<ParamType>(_k));
			const ParamType x_c_3 = static_cast<ParamType>(_bins[0]) * static_cast<ParamType>(_bins[0]) * static_cast<ParamType>(_bins[0]);

			// mass balance
			res[0] -= factor * massDensityShapeFactor * (B_0 * x_c_3 + substrateConversion);

			// update growth flux reconstruction parameters
			_reconstruction->updateParams(k_g_times_s_g, B_0);
		}

		// get pointer to crystal bins
		StateType const* yCrystal = y;
		ResidualType* resCrystal = res;
		if (_mode.hasPBM()) // if we solve the mass balance, then we have the solubility entry (last state entry) and the solute entry (first position), which is why we advance the pointer.
		{
			yCrystal++;
			resCrystal++;
		}

		// define aggregation/fragmentation related local parameters
		ResidualType source = 0.0;
		ResidualType sink = 0.0;
		ParamType source_factor = 0.0;

		for (int i = 0; i < _nBins; ++i)
		{
			if (_mode.hasPBM())
			{
				PBMKernel<StateType, ResidualType, ParamType, FactorType>(yCrystal, resCrystal, factor, i);
			}

			if (_mode.hasAggregation())
			{
				// reset the source and sink terms
				source = 0.0;
				sink = 0.0;

				const ParamType x_i_cube = static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]);

				// aggregation source terms
				for (int p = static_cast<int>(_agg->indexTracker[i]); p < static_cast<int>(_agg->indexTracker[i + 1]); ++p)
				{
					const int j = static_cast<int>(_agg->index.first[p]);
					const int k = static_cast<int>(_agg->index.second[p]);

					const ParamType x_jCubeplusX_kCube = static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) + static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]);

					// source term
					source_factor = static_cast<ParamType>(_aggregationRateConstant) * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i])
						* 1.0 / (2.0 * x_i_cube / x_jCubeplusX_kCube - 1.0);

					if (cadet_unlikely(j == k)) { source_factor *= 0.5; }

					// add different kernels
					switch (_aggregationIndex)
					{
						//constant kernel 0
					case 0:
						break;
						// brownian kernel 1
					case 1:
						source_factor *= (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) / static_cast<ParamType>(_binCenters[k]) / static_cast<ParamType>(_binCenters[j]);
						break;
						// smoluchowski kernel 2
					case 2:
						source_factor *= (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j]));
						break;
						// golovin kernel 3
					case 3:
						source_factor *= static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]);
						break;
						// differential force kernel 4
					case 4:
						source_factor *= (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) - static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]));
						break;
					}

					source += yCrystal[j] * yCrystal[k] * source_factor;
				}

				// aggregation sink terms
				for (int j = 0; j < _nBins; ++j)
				{
					const ParamType sumVolume = pow(x_i_cube + static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]), 1.0 / 3.0);
					const int k = static_cast<int>(_agg->iIndex[i * _nBins + j]);

					ParamType ratio{};
					ParamType sinkCorrectionFactor = 1.0;
					if (k > 0)
					{
						ratio = sumVolume / static_cast<ParamType>(_binCenters[k]);
						sinkCorrectionFactor = 1.0 / (2.0 - ratio * ratio * ratio);
					}

					switch (_aggregationIndex)
					{
						//constant kernel 0
					case 0:
						sink += yCrystal[j] * static_cast<ParamType>(_binSizes[j]) * sinkCorrectionFactor;
						break;
						// brownian kernel 1
					case 1:
						sink += yCrystal[j] * static_cast<ParamType>(_binSizes[j]) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) / static_cast<ParamType>(_binCenters[i]) / static_cast<ParamType>(_binCenters[j]) * sinkCorrectionFactor;
						break;
						// smoluchowski kernel 2
					case 2:
						sink += yCrystal[j] * static_cast<ParamType>(_binSizes[j]) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * sinkCorrectionFactor;
						break;
						// golovin kernel 3
					case 3:
						sink += yCrystal[j] * static_cast<ParamType>(_binSizes[j]) * (static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * sinkCorrectionFactor;
						break;
						// differential force kernel 4
					case 4:
						sink += yCrystal[j] * static_cast<ParamType>(_binSizes[j]) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) - static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * sinkCorrectionFactor;
						break;
					}
				}
				sink *= yCrystal[i] * static_cast<ParamType>(_aggregationRateConstant);

				// add to residual
				resCrystal[i] += factor * source - factor * sink;
			}

					if (_mode.hasFragmentation())
			{
				const ParamType N_j = static_cast<ParamType>(_fragKernelGamma) / (static_cast<ParamType>(_fragKernelGamma) - 1.0);

				// ode input
				ParamType bIntegral = 0.0;
				ParamType selectionFunction = 0.0;
				ResidualType fragSource = 0.0;
				ResidualType fragSink = 0.0;

				// source term
				fragSource = 0.0;
				for (int j = i; j < _nBins; ++j)
				{
					// selection function
					selectionFunction = static_cast<ParamType>(_fragRateConstant) * pow(static_cast<ParamType>(_binCenters[j]), 3.0 * static_cast<ParamType>(_fragSelectionFunctionAlpha));
					if (cadet_likely(i != j))
					{
						bIntegral = N_j * (pow(static_cast<ParamType>(_bins[i + 1]), 3.0 * static_cast<ParamType>(_fragKernelGamma) - 3.0) - pow(static_cast<ParamType>(_bins[i]), 3.0 * static_cast<ParamType>(_fragKernelGamma) - 3.0)) / pow(static_cast<ParamType>(_binCenters[j]), 3.0 * static_cast<ParamType>(_fragKernelGamma) - 3.0);
					}
					else
					{
						bIntegral = N_j * (pow(static_cast<ParamType>(_binCenters[i]), 3.0 * static_cast<ParamType>(_fragKernelGamma) - 3.0) - pow(static_cast<ParamType>(_bins[i]), 3.0 * static_cast<ParamType>(_fragKernelGamma) - 3.0)) / pow(static_cast<ParamType>(_binCenters[j]), 3.0 * static_cast<ParamType>(_fragKernelGamma) - 3.0);
					}
					fragSource += selectionFunction * yCrystal[j] * bIntegral * static_cast<ParamType>(_frag->upsilonSource[j]) * static_cast<ParamType>(_binSizes[j]) / static_cast<ParamType>(_binSizes[i]);
				}

				// sink term
				selectionFunction = static_cast<ParamType>(_fragRateConstant) * pow(static_cast<ParamType>(_binCenters[i]), 3.0 * static_cast<ParamType>(_fragSelectionFunctionAlpha));
				fragSink = yCrystal[i] * selectionFunction * static_cast<ParamType>(_frag->upsilonSink[i]);

				// add to residual
				resCrystal[i] += factor * fragSource - factor * fragSink;
			}
		}
		return 0;
	}

	template <typename StateType, typename ResidualType, typename ParamType>
	int residualCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* yLiquid, StateType const* ySolid, ResidualType* resLiquid, ResidualType* resSolid, double factor, LinearBufferAllocator workSpace) const
	{
		return residualFluxImpl<StateType, ResidualType, ParamType, double>(t, secIdx, colPos, 0,  yLiquid, resLiquid, factor, workSpace);
	}

	template <typename RowIterator>
	void jacobianUpwindKernel(double const* yCrystal, double factor, RowIterator& jac, JacobianParamsBase& jacParams, const int i) const
	{
		int binIdx_i = 0;
		int binIdx_j = 0;

		if (cadet_likely((i > 1) && (i < _nBins)))
		{
			// intermediate cells
			binIdx_i = i - 1;
			// dQ_i/dc
			jac[0 - i] -= factor * (jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) * yCrystal[binIdx_i] - jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) * yCrystal[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]);

			// dQ_i/dn_{i-1}
			jac[-1] += factor * jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]);

			// dQ_i/dn_i
			jac[0] -= factor * jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binSizes[binIdx_i]) * (1.0 / static_cast<double>(_binCenterDists[binIdx_i]) + 1.0 / static_cast<double>(_binCenterDists[binIdx_i - 1]));

			// dQ_i/dn_{i+1}
			jac[1] += factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]);

			//dQ_i/dc_{eq}
			jac[_nComp - 1 - i] -= factor * (jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) * yCrystal[binIdx_i] - jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) * yCrystal[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]);
		}
		else if (cadet_unlikely(i == 1))
		{
			binIdx_i = i - 1;
			// Q_0, left BC
			for (int j = 0; j < _nComp; ++j)
			{
				binIdx_j = j - 1;
				if (cadet_likely((j > 2) && (j + 1 < _nComp)))
				{
					// dQ_0/dni
					jac[j - i] += factor * jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 0))
				{
					// dQ_0/dc
					jac[j - i] -= factor * (yCrystal[0] * jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) - jacParams.dBp_dc - jacParams.dBs_dc) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 1))
				{
					// dQ_0/dn0
					jac[j - i] -= factor * (jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_j + 1]), static_cast<double>(_p))) - jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j])) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 2))
				{
					// dQ_0/dn1
					jac[j - i] += factor * jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j]) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else
				{
					// dQ_0/dceq
					jac[j - i] -= factor * (yCrystal[0] * jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) - jacParams.dBp_dceq - jacParams.dBs_dceq) / static_cast<double>(_binSizes[binIdx_i]);
				}
			}
		}
		else if (cadet_unlikely(i == _nBins))
		{
			binIdx_i = i - 1;
			// Q_{N_x-1}, right BC

			// dQ_{N_x-1}/dc
			jac[0 - i] += factor * yCrystal[binIdx_i - 1] * jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) / static_cast<double>(_binSizes[binIdx_i]);

			// dQ_{N_x-1}/dn_{N_x-2}
			jac[_nComp - 3 - i] += factor * jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]);

			// dQ_{N_x-1}/dn_{N_x-1}
			jac[_nComp - 2 - i] -= factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]);

			// dQ_{N_x-1}/dc_eq
			jac[_nComp - 1 - i] += factor * yCrystal[binIdx_i - 1] * jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) / static_cast<double>(_binSizes[binIdx_i]);
		}
	}

	template <typename RowIterator>
	void jacobianHRKorenKernel(double const* yCrystal, double factor, RowIterator& jac, JacobianParamsHRKoren& jacParams, HRKorenParams<active>& HRKoren, const int i) const
	{
		int binIdx_i = 0;
		int binIdx_j = 0;

		// Q_i, intermediate cells
		if (cadet_likely((i > 2) && (i < _nBins)))
		{
			binIdx_i = i - 1;
			// jacobian, left face, which is the right face of a previous cell
			// dQ_i/dc, j=0
			jac[-i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * jacParams.F_i;

			// dQ_i/dn_{i-2}, j=i-2
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_nim1; // the second term in the product rule is zero

			// dQ_i/dn_{i-1}, j=i-1
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_ni; // full four terms in the product rule

			// dQ_i/dn_i, j=i
			jac[0] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_nip1; // the fourth term in the product rule is zero

			//dQ_i/dc_{eq}, j=_nComp - 1
			jac[_nComp - 1 - i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * jacParams.F_i;

			// update all coefficients
			// HR-related coefficients, right face
			jacParams.r_x_i = static_cast<double>(HRKoren.A_coeff[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1] + jacParams.epsilon) / (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i] + jacParams.epsilon);
			if (cadet_likely(jacParams.r_x_i > 0))
			{
				jacParams.phi = jacParams.r_x_i / (static_cast<double>(HRKoren.R_coeff[binIdx_i]) - 1.0 + jacParams.r_x_i);
				jacParams.F_i = yCrystal[binIdx_i] + jacParams.phi * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);

				// Jacobian related coefficients, right face
				jacParams.ni_difference = 1.0 - jacParams.epsilon / (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i] + jacParams.epsilon);
				jacParams.r_cluster = jacParams.phi;
				jacParams.A_cluster = static_cast<double>(HRKoren.A_coeff[binIdx_i]) / (static_cast<double>(HRKoren.R_coeff[binIdx_i]) - 1.0 + jacParams.r_x_i);
				jacParams.r_square_cluster = jacParams.phi * jacParams.phi;
				jacParams.Ar_cluster = jacParams.phi * jacParams.A_cluster;

				jacParams.dFi_wrt_nim1 = jacParams.Ar_cluster * jacParams.ni_difference - jacParams.ni_difference * jacParams.A_cluster;
				jacParams.dFi_wrt_ni = 1.0 - jacParams.ni_difference * (jacParams.r_square_cluster + jacParams.Ar_cluster) + jacParams.r_cluster * jacParams.ni_difference - jacParams.r_cluster + jacParams.ni_difference * jacParams.A_cluster;
				jacParams.dFi_wrt_nip1 = jacParams.r_square_cluster * jacParams.ni_difference - jacParams.r_cluster * jacParams.ni_difference + jacParams.r_cluster;
			}
			else
			{
				jacParams.phi = 0.0;
				jacParams.F_i = yCrystal[binIdx_i];

				// Jacobian related coefficients, right face
				jacParams.dFi_wrt_nim1 = 0.0;
				jacParams.dFi_wrt_ni = 1.0; // becomes upwind
				jacParams.dFi_wrt_nip1 = 0.0;
			}
			jacParams.vg_right = jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			jacParams.dvG_dc_right = jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.dvG_dceq_right = jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			// jacobian, right face, and diffusion
			// dQ_i/dc, j=0
			jac[-i] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * jacParams.F_i;

			// dQ_i/dn_{i-1}, j=i-1
			jac[-1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_nim1 - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // the second term in the product rule is zero

			// dQ_i/dn_i, j=i
			jac[0] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_ni + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binSizes[binIdx_i]) * (1.0 / static_cast<double>(_binCenterDists[binIdx_i]) + 1.0 / static_cast<double>(_binCenterDists[binIdx_i - 1])); // full four terms in the product rule

			// dQ_i/dn_{i+1}, j=i+1
			jac[1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_nip1 - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]); // the fourth term in the product rule is zero

			//dQ_i/dc_{eq}, j=_nComp - 1
			jac[_nComp - 1 - i] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * jacParams.F_i;
		}

		// Q_0, left BC, also independent on the scheme. Do not change.
		else if (cadet_unlikely(i == 1))
		{
			binIdx_i = i - 1;
			// Q_0, left BC
			for (int j = 0; j < _nComp; ++j)
			{
				binIdx_j = j - 1;
				if (cadet_likely((j > 2) && (j + 1 < _nComp)))
				{
					// dQ_0/dni
					jac[j - i] += factor * jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 0))
				{
					// dQ_0/dc
					jac[j - i] -= factor * (yCrystal[0] * jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) - jacParams.dBp_dc - jacParams.dBs_dc) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 1))
				{
					// dQ_0/dn0
					jac[j - i] -= factor * (jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_j + 1]), static_cast<double>(_p))) - jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j])) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 2))
				{
					// dQ_0/dn1
					jac[j - i] += factor * jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j]) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else
				{
					// dQ_0/dceq
					jac[j - i] -= factor * (yCrystal[0] * jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) - jacParams.dBp_dceq - jacParams.dBs_dceq) / static_cast<double>(_binSizes[binIdx_i]);
				}
			}
		}

		// Q_1, special case for HR.
		else if (cadet_unlikely(i == 2))
		{
			binIdx_i = i - 1;
			// HR-related coefficients, right face
			jacParams.r_x_i = static_cast<double>(HRKoren.A_coeff[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1] + jacParams.epsilon) / (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i] + jacParams.epsilon);
			if (cadet_likely(jacParams.r_x_i > 0))
			{
				jacParams.phi = jacParams.r_x_i / (static_cast<double>(HRKoren.R_coeff[binIdx_i]) - 1.0 + jacParams.r_x_i);
				jacParams.F_i = yCrystal[binIdx_i] + jacParams.phi * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);

				// Jacobian related coefficients, right face
				jacParams.ni_difference = 1.0 - jacParams.epsilon / (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i] + jacParams.epsilon);
				jacParams.r_cluster = jacParams.phi;
				jacParams.A_cluster = static_cast<double>(HRKoren.A_coeff[binIdx_i]) / (static_cast<double>(HRKoren.R_coeff[binIdx_i]) - 1.0 + jacParams.r_x_i);
				jacParams.r_square_cluster = jacParams.phi * jacParams.phi;
				jacParams.Ar_cluster = jacParams.phi * jacParams.A_cluster;

				jacParams.dFi_wrt_nim1 = jacParams.Ar_cluster * jacParams.ni_difference - jacParams.ni_difference * jacParams.A_cluster;
				jacParams.dFi_wrt_ni = 1.0 - jacParams.ni_difference * (jacParams.r_square_cluster + jacParams.Ar_cluster) + jacParams.r_cluster * jacParams.ni_difference - jacParams.r_cluster + jacParams.ni_difference * jacParams.A_cluster;
				jacParams.dFi_wrt_nip1 = jacParams.r_square_cluster * jacParams.ni_difference - jacParams.r_cluster * jacParams.ni_difference + jacParams.r_cluster;
			}
			else
			{
				jacParams.phi = 0.0;
				jacParams.F_i = yCrystal[binIdx_i];

				// Jacobian related coefficients, right face
				jacParams.dFi_wrt_nim1 = 0.0;
				jacParams.dFi_wrt_ni = 1.0; // becomes upwind
				jacParams.dFi_wrt_nip1 = 0.0;
			}
			jacParams.vg_right = jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			jacParams.dvG_dc_right = jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.dvG_dceq_right = jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			// dQ_1/dc, j=0
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) * yCrystal[0]; // left face, upwind scheme
			jac[-2] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * jacParams.F_i; // right face, HR

			// dQ_1/dn0, j=1
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // coming from the upwind scheme
			jac[-1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_nim1; // the second term in the product rule is zero

			// dQ_1/dn1, j=2
			jac[0] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_ni + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binSizes[binIdx_i]) * (1.0 / static_cast<double>(_binCenterDists[binIdx_i]) + 1.0 / static_cast<double>(_binCenterDists[binIdx_i - 1])); // full four terms in the product rule

			// dQ_1/dn2, j=3
			jac[1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_nip1 - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]); // the fourth term in the product rule is zero

			// dQ_1/dceq, j=_nComp -1
			jac[_nComp - 3] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) * yCrystal[0];
			jac[_nComp - 3] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * jacParams.F_i;
		}

		// Q_{N_x-1}, right BC
		else if (cadet_unlikely(i == _nBins))
		{
			binIdx_i = i - 1;
			// dQ_{N_x-1}/dc
			jac[-i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * jacParams.F_i;

			// dQ_{N_x-1}/dn_{N_x-3}
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_nim1; // the second term in the product rule is zero

			// dQ_{N_x-1}/dn_{N_x-2}
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_ni + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // full four terms in the product rule

			// dQ_{N_x-1}/dn_{N_x-1}
			jac[0] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * jacParams.dFi_wrt_nip1 - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // the fourth term in the product rule is zero

			//dQ_i/dc_{eq}, j=_nComp - 1
			jac[1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * jacParams.F_i;
		}
	}

	template <typename RowIterator>
	void jacobianWENO23Kernel(double const* yCrystal, double factor, RowIterator& jac, JacobianParamsWENO23& jacParams, WENO23Params<active>& WENO23, const int i) const
	{
		int binIdx_i = 0;
		int binIdx_j = 0; // note: binIdx_j is used only when B_s is involved.

		// Q_i, intermediate cells
		if (cadet_likely((i > 2) && (i < _nBins)))
		{
			binIdx_i = i - 1;

			// jacobian, left face, which is the right face of a previous cell
			// dQ_i/dc, j=0
			jac[-i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0_right * jacParams.q_0_right + jacParams.W_1_right * jacParams.q_1_right);

			// dQ_i/dn_{i-2}, j=i-2
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_m1 * jacParams.q_0_right + jacParams.dw1_right_dni_m1 * jacParams.q_1_right + jacParams.W_1_right * jacParams.dq1_right_dni_m1); // the second term in the product rule is zero

			// dQ_i/dn_{i-1}, j=i-1
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni * jacParams.q_0_right + jacParams.W_0_right * jacParams.dq0_right_dni + jacParams.dw1_right_dni * jacParams.q_1_right + jacParams.W_1_right * jacParams.dq1_right_dni); // full four terms in the product rule

			// dQ_i/dn_i, j=i
			jac[0] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_p1 * jacParams.q_0_right + jacParams.W_0_right * jacParams.dq0_right_dni_p1 + jacParams.dw1_right_dni_p1 * jacParams.q_1_right); // the fourth term in the product rule is zero

			//dQ_i/dc_{eq}, j=_nComp - 1
			jac[_nComp - 1 - i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0_right * jacParams.q_0_right + jacParams.W_1_right * jacParams.q_1_right);

			// update all coefficients
			// WENO23-related coefficients, right face
			jacParams.IS_0_right = static_cast<double>(WENO23.IS_0_coeff[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);
			jacParams.IS_1_right = static_cast<double>(WENO23.IS_1_coeff[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);
			jacParams.alpha_0_right = static_cast<double>(WENO23.C_right_coeff[binIdx_i]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0_right) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0_right);
			jacParams.alpha_1_right = (1.0 - static_cast<double>(WENO23.C_right_coeff[binIdx_i])) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1_right) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1_right);
			jacParams.W_0_right = jacParams.alpha_0_right / (jacParams.alpha_0_right + jacParams.alpha_1_right);
			jacParams.W_1_right = 1.0 - jacParams.W_0_right;
			jacParams.q_0_right = static_cast<double>(WENO23.q_0_right_coeff[binIdx_i]) * yCrystal[binIdx_i] + (1.0 - static_cast<double>(WENO23.q_0_right_coeff[binIdx_i])) * yCrystal[binIdx_i + 1];
			jacParams.q_1_right = static_cast<double>(WENO23.q_1_right_coeff[binIdx_i]) * yCrystal[binIdx_i] + (1.0 - static_cast<double>(WENO23.q_1_right_coeff[binIdx_i])) * yCrystal[binIdx_i - 1];

			jacParams.dvG_dc_right = jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.dvG_dceq_right = jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			jacParams.vg_right = jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			// Jacobian related coefficients, right face
			jacParams.four_w0_B0_right = 4.0 * jacParams.W_0_right * static_cast<double>(WENO23.IS_0_coeff[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0_right);
			jacParams.four_w1_B1_right = 4.0 * jacParams.W_1_right * static_cast<double>(WENO23.IS_1_coeff[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1_right);
			jacParams.four_w0_w0_B0_right = jacParams.four_w0_B0_right * jacParams.W_0_right;
			jacParams.four_w0_w1_B1_right = jacParams.four_w1_B1_right * jacParams.W_0_right;
			jacParams.four_w0_w1_B0_right = jacParams.four_w0_B0_right * jacParams.W_1_right;
			jacParams.four_w1_w1_B1_right = jacParams.four_w1_B1_right * jacParams.W_1_right;

			jacParams.dw0_right_dni_m1 = -jacParams.four_w0_w1_B1_right;
			jacParams.dw0_right_dni = jacParams.four_w0_B0_right - jacParams.four_w0_w0_B0_right + jacParams.four_w0_w1_B1_right;
			jacParams.dw0_right_dni_p1 = jacParams.four_w0_w0_B0_right - jacParams.four_w0_B0_right;

			jacParams.dw1_right_dni_m1 = jacParams.four_w1_B1_right - jacParams.four_w1_w1_B1_right;
			jacParams.dw1_right_dni = -jacParams.four_w1_B1_right - jacParams.four_w0_w1_B0_right + jacParams.four_w1_w1_B1_right;
			jacParams.dw1_right_dni_p1 = jacParams.four_w0_w1_B0_right;

			jacParams.dq0_right_dni = static_cast<double>(WENO23.q_0_right_coeff[binIdx_i]);
			jacParams.dq0_right_dni_p1 = 1.0 - static_cast<double>(WENO23.q_0_right_coeff[binIdx_i]);

			jacParams.dq1_right_dni_m1 = 1.0 - static_cast<double>(WENO23.q_1_right_coeff[binIdx_i]);
			jacParams.dq1_right_dni = static_cast<double>(WENO23.q_1_right_coeff[binIdx_i]);

			// jacobian, right face, and diffusion
			// dQ_i/dc, j=0
			jac[-i] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0_right * jacParams.q_0_right + jacParams.W_1_right * jacParams.q_1_right);

			// dQ_i/dn_{i-1}, j=i-1
			jac[-1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_m1 * jacParams.q_0_right + jacParams.dw1_right_dni_m1 * jacParams.q_1_right + jacParams.W_1_right * jacParams.dq1_right_dni_m1) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // the second term in the product rule is zero

			// dQ_i/dn_i, j=i
			jac[0] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni * jacParams.q_0_right + jacParams.W_0_right * jacParams.dq0_right_dni + jacParams.dw1_right_dni * jacParams.q_1_right + jacParams.W_1_right * jacParams.dq1_right_dni) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binSizes[binIdx_i]) * (1.0 / static_cast<double>(_binCenterDists[binIdx_i]) + 1.0 / static_cast<double>(_binCenterDists[binIdx_i - 1])); // full four terms in the product rule

			// dQ_i/dn_{i+1}, j=i+1
			jac[1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_p1 * jacParams.q_0_right + jacParams.W_0_right * jacParams.dq0_right_dni_p1 + jacParams.dw1_right_dni_p1 * jacParams.q_1_right) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]); // the fourth term in the product rule is zero

			//dQ_i/dc_{eq}, j=_nComp - 1
			jac[_nComp - 1 - i] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0_right * jacParams.q_0_right + jacParams.W_1_right * jacParams.q_1_right);
		}

		// Q_0, left BC, also independent on the scheme. Do not change.
		else if (cadet_unlikely(i == 1))
		{
			binIdx_i = i - 1;
			// Q_0, left BC
			for (int j = 0; j < _nComp; ++j)
			{
				binIdx_j = j - 1;
				if (cadet_likely((j > 2) && (j + 1 < _nComp)))
				{
					// dQ_0/dni
					jac[j - i] += factor * jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 0))
				{
					// dQ_0/dc
					jac[j - i] -= factor * (yCrystal[0] * jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) - jacParams.dBp_dc - jacParams.dBs_dc) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 1))
				{
					// dQ_0/dn0
					jac[j - i] -= factor * (jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_j + 1]), static_cast<double>(_p))) - jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j])) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 2))
				{
					// dQ_0/dn1
					jac[j - i] += factor * jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j]) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else
				{
					// dQ_0/dceq
					jac[j - i] -= factor * (yCrystal[0] * jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) - jacParams.dBp_dceq - jacParams.dBs_dceq) / static_cast<double>(_binSizes[binIdx_i]);
				}
			}
		}

		// Q_1, special case for weno3.
		else if (cadet_unlikely(i == 2))
		{
			binIdx_i = i - 1;
			// WENO23-related coefficients, right face
			jacParams.IS_0_right = static_cast<double>(WENO23.IS_0_coeff[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);
			jacParams.IS_1_right = static_cast<double>(WENO23.IS_1_coeff[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);
			jacParams.alpha_0_right = static_cast<double>(WENO23.C_right_coeff[binIdx_i]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0_right) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0_right);
			jacParams.alpha_1_right = (1.0 - static_cast<double>(WENO23.C_right_coeff[binIdx_i])) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1_right) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1_right);
			jacParams.W_0_right = jacParams.alpha_0_right / (jacParams.alpha_0_right + jacParams.alpha_1_right);
			jacParams.W_1_right = 1.0 - jacParams.W_0_right;
			jacParams.q_0_right = static_cast<double>(WENO23.q_0_right_coeff[binIdx_i]) * yCrystal[binIdx_i] + (1.0 - static_cast<double>(WENO23.q_0_right_coeff[binIdx_i])) * yCrystal[binIdx_i + 1];
			jacParams.q_1_right = static_cast<double>(WENO23.q_1_right_coeff[binIdx_i]) * yCrystal[binIdx_i] + (1.0 - static_cast<double>(WENO23.q_1_right_coeff[binIdx_i])) * yCrystal[binIdx_i - 1];

			jacParams.dvG_dc_right = jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.dvG_dceq_right = jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			jacParams.vg_right = jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			// Jacobian related coefficients, right face
			jacParams.four_w0_B0_right = 4.0 * jacParams.W_0_right * static_cast<double>(WENO23.IS_0_coeff[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0_right);
			jacParams.four_w1_B1_right = 4.0 * jacParams.W_1_right * static_cast<double>(WENO23.IS_1_coeff[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1_right);
			jacParams.four_w0_w0_B0_right = jacParams.four_w0_B0_right * jacParams.W_0_right;
			jacParams.four_w0_w1_B1_right = jacParams.four_w1_B1_right * jacParams.W_0_right;
			jacParams.four_w0_w1_B0_right = jacParams.four_w0_B0_right * jacParams.W_1_right;
			jacParams.four_w1_w1_B1_right = jacParams.four_w1_B1_right * jacParams.W_1_right;

			jacParams.dw0_right_dni_m1 = -jacParams.four_w0_w1_B1_right;
			jacParams.dw0_right_dni = jacParams.four_w0_B0_right - jacParams.four_w0_w0_B0_right + jacParams.four_w0_w1_B1_right;
			jacParams.dw0_right_dni_p1 = jacParams.four_w0_w0_B0_right - jacParams.four_w0_B0_right;

			jacParams.dw1_right_dni_m1 = jacParams.four_w1_B1_right - jacParams.four_w1_w1_B1_right;
			jacParams.dw1_right_dni = -jacParams.four_w1_B1_right - jacParams.four_w0_w1_B0_right + jacParams.four_w1_w1_B1_right;
			jacParams.dw1_right_dni_p1 = jacParams.four_w0_w1_B0_right;

			jacParams.dq0_right_dni = static_cast<double>(WENO23.q_0_right_coeff[binIdx_i]);
			jacParams.dq0_right_dni_p1 = 1.0 - static_cast<double>(WENO23.q_0_right_coeff[binIdx_i]);

			jacParams.dq1_right_dni_m1 = 1.0 - static_cast<double>(WENO23.q_1_right_coeff[binIdx_i]);
			jacParams.dq1_right_dni = static_cast<double>(WENO23.q_1_right_coeff[binIdx_i]);

			// dQ_1/dc, j=0
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) * yCrystal[0]; // left face, upwind scheme
			jac[-2] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0_right * jacParams.q_0_right + jacParams.W_1_right * jacParams.q_1_right); // right face, WENO3

			// dQ_1/dn0, j=1
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // coming from the upwind scheme
			jac[-1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_m1 * jacParams.q_0_right + jacParams.dw1_right_dni_m1 * jacParams.q_1_right + jacParams.W_1_right * jacParams.dq1_right_dni_m1); // the second term in the product rule is zero

			// dQ_1/dn1, j=2
			jac[0] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni * jacParams.q_0_right + jacParams.W_0_right * jacParams.dq0_right_dni + jacParams.dw1_right_dni * jacParams.q_1_right + jacParams.W_1_right * jacParams.dq1_right_dni) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binSizes[binIdx_i]) * (1.0 / static_cast<double>(_binCenterDists[binIdx_i]) + 1.0 / static_cast<double>(_binCenterDists[binIdx_i - 1])); // full four terms in the product rule

			// dQ_1/dn2, j=3
			jac[1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_p1 * jacParams.q_0_right + jacParams.W_0_right * jacParams.dq0_right_dni_p1 + jacParams.dw1_right_dni_p1 * jacParams.q_1_right) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]); // the fourth term in the product rule is zero

			// dQ_1/dceq, j=_nComp -1
			jac[_nComp - 3] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) * yCrystal[0];
			jac[_nComp - 3] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0_right * jacParams.q_0_right + jacParams.W_1_right * jacParams.q_1_right);
		}

		// Q_{N_x-1}, right BC
		else if (cadet_unlikely(i == _nBins))
		{
			binIdx_i = i - 1;
			// dQ_{N_x-1}/dc
			jac[-i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0_right * jacParams.q_0_right + jacParams.W_1_right * jacParams.q_1_right);

			// dQ_{N_x-1}/dn_{N_x-3}
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_m1 * jacParams.q_0_right + jacParams.dw1_right_dni_m1 * jacParams.q_1_right + jacParams.W_1_right * jacParams.dq1_right_dni_m1); // the second term in the product rule is zero

			// dQ_{N_x-1}/dn_{N_x-2}
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni * jacParams.q_0_right + jacParams.W_0_right * jacParams.dq0_right_dni + jacParams.dw1_right_dni * jacParams.q_1_right + jacParams.W_1_right * jacParams.dq1_right_dni) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // full four terms in the product rule

			// dQ_{N_x-1}/dn_{N_x-1}
			jac[0] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_p1 * jacParams.q_0_right + jacParams.W_0_right * jacParams.dq0_right_dni_p1 + jacParams.dw1_right_dni_p1 * jacParams.q_1_right) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // the fourth term in the product rule is zero

			//dQ_i/dc_{eq}, j=_nComp - 1
			jac[1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0_right * jacParams.q_0_right + jacParams.W_1_right * jacParams.q_1_right);
		}
	}

	template <typename RowIterator>
	void jacobianWENO35Kernel(double const* yCrystal, double factor, RowIterator& jac, JacobianParamsWENO35& jacParams, WENO35Params<active>& WENO35, const int i) const
	{
		int binIdx_i = 0;
		int binIdx_j = 0; // note: binIdx_j is used only when B_s is involved.

		// Q_i, intermediate cells
		if (cadet_likely((i > 3) && (i + 3 < _nComp)))
		{
			binIdx_i = i - 1;
			// jacobian, left face, which is the right face of a previous cell
			// dQ_i/dc
			jac[-i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0 * jacParams.q_0 + jacParams.W_1 * jacParams.q_1 + jacParams.W_2 * jacParams.q_2);

			// dQ_i/dn_{i-2}
			jac[-3] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_m2 * jacParams.q_0 + jacParams.dw1_wrt_dni_m2 * jacParams.q_1 + jacParams.dw2_wrt_dni_m2 * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni_m2); // the second and fourth terms are zero

			// dQ_i/dn_{i-1}
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_m1 * jacParams.q_0 + jacParams.dw1_wrt_dni_m1 * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni_m1 + jacParams.dw2_wrt_dni_m1 * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni_m1); // the second term is zero

			// dQ_i/dn_{i}
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni + jacParams.dw1_wrt_dni * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni + jacParams.dw2_wrt_dni * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni); // full six terms of the product

			// dQ_i/dn_{i+1}
			jac[0] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_p1 * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni_p1 + jacParams.dw1_wrt_dni_p1 * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni_p1 + jacParams.dw2_wrt_dni_p1 * jacParams.q_2); // the last term is zero

			// dQ_i/dn_{i+2}
			jac[1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_p2 * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni_p2 + jacParams.dw1_wrt_dni_p2 * jacParams.q_1 + jacParams.dw2_wrt_dni_p2 * jacParams.q_2); // the fourth and sixth terms are zero

			// dQ_2/dc_eq
			jac[_nComp - 1 - i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0 * jacParams.q_0 + jacParams.W_1 * jacParams.q_1 + jacParams.W_2 * jacParams.q_2);

			// update all coefficients
			// WENO5, right face
			jacParams.IS_0 = static_cast<double>(WENO35.IS_0_coeff_1[binIdx_i]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) + static_cast<double>(WENO35.IS_0_coeff_2[binIdx_i]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]) + static_cast<double>(WENO35.IS_0_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]);
			jacParams.IS_1 = static_cast<double>(WENO35.IS_1_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) + static_cast<double>(WENO35.IS_1_coeff_2[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) + static_cast<double>(WENO35.IS_1_coeff_3[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);
			jacParams.IS_2 = static_cast<double>(WENO35.IS_2_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + static_cast<double>(WENO35.IS_2_coeff_2[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + static_cast<double>(WENO35.IS_2_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);
			jacParams.alpha_0 = static_cast<double>(WENO35.C_0[binIdx_i]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0);
			jacParams.alpha_1 = static_cast<double>(WENO35.C_1[binIdx_i]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1);
			jacParams.alpha_2 = static_cast<double>(WENO35.C_2[binIdx_i]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2);
			jacParams.W_0 = jacParams.alpha_0 / (jacParams.alpha_0 + jacParams.alpha_1 + jacParams.alpha_2);
			jacParams.W_1 = jacParams.alpha_1 / (jacParams.alpha_0 + jacParams.alpha_1 + jacParams.alpha_2);
			jacParams.W_2 = 1.0 - jacParams.W_0 - jacParams.W_1;
			jacParams.q_0 = yCrystal[binIdx_i + 1] + static_cast<double>(WENO35.q_0_coeff_1[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]) + static_cast<double>(WENO35.q_0_coeff_2[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i + 2]);
			jacParams.q_1 = yCrystal[binIdx_i] + static_cast<double>(WENO35.q_1_coeff_1[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]) + static_cast<double>(WENO35.q_1_coeff_2[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);
			jacParams.q_2 = yCrystal[binIdx_i - 1] + static_cast<double>(WENO35.q_2_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + static_cast<double>(WENO35.q_2_coeff_2[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);

			jacParams.dvG_dc_right = jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.dvG_dceq_right = jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.vg_right = jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			// jacobian
			jacParams.dIS0_wrt_dni_p2 = 2.0 * static_cast<double>(WENO35.IS_0_coeff_1[binIdx_i]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) + static_cast<double>(WENO35.IS_0_coeff_2[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]);
			jacParams.dIS0_wrt_dni_p1 = -2.0 * static_cast<double>(WENO35.IS_0_coeff_1[binIdx_i]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) + static_cast<double>(WENO35.IS_0_coeff_2[binIdx_i]) * (2.0 * yCrystal[binIdx_i + 1] - yCrystal[binIdx_i + 2] - yCrystal[binIdx_i]) - 2.0 * static_cast<double>(WENO35.IS_0_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]);
			jacParams.dIS0_wrt_dni = static_cast<double>(WENO35.IS_0_coeff_2[binIdx_i]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) + 2.0 * static_cast<double>(WENO35.IS_0_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]);

			jacParams.dIS1_wrt_dni_p1 = static_cast<double>(WENO35.IS_1_coeff_2[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) + 2.0 * static_cast<double>(WENO35.IS_1_coeff_3[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);
			jacParams.dIS1_wrt_dni = -2.0 * static_cast<double>(WENO35.IS_1_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) + static_cast<double>(WENO35.IS_1_coeff_2[binIdx_i]) * (2.0 * yCrystal[binIdx_i] - yCrystal[binIdx_i + 1] - yCrystal[binIdx_i - 1]) - 2.0 * static_cast<double>(WENO35.IS_1_coeff_3[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);
			jacParams.dIS1_wrt_dni_m1 = 2.0 * static_cast<double>(WENO35.IS_1_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) + static_cast<double>(WENO35.IS_1_coeff_2[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);

			jacParams.dIS2_wrt_dni = static_cast<double>(WENO35.IS_2_coeff_2[binIdx_i]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + 2.0 * static_cast<double>(WENO35.IS_2_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);
			jacParams.dIS2_wrt_dni_m1 = -2.0 * static_cast<double>(WENO35.IS_2_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + static_cast<double>(WENO35.IS_2_coeff_2[binIdx_i]) * (2.0 * yCrystal[binIdx_i - 1] - yCrystal[binIdx_i] - yCrystal[binIdx_i - 2]) - 2.0 * static_cast<double>(WENO35.IS_2_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);
			jacParams.dIS2_wrt_dni_m2 = 2.0 * static_cast<double>(WENO35.IS_2_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + static_cast<double>(WENO35.IS_2_coeff_2[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);

			jacParams.dq0_wrt_dni_p2 = -static_cast<double>(WENO35.q_0_coeff_2[binIdx_i]);
			jacParams.dq0_wrt_dni_p1 = static_cast<double>(WENO35.q_0_coeff_2[binIdx_i]) - static_cast<double>(WENO35.q_0_coeff_1[binIdx_i]) + 1.0;
			jacParams.dq0_wrt_dni = static_cast<double>(WENO35.q_0_coeff_1[binIdx_i]);

			jacParams.dq1_wrt_dni_p1 = static_cast<double>(WENO35.q_1_coeff_1[binIdx_i]);
			jacParams.dq1_wrt_dni = -static_cast<double>(WENO35.q_1_coeff_1[binIdx_i]) + static_cast<double>(WENO35.q_1_coeff_2[binIdx_i]) + 1.0;
			jacParams.dq1_wrt_dni_m1 = -static_cast<double>(WENO35.q_1_coeff_2[binIdx_i]);

			jacParams.dq2_wrt_dni = static_cast<double>(WENO35.q_2_coeff_2[binIdx_i]);
			jacParams.dq2_wrt_dni_m1 = -static_cast<double>(WENO35.q_2_coeff_1[binIdx_i]) - static_cast<double>(WENO35.q_2_coeff_2[binIdx_i]) + 1.0;
			jacParams.dq2_wrt_dni_m2 = static_cast<double>(WENO35.q_2_coeff_1[binIdx_i]);

			// dw_0/dx
			jacParams.w0_w2_is2 = 2.0 * jacParams.W_0 * jacParams.W_2 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2);
			jacParams.w0_w1_is1 = 2.0 * jacParams.W_0 * jacParams.W_1 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1);
			jacParams.w0_w0_is0 = 2.0 * jacParams.W_0 * jacParams.W_0 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0);
			jacParams.w0_is0 = -2.0 * jacParams.W_0 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0);

			jacParams.dw0_wrt_dni_m2 = jacParams.w0_w2_is2 * jacParams.dIS2_wrt_dni_m2;
			jacParams.dw0_wrt_dni_m1 = jacParams.w0_w1_is1 * jacParams.dIS1_wrt_dni_m1 + jacParams.w0_w2_is2 * jacParams.dIS2_wrt_dni_m1;
			jacParams.dw0_wrt_dni = jacParams.w0_is0 * jacParams.dIS0_wrt_dni + jacParams.w0_w0_is0 * jacParams.dIS0_wrt_dni + jacParams.w0_w1_is1 * jacParams.dIS1_wrt_dni + jacParams.w0_w2_is2 * jacParams.dIS2_wrt_dni;
			jacParams.dw0_wrt_dni_p1 = jacParams.w0_is0 * jacParams.dIS0_wrt_dni_p1 + jacParams.w0_w0_is0 * jacParams.dIS0_wrt_dni_p1 + jacParams.w0_w1_is1 * jacParams.dIS1_wrt_dni_p1;
			jacParams.dw0_wrt_dni_p2 = jacParams.w0_is0 * jacParams.dIS0_wrt_dni_p2 + jacParams.w0_w0_is0 * jacParams.dIS0_wrt_dni_p2;


			// dw_1/dx
			jacParams.w1_w2_is2 = 2.0 * jacParams.W_1 * jacParams.W_2 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2);
			jacParams.w1_w1_is1 = 2.0 * jacParams.W_1 * jacParams.W_1 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1);
			jacParams.w1_w0_is0 = 2.0 * jacParams.W_1 * jacParams.W_0 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0);
			jacParams.w1_is1 = -2.0 * jacParams.W_1 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1);

			jacParams.dw1_wrt_dni_m2 = jacParams.w1_w2_is2 * jacParams.dIS2_wrt_dni_m2;
			jacParams.dw1_wrt_dni_m1 = jacParams.w1_is1 * jacParams.dIS1_wrt_dni_m1 + jacParams.w1_w1_is1 * jacParams.dIS1_wrt_dni_m1 + jacParams.w1_w2_is2 * jacParams.dIS2_wrt_dni_m1;
			jacParams.dw1_wrt_dni = jacParams.w1_is1 * jacParams.dIS1_wrt_dni + jacParams.w1_w0_is0 * jacParams.dIS0_wrt_dni + jacParams.w1_w1_is1 * jacParams.dIS1_wrt_dni + jacParams.w1_w2_is2 * jacParams.dIS2_wrt_dni;
			jacParams.dw1_wrt_dni_p1 = jacParams.w1_is1 * jacParams.dIS1_wrt_dni_p1 + jacParams.w1_w0_is0 * jacParams.dIS0_wrt_dni_p1 + jacParams.w1_w1_is1 * jacParams.dIS1_wrt_dni_p1;
			jacParams.dw1_wrt_dni_p2 = jacParams.w1_w0_is0 * jacParams.dIS0_wrt_dni_p2;

			// dw2_dx
			jacParams.w2_is2 = -2.0 * jacParams.W_2 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2);
			jacParams.w2_w2_is2 = 2.0 * jacParams.W_2 * jacParams.W_2 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2);
			jacParams.w2_w1_is1 = 2.0 * jacParams.W_2 * jacParams.W_1 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1);
			jacParams.w2_w0_is0 = 2.0 * jacParams.W_2 * jacParams.W_0 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0);

			jacParams.dw2_wrt_dni_m2 = jacParams.w2_is2 * jacParams.dIS2_wrt_dni_m2 + jacParams.w2_w2_is2 * jacParams.dIS2_wrt_dni_m2;
			jacParams.dw2_wrt_dni_m1 = jacParams.w2_is2 * jacParams.dIS2_wrt_dni_m1 + jacParams.w2_w1_is1 * jacParams.dIS1_wrt_dni_m1 + jacParams.w2_w2_is2 * jacParams.dIS2_wrt_dni_m1;
			jacParams.dw2_wrt_dni = jacParams.w2_is2 * jacParams.dIS2_wrt_dni + jacParams.w2_w0_is0 * jacParams.dIS0_wrt_dni + jacParams.w2_w1_is1 * jacParams.dIS1_wrt_dni + jacParams.w2_w2_is2 * jacParams.dIS2_wrt_dni;
			jacParams.dw2_wrt_dni_p1 = jacParams.w2_w0_is0 * jacParams.dIS0_wrt_dni_p1 + jacParams.w2_w1_is1 * jacParams.dIS1_wrt_dni_p1;
			jacParams.dw2_wrt_dni_p2 = jacParams.w2_w0_is0 * jacParams.dIS0_wrt_dni_p2;

			// jacobian, right face, and diffusion
			// dQ_i/dc
			jac[-i] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0 * jacParams.q_0 + jacParams.W_1 * jacParams.q_1 + jacParams.W_2 * jacParams.q_2);

			// dQ_i/dn_{i-2}
			jac[-2] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_m2 * jacParams.q_0 + jacParams.dw1_wrt_dni_m2 * jacParams.q_1 + jacParams.dw2_wrt_dni_m2 * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni_m2); // the second and fourth terms are zero

			// dQ_i/dn_{i-1}
			jac[-1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_m1 * jacParams.q_0 + jacParams.dw1_wrt_dni_m1 * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni_m1 + jacParams.dw2_wrt_dni_m1 * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni_m1) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // the second term is zero

			// dQ_i/dn_{i}
			jac[0] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni + jacParams.dw1_wrt_dni * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni + jacParams.dw2_wrt_dni * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binSizes[binIdx_i]) * (1.0 / static_cast<double>(_binCenterDists[binIdx_i]) + 1.0 / static_cast<double>(_binCenterDists[binIdx_i - 1])); // full six terms of the product

			// dQ_i/dn_{i+1}
			jac[1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_p1 * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni_p1 + jacParams.dw1_wrt_dni_p1 * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni_p1 + jacParams.dw2_wrt_dni_p1 * jacParams.q_2) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]); // the last term is zero

			// dQ_i/dn_{i+2}
			jac[2] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_p2 * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni_p2 + jacParams.dw1_wrt_dni_p2 * jacParams.q_1 + jacParams.dw2_wrt_dni_p2 * jacParams.q_2); // the fourth and sixth terms are zero

			// dQ_2/dc_eq
			jac[_nComp - 1 - i] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0 * jacParams.q_0 + jacParams.W_1 * jacParams.q_1 + jacParams.W_2 * jacParams.q_2);

		}

		// Q_0, left BC, also independent on the scheme. Do not change.
		else if (cadet_unlikely(i == 1))
		{
			binIdx_i = i - 1;
			// Q_0, left BC
			for (int j = 0; j < _nComp; ++j)
			{
				binIdx_j = j - 1;
				if (cadet_likely((j > 2) && (j + 1 < _nComp)))
				{
					// dQ_0/dni
					jac[j - i] += factor * jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 0))
				{
					// dQ_0/dc
					jac[j - i] -= factor * (yCrystal[0] * jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) - jacParams.dBp_dc - jacParams.dBs_dc) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 1))
				{
					// dQ_0/dn0
					jac[j - i] -= factor * (jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_j + 1]), static_cast<double>(_p))) - jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j])) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else if (cadet_unlikely(j == 2))
				{
					// dQ_0/dn1
					jac[j - i] += factor * jacParams.dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j]) / static_cast<double>(_binSizes[binIdx_i]) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]);
				}
				else
				{
					// dQ_0/dceq
					jac[j - i] -= factor * (yCrystal[0] * jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p))) - jacParams.dBp_dceq - jacParams.dBs_dceq) / static_cast<double>(_binSizes[binIdx_i]);
				}
			}
		}

		// Q_1, special case for weno35.
		else if (cadet_unlikely(i == 2))
		{
			binIdx_i = i - 1;
			jacParams.dvG_dc_right = jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.dvG_dceq_right = jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.vg_right = jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			// dQ_1/dc, j=0
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) * yCrystal[0]; // left face, upwind scheme
			jac[-2] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0_weno3_left * jacParams.q_0_weno3_left + jacParams.W_1_weno3_left * jacParams.q_1_weno3_left); // right face, WENO3

			// dQ_1/dn0, j=1
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // coming from the upwind scheme
			jac[-1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_m1_weno3_left * jacParams.q_0_weno3_left + jacParams.dw1_right_dni_m1_weno3_left * jacParams.q_1_weno3_left + jacParams.W_1_weno3_left * jacParams.dq1_right_dni_m1_weno3_left); // the second term in the product rule is zero

			// dQ_1/dn1, j=2
			jac[0] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_weno3_left * jacParams.q_0_weno3_left + jacParams.W_0_weno3_left * jacParams.dq0_right_dni_weno3_left + jacParams.dw1_right_dni_weno3_left * jacParams.q_1_weno3_left + jacParams.W_1_weno3_left * jacParams.dq1_right_dni_weno3_left) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binSizes[binIdx_i]) * (1.0 / static_cast<double>(_binCenterDists[binIdx_i]) + 1.0 / static_cast<double>(_binCenterDists[binIdx_i - 1])); // full four terms in the product rule

			// dQ_1/dn2, j=3
			jac[1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_p1_weno3_left * jacParams.q_0_weno3_left + jacParams.W_0_weno3_left * jacParams.dq0_right_dni_p1_weno3_left + jacParams.dw1_right_dni_p1_weno3_left * jacParams.q_1_weno3_left) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]); // the fourth term in the product rule is zero

			// dQ_1/dceq, j=_nComp -1
			jac[_nComp - 3] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i]), static_cast<double>(_p))) * yCrystal[0];
			jac[_nComp - 3] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0_weno3_left * jacParams.q_0_weno3_left + jacParams.W_1_weno3_left * jacParams.q_1_weno3_left);
		}

		// Q_2, special case for weno35.
		else if (cadet_unlikely(i == 3))
		{
			binIdx_i = i - 1;
			// left face uses WENO3
			// dQ_2/dc
			jac[-3] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0_weno3_left * jacParams.q_0_weno3_left + jacParams.W_1_weno3_left * jacParams.q_1_weno3_left);

			// dQ_2/dn0
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_m1_weno3_left * jacParams.q_0_weno3_left + jacParams.dw1_right_dni_m1_weno3_left * jacParams.q_1_weno3_left + jacParams.W_1_weno3_left * jacParams.dq1_right_dni_m1_weno3_left); // the second term in the product rule is zero

			// dQ_2/dn1
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_weno3_left * jacParams.q_0_weno3_left + jacParams.W_0_weno3_left * jacParams.dq0_right_dni_weno3_left + jacParams.dw1_right_dni_weno3_left * jacParams.q_1_weno3_left + jacParams.W_1_weno3_left * jacParams.dq1_right_dni_weno3_left); // full four terms in the product rule

			// dQ_2/dn2
			jac[0] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_p1_weno3_left * jacParams.q_0_weno3_left + jacParams.W_0_weno3_left * jacParams.dq0_right_dni_p1_weno3_left + jacParams.dw1_right_dni_p1_weno3_left * jacParams.q_1_weno3_left); // the fourth term in the product rule is zero

			// dQ_2/dc_eq
			jac[_nComp - 4] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0_weno3_left * jacParams.q_0_weno3_left + jacParams.W_1_weno3_left * jacParams.q_1_weno3_left);

			// WENO5, right face
			jacParams.IS_0 = static_cast<double>(WENO35.IS_0_coeff_1[binIdx_i]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) + static_cast<double>(WENO35.IS_0_coeff_2[binIdx_i]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]) + static_cast<double>(WENO35.IS_0_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]);
			jacParams.IS_1 = static_cast<double>(WENO35.IS_1_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) + static_cast<double>(WENO35.IS_1_coeff_2[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) + static_cast<double>(WENO35.IS_1_coeff_3[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);
			jacParams.IS_2 = static_cast<double>(WENO35.IS_2_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + static_cast<double>(WENO35.IS_2_coeff_2[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + static_cast<double>(WENO35.IS_2_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);
			jacParams.alpha_0 = static_cast<double>(WENO35.C_0[binIdx_i]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0);
			jacParams.alpha_1 = static_cast<double>(WENO35.C_1[binIdx_i]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1);
			jacParams.alpha_2 = static_cast<double>(WENO35.C_2[binIdx_i]) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2) / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2);
			jacParams.W_0 = jacParams.alpha_0 / (jacParams.alpha_0 + jacParams.alpha_1 + jacParams.alpha_2);
			jacParams.W_1 = jacParams.alpha_1 / (jacParams.alpha_0 + jacParams.alpha_1 + jacParams.alpha_2);
			jacParams.W_2 = 1.0 - jacParams.W_0 - jacParams.W_1;
			jacParams.q_0 = yCrystal[binIdx_i + 1] + static_cast<double>(WENO35.q_0_coeff_1[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]) + static_cast<double>(WENO35.q_0_coeff_2[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i + 2]);
			jacParams.q_1 = yCrystal[binIdx_i] + static_cast<double>(WENO35.q_1_coeff_1[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]) + static_cast<double>(WENO35.q_1_coeff_2[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);
			jacParams.q_2 = yCrystal[binIdx_i - 1] + static_cast<double>(WENO35.q_2_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + static_cast<double>(WENO35.q_2_coeff_2[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);

			jacParams.dvG_dc_right = jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.dvG_dceq_right = jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.vg_right = jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			// jacobian
			jacParams.dIS0_wrt_dni_p2 = 2.0 * static_cast<double>(WENO35.IS_0_coeff_1[binIdx_i]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) + static_cast<double>(WENO35.IS_0_coeff_2[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]);
			jacParams.dIS0_wrt_dni_p1 = -2.0 * static_cast<double>(WENO35.IS_0_coeff_1[binIdx_i]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) + static_cast<double>(WENO35.IS_0_coeff_2[binIdx_i]) * (2.0 * yCrystal[binIdx_i + 1] - yCrystal[binIdx_i + 2] - yCrystal[binIdx_i]) - 2.0 * static_cast<double>(WENO35.IS_0_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]);
			jacParams.dIS0_wrt_dni = static_cast<double>(WENO35.IS_0_coeff_2[binIdx_i]) * (yCrystal[binIdx_i + 2] - yCrystal[binIdx_i + 1]) + 2.0 * static_cast<double>(WENO35.IS_0_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i + 1]);

			jacParams.dIS1_wrt_dni_p1 = static_cast<double>(WENO35.IS_1_coeff_2[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) + 2.0 * static_cast<double>(WENO35.IS_1_coeff_3[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);
			jacParams.dIS1_wrt_dni = -2.0 * static_cast<double>(WENO35.IS_1_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) + static_cast<double>(WENO35.IS_1_coeff_2[binIdx_i]) * (2.0 * yCrystal[binIdx_i] - yCrystal[binIdx_i + 1] - yCrystal[binIdx_i - 1]) - 2.0 * static_cast<double>(WENO35.IS_1_coeff_3[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);
			jacParams.dIS1_wrt_dni_m1 = 2.0 * static_cast<double>(WENO35.IS_1_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 1] - yCrystal[binIdx_i]) + static_cast<double>(WENO35.IS_1_coeff_2[binIdx_i]) * (yCrystal[binIdx_i + 1] - yCrystal[binIdx_i]);

			jacParams.dIS2_wrt_dni = static_cast<double>(WENO35.IS_2_coeff_2[binIdx_i]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + 2.0 * static_cast<double>(WENO35.IS_2_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);
			jacParams.dIS2_wrt_dni_m1 = -2.0 * static_cast<double>(WENO35.IS_2_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + static_cast<double>(WENO35.IS_2_coeff_2[binIdx_i]) * (2.0 * yCrystal[binIdx_i - 1] - yCrystal[binIdx_i] - yCrystal[binIdx_i - 2]) - 2.0 * static_cast<double>(WENO35.IS_2_coeff_3[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);
			jacParams.dIS2_wrt_dni_m2 = 2.0 * static_cast<double>(WENO35.IS_2_coeff_1[binIdx_i]) * (yCrystal[binIdx_i - 2] - yCrystal[binIdx_i - 1]) + static_cast<double>(WENO35.IS_2_coeff_2[binIdx_i]) * (yCrystal[binIdx_i] - yCrystal[binIdx_i - 1]);

			jacParams.dq0_wrt_dni_p2 = -static_cast<double>(WENO35.q_0_coeff_2[binIdx_i]);
			jacParams.dq0_wrt_dni_p1 = static_cast<double>(WENO35.q_0_coeff_2[binIdx_i]) - static_cast<double>(WENO35.q_0_coeff_1[binIdx_i]) + 1.0;
			jacParams.dq0_wrt_dni = static_cast<double>(WENO35.q_0_coeff_1[binIdx_i]);

			jacParams.dq1_wrt_dni_p1 = static_cast<double>(WENO35.q_1_coeff_1[binIdx_i]);
			jacParams.dq1_wrt_dni = -static_cast<double>(WENO35.q_1_coeff_1[binIdx_i]) + static_cast<double>(WENO35.q_1_coeff_2[binIdx_i]) + 1.0;
			jacParams.dq1_wrt_dni_m1 = -static_cast<double>(WENO35.q_1_coeff_2[binIdx_i]);

			jacParams.dq2_wrt_dni = static_cast<double>(WENO35.q_2_coeff_2[binIdx_i]);
			jacParams.dq2_wrt_dni_m1 = -static_cast<double>(WENO35.q_2_coeff_1[binIdx_i]) - static_cast<double>(WENO35.q_2_coeff_2[binIdx_i]) + 1.0;
			jacParams.dq2_wrt_dni_m2 = static_cast<double>(WENO35.q_2_coeff_1[binIdx_i]);

			// dw_0/dx
			jacParams.w0_w2_is2 = 2.0 * jacParams.W_0 * jacParams.W_2 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2);
			jacParams.w0_w1_is1 = 2.0 * jacParams.W_0 * jacParams.W_1 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1);
			jacParams.w0_w0_is0 = 2.0 * jacParams.W_0 * jacParams.W_0 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0);
			jacParams.w0_is0 = -2.0 * jacParams.W_0 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0);

			jacParams.dw0_wrt_dni_m2 = jacParams.w0_w2_is2 * jacParams.dIS2_wrt_dni_m2;
			jacParams.dw0_wrt_dni_m1 = jacParams.w0_w1_is1 * jacParams.dIS1_wrt_dni_m1 + jacParams.w0_w2_is2 * jacParams.dIS2_wrt_dni_m1;
			jacParams.dw0_wrt_dni = jacParams.w0_is0 * jacParams.dIS0_wrt_dni + jacParams.w0_w0_is0 * jacParams.dIS0_wrt_dni + jacParams.w0_w1_is1 * jacParams.dIS1_wrt_dni + jacParams.w0_w2_is2 * jacParams.dIS2_wrt_dni;
			jacParams.dw0_wrt_dni_p1 = jacParams.w0_is0 * jacParams.dIS0_wrt_dni_p1 + jacParams.w0_w0_is0 * jacParams.dIS0_wrt_dni_p1 + jacParams.w0_w1_is1 * jacParams.dIS1_wrt_dni_p1;
			jacParams.dw0_wrt_dni_p2 = jacParams.w0_is0 * jacParams.dIS0_wrt_dni_p2 + jacParams.w0_w0_is0 * jacParams.dIS0_wrt_dni_p2;

			// dw_1/dx
			jacParams.w1_w2_is2 = 2.0 * jacParams.W_1 * jacParams.W_2 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2);
			jacParams.w1_w1_is1 = 2.0 * jacParams.W_1 * jacParams.W_1 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1);
			jacParams.w1_w0_is0 = 2.0 * jacParams.W_1 * jacParams.W_0 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0);
			jacParams.w1_is1 = -2.0 * jacParams.W_1 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1);

			jacParams.dw1_wrt_dni_m2 = jacParams.w1_w2_is2 * jacParams.dIS2_wrt_dni_m2;
			jacParams.dw1_wrt_dni_m1 = jacParams.w1_is1 * jacParams.dIS1_wrt_dni_m1 + jacParams.w1_w1_is1 * jacParams.dIS1_wrt_dni_m1 + jacParams.w1_w2_is2 * jacParams.dIS2_wrt_dni_m1;
			jacParams.dw1_wrt_dni = jacParams.w1_is1 * jacParams.dIS1_wrt_dni + jacParams.w1_w0_is0 * jacParams.dIS0_wrt_dni + jacParams.w1_w1_is1 * jacParams.dIS1_wrt_dni + jacParams.w1_w2_is2 * jacParams.dIS2_wrt_dni;
			jacParams.dw1_wrt_dni_p1 = jacParams.w1_is1 * jacParams.dIS1_wrt_dni_p1 + jacParams.w1_w0_is0 * jacParams.dIS0_wrt_dni_p1 + jacParams.w1_w1_is1 * jacParams.dIS1_wrt_dni_p1;
			jacParams.dw1_wrt_dni_p2 = jacParams.w1_w0_is0 * jacParams.dIS0_wrt_dni_p2;

			// dw2_dx
			jacParams.w2_is2 = -2.0 * jacParams.W_2 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2);
			jacParams.w2_w2_is2 = 2.0 * jacParams.W_2 * jacParams.W_2 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_2);
			jacParams.w2_w1_is1 = 2.0 * jacParams.W_2 * jacParams.W_1 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_1);
			jacParams.w2_w0_is0 = 2.0 * jacParams.W_2 * jacParams.W_0 / (static_cast<double>(_binSizes[binIdx_i]) + jacParams.IS_0);

			jacParams.dw2_wrt_dni_m2 = jacParams.w2_is2 * jacParams.dIS2_wrt_dni_m2 + jacParams.w2_w2_is2 * jacParams.dIS2_wrt_dni_m2;
			jacParams.dw2_wrt_dni_m1 = jacParams.w2_is2 * jacParams.dIS2_wrt_dni_m1 + jacParams.w2_w1_is1 * jacParams.dIS1_wrt_dni_m1 + jacParams.w2_w2_is2 * jacParams.dIS2_wrt_dni_m1;
			jacParams.dw2_wrt_dni = jacParams.w2_is2 * jacParams.dIS2_wrt_dni + jacParams.w2_w0_is0 * jacParams.dIS0_wrt_dni + jacParams.w2_w1_is1 * jacParams.dIS1_wrt_dni + jacParams.w2_w2_is2 * jacParams.dIS2_wrt_dni;
			jacParams.dw2_wrt_dni_p1 = jacParams.w2_w0_is0 * jacParams.dIS0_wrt_dni_p1 + jacParams.w2_w1_is1 * jacParams.dIS1_wrt_dni_p1;
			jacParams.dw2_wrt_dni_p2 = jacParams.w2_w0_is0 * jacParams.dIS0_wrt_dni_p2;

			// dQ_2/dc
			jac[-3] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0 * jacParams.q_0 + jacParams.W_1 * jacParams.q_1 + jacParams.W_2 * jacParams.q_2);

			// dQ_2/dn0
			jac[-2] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_m2 * jacParams.q_0 + jacParams.dw1_wrt_dni_m2 * jacParams.q_1 + jacParams.dw2_wrt_dni_m2 * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni_m2); // the second and fourth terms are zero

			// dQ_2/dn1
			jac[-1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_m1 * jacParams.q_0 + jacParams.dw1_wrt_dni_m1 * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni_m1 + jacParams.dw2_wrt_dni_m1 * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni_m1) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // the second term is zero

			// dQ_2/dn2
			jac[0] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni + jacParams.dw1_wrt_dni * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni + jacParams.dw2_wrt_dni * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binSizes[binIdx_i]) * (1.0 / static_cast<double>(_binCenterDists[binIdx_i]) + 1.0 / static_cast<double>(_binCenterDists[binIdx_i - 1])); // full six terms of the product

			// dQ_2/dn3
			jac[1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_p1 * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni_p1 + jacParams.dw1_wrt_dni_p1 * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni_p1 + jacParams.dw2_wrt_dni_p1 * jacParams.q_2) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]); // the last term is zero

			// dQ_2/dn4
			jac[2] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_p2 * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni_p2 + jacParams.dw1_wrt_dni_p2 * jacParams.q_1 + jacParams.dw2_wrt_dni_p2 * jacParams.q_2); // the fourth and sixth terms are zero

			// dQ_2/dc_eq
			jac[_nComp - 4] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0 * jacParams.q_0 + jacParams.W_1 * jacParams.q_1 + jacParams.W_2 * jacParams.q_2);


		}

		// Q_{N_x-2}, special case for WENO35
		else if (cadet_unlikely(i == _nComp - 3))
		{
			binIdx_i = i - 1;
			// left face is weno5
			// dQ_{N_x-2}/dc
			jac[-i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0 * jacParams.q_0 + jacParams.W_1 * jacParams.q_1 + jacParams.W_2 * jacParams.q_2);

			// dQ_{N_x-1}/dn_{N_x-4}
			jac[-3] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_m2 * jacParams.q_0 + jacParams.dw1_wrt_dni_m2 * jacParams.q_1 + jacParams.dw2_wrt_dni_m2 * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni_m2); // the second and fourth terms are zero

			// dQ_{N_x-1}/dn_{N_x-3}
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_m1 * jacParams.q_0 + jacParams.dw1_wrt_dni_m1 * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni_m1 + jacParams.dw2_wrt_dni_m1 * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni_m1); // the second term is zero

			// dQ_{N_x-1}/dn_{N_x-2}
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni + jacParams.dw1_wrt_dni * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni + jacParams.dw2_wrt_dni * jacParams.q_2 + jacParams.W_2 * jacParams.dq2_wrt_dni); // full six terms of the product

			// dQ_{N_x-1}/dn_{N_x-1}
			jac[0] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_p1 * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni_p1 + jacParams.dw1_wrt_dni_p1 * jacParams.q_1 + jacParams.W_1 * jacParams.dq1_wrt_dni_p1 + jacParams.dw2_wrt_dni_p1 * jacParams.q_2); // the last term is zero

			// dQ_{N_x-1}/dn_{N_x}
			jac[1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_wrt_dni_p2 * jacParams.q_0 + jacParams.W_0 * jacParams.dq0_wrt_dni_p2 + jacParams.dw1_wrt_dni_p2 * jacParams.q_1 + jacParams.dw2_wrt_dni_p2 * jacParams.q_2); // the fourth and sixth terms are zero

			// dQ_2/dc_eq
			jac[_nComp - 1 - i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0 * jacParams.q_0 + jacParams.W_1 * jacParams.q_1 + jacParams.W_2 * jacParams.q_2);

			// right face is weno3
			jacParams.dvG_dc_right = jacParams.dvG_dc_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.dvG_dceq_right = jacParams.dvG_dceq_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));
			jacParams.vg_right = jacParams.vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_bins[binIdx_i + 1]), static_cast<double>(_p)));

			// dQ_{N_x-2}/dc
			jac[-i] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0_weno3_right * jacParams.q_0_weno3_right + jacParams.W_1_weno3_right * jacParams.q_1_weno3_right);

			// dQ_{N_x-2}/dn_{N_x-3}
			jac[-1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_m1_weno3_right * jacParams.q_0_weno3_right + jacParams.dw1_right_dni_m1_weno3_right * jacParams.q_1_weno3_right + jacParams.W_1_weno3_right * jacParams.dq1_right_dni_m1_weno3_right) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // the second term in the product rule is zero

			// dQ_{N_x-2}/dn_{N_x-2}
			jac[0] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_weno3_right * jacParams.q_0_weno3_right + jacParams.W_0_weno3_right * jacParams.dq0_right_dni_weno3_right + jacParams.dw1_right_dni_weno3_right * jacParams.q_1_weno3_right + jacParams.W_1_weno3_right * jacParams.dq1_right_dni_weno3_right) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binSizes[binIdx_i]) * (1.0 / static_cast<double>(_binCenterDists[binIdx_i]) + 1.0 / static_cast<double>(_binCenterDists[binIdx_i - 1])); // full four terms in the product rule

			// dQ_{N_x-2}/dn_{N_x-1}
			jac[1] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_p1_weno3_right * jacParams.q_0_weno3_right + jacParams.W_0_weno3_right * jacParams.dq0_right_dni_p1_weno3_right + jacParams.dw1_right_dni_p1_weno3_right * jacParams.q_1_weno3_right) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i]) / static_cast<double>(_binSizes[binIdx_i]); // the fourth term in the product rule is zero

			// dQ_{N_x-2}/dc_eq
			jac[_nComp - 1 - i] -= factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0_weno3_right * jacParams.q_0_weno3_right + jacParams.W_1_weno3_right * jacParams.q_1_weno3_right);
		}

		// Q_{N_x-1}, right BC
		else if (cadet_unlikely(i == _nComp - 2))
		{
			binIdx_i = i - 1;
			// left face is weno3
			// dQ_{N_x-1}/dc
			jac[-i] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dc_right * (jacParams.W_0_weno3_right * jacParams.q_0_weno3_right + jacParams.W_1_weno3_right * jacParams.q_1_weno3_right);

			// dQ_{N_x-1}/dn_{N_x-3}
			jac[-2] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_m1_weno3_right * jacParams.q_0_weno3_right + jacParams.dw1_right_dni_m1_weno3_right * jacParams.q_1_weno3_right + jacParams.W_1_weno3_right * jacParams.dq1_right_dni_m1_weno3_right); // the second term in the product rule is zero

			// dQ_{N_x-1}/dn_{N_x-2}
			jac[-1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_weno3_right * jacParams.q_0_weno3_right + jacParams.W_0_weno3_right * jacParams.dq0_right_dni_weno3_right + jacParams.dw1_right_dni_weno3_right * jacParams.q_1_weno3_right + jacParams.W_1_weno3_right * jacParams.dq1_right_dni_weno3_right) + factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // full four terms in the product rule

			// dQ_{N_x-1}/dn_{N_x-1}
			jac[0] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.vg_right * (jacParams.dw0_right_dni_p1_weno3_right * jacParams.q_0_weno3_right + jacParams.W_0_weno3_right * jacParams.dq0_right_dni_p1_weno3_right + jacParams.dw1_right_dni_p1_weno3_right * jacParams.q_1_weno3_right) - factor * static_cast<double>(_growthDispersionRate) / static_cast<double>(_binCenterDists[binIdx_i - 1]) / static_cast<double>(_binSizes[binIdx_i]); // the fourth term in the product rule is zero

			// dQ_i/dc_{eq}
			jac[1] += factor / static_cast<double>(_binSizes[binIdx_i]) * jacParams.dvG_dceq_right * (jacParams.W_0_weno3_right * jacParams.q_0_weno3_right + jacParams.W_1_weno3_right * jacParams.q_1_weno3_right);
		}
	}

	template <typename RowIterator>
	void jacobianFluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nStates, double const* y, double factor, RowIterator& jac, LinearBufferAllocator workSpace) const
	{
		if (_mode.hasPBM())
		{
			double const* const yCrystal = y + 1; // pointer to crystal bins. Jump over the solute entry, which only exists if we solve the mass balance
			// supersaturation s, rewrite it to 0.0 if it drops below 0.0
			double const s = (cadet_likely(y[0] / y[_nComp - 1] - 1.0 > 0)) ? y[0] / y[_nComp - 1] - 1.0 : 0.0;

			// compute the summations 
			double M = 0.0;
			double substrateConversion = 0.0;
			for (int i = 0; i < _nBins; ++i)
			{
				M += yCrystal[i] * static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) * static_cast<double>(_binSizes[i]);
				substrateConversion += yCrystal[i] * static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_binCenters[i]), static_cast<double>(_p))) * static_cast<double>(_binSizes[i]);
			}

			// define variables to aid jacobian implementation
			double const vG_factor = static_cast<double>(_growthRateConstant) * pow(s, static_cast<double>(_g));

			double const rho_kv = static_cast<double>(_nucleiMassDensity) * static_cast<double>(_volShapeFactor);
			double const x_c_3 = static_cast<double>(_bins[0]) * static_cast<double>(_bins[0]) * static_cast<double>(_bins[0]);

			// rewrite the derivatives to 0.0 if s = 0.0
			double const dBp_dc = cadet_unlikely(s > 0.0) ? static_cast<double>(_primaryNucleationRate) * static_cast<double>(_u) * pow(s, static_cast<double>(_u) - 1.0) / y[_nComp - 1] : 0.0;
			double const dBp_dceq = cadet_unlikely(s > 0.0) ? -static_cast<double>(_primaryNucleationRate) * static_cast<double>(_u) * pow(s, static_cast<double>(_u) - 1.0) * y[0] / y[_nComp - 1] / y[_nComp - 1] : 0.0;
			double const dBs_dc = cadet_unlikely(s > 0.0) ? rho_kv * static_cast<double>(_secondaryNucleationRate) * static_cast<double>(_b) * pow(s, static_cast<double>(_b) - 1.0) * M / y[_nComp - 1] : 0.0;
			double const dBs_dceq = cadet_unlikely(s > 0.0) ? -rho_kv * static_cast<double>(_secondaryNucleationRate) * static_cast<double>(_b) * pow(s, static_cast<double>(_b) - 1.0) * M * y[0] / y[_nComp - 1] / y[_nComp - 1] : 0.0;
			double const dvG_dc_factor = cadet_unlikely(s > 0.0) ? static_cast<double>(_growthRateConstant) * static_cast<double>(_g) * pow(s, static_cast<double>(_g) - 1.0) / y[_nComp - 1] : 0.0;
			double const dvG_dceq_factor = cadet_unlikely(s > 0.0) ? -static_cast<double>(_growthRateConstant) * static_cast<double>(_g) * pow(s, static_cast<double>(_g) - 1.0) * y[0] / y[_nComp - 1] / y[_nComp - 1] : 0.0;
			double const dBs_dni_factor = cadet_unlikely(s > 0.0) ? rho_kv * static_cast<double>(_secondaryNucleationRate) * pow(s, static_cast<double>(_b)) : 0.0;

			// the jacobian has a shape: (_nComp) x (_nComp), undefined ones are 0.0.
			// the first loop i iterates over rows, the second loop j iterates over columns. The offset i is used to move the jac index to 0 at the beginning of iterating over j.
			int binIdx_i = 0;
			int binIdx_j = 0;

			// Q_c, mass balance, independent of the discretization method
			for (int j = 0; j < _nComp; ++j)
			{
				binIdx_j = j - 1;
				if (cadet_likely((j > 0) && (j + 1 < _nComp)))
				{
					// dQ_c/dn_i
					jac[j] -= factor * rho_kv * (x_c_3 * dBs_dni_factor * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j]) + 3.0 * vG_factor * (static_cast<double>(_a) + static_cast<double>(_growthConstant) * pow(static_cast<double>(_binCenters[binIdx_j]), static_cast<double>(_p))) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binCenters[binIdx_j]) * static_cast<double>(_binSizes[binIdx_j]));
				}
				else if (cadet_unlikely(j == 0))
				{
					// dQ_c/dc
					jac[j] -= factor * rho_kv * ((dBp_dc + dBs_dc) * x_c_3 + 3.0 * dvG_dc_factor * substrateConversion);
				}
				else
				{
					// dQ_c/dc_eq
					jac[j] -= factor * rho_kv * ((dBp_dceq + dBs_dceq) * x_c_3 + 3.0 * dvG_dceq_factor * substrateConversion);
				}
			}
			++jac;

			if (_growthSchemeOrder == detail::PBMReconstruction::WENO35)
			{
				static_cast<JacobianParamsWENO35*>(_jacParams)->updateParams(vG_factor, dBp_dc, dBp_dceq, dBs_dc, dBs_dceq, dvG_dc_factor, dvG_dceq_factor, dBs_dni_factor, yCrystal, _nComp, _binSizes, *static_cast<WENO35Params<active>*>(_reconstruction));
			}
			else
				_jacParams->updateParams(vG_factor, dBp_dc, dBp_dceq, dBs_dc, dBs_dceq, dvG_dc_factor, dvG_dceq_factor, dBs_dni_factor);
		}

		double const* const yCrystal = y; // pointer to crystal bins

		// main loop iterating over Jacobian rows
		for (int l = 1; l <= _nBins; ++l)
		{
			const int i = l - 1; // index for fragmetation and aggregation

			if (_mode.hasPBM())
			{
				switch (_growthSchemeOrder)
				{
				case detail::PBMReconstruction::Upwind:
				{
					jacobianUpwindKernel<RowIterator>(yCrystal + 1, factor, jac, *_jacParams, l);
				}
				break;
				case detail::PBMReconstruction::HRKoren:
				{
					jacobianHRKorenKernel<RowIterator>(yCrystal + 1, factor, jac, *static_cast<JacobianParamsHRKoren*>(_jacParams), *static_cast<HRKorenParams<active>*>(_reconstruction), l);
				}
				break;
				case detail::PBMReconstruction::WENO23:
				{
					jacobianWENO23Kernel<RowIterator>(yCrystal + 1, factor, jac, *static_cast<JacobianParamsWENO23*>(_jacParams), *static_cast<WENO23Params<active>*>(_reconstruction), l);
				}
				break;
				case detail::PBMReconstruction::WENO35:
				{
					jacobianWENO35Kernel<RowIterator>(yCrystal + 1, factor, jac, *static_cast<JacobianParamsWENO35*>(_jacParams), *static_cast<WENO35Params<active>*>(_reconstruction), l);
				}
				break;
				}
			}

			if (_mode.hasAggregation())
			{
				// aggregation source term
				const double x_i_cube = static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]);

				for (int p = static_cast<int>(_agg->indexTracker[i]); p < static_cast<int>(_agg->indexTracker[i + 1]); ++p)
				{
					const int j = static_cast<int>(_agg->index.first[p]);
					const int k = static_cast<int>(_agg->index.second[p]);

					const double x_jCubeplusX_kCube = static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]) + static_cast<double>(_binCenters[k]) * static_cast<double>(_binCenters[k]) * static_cast<double>(_binCenters[k]);

					// source term
					double aggSourceFactor = static_cast<double>(_aggregationRateConstant) * static_cast<double>(_binSizes[j]) * static_cast<double>(_binSizes[k]) / static_cast<double>(_binSizes[i])
						* 1.0 / (2.0 * x_i_cube / x_jCubeplusX_kCube - 1.0);

					if (cadet_unlikely(j == k)) { aggSourceFactor *= 0.5; }

					// add different kernels
					switch (_aggregationIndex)
					{
						//constant kernel 0
					case 0:
						break;
						// brownian kernel 1
					case 1:
						aggSourceFactor *= (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) / static_cast<double>(_binCenters[k]) / static_cast<double>(_binCenters[j]);
						break;
						// smoluchowski kernel 2
					case 2:
						aggSourceFactor *= (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j]));
						break;
						// golovin kernel 3
					case 3:
						aggSourceFactor *= static_cast<double>(_binCenters[k]) * static_cast<double>(_binCenters[k]) * static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]);
						break;
						// differential force kernel 4
					case 4:
						aggSourceFactor *= (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[k]) * static_cast<double>(_binCenters[k]) - static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]));
						break;
					}

					// add to the jacobian
					//j=k case is covered: source_factor is already multiplied by 0.5
					jac[j - i] += factor * yCrystal[k] * aggSourceFactor;
					jac[k - i] += factor * yCrystal[j] * aggSourceFactor;
				}

				// aggregation sink terms
				for (int j = 0; j < _nBins; ++j)
				{
					const double sumVolume = pow(x_i_cube + static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]), 1.0 / 3.0);
					const int k = static_cast<int>(_agg->iIndex[i * _nBins + j]);
					double ratio{};
					double sinkCorrectionFactor = 1.0;
					if (k > 0)
					{
						ratio = sumVolume / static_cast<double>(_binCenters[k]);
						sinkCorrectionFactor = 1.0 / (2.0 - ratio * ratio * ratio);
					}

					double aggregationSinkFactor = static_cast<double>(_aggregationRateConstant);

					switch (_aggregationIndex)
					{
						//constant kernel 0
					case 0:
						aggregationSinkFactor *= static_cast<double>(_binSizes[j]) * sinkCorrectionFactor;
						break;
						// brownian kernel 1
					case 1:
						aggregationSinkFactor *= static_cast<double>(_binSizes[j]) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) / static_cast<double>(_binCenters[i]) / static_cast<double>(_binCenters[j]) * sinkCorrectionFactor;
						break;
						// smoluchowski kernel 2
					case 2:
						aggregationSinkFactor *= static_cast<double>(_binSizes[j]) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * sinkCorrectionFactor;
						break;
						// golovin kernel 3
					case 3:
						aggregationSinkFactor *= static_cast<double>(_binSizes[j]) * (static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j])) * sinkCorrectionFactor;
						break;
						// differential force kernel 4
					case 4:
						aggregationSinkFactor *= static_cast<double>(_binSizes[j]) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) - static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j])) * sinkCorrectionFactor;
						break;
					}

					jac[j - i] -= factor * aggregationSinkFactor * yCrystal[i];  // wrt j
					jac[0] -= factor * aggregationSinkFactor * yCrystal[j];  // wrt i
				}
			}

			if (_mode.hasFragmentation())
			{
				// jacobian, when adding to growth terms, remember to change the index to binIdx_i! Q_ceq is not considered.
				double selectionFunction = 0.0;
				double bIntegral = 0.0;

				const double N_j = static_cast<double>(_fragKernelGamma) / (static_cast<double>(_fragKernelGamma) - 1.0);

				// add fragmentation source to jac
				for (int j = i; j < _nBins; ++j)
				{
					// update selection function
					selectionFunction = static_cast<double>(_fragRateConstant) * pow(static_cast<double>(_binCenters[j]), 3.0 * static_cast<double>(_fragSelectionFunctionAlpha));
					if (cadet_likely(i != j))
					{
						bIntegral = N_j * (pow(static_cast<double>(_bins[i + 1]), 3.0 * static_cast<double>(_fragKernelGamma) - 3.0) - pow(static_cast<double>(_bins[i]), 3.0 * static_cast<double>(_fragKernelGamma) - 3.0)) / pow(static_cast<double>(_binCenters[j]), 3.0 * static_cast<double>(_fragKernelGamma) - 3.0);
					}
					else
					{
						bIntegral = N_j * (pow(static_cast<double>(_binCenters[i]), 3.0 * static_cast<double>(_fragKernelGamma) - 3.0) - pow(static_cast<double>(_bins[i]), 3.0 * static_cast<double>(_fragKernelGamma) - 3.0)) / pow(static_cast<double>(_binCenters[j]), 3.0 * static_cast<double>(_fragKernelGamma) - 3.0);
					}
					jac[-i + j] += factor * selectionFunction * bIntegral * static_cast<double>(_frag->upsilonSource[j]) * static_cast<double>(_binSizes[j]) / static_cast<double>(_binSizes[i]);
				}

				// add fragmentation sink to jac
				selectionFunction = static_cast<double>(_fragRateConstant) * pow(static_cast<double>(_binCenters[i]), 3.0 * static_cast<double>(_fragSelectionFunctionAlpha));
				jac[0] -= factor * selectionFunction * static_cast<double>(_frag->upsilonSink[i]);
			}

			// go to the next row
			++jac;
		}

		if (_mode.hasPBM())
			++jac; // the last row corresponding to Q_{ceq} is zero, and the interface expects the iterator to be moved to the end.

	}

	template <typename RowIteratorLiquid, typename RowIteratorSolid>
	void jacobianCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid, double factor, RowIteratorLiquid& jacLiquid, RowIteratorSolid& jacSolid, LinearBufferAllocator workSpace) const
	{
		jacobianFluxImpl<RowIteratorLiquid>(t, secIdx, colPos, 0, yLiquid, factor, jacLiquid, workSpace);
	}
};

namespace reaction
{
	void registerCrystallizationReaction(std::unordered_map<std::string, std::function<model::IDynamicReactionModel* ()>>& reactions)
	{
		reactions[CrystallizationReaction::identifier()] = []() { return new CrystallizationReaction(); };
	}
}  // namespace reaction

}  // namespace model

}  // namespace cadet

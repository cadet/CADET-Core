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

namespace cadet
{
namespace model
{
namespace detail
{
	enum class CrystallizationMode : int
	{
		/**
		 * Population mass balance equation
		 */
		PurePBM,

		/**
		 * Aggregation equation
		 */
		PureAggregation,

		/**
		 * Breakage equation
		 */
		PureBreakage,

		/**
		 * Combined aggregation and breakage equation
		 */
		AggregationBreakage
	};

	struct ModeFlags
	{
		bool hasMassBalance = false;
		bool hasAggregation = false;
		bool hasFragmentation = false;
	};

	std::unordered_map<CrystallizationMode, ModeFlags> modeFlagMap = {
		{CrystallizationMode::PurePBM, {true, false, false}},
		{CrystallizationMode::PureAggregation, {true, true, false}},
		{CrystallizationMode::PureBreakage, {false, false, true}},
		{CrystallizationMode::AggregationBreakage, {true, true, true}}
	};

	struct AggCoefficients
	{
		std::pair<std::vector<unsigned short int>, std::vector<unsigned short int>> index;  // a pair of vectors to store combínation of j and k
		std::vector<unsigned int> index_tracker;           // the idx can be a large number  
		std::vector<short int> i_index;

		active sum_volume_cube_root = 0.0;
		unsigned int count = 0;

		// the sizes of j_index and k_index are not known a priori
		AggCoefficients(const std::vector<active>& binCenters, const std::vector<active>& bins, const std::vector<active>& binSizes)
			: index(), index_tracker(binCenters.size() + 1), i_index(binCenters.size()* binCenters.size())
		{
			// make the first element of index_tracker 0
			index_tracker[0] = 0;

			// loop over daughter particle i
			for (unsigned int i = 0; i < binCenters.size(); ++i)
			{
				// source term, find lij in set A
				for (int j = 0; j < i; ++j)
				{
					for (int k = j; k < i + 1; ++k)
					{
						sum_volume_cube_root = pow(binCenters[k] * binCenters[k] * binCenters[k] + binCenters[j] * binCenters[j] * binCenters[j], 1.0 / 3.0);
						if ((sum_volume_cube_root > bins[i]) && (sum_volume_cube_root <= bins[i + 1]))
						{
							index.first.emplace_back(static_cast<short int>(j));
							index.second.emplace_back(static_cast<short int>(k));
							++count;
						}
					}
				}
				index_tracker[i + 1] = count;

				// sink term, loop over another mother particle j
				for (unsigned int j = 0; j < binCenters.size(); ++j)
				{
					sum_volume_cube_root = pow(binCenters[i] * binCenters[i] * binCenters[i] + binCenters[j] * binCenters[j] * binCenters[j], 1.0 / 3.0);
					int k = (i > j) ? i + 1 : j + 1;

					i_index[i * binCenters.size() + j] = -1;        // -1 indicates there is no match
					for (k; k < binCenters.size(); ++k)
					{
						if ((sum_volume_cube_root > bins[k]) && (sum_volume_cube_root <= bins[k + 1]))
						{
							i_index[i * binCenters.size() + j] = k;  // if there is a match, overwrite it
							break;
						}
					}
				}
			}
		}
	};

	struct FragCoefficients
	{
		// define breakage-related local parameters
		std::vector<active> Upsilon_source;
		std::vector<active> Upsilon_sink;

		// constructor
		FragCoefficients(const std::vector<active>& binCenters, const std::vector<active>& bins, const active& breakageKernelGamma)
			: Upsilon_source(binCenters.size()), Upsilon_sink(binCenters.size())
		{
			const active N_j = breakageKernelGamma / (breakageKernelGamma - 1.0);

			active Upsilon_birth_sum = 0.0;
			active Upsilon_death_sum = 0.0;
			active b_integral_birth = 0.0;

			// calculate upsilon and store it
			for (int i = 0; i < binCenters.size(); ++i)
			{
				// reset the sum for each i
				Upsilon_birth_sum = 0.0;
				Upsilon_death_sum = 0.0;
				const active x_i_3 = binCenters[i] * binCenters[i] * binCenters[i];
				for (int j = 0; j < i + 1; ++j)
				{
					const active x_j_3 = binCenters[j] * binCenters[j] * binCenters[j];
					if (cadet_likely(i != j))
					{
						b_integral_birth = N_j * (pow(bins[j + 1], 3.0 * breakageKernelGamma - 3.0) - pow(bins[j], 3.0 * breakageKernelGamma - 3.0)) / pow(binCenters[i], 3.0 * breakageKernelGamma - 3.0);
						Upsilon_birth_sum += b_integral_birth * (x_i_3 - x_j_3);
						Upsilon_death_sum += b_integral_birth * x_j_3;
					}
					else
					{
						b_integral_birth = N_j * (pow(binCenters[j], 3.0 * breakageKernelGamma - 3.0) - pow(bins[j], 3.0 * breakageKernelGamma - 3.0)) / pow(binCenters[i], 3.0 * breakageKernelGamma - 3.0);
						Upsilon_death_sum += b_integral_birth * x_j_3;
					}
				}
				if (cadet_likely(i > 0))
				{
					Upsilon_source[i] = (N_j - 1.0) * x_i_3 / Upsilon_birth_sum;
				}
				else
				{
					Upsilon_source[i] = (N_j - 1.0) * x_i_3;
				}
				Upsilon_sink[i] = Upsilon_death_sum * Upsilon_source[i] / x_i_3;
			}
		}
	};

} // namespace detail

/**
* @brief Defines the crystallization reaction model
*/
class CrystallizationReaction : public IDynamicReactionModel
{
public:

	CrystallizationReaction() : _nComp(0), _nBins(0), _bins(0), _binCenters(0), _binSizes(0), _agg(nullptr), _frag(nullptr) { }
	virtual ~CrystallizationReaction() CADET_NOEXCEPT
	{
		clearSchemeCoefficients();
	}

	static const char* identifier() { return "CRYSTALLIZATION"; }
	virtual const char* name() const CADET_NOEXCEPT { return identifier(); }

	virtual bool requiresConfiguration() const CADET_NOEXCEPT { return true; }
	virtual bool usesParamProviderInDiscretizationConfig() const CADET_NOEXCEPT { return false; }

	virtual bool configure(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		readScalarParameterOrArray(_bins, paramProvider, "CRY_BINS", 1);

		if (paramProvider.exists("CRYSTALLIZATION_MODE"))
			_mode = static_cast<detail::CrystallizationMode>(paramProvider.getInt("CRYSTALLIZATION_MODE"));
		else
			_mode = detail::CrystallizationMode::PurePBM; // For backwards compatibility

		_usePBM = detail::modeFlagMap[_mode].hasMassBalance;
		_useAgg = detail::modeFlagMap[_mode].hasAggregation;
		_useFrag = detail::modeFlagMap[_mode].hasFragmentation;

		if (_bins.size() != _nBins + 1)
			throw InvalidParameterException("Expected CRY_BINS to have " + std::to_string(_nBins + 1) + " elements (got " + std::to_string(_bins.size()) + ")");

		registerParam1DArray(_parameters, _bins, [=](bool multi, unsigned int idx) { return makeParamId(hashString("CRY_BINS"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, idx, SectionIndep); });

		_binCenters.resize(_nBins);
		_binSizes.resize(_nBins);
		updateBinCoords();

		_aggregationRateConstant = paramProvider.getDouble("CRY_AGGREGATION_RATE_CONSTANT");
		_parameters[makeParamId(hashString("CRY_AGGREGATION_RATE_CONSTANT"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_aggregationRateConstant;

		_aggregationIndex = paramProvider.getInt("CRY_AGGREGATION_INDEX");

		_breakageRateConstant = paramProvider.getDouble("CRY_BREAKAGE_RATE_CONSTANT");
		_parameters[makeParamId(hashString("CRY_BREAKAGE_RATE_CONSTANT"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_breakageRateConstant;

		_breakageKernelGamma = paramProvider.getDouble("CRY_BREAKAGE_KERNEL_GAMMA");
		if (_breakageKernelGamma <= 1.0)
			throw InvalidParameterException("CRY_BREAKAGE_KERNEL_GAMMA needs to be larger than 1.0");
		_parameters[makeParamId(hashString("CRY_BREAKAGE_KERNEL_GAMMA"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_breakageKernelGamma;

		_breakageSelectionFunctionAlpha = paramProvider.getDouble("CRY_BREAKAGE_SELECTION_FUNCTION_ALPHA");
		_parameters[makeParamId(hashString("CRY_BREAKAGE_SELECTION_FUNCTION_ALPHA"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_breakageSelectionFunctionAlpha;

		clearSchemeCoefficients();

		if (_breakageRateConstant != 0.0)
			_frag = new detail::FragCoefficients(_binCenters, _bins, _breakageKernelGamma);

		if (_aggregationRateConstant != 0.0)
			_agg = new detail::AggCoefficients(_binCenters, _bins, _binSizes);

		return true;
	}

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		_nComp = nComp;
		_nBins = _nComp;

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

	std::unordered_map<ParameterId, active*> _parameters; //!< Map used to translate ParameterIds to actual variables
	cadet::model::detail::CrystallizationMode _mode; //!< Crystallization mode, i.e. specification of considered effects
	bool _usePBM; //!< Apply population mass balance
	bool _useFrag; //!< Apply fragmentation term
	bool _useAgg; //!< Apply aggregation term
	int _nComp; //!< Number of components
	int _nBins; //!< Number of crystal size bins

	std::vector<active> _bins;
	std::vector<active> _binCenters;
	std::vector<active> _binSizes;

	active _breakageRateConstant; // constant breakage rate constant
	active _breakageKernelGamma; // gamma in the breakage kernel
	active _breakageSelectionFunctionAlpha; // alpha in the selection function

	detail::FragCoefficients* _frag;

	int _aggregationIndex; // determines which kernel to use
	active _aggregationRateConstant; // aggregation rate constant

	detail::AggCoefficients* _agg;

	void updateBinCoords() CADET_NOEXCEPT
	{
		for (int i = 0; i < _nBins; ++i)
		{
			_binCenters[i] = 0.5 * (_bins[i] + _bins[i + 1]);
			_binSizes[i] = _bins[i + 1] - _bins[i];
		}
	}

	void clearSchemeCoefficients() CADET_NOEXCEPT
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
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	int residualLiquidImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* y, ResidualType* res, const FactorType& factor, LinearBufferAllocator workSpace) const
	{
		typedef typename DoubleActivePromoter<StateType, ParamType>::type StateParam;
		// Pointer to crystal bins
		StateType const* const yCrystal = y;
		ResidualType* const resCrystal = res;

		// define aggregation-related local parameters
		StateParam source = 0.0;
		StateParam sink = 0.0;
		ParamType source_factor = 0.0;

		for (int i = 0; i < _nBins; ++i)
		{
			if (_useAgg)
			{
				// reset the source and sink terms
				source = 0.0;
				sink = 0.0;

				const ParamType x_i_cube = static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]);

				// aggregation source terms
				for (int p = static_cast<int>(_agg->index_tracker[i]); p < static_cast<int>(_agg->index_tracker[i + 1]); ++p)
				{
					const int j = static_cast<int>(_agg->index.first[p]);
					const int k = static_cast<int>(_agg->index.second[p]);

					const ParamType x_j_cube_plus_x_k_cube = static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) + static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]);

					// source term
					source_factor = static_cast<ParamType>(_aggregationRateConstant) * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i])
						* 1.0 / (2.0 * x_i_cube / x_j_cube_plus_x_k_cube - 1.0);

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
					const ParamType sum_volume = pow(x_i_cube + static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]), 1.0 / 3.0);
					const int k = static_cast<int>(_agg->i_index[i * _nBins + j]);

					ParamType ratio{};
					ParamType sink_correction_factor = 1.0;
					if (k > 0)
					{
						ratio = sum_volume / static_cast<ParamType>(_binCenters[k]);
						sink_correction_factor = 1.0 / (2.0 - ratio * ratio * ratio);
					}

					switch (_aggregationIndex)
					{
						//constant kernel 0
					case 0:
						sink += yCrystal[j] * static_cast<ParamType>(_binSizes[j]) * sink_correction_factor;
						break;
						// brownian kernel 1
					case 1:
						sink += yCrystal[j] * static_cast<ParamType>(_binSizes[j]) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) / static_cast<ParamType>(_binCenters[i]) / static_cast<ParamType>(_binCenters[j]) * sink_correction_factor;
						break;
						// smoluchowski kernel 2
					case 2:
						sink += yCrystal[j] * static_cast<ParamType>(_binSizes[j]) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * sink_correction_factor;
						break;
						// golovin kernel 3
					case 3:
						sink += yCrystal[j] * static_cast<ParamType>(_binSizes[j]) * (static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * sink_correction_factor;
						break;
						// differential force kernel 4
					case 4:
						sink += yCrystal[j] * static_cast<ParamType>(_binSizes[j]) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) - static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * sink_correction_factor;
						break;
					}
				}
				sink *= yCrystal[i] * static_cast<ParamType>(_aggregationRateConstant);

				// add to residual
				resCrystal[i] += factor * source - factor * sink;
			}

			if (_useFrag)
			{
				const ParamType N_j = static_cast<ParamType>(_breakageKernelGamma) / (static_cast<ParamType>(_breakageKernelGamma) - 1.0);

				// ode input
				ParamType b_integral = 0.0;
				ParamType selection_function = 0.0;
				StateParam breakage_source = 0.0;
				StateParam breakage_sink = 0.0;

				// source term
				breakage_source = 0.0;
				for (int j = i; j < _nBins; ++j)
				{
					// selection function
					selection_function = static_cast<ParamType>(_breakageRateConstant) * pow(static_cast<ParamType>(_binCenters[j]), 3.0 * static_cast<ParamType>(_breakageSelectionFunctionAlpha));
					if (cadet_likely(i != j))
					{
						b_integral = N_j * (pow(static_cast<ParamType>(_bins[i + 1]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0) - pow(static_cast<ParamType>(_bins[i]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0)) / pow(static_cast<ParamType>(_binCenters[j]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0);
					}
					else
					{
						b_integral = N_j * (pow(static_cast<ParamType>(_binCenters[i]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0) - pow(static_cast<ParamType>(_bins[i]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0)) / pow(static_cast<ParamType>(_binCenters[j]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0);
					}
					breakage_source += selection_function * yCrystal[j] * b_integral * static_cast<ParamType>(_frag->Upsilon_source[j]) * static_cast<ParamType>(_binSizes[j]) / static_cast<ParamType>(_binSizes[i]);
				}

				// sink term
				selection_function = static_cast<ParamType>(_breakageRateConstant) * pow(static_cast<ParamType>(_binCenters[i]), 3.0 * static_cast<ParamType>(_breakageSelectionFunctionAlpha));
				breakage_sink = yCrystal[i] * selection_function * static_cast<ParamType>(_frag->Upsilon_sink[i]);

				// add to residual
				resCrystal[i] += factor * breakage_source - factor * breakage_sink;
			}
		}
		return 0;
	}

	template <typename StateType, typename ResidualType, typename ParamType>
	int residualCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* yLiquid, StateType const* ySolid, ResidualType* resLiquid, ResidualType* resSolid, double factor, LinearBufferAllocator workSpace) const
	{
		return residualLiquidImpl<StateType, ResidualType, ParamType, double>(t, secIdx, colPos, yLiquid, resLiquid, factor, workSpace);
	}

	template <typename RowIterator>
	void jacobianLiquidImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double factor, RowIterator& jac, LinearBufferAllocator workSpace) const
	{
		// Pointer to crystal bins
		double const* const yCrystal = y;

		// jacobian, when adding to growth terms
		for (int i = 0; i < _nBins; ++i)
		{
			if (_useAgg)
			{
				// aggregation source term
				const double x_i_cube = static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]);

				for (int p = static_cast<int>(_agg->index_tracker[i]); p < static_cast<int>(_agg->index_tracker[i + 1]); ++p)
				{
					const int j = static_cast<int>(_agg->index.first[p]);
					const int k = static_cast<int>(_agg->index.second[p]);

					const double x_j_cube_plus_x_k_cube = static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]) + static_cast<double>(_binCenters[k]) * static_cast<double>(_binCenters[k]) * static_cast<double>(_binCenters[k]);

					// source term
					double agg_source_factor = static_cast<double>(_aggregationRateConstant) * static_cast<double>(_binSizes[j]) * static_cast<double>(_binSizes[k]) / static_cast<double>(_binSizes[i])
						* 1.0 / (2.0 * x_i_cube / x_j_cube_plus_x_k_cube - 1.0);

					if (cadet_unlikely(j == k)) { agg_source_factor *= 0.5; }

					// add different kernels
					switch (_aggregationIndex)
					{
						//constant kernel 0
					case 0:
						break;
						// brownian kernel 1
					case 1:
						agg_source_factor *= (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) / static_cast<double>(_binCenters[k]) / static_cast<double>(_binCenters[j]);
						break;
						// smoluchowski kernel 2
					case 2:
						agg_source_factor *= (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j]));
						break;
						// golovin kernel 3
					case 3:
						agg_source_factor *= static_cast<double>(_binCenters[k]) * static_cast<double>(_binCenters[k]) * static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]);
						break;
						// differential force kernel 4
					case 4:
						agg_source_factor *= (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[k]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[k]) * static_cast<double>(_binCenters[k]) - static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]));
						break;
					}

					// add to the jacobian
					//j=k case is covered: source_factor is already multiplied by 0.5
					jac[j - i] += factor * yCrystal[k] * agg_source_factor;
					jac[k - i] += factor * yCrystal[j] * agg_source_factor;
				}

				// aggregation sink terms
				for (int j = 0; j < _nBins; ++j)
				{
					const double sum_volume = pow(x_i_cube + static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]), 1.0 / 3.0);
					const int k = static_cast<int>(_agg->i_index[i * _nBins + j]);
					double ratio{};
					double sink_correction_factor = 1.0;
					if (k > 0)
					{
						ratio = sum_volume / static_cast<double>(_binCenters[k]);
						sink_correction_factor = 1.0 / (2.0 - ratio * ratio * ratio);
					}

					double aggregation_sink_factor = static_cast<double>(_aggregationRateConstant);

					switch (_aggregationIndex)
					{
						//constant kernel 0
					case 0:
						aggregation_sink_factor *= static_cast<double>(_binSizes[j]) * sink_correction_factor;
						break;
						// brownian kernel 1
					case 1:
						aggregation_sink_factor *= static_cast<double>(_binSizes[j]) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) / static_cast<double>(_binCenters[i]) / static_cast<double>(_binCenters[j]) * sink_correction_factor;
						break;
						// smoluchowski kernel 2
					case 2:
						aggregation_sink_factor *= static_cast<double>(_binSizes[j]) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * sink_correction_factor;
						break;
						// golovin kernel 3
					case 3:
						aggregation_sink_factor *= static_cast<double>(_binSizes[j]) * (static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j])) * sink_correction_factor;
						break;
						// differential force kernel 4
					case 4:
						aggregation_sink_factor *= static_cast<double>(_binSizes[j]) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[i]) + static_cast<double>(_binCenters[j])) * (static_cast<double>(_binCenters[i]) * static_cast<double>(_binCenters[i]) - static_cast<double>(_binCenters[j]) * static_cast<double>(_binCenters[j])) * sink_correction_factor;
						break;
					}

					jac[j - i] -= factor * aggregation_sink_factor * yCrystal[i];  // wrt j
					jac[0] -= factor * aggregation_sink_factor * yCrystal[j];  // wrt i
				}
			}

			if (_useFrag)
			{
				// jacobian, when adding to growth terms, remember to change the index to binIdx_i! Q_ceq is not considered.
				double selection_function = 0.0;
				double b_integral = 0.0;

				const double N_j = static_cast<double>(_breakageKernelGamma) / (static_cast<double>(_breakageKernelGamma) - 1.0);

				for (int i = 0; i < _nBins; ++i)
				{
					// add breakage source to jac
					for (int j = i; j < _nBins; ++j)
					{
						// update selection function
						selection_function = static_cast<double>(_breakageRateConstant) * pow(static_cast<double>(_binCenters[j]), 3.0 * static_cast<double>(_breakageSelectionFunctionAlpha));
						if (cadet_likely(i != j))
						{
							b_integral = N_j * (pow(static_cast<double>(_bins[i + 1]), 3.0 * static_cast<double>(_breakageKernelGamma) - 3.0) - pow(static_cast<double>(_bins[i]), 3.0 * static_cast<double>(_breakageKernelGamma) - 3.0)) / pow(static_cast<double>(_binCenters[j]), 3.0 * static_cast<double>(_breakageKernelGamma) - 3.0);
						}
						else
						{
							b_integral = N_j * (pow(static_cast<double>(_binCenters[i]), 3.0 * static_cast<double>(_breakageKernelGamma) - 3.0) - pow(static_cast<double>(_bins[i]), 3.0 * static_cast<double>(_breakageKernelGamma) - 3.0)) / pow(static_cast<double>(_binCenters[j]), 3.0 * static_cast<double>(_breakageKernelGamma) - 3.0);
						}
						jac[-i + j] += factor * selection_function * b_integral * static_cast<double>(_frag->Upsilon_source[j]) * static_cast<double>(_binSizes[j]) / static_cast<double>(_binSizes[i]);
					}

					// add breakage sink to jac
					selection_function = static_cast<double>(_breakageRateConstant) * pow(static_cast<double>(_binCenters[i]), 3.0 * static_cast<double>(_breakageSelectionFunctionAlpha));
					jac[0] -= factor * selection_function * static_cast<double>(_frag->Upsilon_sink[i]);
				}
			}

			// go to the next row
			++jac;
		}
	}

	template <typename RowIteratorLiquid, typename RowIteratorSolid>
	void jacobianCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid, double factor, RowIteratorLiquid& jacLiquid, RowIteratorSolid& jacSolid, LinearBufferAllocator workSpace) const
	{
		jacobianLiquidImpl<RowIteratorLiquid>(t, secIdx, colPos, yLiquid, factor, jacLiquid, workSpace);
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

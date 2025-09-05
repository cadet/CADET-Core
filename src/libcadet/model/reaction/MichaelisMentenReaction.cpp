// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
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
#include "model/ExternalFunctionSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "linalg/ActiveDenseMatrix.hpp"
#include "MathUtil.hpp"
#include "Memory.hpp"

#include <algorithm>
#include <unordered_map>
#include <string>
#include <vector>

/*<codegen>
{
    "name": "MichaelisMentenParamHandler",
    "externalName": "ExtMichaelisMentenParamHandler",
    "parameters":
        [
            { "type": "ScalarReactionDependentParameter", "varName": "vMax", "confName": "MM_VMAX"},
            { "type": "ComponentDependentReactionDependentParameter", "varName": "kMM", "confName": "MM_KMM"},
            
            { "type": "ComponentDependentReactionDependentParameter", "varName": "kInhibitComp", "confName": "MM_KI_C"},
            { "type": "ComponentDependentReactionDependentParameter", "varName": "kInhibitUnComp", "confName": "MM_KI_UC"},
            { "type": "ComponentDependentReactionDependentParameter", "varName": "kInhibit", "confName": "MM_KI"}
        ]
}
</codegen>*/

/* Parameter description
 ------------------------
 MM_VMAX    - Maximum reaction rate (dim: number of reactions)
 MM_KMM     - Michaelis-Menten constant for each component (dim: number of reactions x number of components)
 MM_KI_C    - Competitive inhibition constants (\tilde{K}_{i,j,k}) ( dim: number of reactions x number of components x number of componetns)
 MM_KI_UC   - Uncompetitive inhibition constants (K_{i,j,k}) (dim: number of reactions x number of components x number of componetns)
 MM_KI      - Old interface for non competitive inhibition constants (K_{i,k}) (dim: number of reactions x number of components) 
                In the new interface non competative inhibition is initialized with \tilde{K}_{i,j,k} = K_{i,j,k}
*/


namespace cadet
{

namespace model
{

inline const char* MichaelisMentenParamHandler::identifier() CADET_NOEXCEPT { return "MICHAELIS_MENTEN"; }

inline bool MichaelisMentenParamHandler::validateConfig(unsigned int nReactions, unsigned int nComp, unsigned int const* nBoundStates)
{
	return true;
}

inline const char* ExtMichaelisMentenParamHandler::identifier() CADET_NOEXCEPT { return "EXT_MICHAELIS_MENTEN"; }

inline bool ExtMichaelisMentenParamHandler::validateConfig(unsigned int nReactions, unsigned int nComp, unsigned int const* nBoundStates)
{
	return true;
}

namespace
{
    /**
        * @brief Registers a matrix-valued parameter (row-major storage) with components as rows
        * @details The matrix-valued parameter has as many rows as there are components in the system.
        * @param [in,out] parameters Parameter map
        * @param [in] unitOpIdx Unit operation id
        * @param [in] parTypeIdx Particle type index
        * @param [in] paramName Name of the parameter
        * @param [in] mat Matrix to register
        */
    inline void registerCompRowMatrix(std::unordered_map<ParameterId, active*>& parameters, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx, const std::string& paramName, cadet::linalg::ActiveDenseMatrix& mat)
    {
        const cadet::StringHash hashName = cadet::hashStringRuntime(paramName);
        cadet::registerParam2DArray(parameters, mat.data(), mat.elements(), [=](bool multi, unsigned int row, unsigned int col)
            {
                return cadet::makeParamId(hashName, unitOpIdx, row, parTypeIdx, cadet::BoundStateIndep, col, cadet::SectionIndep);
            },
            mat.columns()
        );
    }
}

/**
    * @brief Defines a Multi-Substrate Michaelis-Menten reaction with complex inhibition
    * @details Implements the extended Michaelis-Menten kinetics with multiple types of inhibition:
    *          \f[ \begin{align}
    *              S \nu,
    *          \end{align} \f]
    *          where \f$ S \f$ is the stoichiometric matrix and the fluxes are given by
    *          \f[ \begin{align}
    *              \nu_i = v_{\max} \prod_{j=1}^{N} \left( \frac{c_j}{K_{m,j} (1 + \sum_{k=1}^{N} \tilde{K}_{i,j,k} c_k) + c_j \cdot (1 + \sum_{k=1}^{N} K_{i,j,k} c_k)} \right)
    *          \end{align} \f]
    *          Substrate components \f$ c_j \f$ are identified by negative entries in the stoichiometry matrix.
    *          This model supports:
    *          - Competitive inhibition via \f$ \tilde{K}_{i,j,k} \f$ (MM_KIC)
    *          - Uncompetitive inhibition via \f$ K_{i,j,k} \f$ (MM_KIUC)
    *          - Non-competitive inhibition through combination of both parameters
    * @tparam ParamHandler_t Type that can add support for external function dependence
    */
template <class ParamHandler_t>
class MichaelisMentenBase : public DynamicReactionModelBase
{
public:
    MichaelisMentenBase() { }
    virtual ~MichaelisMentenBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _paramHandler.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }
	virtual bool requiresWorkspace() const CADET_NOEXCEPT { return true; }
    virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
    {
        return _paramHandler.cacheSize(_stoichiometryBulk.columns(), nComp, totalNumBoundStates) +
            std::max(_stoichiometryBulk.columns() * sizeof(active),
                2 * (_nComp + totalNumBoundStates) * sizeof(double));
    }

    virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
    {
        DynamicReactionModelBase::configureModelDiscretization(paramProvider, nComp, nBound, boundOffset);

        if (paramProvider.exists("MM_STOICHIOMETRY_BULK"))
        {
            const std::size_t numElements = paramProvider.numElements("MM_STOICHIOMETRY_BULK");
            if (numElements % nComp != 0)
                throw InvalidParameterException("Size of field MM_STOICHIOMETRY_BULK must be a positive multiple of NCOMP (" + std::to_string(nComp) + ")");

            const unsigned int nReactions = numElements / nComp;

            _stoichiometryBulk.resize(nComp, nReactions);
            _idxSubstrate.resize(nReactions);
        }

        return true;
    }

    virtual unsigned int numReactionsLiquid() const CADET_NOEXCEPT { return _stoichiometryBulk.columns(); }
    virtual unsigned int numReactionsCombined() const CADET_NOEXCEPT { return 0; }

    CADET_DYNAMICREACTIONMODEL_BOILERPLATE

protected:
    ParamHandler_t _paramHandler; //!< Handles parameters and their dependence on external functions

    linalg::ActiveDenseMatrix _stoichiometryBulk;
    std::vector<std::vector<int>> _idxSubstrate; //!< Indices of substrate components for each reaction [reaction][substrate indices]
    std::vector<std::vector<std::unordered_set<int>>> _idxCompInhibitors; //!< Indices of competitive inhibitors [reaction][substrate][inhibitor indices]
    std::vector<std::vector<std::unordered_set<int>>> _idxUncompInhibitors; //!< Indices of uncompetitive inhibitors [reaction][substrate][inhibitor indices]
    bool _oldInterface;

    // Helper function to calculate the parameter index for inhibition
    inline unsigned int getInhibitionParamIndex(unsigned int reaction, int substrate, unsigned int inhibitor, bool oldInterface) const
    {
        // Calculate 3D parameter index: (reaction, substrate, inhibitor)
        if (!oldInterface)
            return (reaction * _nComp + static_cast<unsigned int>(substrate)) * _nComp + inhibitor;
        else
            return reaction * _nComp + inhibitor;
    }

    // Helper function to calculate the parameter index for km value
    inline unsigned int getKmmParamIndex(unsigned int reaction, int substrate, bool oldInterface) const
    {
        // Calculate 2D parameter index: (reaction, substrate)
        if (!oldInterface)
            return (reaction * _nComp + static_cast<unsigned int>(substrate));
        else
            return reaction;
    }

    virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
    {
        _paramHandler.configure(paramProvider, _stoichiometryBulk.columns(), _nComp, _nBoundStates);
        _paramHandler.registerParameters(_parameters, unitOpIdx, parTypeIdx, _nComp, _nBoundStates);
        _oldInterface = false;

        // handle old interface
        // - kmm is the number of reactions 
        // - kI refers to the non competative inhibition konstants
        if ((_stoichiometryBulk.columns() > 0) && (_paramHandler.kMM().size() == _stoichiometryBulk.columns()) || (paramProvider.exists("MM_KI")))
        {
            _oldInterface = true;
            LOG(Warning) << "MM_KI is only is only supported for backwards compatibility, but the implementation of Michaelis Menten kinetics has changed, please refer to the documentation for CADET version >= 5.2.0";
        }

        // Parameter validations
        if ((_stoichiometryBulk.columns() > 0) && ((_paramHandler.vMax().size() < _stoichiometryBulk.columns())))
            throw InvalidParameterException("MM_VMAX have to have the size (number of reactions)");

        if (_oldInterface && (_stoichiometryBulk.columns() > 0) && (_paramHandler.kMM().size() < _stoichiometryBulk.columns()))
            throw InvalidParameterException("MM_KMM have to have the size (number of reactions)");

        if (!_oldInterface && (_stoichiometryBulk.columns() > 0) && (_paramHandler.kMM().size() < _stoichiometryBulk.columns() * _nComp))
            throw InvalidParameterException("MM_KMM must have the size (number of reactions) x (number of components)");

        // kInhibitComp and kInhibitUnComp are 3D parameters with size (number of reactions) x (number of substrates) x (number of components)
        if (paramProvider.exists("MM_STOICHIOMETRY_BULK"))
        {
            const std::vector<double> s = paramProvider.getDoubleArray("MM_STOICHIOMETRY_BULK");
            std::vector<double> KIC(_stoichiometryBulk.columns() * _nComp * _nComp); 
            std::vector<double> KIUC(_stoichiometryBulk.columns() * _nComp * _nComp);
            bool hasCompetiveInhibition = false;
            
            if (paramProvider.exists("MM_KI_C"))
            {
                KIC = paramProvider.getDoubleArray("MM_KI_C");
                if (KIC.size() != _stoichiometryBulk.columns() * _nComp * _nComp)
                    throw InvalidParameterException("MM_KI_C must have the size (number of reactions) x (number of components) x (number of components)");
            }

            if (paramProvider.exists("MM_KI_UC"))
            {
                KIUC = paramProvider.getDoubleArray("MM_KI_UC");
                if (KIUC.size() != _stoichiometryBulk.columns() * _nComp * _nComp)
                    throw InvalidParameterException("MM_KI_UC must have the size (number of reactions) x (number of components) x (number of components)");
            }
            if (paramProvider.exists("MM_KI"))
            {
                KIUC = paramProvider.getDoubleArray("MM_KI");
                KIC = paramProvider.getDoubleArray("MM_KI");
            }


            if (s.size() != _stoichiometryBulk.elements())
                throw InvalidParameterException("MM_STOICHIOMETRY_BULK size mismatch: Expected " +
                    std::to_string(_stoichiometryBulk.elements()) + " elements but got " + std::to_string(s.size()));

            std::copy(s.begin(), s.end(), _stoichiometryBulk.data());

            // Identify substrates (negative entries in stoichiometry matrix)
            const unsigned int nReactions = static_cast<unsigned int>(_stoichiometryBulk.columns());
            _idxSubstrate.clear();
            _idxCompInhibitors.resize(nReactions);
            _idxUncompInhibitors.resize(nReactions);

            for (unsigned int r = 0; r < nReactions; ++r)
            {
                std::vector<int> idxSubstrateReaction_r;
                for (unsigned int c = 0; c < _nComp; ++c)
                {
                    if (_stoichiometryBulk.native(c, r) < 0.0)
                        idxSubstrateReaction_r.push_back(static_cast<int>(c));
                }

                if (idxSubstrateReaction_r.empty())
                    throw InvalidParameterException("Michaelis Menten: No substrates found for reaction " + std::to_string(r));
                if (idxSubstrateReaction_r.size() > 1 && _oldInterface)
                    throw InvalidParameterException("Michaelis Menten: The old interface is used which does not support multiple substrats. Please refer to the documentation of CADET version >= 5.2.0");

                _idxSubstrate.push_back(idxSubstrateReaction_r);

                // Pre-identify inhibitors for each substrate in this reaction
                const size_t numSubstrates = idxSubstrateReaction_r.size();


                _idxCompInhibitors[r].resize(numSubstrates);
                _idxUncompInhibitors[r].resize(numSubstrates);

                for (size_t subIdx = 0; subIdx < numSubstrates; ++subIdx)
                {
                    int j = idxSubstrateReaction_r[subIdx]; // Substrate index
                    std::unordered_set<int> compInhibitors;
                    std::unordered_set<int> uncompInhibitors;

                    for (unsigned int k = 0; k < _nComp; ++k)
                    {
                        const unsigned int paramIdx = getInhibitionParamIndex(r, j, k, _oldInterface);
                        double kIC = static_cast<double>(KIC[paramIdx]);
                        double kIUC = static_cast<double>(KIUC[paramIdx]);

                        if (kIC > 0.0)
                            compInhibitors.insert(static_cast<int>(k));

                        if (kIUC > 0.0)
                            uncompInhibitors.insert(static_cast<int>(k));
                    }

                    _idxCompInhibitors[r][subIdx] = std::move(compInhibitors);
                    _idxUncompInhibitors[r][subIdx] = std::move(uncompInhibitors);
                }

            }
        }

        registerCompRowMatrix(_parameters, unitOpIdx, parTypeIdx, "MM_STOICHIOMETRY_BULK", _stoichiometryBulk);
        return true;
    }

    template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
    int residualLiquidImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
        StateType const* y, ResidualType* res, const FactorType& factor, LinearBufferAllocator workSpace) const
    {
        typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

        // Calculate fluxes
        typedef typename DoubleActivePromoter<StateType, ParamType>::type flux_t;
        BufferedArray<flux_t> fluxes = workSpace.array<flux_t>(_stoichiometryBulk.columns());

        for (unsigned int r = 0; r < static_cast<unsigned int>(_stoichiometryBulk.columns()); ++r)
        {
            unsigned int nSub = _idxSubstrate[r].size();
            flux_t vProd = 1.0;

            // Product over all substrates
            for (unsigned int subIdx = 0; subIdx < nSub; ++subIdx)
            {
                int j = _idxSubstrate[r][subIdx]; // Substrate index j

                // Calculate competitive inhibition sum (K_tilde) using pre-identified inhibitors
                flux_t compInhSum = 0.0;
                for (int k : _idxCompInhibitors[r][subIdx])
                {
                    const unsigned int paramIdx = getInhibitionParamIndex(r, j, static_cast<unsigned int>(k), _oldInterface);


                    flux_t kIC = 0.0;
                    if (!_oldInterface)
                        kIC = static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->kInhibitComp[paramIdx]);
                    else
                        kIC = static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->kInhibit[paramIdx]);


                    compInhSum += y[k] / kIC;
                }

                // Uncompetitive inhibition using pre-identified inhibitors
                flux_t uncompInhSum = 0.0;
                for (int k : _idxUncompInhibitors[r][subIdx])
                {
                    const unsigned int paramIdx = getInhibitionParamIndex(r, j, static_cast<unsigned int>(k), _oldInterface);

                    flux_t kIU = 0.0;
                    if (!_oldInterface)
                        kIU = static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->kInhibitUnComp[paramIdx]);
                    else
                        kIU = static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->kInhibit[paramIdx]);

                    uncompInhSum += y[k] / kIU;
                }

                // KMM for this substrate

                const unsigned int kmmIdx = getKmmParamIndex(r, static_cast<unsigned int>(j), _oldInterface);
                const flux_t kMM_j = static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->kMM[kmmIdx]);



                // Calculate substrate contribution to production rate
                const flux_t numerator = y[j];
                const flux_t denominator = (kMM_j * (1.0 + compInhSum)) + y[j] * (1.0 + uncompInhSum);

                vProd *= numerator / denominator;
            }

            // Multiplication with vMax
            const flux_t vMax = static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->vMax[r]);
            fluxes[r] = vMax * vProd;
        }

        // Add reaction terms to residual
        _stoichiometryBulk.multiplyVector(static_cast<flux_t*>(fluxes), factor, res);

        return 0;
    }

    template <typename StateType, typename ResidualType, typename ParamType>
    int residualCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
        StateType const* yLiquid, StateType const* ySolid, ResidualType* resLiquid, ResidualType* resSolid, double factor, LinearBufferAllocator workSpace) const
    {
        std::fill_n(resLiquid, _nComp, 0.0);

        if (_nTotalBoundStates == 0)
            return 0;

        std::fill_n(resSolid, _nTotalBoundStates, 0.0);

        return 0;
    }

    template <typename RowIterator>
    void jacobianLiquidImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double factor, const RowIterator& jac, LinearBufferAllocator workSpace) const
    {
        typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

        for (unsigned int r = 0; r < static_cast<unsigned int>(_stoichiometryBulk.columns()); ++r)
        {
            unsigned int nSub = _idxSubstrate[r].size();

            // Calculate basic values for the flux
            std::vector<double> substratValues(nSub);    // Values for each substrate term
            std::vector<double> compInhSums(nSub, 0.0);  // Competitive inhibition sums
            std::vector<double> uncompInhSums(nSub, 0.0); // Uncompetitive inhibition sums
            double flux = 1.0; // Total flux

            // Calculate inhibition sums and substrate values
            for (unsigned int subIdx = 0; subIdx < nSub; ++subIdx)
            {
                int j = _idxSubstrate[r][subIdx]; // Substrate index j

                // Competitive inhibition - using pre-identified inhibitors
                double compInhSum = 0.0;
                for (int k : _idxCompInhibitors[r][subIdx])
                {
                    const unsigned int paramIdx = getInhibitionParamIndex(r, j, static_cast<unsigned int>(k), _oldInterface);

                    double kIC = 0.0;
                    if (!_oldInterface)
                        kIC = static_cast<double>(p->kInhibitComp[paramIdx]);
                    else
                        kIC = static_cast<double>(p->kInhibit[paramIdx]);

                    compInhSum += y[k] / kIC;
                }
                compInhSums[subIdx] = compInhSum;

                // Uncompetitive inhibition - using pre-identified inhibitors
                double uncompInhSum = 0.0;
                for (int k : _idxUncompInhibitors[r][subIdx])
                {
                    const unsigned int paramIdx = getInhibitionParamIndex(r, j, static_cast<unsigned int>(k), _oldInterface);

                    double kIU = 0.0;
                    if (!_oldInterface)
                        kIU = static_cast<double>(p->kInhibitUnComp[paramIdx]);
                    else
                        kIU = static_cast<double>(p->kInhibit[paramIdx]);

                    uncompInhSum += y[k] / kIU;
                }
                uncompInhSums[subIdx] = uncompInhSum;

                // KMM for this substrate
                const unsigned int kmmIdx = getKmmParamIndex(r, static_cast<unsigned int>(j), _oldInterface);
                const double kMM_j = static_cast<double>(p->kMM[kmmIdx]);

                // Calculate substrate term
                const double numerator = y[j];
                const double denominator = (kMM_j * (1.0 + compInhSum)) + y[j] * (1.0 + uncompInhSum);
                const double subValue = numerator / denominator;

                substratValues[subIdx] = subValue;
                flux *= subValue;
            }

            // Multiplication with vMax
            const double vMax = static_cast<double>(p->vMax[r]);
            flux *= vMax;

            // Calculate Jacobian derivatives for all components
            for (unsigned int comp = 0; comp < _nComp; ++comp)
            {
                double dvdy = 0.0;

                // Case 1: Component is a substrate
                bool isSubstrate = false;
                int substrateIdx = -1;
                for (unsigned int subIdx = 0; subIdx < nSub; ++subIdx)
                {
                    if (static_cast<int>(comp) == _idxSubstrate[r][subIdx])
                    {
                        isSubstrate = true;
                        substrateIdx = static_cast<int>(subIdx);
                        break;
                    }
                }

                if (isSubstrate)
                {
                    int j = static_cast<int>(comp); // Component is substrate j

                    // KMM for this substrate
                    const unsigned int kmmIdx = getKmmParamIndex(r, comp, _oldInterface);
                    const double kMM_j = static_cast<double>(p->kMM[kmmIdx]);

                    // Inhibition sums
                    const double compInhSum = compInhSums[substrateIdx];
                    const double uncompInhSum = uncompInhSums[substrateIdx];

                    // Denominator of substrate rate
                    const double denominator = (kMM_j * (1.0 + compInhSum)) + y[j] * (1.0 + uncompInhSum);

                    // Derivative dv/dj for a substrate j
                    const double dsubs_dj = (kMM_j * (1.0 + compInhSum)) / (denominator * denominator);

                    // Calculate total derivative for the substrate
                    const double factor_j = flux / substratValues[substrateIdx];
                    dvdy += factor_j * dsubs_dj;
                }

                // Case 2: Check for each substrate if this component is an inhibitor
                for (unsigned int subIdx = 0; subIdx < nSub; ++subIdx)
                {
                    int j = _idxSubstrate[r][subIdx]; // Substrate index j

                    // Check if comp is a competitive inhibitor for substrate j
                    if (_idxCompInhibitors[r][subIdx].count(static_cast<int>(comp)) > 0)
                    {
                        const unsigned int paramIdx = getInhibitionParamIndex(r, j, comp, _oldInterface);

                        double kIC = 0.0;
                        if (!_oldInterface)
                            kIC = static_cast<double>(p->kInhibitComp[paramIdx]);
                        else
                            kIC = static_cast<double>(p->kInhibit[paramIdx]);

                        const double kICinv = 1.0 / kIC;

                        // KMM for this substrate
                        const unsigned int kmmIdx = getKmmParamIndex(r, static_cast<unsigned int>(j), _oldInterface);
                        const double kMM_j = static_cast<double>(p->kMM[kmmIdx]);

                        // Inhibition sums
                        const double compInhSum = compInhSums[subIdx];
                        const double uncompInhSum = uncompInhSums[subIdx];

                        // Denominator of substrate rate
                        const double denominator = (kMM_j * (1.0 + compInhSum)) + y[j] * (1.0 + uncompInhSum);

                        // Derivative for competitive inhibition
                        const double dsubs_di_comp = -y[j] * kMM_j * kICinv / (denominator * denominator);

                        // Calculate total derivative for the competitive inhibitor
                        const double factor_i = flux / substratValues[subIdx];
                        dvdy += factor_i * dsubs_di_comp;
                    }

                    // Check if comp is an uncompetitive inhibitor for substrate j
                    if (_idxUncompInhibitors[r][subIdx].count(static_cast<int>(comp)) > 0)
                    {
                        const unsigned int paramIdx = getInhibitionParamIndex(r, j, comp, _oldInterface);


                        double kIU = 0.0;
                        if (!_oldInterface)
                            kIU = static_cast<double>(p->kInhibitUnComp[paramIdx]);
                        else
                            kIU = static_cast<double>(p->kInhibit[paramIdx]);

                        const double kIUinv = 1.0 / kIU;

                        // KMM for this substrate
                        const unsigned int kmmIdx = getKmmParamIndex(r, static_cast<unsigned int>(j), _oldInterface);
                        const double kMM_j = static_cast<double>(p->kMM[kmmIdx]);

                        // Inhibition sums
                        const double compInhSum = compInhSums[subIdx];
                        const double uncompInhSum = uncompInhSums[subIdx];

                        // Denominator of substrate rate
                        const double denominator = (kMM_j * (1.0 + compInhSum)) + y[j] * (1.0 + uncompInhSum);

                        // Derivative for uncompetitive inhibition
                        const double dsubs_di_uncomp = -y[j] * y[j] * kIUinv / (denominator * denominator);

                        // Calculate total derivative for the uncompetitive inhibitor
                        const double factor_i = flux / substratValues[subIdx];
                        dvdy += factor_i * dsubs_di_uncomp;
                    }
                }

                if (std::abs(dvdy) < 1e-12)
                    continue;

                // Gradients to the Jacobian matrix
                RowIterator curJac = jac;
                for (unsigned int row = 0; row < _nComp; ++row, ++curJac)
                {
                    const double colFactor = static_cast<double>(_stoichiometryBulk.native(row, r)) * factor;
                    curJac[static_cast<int>(comp) - static_cast<int>(row)] += colFactor * dvdy;
                }
            }
        }
    }

    template <typename RowIteratorLiquid, typename RowIteratorSolid>
    void jacobianCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid, double factor, const RowIteratorLiquid& jacLiquid, const RowIteratorSolid& jacSolid, LinearBufferAllocator workSpace) const
    {
        // No combined phase reaction implemented
    }
};

typedef MichaelisMentenBase<MichaelisMentenParamHandler> MichaelisMentenReaction;
typedef MichaelisMentenBase<ExtMichaelisMentenParamHandler> ExternalMichaelisMentenReaction;

namespace reaction
{
    void registerMichaelisMentenReaction(std::unordered_map<std::string, std::function<model::IDynamicReactionModel* ()>>& reactions)
    {
        reactions[MichaelisMentenReaction::identifier()] = []() { return new MichaelisMentenReaction(); };
        reactions[ExternalMichaelisMentenReaction::identifier()] = []() { return new ExternalMichaelisMentenReaction(); };
    }
}  // namespace reaction

}  // namespace model

}  // namespace cadet

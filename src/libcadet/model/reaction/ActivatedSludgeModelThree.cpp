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
	"name": "ActivatedSludgeModelThreeParamHandler",
	"externalName": "ExtActivatedSludgeModelThreeParamHandler",
	"parameters":
		[
			{ "type": "ScalarParameter", "varName": "kh20", "confName": "ASM3_KH20"},
			{ "type": "ScalarParameter", "varName": "T", "confName": "ASM3_T"},
			{ "type": "ScalarParameter", "varName": "io2", "confName": "ASM3_IO2"},
			{ "type": "ScalarParameter", "varName": "V", "confName": "ASM3_V"},
			{ "type": "ScalarParameter", "varName": "k_sto20", "confName": "ASM3_KSTO20"},
			{ "type": "ScalarParameter", "varName": "kx", "confName": "ASM3_KX"},
			{ "type": "ScalarParameter", "varName": "kho2", "confName": "ASM3_KHO2"},
			{ "type": "ScalarParameter", "varName": "khss", "confName": "ASM3_KHSS"},
			{ "type": "ScalarParameter", "varName": "khn03", "confName": "ASM3_KHNO3"},
			{ "type": "ScalarParameter", "varName": "etahno3", "confName": "ASM3_ETA_HNO3"},
			{ "type": "ScalarParameter", "varName": "khnh4", "confName": "ASM3_KHNH4"},
			{ "type": "ScalarParameter", "varName": "khalk", "confName": "ASM3_KHALK"},
			{ "type": "ScalarParameter", "varName": "khsto", "confName": "ASM3_KHSTO"},
			{ "type": "ScalarParameter", "varName": "muh2o", "confName": "ASM3_MU_H20"},
			{ "type": "ScalarParameter", "varName": "etahend", "confName": "ASM3_ETAH_END"},
			{ "type": "ScalarParameter", "varName": "bh20", "confName": "ASM3_BH20"},
			{ "type": "ScalarParameter", "varName": "muAUT20", "confName": "ASM3_MU_AUT20"},
			{ "type": "ScalarParameter", "varName": "kno2", "confName": "ASM3_KNO2"},
			{ "type": "ScalarParameter", "varName": "knnh4", "confName": "ASM3_KNNH4"},
			{ "type": "ScalarParameter", "varName": "knalk", "confName": "ASM3_KNALK"},
			{ "type": "ScalarParameter", "varName": "baut20", "confName": "ASM3_BAUT20"},
			{ "type": "ScalarParameter", "varName": "etanend", "confName": "ASM3_ETAN_END"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------


*/


namespace cadet
{

namespace model
{

inline const char* ActivatedSludgeModelThreeParamHandler::identifier() CADET_NOEXCEPT { return "ACTIVATED_SLUDGE_MODEL3"; }

inline bool ActivatedSludgeModelThreeParamHandler::validateConfig(unsigned int nReactions, unsigned int nComp, unsigned int const* nBoundStates)
{
	return true;
}

inline const char* ExtActivatedSludgeModelThreeParamHandler::identifier() CADET_NOEXCEPT { return "EXT_ACTIVATED_SLUDGE_MODEL3"; }

inline bool ExtActivatedSludgeModelThreeParamHandler::validateConfig(unsigned int nReactions, unsigned int nComp, unsigned int const* nBoundStates)
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
 * @brief Defines a Michaelis-Menten reaction kinetic with simple inhibition
 * @details Implements the Michaelis-Menten kinetics: \f[ \begin{align}
 *              S \nu,
 *          \end{align} \f]
 *          where \f$ S \f$ is the stoichiometric matrix and the fluxes are given by
 *          \f[ \begin{align}
 *              \nu_i = \frac{\mu_{\mathrm{max},i} c_S}{k_{\mathrm{MM},i} + c_S}.
 *          \end{align} \f]
 *          The substrate component \f$ c_S \f$ is identified by the index of the
 *          first negative entry in the stoichiometry of this reaction.
 *			In addition, the reaction might be inhibited by other components. In this
 *          case, the flux has the form
 *          \f[ \begin{align}
 *              \nu_i = \frac{\mu_{\mathrm{max},i} c_S}{k_{\mathrm{MM},i} + c_S} \cdot \frac{1}{1 + \sum_i \frac{1+ k_{\mathrm{I},i,j}}{k_{\mathrm{I},i,j} + c_{\mathrm{I},j}}}.
 *          \end{align} \f]
 *          The value of \f$ k_{\mathrm{I},i,j} \f$ decides whether component \f$ j \f$
 *          inhibits reaction \f$ i \f$. If \f$ k_{\mathrm{I},i,j} \leq 0 \f$, the component
 *          does not inhibit the reaction.
 *          Only reactions in liquid phase are supported (no solid phase or cross-phase reactions).
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class ActivatedSludgeModelThreeBase : public DynamicReactionModelBase
{
public:

	ActivatedSludgeModelThreeBase() CADET_NOEXCEPT { }
	virtual ~ActivatedSludgeModelThreeBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _paramHandler.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }
	virtual bool requiresWorkspace() const CADET_NOEXCEPT { return true; }
	virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
	{
		return _paramHandler.cacheSize(_stoichiometry.columns(), nComp, totalNumBoundStates) + std::max(_stoichiometry.columns() * sizeof(active), 2 * (_nComp + totalNumBoundStates) * sizeof(double));
	}

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		DynamicReactionModelBase::configureModelDiscretization(paramProvider, nComp, nBound, boundOffset);
		_stoichiometry.resize(nComp, 13);
		return true;
	}

	virtual unsigned int numReactionsLiquid() const CADET_NOEXCEPT { return _stoichiometry.columns(); }
	virtual unsigned int numReactionsCombined() const CADET_NOEXCEPT { return 0; }

	CADET_DYNAMICREACTIONMODEL_BOILERPLATE

protected:
	ParamHandler_t _paramHandler; //!< Handles parameters and their dependence on external functions

	linalg::ActiveDenseMatrix _stoichiometry;

	unsigned int _idxSO = 0; //!< SO component index, default 0
	unsigned int _idxSS_nad = 1; //!< SS / SS_ad component index, default 1
	unsigned int _idxSNH = 2; //!< SNH component index, default 2
	unsigned int _idxSNO = 3; //!< SNO component index, default 3
	unsigned int _idxSN2 = 4; //!< SN2 component index, default 4
	unsigned int _idxSALK = 5; //!< SALK component index, default 5
	unsigned int _idxSI_nad = 6; //!< SI component index, default 6
	unsigned int _idxXI = 7; //!< XI component index, default 7
	unsigned int _idxXS = 8; //!< XS component index, default 8
	unsigned int _idxXH = 9; //!< XH component index, default 9
	unsigned int _idxXSTO = 10; //!< XSTO component index, default 10
	unsigned int _idxXA = 11; //!< XA component index, default 11
	unsigned int _idxXMI = 12; //!< XMI component index, default 12
	unsigned int _idxSI_ad = 13;
	unsigned int _idxSS_ad = 14;

	bool  _fractionate = false;

	

	virtual bool configureStoich(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		return true;
	}


	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		_paramHandler.configure(paramProvider, _stoichiometry.columns(), _nComp, _nBoundStates);
		_paramHandler.registerParameters(_parameters, unitOpIdx, parTypeIdx, _nComp, _nBoundStates);


		_stoichiometry.resize(_nComp, 13);
		_stoichiometry.setAll(0);

		if (paramProvider.exists("ASM3_FRACTIONATE"))
			_fractionate = paramProvider.getBool("ASM3_FRACTIONATE");

		if (paramProvider.exists("ASM3_COMP_IDX")) {
			const std::vector<uint64_t> compIdx = paramProvider.getUint64Array("ASM3_COMP_IDX");
			if (compIdx.size() != 13) {
				throw InvalidParameterException("ASM3 configuration: ASM3_COMP_IDX must have 13 elements");
			}
			else {
				LOG(Debug) << "ASM3_COMP_IDX set: " << compIdx;
				_idxSO = compIdx[0];
				_idxSS_nad = compIdx[1];
				_idxSNH = compIdx[2];
				_idxSNO = compIdx[3];
				_idxSN2 = compIdx[4];
				_idxSALK = compIdx[5];
				_idxSI_nad = compIdx[6];
				_idxXI = compIdx[7];
				_idxXS = compIdx[8];
				_idxXH = compIdx[9];
				_idxXSTO = compIdx[10];
				_idxXA = compIdx[11];
				_idxXMI = compIdx[12];
				if (_fractionate)
				{
					//todo warníng 
					_idxSS_ad = compIdx[13];
					_idxSI_ad = compIdx[14];
				}
			}
		}
		else {
			LOG(Debug) << "ASM3_COMP_IDX not set, using defaults";
		}

		// parameter set ASM3h
		const double iNSI = paramProvider.getDouble("ASM3_INSI");
		const double iNSS = paramProvider.getDouble("ASM3_INSS");
		const double iNXI = paramProvider.getDouble("ASM3_INXI");
		const double iNXS = paramProvider.getDouble("ASM3_INXS");
		const double iNBM = paramProvider.getDouble("ASM3_INBM");

		const double fSI = paramProvider.getDouble("ASM3_FSI");
		const double fXI = paramProvider.getDouble("ASM3_FXI");

		const double YH_aer = paramProvider.getDouble("ASM3_YH_AER");
		if (YH_aer < 1e-16)
			throw InvalidParameterException("ASM3 configuration: YH_aer must be bigger than zero");
		const double YH_anox = paramProvider.getDouble("ASM3_YH_ANOX");
		const double YSTO_aer = paramProvider.getDouble("ASM3_YSTO_AER");
		if (YSTO_aer < 1e-16)
			throw InvalidParameterException("ASM3 configuration: YSTO_aer must be bigger than zero");
		const double YSTO_anox = paramProvider.getDouble("ASM3_YSTO_ANOX");
		const double YA = paramProvider.getDouble("ASM3_YA");
		if (YA < 1e-16)
			throw InvalidParameterException("ASM3 configuration: YA must be bigger than zero");

		const double fiSS_BM_prod = paramProvider.getDouble("ASM3_FISS_BM_PROD");
		const double iVSS_BM = paramProvider.getDouble("ASM3_IVSS_BM");
		const double iTSS_VSS_BM = paramProvider.getDouble("ASM3_ITSS_VSS_BM");

		// internal variables
		const double fXMI_BM = fiSS_BM_prod * fXI * iVSS_BM * (iTSS_VSS_BM - 1);
		const double c1n = iNXS - iNSI * fSI - (1 - fSI) * iNSS;
		const double c2n = iNSS;
		const double c3n = iNSS;
		const double c4n = -iNBM;
		const double c5n = c4n;
		const double c6n = -fXI * iNXI + iNBM;
		const double c7n = c6n;
		const double c10n = -1 / YA - iNBM;
		const double c11n = -fXI * iNXI + iNBM;
		const double c12n = c11n;
		const double c3no = (YSTO_anox - 1.0) / (40.0 / 14.0);
		const double c5no = (1.0 - 1 / YH_anox) / (40.0 / 14.0);
		const double c7no = (fXI - 1) / (40.0 / 14.0);
		const double c9no = -14.0 / 40.0;
		const double c10no = 1 / YA;
		const double c12no = c7no;
		const double c1a = c1n / 14.0;
		const double c2a = c2n / 14.0;
		const double c3a = (c3n - c3no) / 14.0;
		const double c4a = c4n / 14.0;
		const double c5a = (c5n - c5no) / 14.0;
		const double c6a = c6n / 14.0;
		const double c7a = (c7n - c7no) / 14.0;
		const double c9a = 1 / 40.0;
		const double c10a = (c10n - c10no) / 14.0;
		const double c11a = c11n / 14.0;
		const double c12a = (c12n - c12no) / 14.0;

		// SO
		_stoichiometry.native(_idxSO, 1) = YSTO_aer - 1;
		_stoichiometry.native(_idxSO, 3) = 1 - 1 / YH_aer;
		_stoichiometry.native(_idxSO, 5) = -1 * (1 - fXI);
		_stoichiometry.native(_idxSO, 7) = -1;
		_stoichiometry.native(_idxSO, 9) = -(64.0 / 14.0) * 1 / YA + 1;
		_stoichiometry.native(_idxSO, 10) = -1 * (1 - fXI);
		_stoichiometry.native(_idxSO, 12) = 1;

		// SS_nad
		_stoichiometry.native(_idxSS_nad, 0) = (1 - fSI);
		_stoichiometry.native(_idxSS_nad, 1) = -1;
		_stoichiometry.native(_idxSS_nad, 2) = -1;

		// SS_ad
		if (_fractionate)
		{
			_stoichiometry.native(_idxSS_ad, 0) = (1 - fSI);
			_stoichiometry.native(_idxSS_ad, 1) = -1;
			_stoichiometry.native(_idxSS_ad, 2) = -1;
		}
		// SNH
		_stoichiometry.native(_idxSNH, 0) = c1n;
		_stoichiometry.native(_idxSNH, 1) = c2n;
		_stoichiometry.native(_idxSNH, 2) = c3n;
		_stoichiometry.native(_idxSNH, 3) = c4n;
		_stoichiometry.native(_idxSNH, 4) = c5n;
		_stoichiometry.native(_idxSNH, 5) = c6n;
		_stoichiometry.native(_idxSNH, 6) = c7n;
		_stoichiometry.native(_idxSNH, 9) = c10n;
		_stoichiometry.native(_idxSNH, 10) = c11n;
		_stoichiometry.native(_idxSNH, 11) = c12n;

		// SNO
		_stoichiometry.native(_idxSNO, 2) = c3no;
		_stoichiometry.native(_idxSNO, 4) = c5no;
		_stoichiometry.native(_idxSNO, 6) = c7no;
		_stoichiometry.native(_idxSNO, 8) = c9no;
		_stoichiometry.native(_idxSNO, 9) = c10no;
		_stoichiometry.native(_idxSNO, 11) = c12no;

		// SN2
		_stoichiometry.native(_idxSN2, 2) = -c3no;
		_stoichiometry.native(_idxSN2, 4) = -c5no;
		_stoichiometry.native(_idxSN2, 6) = -c7no;
		_stoichiometry.native(_idxSN2, 8) = -c9no;
		_stoichiometry.native(_idxSN2, 11) = -c12no;

		// SALK
		_stoichiometry.native(_idxSALK, 0) = c1a;
		_stoichiometry.native(_idxSALK, 1) = c2a;
		_stoichiometry.native(_idxSALK, 2) = c3a;
		_stoichiometry.native(_idxSALK, 3) = c4a;
		_stoichiometry.native(_idxSALK, 4) = c5a;
		_stoichiometry.native(_idxSALK, 5) = c6a;
		_stoichiometry.native(_idxSALK, 6) = c7a;
		_stoichiometry.native(_idxSALK, 8) = c9a;
		_stoichiometry.native(_idxSALK, 9) = c10a;
		_stoichiometry.native(_idxSALK, 10) = c11a;
		_stoichiometry.native(_idxSALK, 11) = c12a;

		// SI_ad
		_stoichiometry.native(_idxSI_nad, 0) = fSI;

		//SI_nad
		if(_fractionate)
			_stoichiometry.native(_idxSI_ad, 0) = fSI;

		// XI
		_stoichiometry.native(_idxXI, 5) = fXI;
		_stoichiometry.native(_idxXI, 6) = fXI;
		_stoichiometry.native(_idxXI, 10) = fXI;
		_stoichiometry.native(_idxXI, 11) = fXI;

		// XS
		_stoichiometry.native(_idxXS, 0) = -1;

		// XH
		_stoichiometry.native(_idxXH, 3) = 1;
		_stoichiometry.native(_idxXH, 4) = 1;
		_stoichiometry.native(_idxXH, 5) = -1;
		_stoichiometry.native(_idxXH, 6) = -1;

		// XSTO
		_stoichiometry.native(_idxXSTO, 1) = YSTO_aer;
		_stoichiometry.native(_idxXSTO, 2) = YSTO_anox;
		_stoichiometry.native(_idxXSTO, 3) = -1 / YH_aer;
		_stoichiometry.native(_idxXSTO, 4) = -1 / YH_anox;
		_stoichiometry.native(_idxXSTO, 7) = -1;
		_stoichiometry.native(_idxXSTO, 8) = -1;

		// XA
		_stoichiometry.native(_idxXA, 9) = 1;
		_stoichiometry.native(_idxXA, 10) = -1;
		_stoichiometry.native(_idxXA, 11) = -1;

		// XMI
		_stoichiometry.native(_idxXMI, 5) = fXMI_BM;
		_stoichiometry.native(_idxXMI, 6) = fXMI_BM;
		_stoichiometry.native(_idxXMI, 10) = fXMI_BM;
		_stoichiometry.native(_idxXMI, 11) = fXMI_BM;


		//registerCompRowMatrix(_parameters, unitOpIdx, parTypeIdx, "MM_STOICHIOMETRY_BULK", _stoichiometryBulk); todo

		return true;
	}

		template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	int residualLiquidImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* y, ResidualType* res, const FactorType& factor, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Calculate fluxes
		typedef typename DoubleActivePromoter<StateType, ParamType>::type flux_t;
		BufferedArray<flux_t> fluxes = workSpace.array<flux_t>(_stoichiometry.columns());
		
		const flux_t kh20		= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->kh20);
		const flux_t T		= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->T);
		const flux_t io2		= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->io2);
		const flux_t V		= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->V);
		const flux_t k_sto20	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->k_sto20);
		const flux_t kx		= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->kx);
		const flux_t kho2		= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->kho2);
		const flux_t khss		= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->khss);
		const flux_t khn03	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->khn03);
		const flux_t etahno3	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->etahno3);
		const flux_t khnh4	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->khnh4);
		const flux_t khalk	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->khalk);
		const flux_t khsto	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->khsto);
		const flux_t muh2o	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->muh2o);
		const flux_t etahend	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->etahend);
		const flux_t bh20		= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->bh20);
		const flux_t muAUT20	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->muAUT20);
		const flux_t kno2		= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->kno2);
		const flux_t knnh4	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->knnh4);
		const flux_t knalk	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->knalk);
		const flux_t baut20	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->baut20);
		const flux_t etanend	= static_cast<typename DoubleActiveDemoter<flux_t, ParamType>::type>(p->etanend);

		// derived parameters
		const double ft04 = exp(-0.04 * (20.0 - static_cast<double>(T)));
		const double ft07 = exp(-0.06952 * (20.0 - static_cast<double>(T)));
		const double ft105 = exp(-0.105 * (20.0 - static_cast<double>(T)));
		const double k_sto = static_cast<double>(k_sto20) * ft07;
		const double muH = static_cast<double>(muh2o) * ft07;
		const double bH = static_cast<double>(bh20) * ft07;
		const double muAUT = static_cast<double>(muAUT20) * ft105;
		const double bAUT = static_cast<double>(baut20) * ft105;

		StateType SO = y[_idxSO];
		StateType SS_nad = y[_idxSS_nad];
		StateType SS_ad = 0.0;
		if (_fractionate)
			SS_ad = y[_idxSS_ad];
		StateType SNH = y[_idxSNH];
		StateType SNO = y[_idxSNO];
		StateType SN2 = y[_idxSN2];
		StateType SALK = y[_idxSALK];
		// StateType SI = y[_idxSI]; // unused
		// StateType XI = y[_idxXI]; // unused
		StateType XS = y[_idxXS];
		StateType XH = y[_idxXH];
		StateType XSTO = y[_idxXSTO];
		StateType XA = y[_idxXA];
		// StateType XMI = y[_idxXMI]; // unused

		
		// p1: Hydrolysis of organic structures
		fluxes[0] = kh20 * ft04 * (XS/XH) / ((XS/XH) + kx) * XH;
		if (XH < 0.1)
			fluxes[0] = kh20 * ft04 * (XS/0.1) / ((XS/0.1) + kx) * XH;

		
		// p2: Aerobic storage of SS
		fluxes[1] = k_sto * SO / (SO + kho2) * ( SS_ad + SS_nad ) / ( ( SS_ad + SS_nad ) + khss )  * XH;

		// p3: Anoxic storage of SS
		fluxes[2] = k_sto * etahno3 * kho2 / (SO + kho2) * ( SS_ad + SS_nad ) / ( ( SS_ad + SS_nad ) + khss )  * SNO / (SNO + khn03) * XH;

		// p4: Aerobic growth of heterotrophic biomass (XH)
		fluxes[3] = muH * SO / (SO + kho2) * SNH / (SNH + khnh4) * SALK / (SALK + khalk) * ((XSTO/XH)) / (((XSTO/XH)) + khsto) * XH;
		if (XH < 0.1)
			fluxes[3] = muH * SO / (SO + kho2) * SNH / (SNH + khnh4) * SALK / (SALK + khalk) * (XSTO / 0.1) / ((XSTO / 0.1) + khsto) * XH;

		// p5: Anoxic growth of heterotrophic biomass (XH, denitrification)
		fluxes[4] = muH * etahno3 * kho2 / (kho2 + SO) * SNH / (khnh4 + SNH) * SALK / (khalk + SALK) * (XSTO / XH) * 1 / (khsto + (XSTO / XH)) * SNO / (khn03 + SNO) * XH;
		if (XH < 0.1)
			fluxes[4] = muH * etahno3 * kho2 / (kho2 + SO) * SNH / (khnh4 + SNH) * SALK / (khalk + SALK) * (XSTO / 0.1) * 1 / (khsto + (XSTO/0.1)) * SNO / (khn03 + SNO) * XH;

		// r6: Aerobic endogenous respiration of heterotroph microorganisms (XH)
		fluxes[5] = bH * SO / (SO + kho2) * XH;

		// r7: Anoxic endogenous respiration of heterotroph microorganisms (XH)
		fluxes[6] = bH * etahend * kho2 / (SO + kho2) * SNO / (SNO + khn03) * XH;

		// r8: Aerobic respiration of internal cell storage products
		fluxes[7] = bH * SO / (SO + kho2) * XSTO;

		// r9: Anoxic respiration of internal cell storage products 
		fluxes[8] = bH * etahend * kho2 / (SO + kho2) * SNO / (SNO + khn03) * XSTO;

		// r10: Aerobic growth of autotrophic biomass (XAUT, nitrification)
		fluxes[9] = muAUT * SO / (SO + kno2) * SNH / (SNH + knnh4) * SALK / (SALK + knalk) * XA;

		// r11: Aerobic endogenous respiration of autotrophic biomass (XAUT)
		fluxes[10] = bAUT * SO / (SO + kho2) * XA;

		// r12: Anoxic endogenous respiration of autotrophic biomass (XAUT)
		fluxes[11] = bAUT * etanend * SNO / (SNO + khn03) * kho2 / (SO + kho2) * XA;

		// r13: Aeration
		// TODO: is V in litres?
		fluxes[12] = io2 / V;
			
		// Add reaction terms to residual
		_stoichiometry.multiplyVector(static_cast<flux_t*>(fluxes), factor, res);

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

		//parameters
		const double kh20 = static_cast<double>(p->kh20);
		const double T = static_cast<double>(p->T);
		const double io2 = static_cast<double>(p->io2);
		const double V = static_cast<double>(p->V);
		const double k_sto20 = static_cast<double>(p->k_sto20);
		const double kx = static_cast<double>(p->kx);
		const double kho2 = static_cast<double>(p->kho2);
		const double khss = static_cast<double>(p->khss);
		const double khn03 = static_cast<double>(p->khn03);
		const double etahno3 = static_cast<double>(p->etahno3);
		const double khnh4 = static_cast<double>(p->khnh4);
		const double khalk = static_cast<double>(p->khalk);
		const double khsto = static_cast<double>(p->khsto);
		const double muh2o = static_cast<double>(p->muh2o);
		const double etahend = static_cast<double>(p->etahend);
		const double bh20 = static_cast<double>(p->bh20);
		const double muAUT20 = static_cast<double>(p->muAUT20);
		const double kno2 = static_cast<double>(p->kno2);
		const double knnh4 = static_cast<double>(p->knnh4);
		const double knalk = static_cast<double>(p->knalk);
		const double baut20 = static_cast<double>(p->baut20);
		const double etanend = static_cast<double>(p->etanend);

		// derived parameters
		const double ft04 = exp(-0.04 * (20.0 - static_cast<double>(T)));
		const double ft07 = exp(-0.06952 * (20.0 - static_cast<double>(T)));
		const double ft105 = exp(-0.105 * (20.0 - static_cast<double>(T)));
		const double k_sto = static_cast<double>(k_sto20) * ft07;
		const double muH = static_cast<double>(muh2o) * ft07;
		const double bH = static_cast<double>(bh20) * ft07;
		const double muAUT = static_cast<double>(muAUT20) * ft105;
		const double bAUT = static_cast<double>(baut20) * ft105;

		double SO = y[_idxSO];
		double SS_nad = y[_idxSS_nad];
		double SS_ad = 0.0;
		double SNH = y[_idxSNH];
		double SNO = y[_idxSNO];
		double SN2 = y[_idxSN2];
		double SALK = y[_idxSALK];
		double SI_nad = y[_idxSI_nad];
		double SI_ad = 0.0;
		double XS = y[_idxXS];
		double XH = y[_idxXH];
		double XSTO = y[_idxXSTO];
		double XA = y[_idxXA];

		if (_fractionate)
		{
			SS_ad = y[_idxSS_ad];
			SI_ad = y[_idxSI_ad];
		}


		double d[13][15] = {};


		// p1: Hydrolysis: kh20 * ft04 * XS/XH_S / (XS/XH_S + kx) * XH;
		d[0][_idxXS] = kh20 * ft04
			* XH / ((XS + XH * kx)
				* (XS + XH * kx)) * XH;
		d[0][_idxXH] = kh20 * ft04
			* (XS * XS) / ((XS + kx * XH)
				* (XS + kx * XH));
		if (XH < 0.1)
		{
			d[0][_idxXS] = kh20 * ft04
				* 0.1 / ((XS + 0.1 * kx) * (XS + 0.1 * kx)) * XH;
			d[0][_idxXH] = 0.0;
		}

		// p2: Aerobic storage of SS: k_sto * SO / (SO + kho2) * ( SS_ad + SS_nad ) / ( ( SS_ad + SS_nad ) + khss )  * XH;
		d[1][_idxSO] = k_sto
			* (SS_ad + SS_nad) / ((SS_ad + SS_nad) + khss)
			* kho2 / ((SO + kho2) * (SO + kho2)) * XH;
		if (_fractionate)
			d[1][_idxSS_ad] = k_sto
				* SO / (SO + kho2)
				* khss / ((SS_ad + SS_nad + khss) * (SS_ad + SS_nad + khss)) * XH;
		d[1][_idxSS_nad] = k_sto
			* SO / (SO + kho2)
			* khss / ((SS_ad + SS_nad + khss) * (SS_ad + SS_nad + khss)) * XH;
		d[1][_idxXH] = k_sto
			* SO / (SO + kho2)
			* (SS_ad + SS_nad) / ((SS_ad + SS_nad) + khss);

		// p3: Anoxic storage of SS: k_sto * etahno3 * kho2 / (SO + kho2) * ( SS_ad + SS_nad ) / ( ( SS_ad + SS_nad ) + khss )  * SNO / (SNO + khn03) * XH;
		d[2][_idxSO] = k_sto * etahno3
			* -kho2 / ((SO + kho2) * (SO + kho2))
			* (SS_ad + SS_nad) / ((SS_ad + SS_nad) + khss)
			* SNO / (SNO + khn03) * XH;
		if(_fractionate)
			d[2][_idxSS_ad] = k_sto * etahno3
				* kho2 / (SO + kho2)
				* khss / ((SS_ad + SS_nad + khss) * (SS_ad + SS_nad + khss))
				* SNO / (SNO + khn03) * XH;
		d[2][_idxSS_nad] = k_sto * etahno3
			* kho2 / (SO + kho2)
			* khss / ((SS_ad + SS_nad + khss) * (SS_ad + SS_nad + khss))
			* SNO / (SNO + khn03) * XH;
		d[2][_idxSNO] = k_sto * etahno3
			* kho2 / (SO + kho2)
			* (SS_ad + SS_nad) / ((SS_ad + SS_nad) + khss)
			* khn03 / ((SNO + khn03) * (SNO + khn03)) * XH;
		d[2][_idxXH] = k_sto * etahno3
			* kho2 / (SO + kho2)
			* (SS_ad + SS_nad) / ((SS_ad + SS_nad) + khss)
			* SNO / (SNO + khn03);

		// p4: Aerobic growth: muH * SO / (SO + kho2) * SNH / (SNH + khnh4) * SALK / (SALK + khalk) * (XSTO/XH_S) / ((XSTO/XH_S) + khsto) * XH;
		d[3][_idxSO] = muH
			* kho2 / ((SO + kho2) * (SO + kho2))
			* SNH / (SNH + khnh4)
			* SALK / (SALK + khalk)
			* (XSTO / XH) / ((XSTO / XH) + khsto) * XH;
		d[3][_idxSNH] = muH
			* SO / (SO + kho2)
			* khnh4 / ((SNH + khnh4) * (SNH + khnh4))
			* SALK / (SALK + khalk)
			* (XSTO / XH) / ((XSTO / XH) + khsto) * XH;
		d[3][_idxSALK] = muH
			* SO / (SO + kho2)
			* SNH / (SNH + khnh4)
			* khalk / ((SALK + khalk) * (SALK + khalk))
			* (XSTO / XH) / ((XSTO / XH) + khsto) * XH;
		d[3][_idxXSTO] = muH
			* SO / (SO + kho2)
			* SNH / (SNH + khnh4)
			* SALK / (SALK + khalk)
			* (khsto * XH) / ((XSTO + khsto * XH) * (XSTO + khsto * XH)) * XH;
		d[3][_idxXH] = muH
			* SO / (SO + kho2)
			* SNH / (SNH + khnh4)
			* SALK / (SALK + khalk)
			* (XSTO * XSTO) / ((XSTO + khsto * XH) * (XSTO + khsto * XH));

		if (XH < 0.1)
		{
			d[3][_idxSO] = muH
				* -kho2 / ((SO + kho2) * (SO + kho2))
				* SNH / (SNH + khnh4)
				* SALK / (SALK + khalk)
				* (XSTO / 0.1) / ((XSTO / 0.1) + khsto) * 0.1;
			d[3][_idxSNH] = muH
				* SO / (SO + kho2)
				* khnh4 / ((SNH + khnh4) * (SNH + khnh4))
				* khnh4 / ((SNH + khnh4) * (SNH + khnh4))
				* SALK / (SALK + khalk)
				* (XSTO / 0.1) / ((XSTO / 0.1) + khsto) * 0.1;
			d[3][_idxSALK] = muH
				* SO / (SO + kho2)
				* SNH / (SNH + khnh4)
				* khalk / ((SALK + khalk) * (SALK + khalk))
				* (XSTO / 0.1) / ((XSTO / 0.1) + khsto) * 0.1;
			d[3][_idxXSTO] = muH
				* SO / (SO + kho2)
				* SNH / (SNH + khnh4)
				* SALK / (SALK + khalk)
				* (khsto * 0.1) / ((XSTO + khsto * 0.1) * (XSTO + khsto * 0.1)) * 0.1;
			d[3][_idxXH] = 0.0;
		}

		// p5: Anoxic growth: muH * etahno3 * kho2 / (kho2 + SO) * SNH / (khnh4 + SNH) * SALK / (khalk + SALK)  * SNO / (khn03 + SNO)* (XSTO / XH_S) / (khsto + (XSTO/XH)_S) * XH;
		d[4][_idxSO] = muH * etahno3
			* -kho2 / ((kho2 + SO) * (kho2 + SO))
			* SNH / (khnh4 + SNH)
			* SALK / (khalk + SALK)
			* (XSTO / XH) / (khsto + (XSTO / XH))
			* SNO / (khn03 + SNO) * XH;
		d[4][_idxSNH] = muH * etahno3
			* kho2 / (kho2 + SO)
			* khnh4 / ((khnh4 + SNH) * (khnh4 + SNH))
			* SALK / (khalk + SALK)
			* (XSTO / XH) / (khsto + (XSTO / XH))
			* SNO / (khn03 + SNO) * XH;
		d[4][_idxSALK] = muH * etahno3
			* kho2 / (kho2 + SO)
			* SNH / (khnh4 + SNH)
			* khalk / ((khalk + SALK) * (khalk + SALK))
			* (XSTO / XH) / (khsto + (XSTO / XH))
			* SNO / (khn03 + SNO) * XH;
		d[4][_idxSNO] = muH * etahno3
			* kho2 / (kho2 + SO)
			* SNH / (khnh4 + SNH)
			* SALK / (khalk + SALK)
			* (XSTO / XH) / (khsto + (XSTO / XH))
			* khn03 / ((khn03 + SNO) * (khn03 + SNO)) * XH;
		d[4][_idxXSTO] = muH * etahno3
			* kho2 / (kho2 + SO)
			* SNH / (khnh4 + SNH)
			* SALK / (khalk + SALK)
			* (1.0 / XH)
			* SNO / (khn03 + SNO)
			* (khsto * XH * XH) / ((XSTO + khsto * XH) * (XSTO + khsto * XH)) * XH;
		d[4][_idxXH] = muH * etahno3
			* kho2 / (kho2 + SO)
			* SNH / (khnh4 + SNH)
			* SALK / (khalk + SALK)
			* SNO / (khn03 + SNO)
			* (XSTO * XSTO) / ((XSTO + khsto * XH) * (XSTO + khsto * XH));


		if (XH < 0.1)
		{
			d[4][_idxSO] = muH * etahno3
				* -kho2 / ((kho2 + SO) * (kho2 + SO))
				* SNH / (khnh4 + SNH)
				* SALK / (khalk + SALK)
				* (XSTO / 0.1) / (khsto + (XSTO / 0.1))
				* SNO / (khn03 + SNO) * 0.1;
			d[4][_idxSNH] = muH * etahno3
				* kho2 / (kho2 + SO)
				* khnh4 / ((khnh4 + SNH) * (khnh4 + SNH))
				* SALK / (khalk + SALK)
				* (XSTO / 0.1) / (khsto + (XSTO / 0.1))
				* SNO / (khn03 + SNO) * 0.1;
			d[4][_idxSALK] = muH * etahno3
				* kho2 / (kho2 + SO)
				* SNH / (khnh4 + SNH)
				* khalk / ((khalk + SALK) * (khalk + SALK))
				* (XSTO / 0.1) / (khsto + (XSTO / 0.1))
				* SNO / (khn03 + SNO) * 0.1;
			d[4][_idxSNO] = muH * etahno3
				* kho2 / (kho2 + SO)
				* SNH / (khnh4 + SNH)
				* SALK / (khalk + SALK)
				* (XSTO / 0.1) / (khsto + (XSTO / 0.1))
				* khn03 / ((khn03 + SNO) * (khn03 + SNO)) * 0.1;
			d[4][_idxXSTO] = muH * etahno3
				* kho2 / (kho2 + SO)
				* SNH / (khnh4 + SNH)
				* SALK / (khalk + SALK)
				* (1.0 / 0.1)
				* SNO / (khn03 + SNO)
				* (khsto * 0.1 * 0.1) / ((XSTO + khsto * 0.1) * (XSTO + khsto * 0.1)) * 0.1;
			d[4][_idxXH] = 0.0;
		}

		//reaction6: bH * SO / (SO + kho2) * XH;
		d[5][_idxSO] = bH
			* kho2 / ((SO + kho2) * (SO + kho2)) * XH;
		d[5][_idxXH] = bH * SO / (SO + kho2);

		//reaction7: bH * etahend * kho2 / (SO + kho2) * SNO / (SNO + khn03) * XH;
		d[6][_idxSO] = bH * etahend
			* -kho2 / ((SO + kho2) * (SO + kho2))
			* SNO / (SNO + khn03) * XH;
		d[6][_idxSNO] = bH * etahend
			* kho2 / (SO + kho2)
			* khn03 / ((SNO + khn03) * (SNO + khn03)) * XH;
		d[6][_idxXH] = bH * etahend
			* kho2 / (SO + kho2)
			* SNO / (SNO + khn03);

		//reaction8: bH * SO / (SO + kho2) * XSTO;
		d[7][_idxSO] = bH
			* kho2 / ((SO + kho2) * (SO + kho2)) * XSTO;
		d[7][_idxXSTO] = bH * SO / (SO + kho2);

		//reaction9: bH * etahend * kho2 / (SO + kho2) * SNO / (SNO + khn03) * XSTO;
		d[8][_idxSO] = bH * etahend
			* -kho2 / ((SO + kho2) * (SO + kho2))
			* SNO / (SNO + khn03) * XSTO;
		d[8][_idxSNO] = bH * etahend
			* kho2 / (SO + kho2)
			* khn03 / ((SNO + khn03) * (SNO + khn03)) * XSTO;
		d[8][_idxXSTO] = bH * etahend
			* kho2 / (SO + kho2)
			* SNO / (SNO + khn03);

		//reaction10: muAUT * SO / (SO + kno2) * SNH / (SNH + knnh4) * SALK / (SALK + knalk) * XA;
		d[9][_idxSO] = muAUT
			* kno2 / ((SO + kno2) * (SO + kno2))
			* SNH / (SNH + knnh4)
			* SALK / (SALK + knalk) * XA;
		d[9][_idxSALK] = muAUT
			* SO / (SO + kno2)
			* SNH / (SNH + knnh4)
			* knalk / ((SALK + knalk) * (SALK + knalk)) * XA;
		d[9][_idxSNH] = muAUT
			* SO / (SO + kno2)
			* SALK / (SALK + knalk)
			* knnh4 / ((SNH + knnh4) * (SNH + knnh4)) * XA;
		d[9][_idxXA] = muAUT
			* SO / (SO + kno2)
			* SNH / (SNH + knnh4)
			* SALK / (SALK + knalk);

		//reaction11: bAUT * SO / (SO + kho2) * XA;
		d[10][_idxSO] = bAUT
			* kho2 / ((SO + kho2) * (SO + kho2)) * XA;
		d[10][_idxXA] = bAUT
			* SO / (SO + kho2);

		//reaction12: bAUT * etanend * SNO / (SNO + khn03) * kho2 / (SO + kho2) * XA;
		d[11][_idxSO] = bAUT * etanend
			* SNO / (SNO + khn03)
			* -kho2 / ((SO + kho2) * (SO + kho2)) * XA;
		d[11][_idxSNO] = bAUT * etanend
			* kho2 / (SO + kho2)
			* khn03 / ((SNO + khn03) * (SNO + khn03)) * XA;
		d[11][_idxXA] = bAUT * etanend
			* SNO / (SNO + khn03)
			* kho2 / (SO + kho2);


		RowIterator curJac = jac;
		for (size_t rIdx = 0; rIdx < _stoichiometry.columns(); rIdx++)
		{
			RowIterator curJac = jac;
			for (int row = 0; row < _stoichiometry.rows(); ++row, ++curJac)
			{
				double colFactor = static_cast<double>(_stoichiometry.native(row, rIdx));
				for (size_t compIdx = 0; compIdx < _stoichiometry.rows(); compIdx++)
				{
					if (_fractionate)
					{
						if (compIdx == _idxSI_ad)
							colFactor *= SI_ad / (SI_ad + SI_nad);
						if (compIdx == _idxSI_nad)
							colFactor *= SI_nad / (SI_ad + SI_nad);
						if (compIdx == _idxSS_ad)
							colFactor *= SS_ad / (SS_ad + SS_nad);
						if (compIdx == _idxSS_nad)
							colFactor *= SS_nad / (SS_ad + SS_nad);
					}
					curJac[compIdx - static_cast<int>(row)] += colFactor * d[rIdx][compIdx];
				}
			}
		}

	}

	template <typename RowIteratorLiquid, typename RowIteratorSolid>
	void jacobianCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid, double factor, const RowIteratorLiquid& jacLiquid, const RowIteratorSolid& jacSolid, LinearBufferAllocator workSpace) const
	{
	}
};

typedef ActivatedSludgeModelThreeBase<ActivatedSludgeModelThreeParamHandler> ActivatedSludgeModelThreeReaction;
typedef ActivatedSludgeModelThreeBase<ExtActivatedSludgeModelThreeParamHandler> ExternalActivatedSludgeModelThreeReaction;

namespace reaction
{
	void registerActivatedSludgeModelThreeReaction(std::unordered_map<std::string, std::function<model::IDynamicReactionModel*()>>& reactions)
	{
		reactions[ActivatedSludgeModelThreeReaction::identifier()] = []() { return new ActivatedSludgeModelThreeReaction(); };
		reactions[ExternalActivatedSludgeModelThreeReaction::identifier()] = []() { return new ExternalActivatedSludgeModelThreeReaction(); };
	}
}  // namespace reaction

}  // namespace model

}  // namespace cadet

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

#include <catch.hpp>
#include "Approx.hpp"
#include "cadet/ExternalFunction.hpp"
#include "cadet/ParameterProvider.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "JacobianHelper.hpp"
#include "ReactionModelTests.hpp"
#include "UnitOperationTests.hpp"
#include "JsonTestModels.hpp"
#include "Utils.hpp"
#include "common/Driver.hpp"

#include "common/JsonParameterProvider.hpp"

#include "ReactionModelFactory.hpp"
#include "model/ReactionModel.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "AdUtils.hpp"
#include "AutoDiff.hpp"
#include "Memory.hpp"

#include <algorithm>
#include <numeric>
#include <cstring>

/**
 * @brief Returns the absolute path to the test/ folder of the project
 * @details Absolute path to the test/ folder of the project without trailing slash
 * @return Absolute path to the test/ folder
 */
const char* getTestDirectory();

namespace
{
	inline cadet::model::IDynamicReactionModel* createDynamicReactionModel(const char* name)
	{
		cadet::ReactionModelFactory rmf;
		cadet::model::IDynamicReactionModel* const rm = rmf.createDynamic(name);
		
		REQUIRE(nullptr != rm);
		return rm;
	}

	class ConstExternalFunction : public cadet::IExternalFunction
	{
	public:
		virtual bool configure(cadet::IParameterProvider* paramProvider) { return true; }
		virtual const char* name() const CADET_NOEXCEPT { return "CONSTFUN"; }
		virtual double externalProfile(double t, double z, double rho, double r, unsigned int sec) { return 1.0; }
		virtual double timeDerivative(double t, double z, double rho, double r, unsigned int sec) { return 0.0; }
		virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections) { }
	};

	class LinearExternalFunction : public cadet::IExternalFunction
	{
	public:
		virtual bool configure(cadet::IParameterProvider* paramProvider) { return true; }
		virtual const char* name() const CADET_NOEXCEPT { return "LINFUN"; }
		virtual double externalProfile(double t, double z, double rho, double r, unsigned int sec) { return t; }
		virtual double timeDerivative(double t, double z, double rho, double r, unsigned int sec) { return 1.0; }
		virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections) { }
	};
}

namespace cadet
{

namespace test
{

namespace reaction
{

	ConfiguredDynamicReactionModel::~ConfiguredDynamicReactionModel()
	{
		::operator delete(_bufferMemory);

		delete[] _extFuns;
		delete[] _boundOffset;
		delete _reaction;
	}

	ConfiguredDynamicReactionModel ConfiguredDynamicReactionModel::create(const char* name, unsigned int nComp, unsigned int const* nBound, const char* config)
	{
		cadet::model::IDynamicReactionModel* const rm = createDynamicReactionModel(name);

		// Calculate offset of bound states
		unsigned int* boundOffset = new unsigned int[nComp];
		boundOffset[0] = 0;
		for (unsigned int i = 1; i < nComp; ++i)
		{
			boundOffset[i] = boundOffset[i-1] + nBound[i-1];
		}
		const unsigned int totalBoundStates = boundOffset[nComp - 1] + nBound[nComp - 1];

		// Configure
		cadet::JsonParameterProvider jpp(config);
		rm->configureModelDiscretization(jpp, nComp, nBound, boundOffset);
		if (rm->requiresConfiguration())
		{
			jpp.set("EXTFUN", std::vector<int>(1, 0));
			REQUIRE(rm->configure(jpp, 0, 0));
		}

		// Assign external functions
		cadet::IExternalFunction* extFuns = new LinearExternalFunction[50];
		rm->setExternalFunctions(&extFuns, 50);

		// Allocate memory buffer
		unsigned int requiredMem = 0;
		if (rm->requiresWorkspace())
			requiredMem = rm->workspaceSize(nComp, totalBoundStates, boundOffset);

		void* buffer = nullptr;
		void* bufferEnd = nullptr;
		if (requiredMem > 0)
		{
			buffer = ::operator new(requiredMem);
			bufferEnd = static_cast<char*>(buffer) + requiredMem;
			std::memset(buffer, 0, requiredMem);
		}

		return ConfiguredDynamicReactionModel(rm, nComp, nBound, boundOffset, buffer, bufferEnd, extFuns);
	}

	void ConfiguredDynamicReactionModel::increaseBufferSize(int inc)
	{
		const int bufSize = requiredBufferSize() + inc;

		::operator delete(_bufferMemory);
		_bufferMemory = ::operator new(bufSize);
		std::memset(_bufferMemory, 0, bufSize);

#ifdef CADET_DEBUG
		_buffer.setBuffer(_bufferMemory, static_cast<char*>(_bufferMemory) + bufSize);
#else
		_buffer.setBuffer(_bufferMemory);
#endif
	}

	int ConfiguredDynamicReactionModel::requiredBufferSize() CADET_NOEXCEPT
	{
		if (_reaction->requiresWorkspace())
		{
			unsigned int totalBoundStates = 0;
			for (unsigned int i = 0; i < _nComp; ++i)
				totalBoundStates += _nBound[i];

			return _reaction->workspaceSize(_nComp, totalBoundStates, _boundOffset);
		}

		return 0;
	}

	void extendModelWithDynamicReactions(cadet::JsonParameterProvider& jpp, UnitOpIdx unit, bool bulk, bool particle, bool particleModifiers)
	{
		auto gs = util::makeModelGroupScope(jpp, unit);
		const int nComp = jpp.getInt("NCOMP");

		const std::string uoType = jpp.getString("UNIT_TYPE");
		const bool isLRM = (uoType == "LUMPED_RATE_MODEL_WITHOUT_PORES");
		const bool isCSTR = (uoType == "CSTR");

		if (bulk)
		{
			const int reactionTypes = 1;
			jpp.set("NREAC_LIQUID", reactionTypes);

			char const* const scopeName = "liquid_reaction_000";
			jpp.addScope(scopeName);
			auto gs2 = util::makeGroupScope(jpp, scopeName);

			jpp.set("TYPE", "MASS_ACTION_LAW");

			const int nReactions = 1;
			std::vector<double> stoichMat(nReactions * nComp, 0.0);
			std::vector<double> expFwd(nReactions * nComp, 0.0);
			std::vector<double> expBwd(nReactions * nComp, 0.0);
			std::vector<double> rateFwd(nReactions, 0.0);
			std::vector<double> rateBwd(nReactions, 0.0);

			util::populate(stoichMat.data(), [](unsigned int idx) { return 2.0 * std::sin(0.5 + idx * 0.7); }, stoichMat.size());
			util::populate(expFwd.data(), [](unsigned int idx) { return std::sin(0.1 + idx * 0.7) + 1.4; }, expFwd.size());
			util::populate(expBwd.data(), [](unsigned int idx) { return std::sin(0.2 + idx * 0.3) + 1.1; }, expBwd.size());
			util::populate(rateFwd.data(), [](unsigned int idx) { return std::sin(0.3 + idx * 0.2) * 0.5 + 1.5; }, rateFwd.size());
			util::populate(rateBwd.data(), [](unsigned int idx) { return std::sin(0.4 + idx * 0.4) * 0.5 + 1.6; }, rateBwd.size());

			jpp.set("MAL_STOICHIOMETRY", stoichMat);
			jpp.set("MAL_EXPONENTS_FWD", expFwd);
			jpp.set("MAL_EXPONENTS_BWD", expBwd);
			jpp.set("MAL_KFWD", rateFwd);
			jpp.set("MAL_KBWD", rateBwd);
		}

		if (!jpp.exists("NPARTYPE") && !isCSTR)
			return;

		int nParType = 0;
		if (!isCSTR)
			nParType = jpp.getInt("NPARTYPE");
		else
			nParType = jpp.getIntArray("NBOUND").size() / jpp.getInt("NCOMP");

		if ((!isLRM && particle))
		{
			const int nReactions = 1;

			for (int i = 0; i < nParType; ++i)
			{
				std::vector<int> nBound;
				int nTotalBound;
				if (isCSTR)
				{
					nBound = jpp.getIntArray("NBOUND");

					int start = i * nComp;
					int end = start + nComp;

					nTotalBound = std::accumulate(
						nBound.begin() + start,
						nBound.begin() + end,
						0
					);

					jpp.addScope("particle_type_" + std::string(3 - std::to_string(i).length(), '0') + std::to_string(i));
					jpp.pushScope("particle_type_" + std::string(3 - std::to_string(i).length(), '0') + std::to_string(i)); // particle_type_xxx
				}
				else
				{
					jpp.pushScope("particle_type_" + std::string(3 - std::to_string(i).length(), '0') + std::to_string(i)); // particle_type_xxx
					nBound = jpp.getIntArray("NBOUND");
					nTotalBound = std::accumulate(nBound.begin(), nBound.end(), 0);
				}

				int nReac = 1;
				jpp.set("NREAC_CROSS_PHASE", nReac);

				jpp.addScope("cross_phase_reaction_000"); //particle_type_xxx/cross_phase_reaction_000
				auto gs3 = util::makeGroupScope(jpp, "cross_phase_reaction_000");

				jpp.set("TYPE", "MASS_ACTION_LAW_CROSS_PHASE");

				std::vector<double> stoichMatLiquid(nReactions * nComp, 0.0);
				std::vector<double> expFwdLiquid(nReactions * nComp, 0.0);
				std::vector<double> expBwdLiquid(nReactions * nComp, 0.0);
				std::vector<double> rateFwdLiquid(nReactions, 0.0);
				std::vector<double> rateBwdLiquid(nReactions, 0.0);

				std::vector<double> stoichMatSolid(nReactions * nTotalBound, 0.0);
				std::vector<double> expFwdSolid(nReactions * nTotalBound, 0.0);
				std::vector<double> expBwdSolid(nReactions * nTotalBound, 0.0);
				std::vector<double> rateFwdSolid(nReactions, 0.0);
				std::vector<double> rateBwdSolid(nReactions, 0.0);

				util::populate(stoichMatLiquid.data(), [](unsigned int idx) { return 2.0 * std::sin(0.5 + idx * 0.7); }, stoichMatLiquid.size());
				util::populate(expFwdLiquid.data(), [](unsigned int idx) { return std::sin(0.1 + idx * 0.7) + 1.4; }, expFwdLiquid.size());
				util::populate(expBwdLiquid.data(), [](unsigned int idx) { return std::sin(0.2 + idx * 0.3) + 1.1; }, expBwdLiquid.size());
				util::populate(rateFwdLiquid.data(), [](unsigned int idx) { return std::sin(0.3 + idx * 0.2) * 0.5 + 1.5; }, rateFwdLiquid.size());
				util::populate(rateBwdLiquid.data(), [](unsigned int idx) { return std::sin(0.4 + idx * 0.4) * 0.5 + 1.6; }, rateBwdLiquid.size());

				util::populate(stoichMatSolid.data(), [](unsigned int idx) { return 2.0 * std::sin(0.3 + idx * 0.9); }, stoichMatSolid.size());
				util::populate(expFwdSolid.data(), [](unsigned int idx) { return std::sin(0.2 + idx * 0.6) + 1.3; }, expFwdSolid.size());
				util::populate(expBwdSolid.data(), [](unsigned int idx) { return std::sin(0.1 + idx * 0.4) + 1.2; }, expBwdSolid.size());
				util::populate(rateFwdSolid.data(), [](unsigned int idx) { return std::sin(0.5 + idx * 0.4) * 0.5 + 1.6; }, rateFwdSolid.size());
				util::populate(rateBwdSolid.data(), [](unsigned int idx) { return std::sin(0.2 + idx * 0.3) * 0.5 + 1.5; }, rateBwdSolid.size());

				jpp.set("MAL_STOICHIOMETRY_LIQUID", stoichMatLiquid);
				jpp.set("MAL_EXPONENTS_LIQUID_FWD", expFwdLiquid);
				jpp.set("MAL_EXPONENTS_LIQUID_BWD", expBwdLiquid);
				jpp.set("MAL_KFWD_LIQUID", rateFwdLiquid);
				jpp.set("MAL_KBWD_LIQUID", rateBwdLiquid);

				jpp.set("MAL_STOICHIOMETRY_SOLID", stoichMatSolid);
				jpp.set("MAL_EXPONENTS_SOLID_FWD", expFwdSolid);
				jpp.set("MAL_EXPONENTS_SOLID_BWD", expBwdSolid);
				jpp.set("MAL_KFWD_SOLID", rateFwdSolid);
				jpp.set("MAL_KBWD_SOLID", rateBwdSolid);

				if (particleModifiers)
				{
					std::vector<double> expFwdLiquidModSolid(nReactions * nTotalBound, 0.0);
					std::vector<double> expBwdLiquidModSolid(nReactions * nTotalBound, 0.0);
					std::vector<double> expFwdSolidModLiquid(nReactions * nComp, 0.0);
					std::vector<double> expBwdSolidModLiquid(nReactions * nComp, 0.0);

					util::populate(expFwdLiquidModSolid.data(), [](unsigned int idx) { return std::sin(idx * 0.7 + 0.1) * 0.5 + 1.6; }, expFwdLiquidModSolid.size());
					util::populate(expBwdLiquidModSolid.data(), [](unsigned int idx) { return std::sin(idx * 0.7 + 0.2) * 0.5 + 1.6; }, expBwdLiquidModSolid.size());
					util::populate(expFwdSolidModLiquid.data(), [](unsigned int idx) { return std::sin(idx * 0.7 + 0.3) * 0.5 + 1.6; }, expFwdSolidModLiquid.size());
					util::populate(expBwdSolidModLiquid.data(), [](unsigned int idx) { return std::sin(idx * 0.7 + 0.4) * 0.5 + 1.6; }, expBwdSolidModLiquid.size());

					jpp.set("MAL_EXPONENTS_LIQUID_FWD_MODSOLID", expFwdLiquidModSolid);
					jpp.set("MAL_EXPONENTS_LIQUID_BWD_MODSOLID", expBwdLiquidModSolid);
					jpp.set("MAL_EXPONENTS_SOLID_FWD_MODLIQUID", expFwdSolidModLiquid);
					jpp.set("MAL_EXPONENTS_SOLID_BWD_MODLIQUID", expBwdSolidModLiquid);
				}

					jpp.popScope();
			}
		}
	}
	
	void testMichaelisMentenToSMAMicroKinetic(const std::string configFilePathMM, const std::string configFilePathSMA, const double absTol, const double relTol)
	{
		// read json model setup file
		const std::string setupFileMM = std::string(getTestDirectory()) + configFilePathMM;
		const std::string setupFileSMA = std::string(getTestDirectory()) + configFilePathSMA;
		JsonParameterProvider pp_setup_MM(JsonParameterProvider::fromFile(setupFileMM));
		JsonParameterProvider pp_setup_SMA(JsonParameterProvider::fromFile(setupFileSMA));

		nlohmann::json* setupJsonMM = pp_setup_MM.data();
		nlohmann::json* setupJsonSMA = pp_setup_SMA.data();

		// MM simulation
		cadet::Driver drvMM;
		drvMM.configure(pp_setup_MM);
		drvMM.run();

		// SMA micro-kinetics simulation
		cadet::Driver drvSMA;
		drvSMA.configure(pp_setup_SMA);
		drvSMA.run();

		cadet::InternalStorageUnitOpRecorder const* const MMData = drvMM.solution()->unitOperation(0);
		cadet::InternalStorageUnitOpRecorder const* const SMAData = drvSMA.solution()->unitOperation(0);

		double const* outletMM = MMData->outlet();
		double const* outletSMA = SMAData->outlet();

		const unsigned int nCompMM = MMData->numComponents();
		const unsigned int nCompSMA = SMAData->numComponents();
		for (unsigned int i = 0; i < SMAData->numDataPoints(); ++i, outletMM += nCompMM, outletSMA += nCompSMA)
		{
			CAPTURE(i);
			CHECK((outletMM[0]) == cadet::test::makeApprox(outletSMA[0], relTol, absTol)); // substrate
			CHECK((outletMM[1]) == cadet::test::makeApprox(outletSMA[1], relTol, absTol)); // product
		}
	}

	void testDynamicJacobianAD(const char* modelName, unsigned int nComp, unsigned int const* nBound, const char* config, double const* point, double absTol, double relTol)
	{
		ConfiguredDynamicReactionModel crm = ConfiguredDynamicReactionModel::create(modelName, nComp, nBound, config);

		const unsigned int numDofs = crm.nComp() + crm.numBoundStates();
		std::vector<double> yState(numDofs, 0.0);
		std::copy_n(point, numDofs, yState.data());

		std::vector<double> dir(numDofs, 0.0);
		std::vector<double> colA(numDofs, 0.0);
		std::vector<double> colB(numDofs, 0.0);

		// Enable AD
		cadet::ad::setDirections(cadet::ad::getMaxDirections());
		cadet::active* adRes = new cadet::active[numDofs];
		cadet::active* adY = new cadet::active[numDofs];

		// Combined liquid and solid phase

		// Evaluate with AD
		ad::prepareAdVectorSeedsForDenseMatrix(adY, 0, numDofs);
		ad::copyToAd(yState.data(), adY, numDofs);
		ad::resetAd(adRes, numDofs);
		crm.model().residualCombinedAdd(1.0, 0u, ColumnPosition{ 0.0, 0.0, 0.0 }, adY, adY + crm.nComp(), adRes, adRes + crm.nComp(), 1.0, crm.buffer());

		// Extract Jacobian
		cadet::linalg::DenseMatrix jacAD;
		jacAD.resize(numDofs, numDofs);
		ad::extractDenseJacobianFromAd(adRes, 0, jacAD);

		// Calculate analytic Jacobian
		cadet::linalg::DenseMatrix jacAna;
		jacAna.resize(numDofs, numDofs);
		crm.model().analyticJacobianCombinedAdd(1.0, 0u, ColumnPosition{ 0.0, 0.0, 0.0 }, yState.data(), yState.data() + crm.nComp(), 1.0, jacAna.row(0), jacAna.row(crm.nComp()), crm.buffer());

		cadet::test::checkJacobianPatternFD(
			[&](double const* lDir, double* res) -> void
			{
				std::fill_n(res, numDofs, 0.0);
				crm.model().residualCombinedAdd(1.0, 0u, ColumnPosition{ 0.0, 0.0, 0.0 }, lDir, lDir + crm.nComp(), res, res + crm.nComp(), 1.0, crm.buffer());
			},
			[&](double const* lDir, double* res) -> void
			{
				jacAna.multiplyVector(lDir, res);
			},
			yState.data(), dir.data(), colA.data(), colB.data(), numDofs, numDofs);

		cadet::test::checkJacobianPatternFD(
			[&](double const* lDir, double* res) -> void
			{
				std::fill_n(res, numDofs, 0.0);
				crm.model().residualCombinedAdd(1.0, 0u, ColumnPosition{ 0.0, 0.0, 0.0 }, lDir, lDir + crm.nComp(), res, res + crm.nComp(), 1.0, crm.buffer());
			},
			[&](double const* lDir, double* res) -> void
			{
				jacAD.multiplyVector(lDir, res);
			},
			yState.data(), dir.data(), colA.data(), colB.data(), numDofs, numDofs);

		// Check Jacobians against each other
		for (unsigned int row = 0; row < numDofs; ++row)
		{
			for (unsigned int col = 0; col < numDofs; ++col)
			{
				CAPTURE(row);
				CAPTURE(col);
				CHECK(jacAna.native(row, col) == makeApprox(jacAD.native(row, col), absTol, relTol));
			}
		}

		// Only liquid phase

		// Evaluate with AD
		ad::resetAd(adRes, numDofs);
		crm.model().residualFluxAdd(1.0, 0u, ColumnPosition{0.0, 0.0, 0.0}, crm.nComp(), adY, adRes, 1.0, crm.buffer());

		// Extract Jacobian
		jacAD.setAll(0.0);
		ad::extractDenseJacobianFromAd(adRes, 0, jacAD);

		// Calculate analytic Jacobian
		jacAna.setAll(0.0);

		crm.model().analyticJacobianAdd(1.0, 0u, ColumnPosition{0.0, 0.0, 0.0}, crm.nComp(), yState.data(), 1.0, jacAna.row(0), crm.buffer());

		delete[] adY;
		delete[] adRes;

		cadet::test::checkJacobianPatternFD(
			[&](double const* lDir, double* res) -> void
				{
					std::fill_n(res, nComp, 0.0);
					crm.model().residualFluxAdd(1.0, 0u, ColumnPosition{0.0, 0.0, 0.0}, crm.nComp(), lDir, res, 1.0, crm.buffer());
				},
			[&](double const* lDir, double* res) -> void 
				{
					jacAna.submatrixMultiplyVector(lDir, 0, 0, nComp, nComp, res);
				},
			yState.data(), dir.data(), colA.data(), colB.data(), nComp, nComp);

		cadet::test::checkJacobianPatternFD(
			[&](double const* lDir, double* res) -> void
				{
					std::fill_n(res, nComp, 0.0);
					crm.model().residualFluxAdd(1.0, 0u, ColumnPosition{0.0, 0.0, 0.0}, crm.nComp(), lDir, res, 1.0, crm.buffer());
				},
			[&](double const* lDir, double* res) -> void 
				{
					jacAD.submatrixMultiplyVector(lDir, 0, 0, nComp, nComp, res);
				},
			yState.data(), dir.data(), colA.data(), colB.data(), nComp, nComp);

		// Check Jacobians against each other
		for (unsigned int row = 0; row < nComp; ++row)
		{
			for (unsigned int col = 0; col < nComp; ++col)
			{
				CAPTURE(row);
				CAPTURE(col);
				CHECK(jacAna.native(row, col) == makeApprox(jacAD.native(row, col), absTol, relTol));
			}
		}
	}

	void testUnitJacobianDynamicReactionsAD(cadet::JsonParameterProvider& jpp, bool bulk, bool particle, bool particleModifiers, const double absTolFDpattern)
	{
		extendModelWithDynamicReactions(jpp, 0, bulk, particle, particleModifiers);
		unitoperation::testJacobianAD(jpp, absTolFDpattern);
	}

	void testUnitJacobianDynamicReactionsAD(const std::string& uoType, const std::string& spatialMethod, bool bulk, bool particle, bool particleModifiers, const double absTolFDpattern)
	{
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding(uoType, spatialMethod);
		testUnitJacobianDynamicReactionsAD(jpp, bulk, particle, particleModifiers, absTolFDpattern);
	}

	void testTimeDerivativeJacobianDynamicReactionsFD(cadet::JsonParameterProvider& jpp, bool bulk, bool particle, bool particleModifiers, double h, double absTol, double relTol)
	{
		extendModelWithDynamicReactions(jpp, 0, bulk, particle, particleModifiers);
		unitoperation::testTimeDerivativeJacobianFD(jpp, h, absTol, relTol);
	}

	void testTimeDerivativeJacobianDynamicReactionsFD(const std::string& uoType, const std::string& spatialMethod, bool bulk, bool particle, bool particleModifiers, double h, double absTol, double relTol)
	{
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding(uoType, spatialMethod);
		testTimeDerivativeJacobianDynamicReactionsFD(jpp, bulk, particle, particleModifiers, h, absTol, relTol);
	}

	void testLiquidReactionJacobianAD(const char* modelName, unsigned int nComp, unsigned int const* nBound, const char* config, double const* point, double absTol, double relTol)
	{
		ConfiguredDynamicReactionModel crm = ConfiguredDynamicReactionModel::create(modelName, nComp, nBound, config);

		const unsigned int numDofs = crm.nComp();
		std::vector<double> yState(numDofs, 0.0);
		std::copy_n(point, numDofs, yState.data());

		std::vector<double> dir(numDofs, 0.0);
		std::vector<double> colA(numDofs, 0.0);
		std::vector<double> colB(numDofs, 0.0);

		// Enable AD
		cadet::ad::setDirections(cadet::ad::getMaxDirections());
		cadet::active* adRes = new cadet::active[numDofs];
		cadet::active* adY = new cadet::active[numDofs];

		// Liquid phase only

		// Evaluate with AD
		ad::prepareAdVectorSeedsForDenseMatrix(adY, 0, numDofs);
		ad::copyToAd(yState.data(), adY, numDofs);
		ad::resetAd(adRes, numDofs);
		crm.model().residualFluxAdd(1.0, 0u, ColumnPosition{0.0, 0.0, 0.0}, crm.nComp(), adY, adRes, 1.0, crm.buffer());

		// Extract Jacobian
		cadet::linalg::DenseMatrix jacAD;
		jacAD.resize(numDofs, numDofs);
		ad::extractDenseJacobianFromAd(adRes, 0, jacAD);

		// Calculate analytic Jacobian
		cadet::linalg::DenseMatrix jacAna;
		jacAna.resize(numDofs, numDofs);
		crm.model().analyticJacobianAdd(1.0, 0u, ColumnPosition{0.0, 0.0, 0.0}, crm.nComp(), yState.data(), 1.0, jacAna.row(0), crm.buffer());

		delete[] adY;
		delete[] adRes;

		cadet::test::checkJacobianPatternFD(
			[&](double const* lDir, double* res) -> void
			{
				std::fill_n(res, numDofs, 0.0);
				crm.model().residualFluxAdd(1.0, 0u, ColumnPosition{0.0, 0.0, 0.0}, crm.nComp(), lDir, res, 1.0, crm.buffer());
			},
			[&](double const* lDir, double* res) -> void 
			{
				jacAna.multiplyVector(lDir, res);
			},
			yState.data(), dir.data(), colA.data(), colB.data(), numDofs, numDofs);

		cadet::test::checkJacobianPatternFD(
			[&](double const* lDir, double* res) -> void
			{
				std::fill_n(res, numDofs, 0.0);
				crm.model().residualFluxAdd(1.0, 0u, ColumnPosition{0.0, 0.0, 0.0}, crm.nComp(), lDir, res, 1.0, crm.buffer());
			},
			[&](double const* lDir, double* res) -> void 
			{
				jacAD.multiplyVector(lDir, res);
			},
			yState.data(), dir.data(), colA.data(), colB.data(), numDofs, numDofs);

		// Check Jacobians against each other
		for (unsigned int row = 0; row < numDofs; ++row)
		{
			for (unsigned int col = 0; col < numDofs; ++col)
			{
				CAPTURE(row);
				CAPTURE(col);
				CHECK(jacAna.native(row, col) == makeApprox(jacAD.native(row, col), absTol, relTol));
			}
		}
	}

} // namespace reaction
} // namespace test
} // namespace cadet

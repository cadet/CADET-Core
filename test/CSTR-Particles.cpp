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
#include "cadet/cadet.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "ColumnTests.hpp"
#include "JacobianHelper.hpp"
#include "JsonTestModels.hpp"
#include "ParticleHelper.hpp"
#include "SimHelper.hpp"
#include "UnitOperationTests.hpp"
#include "Utils.hpp"
#include "common/Driver.hpp"
#include "model/UnitOperation.hpp"
#include "ModelBuilderImpl.hpp"
#include "AdUtils.hpp"
#include "SimulationTypes.hpp"
#include "ParallelSupport.hpp"

#include <cmath>
#include <functional>
#include <iomanip>
#include <sstream>
#include <vector>

namespace
{
	/**
	 * @brief Creates and configures a finite bath and sets its flow rates
	 * @param [in] mb ModelBuilder
	 * @param [in,out] jpp Configuration of the unit operation
	 * @param [in] flowRate Volumetric flow rate through the vessel
	 * @return Runnable unit operation
	 */
	cadet::IUnitOperation* createAndConfigureStirredTank(cadet::IModelBuilder& mb, cadet::JsonParameterProvider& jpp, double flowRate)
	{
		cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(jpp, mb);

		const cadet::active in = flowRate;
		const cadet::active out = flowRate;
		unit->setFlowRates(&in, &out);

		return unit;
	}

	/**
	 * @brief Compares the analytic Jacobian of a finite bath against the one obtained by AD
	 * @details In contrast to the generic test, the flow rates are set, which only a configured
	 *          network would do otherwise.
	 * @param [in,out] jpp Configuration of the unit operation
	 * @param [in] flowRate Volumetric flow rate through the vessel
	 */
	void checkJacobianAD(cadet::JsonParameterProvider& jpp, double flowRate)
	{
		cadet::IModelBuilder* const mb = cadet::createModelBuilder();
		REQUIRE(nullptr != mb);

		cadet::ad::setDirections(cadet::ad::getMaxDirections());

		cadet::IUnitOperation* const unitAna = createAndConfigureStirredTank(*mb, jpp, flowRate);
		cadet::IUnitOperation* const unitAD = createAndConfigureStirredTank(*mb, jpp, flowRate);
		unitAD->useAnalyticJacobian(false);

		REQUIRE(unitAD->requiredADdirs() <= cadet::ad::getMaxDirections());

		cadet::active* adRes = new cadet::active[unitAD->numDofs()];
		cadet::active* adY = new cadet::active[unitAD->numDofs()];

		const cadet::AdJacobianParams noParams{ nullptr, nullptr, 0u };
		const cadet::AdJacobianParams adParams{ adRes, adY, 0u };

		unitAD->prepareADvectors(adParams);

		const unsigned int nDof = unitAna->numDofs();
		std::vector<double> y(nDof, 0.0);
		std::vector<double> jacDir(nDof, 0.0);
		std::vector<double> jacCol1(nDof, 0.0);
		std::vector<double> jacCol2(nDof, 0.0);
		cadet::util::ThreadLocalStorage tls;
		tls.resize(unitAna->threadLocalMemorySize());

		cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, nDof);

		unitAna->notifyDiscontinuousSectionTransition(0.0, 0u, { y.data(), nullptr }, noParams);
		unitAD->notifyDiscontinuousSectionTransition(0.0, 0u, { y.data(), nullptr }, adParams);

		unitAna->residualWithJacobian(cadet::SimulationTime{ 0.0, 0u }, cadet::ConstSimulationState{ y.data(), nullptr }, jacDir.data(), noParams, tls);
		unitAD->residualWithJacobian(cadet::SimulationTime{ 0.0, 0u }, cadet::ConstSimulationState{ y.data(), nullptr }, jacDir.data(), adParams, tls);
		std::fill(jacDir.begin(), jacDir.end(), 0.0);

		// The analytic Jacobian has to have at least the pattern of a finite difference Jacobian
		cadet::test::checkJacobianPatternFD(unitAna, unitAD, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data(), tls);
		cadet::test::checkJacobianPatternFD(unitAna, unitAna, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data(), tls);
		cadet::test::compareJacobian(unitAna, unitAD, nullptr, nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());

		delete[] adRes;
		delete[] adY;
		mb->destroyUnitOperation(unitAna);
		mb->destroyUnitOperation(unitAD);
		destroyModelBuilder(mb);
	}

	/**
	 * @brief Sets the volumetric flow rates of the connections of a model system
	 * @details In contrast to cadet::test::setFlowRates(), this does not require the unit operation to
	 *          have a FLOWRATE_FILTER field, which a finite bath does not have.
	 * @param [in,out] jpp Parameter provider of the full model
	 * @param [in] secIdx Index of the section
	 * @param [in] in Flow rate into unit_000
	 * @param [in] out Flow rate out of unit_000
	 */
	void setConnectionFlowRates(cadet::JsonParameterProvider& jpp, unsigned int secIdx, double in, double out)
	{
		jpp.pushScope("model");
		jpp.pushScope("connections");

		std::ostringstream ss;
		ss << "switch_" << std::setfill('0') << std::setw(3) << secIdx;
		jpp.pushScope(ss.str());

		std::vector<double> con = jpp.getDoubleArray("CONNECTIONS");
		con[6] = in;
		con[13] = out;
		jpp.set("CONNECTIONS", con);

		jpp.popScope();
		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Compares the time derivative Jacobian of a finite bath against centered finite differences
	 */
	void checkTimeDerivativeJacobianFD(cadet::JsonParameterProvider& jpp, double flowRate, double h, double absTol, double relTol)
	{
		cadet::IModelBuilder* const mb = cadet::createModelBuilder();
		REQUIRE(nullptr != mb);

		cadet::IUnitOperation* const unit = createAndConfigureStirredTank(*mb, jpp, flowRate);

		const unsigned int nDof = unit->numDofs();
		std::vector<double> y(nDof, 0.0);
		std::vector<double> yDot(nDof, 0.0);
		std::vector<double> jacDir(nDof, 0.0);
		std::vector<double> jacCol1(nDof, 0.0);
		std::vector<double> jacCol2(nDof, 0.0);
		cadet::util::ThreadLocalStorage tls;
		tls.resize(unit->threadLocalMemorySize());

		cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, nDof);
		cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, nDof);

		const cadet::AdJacobianParams noParams{ nullptr, nullptr, 0u };
		unit->notifyDiscontinuousSectionTransition(0.0, 0u, { y.data(), yDot.data() }, noParams);

		cadet::test::compareTimeDerivativeJacobianFD(unit, unit, y.data(), yDot.data(), jacDir.data(), jacCol1.data(), jacCol2.data(), tls, h, absTol, relTol);

		mb->destroyUnitOperation(unit);
		destroyModelBuilder(mb);
	}
}

TEST_CASE("CSTR with homogeneous particles Jacobian vs AD", "[CSTR],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("HOMOGENEOUS_PARTICLE", "FV");
	checkJacobianAD(jpp, 1e-6);
}

TEST_CASE("CSTR with general rate particles and FV discretization Jacobian vs AD", "[CSTR],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");
	checkJacobianAD(jpp, 1e-6);
}

TEST_CASE("CSTR with general rate particles and DG discretization Jacobian vs AD", "[CSTR],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "DG");
	checkJacobianAD(jpp, 1e-6);
}

TEST_CASE("CSTR with two particle types Jacobian vs AD", "[CSTR],[UnitOp],[Residual],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");

	const double volFrac[] = { 0.3, 0.7 };
	const double parFactor[] = { 0.9 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 2, parFactor, volFrac);

	checkJacobianAD(jpp, 1e-6);
}

TEST_CASE("CSTR without flow Jacobian vs AD", "[CSTR],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("HOMOGENEOUS_PARTICLE", "FV");
	checkJacobianAD(jpp, 0.0);
}

TEST_CASE("CSTR with variable volume Jacobian vs AD", "[CSTR],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	// The constant volume parameterization rejects flow rates that do not cancel, so only the variable
	// volume one reaches a non-vanishing \dot{V}, which the liquid volume column of the Jacobian
	// depends on
	for (const std::string particleType : { "HOMOGENEOUS_PARTICLE", "GENERAL_RATE_PARTICLE" })
	{
		SECTION("Particle type " + particleType)
		{
			cadet::JsonParameterProvider jpp = createFiniteBath(particleType, "FV");
			jpp.remove("LIQUID_VOLUME");
			jpp.remove("COL_POROSITY");
			jpp.set("INIT_LIQUID_VOLUME", 1e-3);
			jpp.set("CONST_SOLID_VOLUME", 1.7027027027027027e-3);
			jpp.set("FLOWRATE_FILTER", std::vector<double>{ 1e-5 });

			checkJacobianAD(jpp, 1e-6);
		}
	}
}

TEST_CASE("CSTR with variable volume time derivative Jacobian vs FD", "[CSTR],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");
	jpp.remove("LIQUID_VOLUME");
	jpp.remove("COL_POROSITY");
	jpp.set("INIT_LIQUID_VOLUME", 1e-3);
	jpp.set("CONST_SOLID_VOLUME", 1.7027027027027027e-3);
	jpp.set("FLOWRATE_FILTER", std::vector<double>{ 1e-5 });

	checkTimeDerivativeJacobianFD(jpp, 1e-6, 1e-6, 0.0, 1e-4);
}

TEST_CASE("CSTR with homogeneous particles time derivative Jacobian vs FD", "[CSTR],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("HOMOGENEOUS_PARTICLE", "FV");
	checkTimeDerivativeJacobianFD(jpp, 1e-6, 1e-6, 0.0, 1e-4);
}

TEST_CASE("CSTR with general rate particles time derivative Jacobian vs FD", "[CSTR],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");
	checkTimeDerivativeJacobianFD(jpp, 1e-6, 1e-6, 0.0, 1e-4);
}

TEST_CASE("CSTR inlet DOF Jacobian", "[CSTR],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");
	cadet::IUnitOperation* const unit = createAndConfigureStirredTank(*mb, jpp, 1e-6);

	cadet::test::unitoperation::testInletDofJacobian(unit, false);

	mb->destroyUnitOperation(unit);
	destroyModelBuilder(mb);
}

TEST_CASE("CSTR consistent initialization with linear binding", "[CSTR],[ConsistentInit],[CI]")
{
	for (int adMode = 0; adMode < 2; ++adMode)
	{
		const bool adEnabled = (adMode > 0);

		for (int bindingMode = 0; bindingMode < 2; ++bindingMode)
		{
			const bool isKinetic = (bindingMode == 0);

			SECTION(std::string(isKinetic ? "Kinetic binding" : "Quasi-stationary binding") + " with AD " + (adEnabled ? "enabled" : "disabled"))
			{
				cadet::IModelBuilder* const mb = cadet::createModelBuilder();
				REQUIRE(nullptr != mb);

				cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");
				cadet::test::setBindingMode(jpp, isKinetic);

				if (adEnabled)
					cadet::ad::setDirections(cadet::ad::getMaxDirections());

				cadet::IUnitOperation* const unit = createAndConfigureStirredTank(*mb, jpp, 1e-6);

				// Fill state vector with values
				std::vector<double> y(unit->numDofs(), 0.0);
				cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, unit->numDofs());

				cadet::test::unitoperation::testConsistentInitialization(unit, adEnabled, y.data(), 1e-14, 1e-10);

				mb->destroyUnitOperation(unit);
				destroyModelBuilder(mb);
			}
		}
	}
}

TEST_CASE("CSTR consistent initialization with a non-binding component", "[CSTR],[ConsistentInit],[CI]")
{
	for (int adMode = 0; adMode < 2; ++adMode)
	{
		const bool adEnabled = (adMode > 0);

		SECTION(std::string("Non-binding component with AD ") + (adEnabled ? "enabled" : "disabled"))
		{
			cadet::IModelBuilder* const mb = cadet::createModelBuilder();
			REQUIRE(nullptr != mb);

			cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");
			cadet::test::setBindingMode(jpp, false);
			{
				auto ms = cadet::test::util::makeOptionalGroupScope(jpp, "model");
				auto us = cadet::test::util::makeOptionalGroupScope(jpp, "unit_000");
				auto ps = cadet::test::util::makeOptionalGroupScope(jpp, "particle_type_000");

				// Bound-state indices differ from component indices, and binding is strong enough
				// that an inconsistent Jacobian of the conservation relations prevents convergence
				jpp.set("NBOUND", std::vector<int>{ 0, 1 });
				jpp.set("INIT_CS", std::vector<double>{ 6.0 });
				jpp.set("SURFACE_DIFFUSION", std::vector<double>{ 5e-11 });

				jpp.pushScope("adsorption");
				jpp.set("LIN_KA", std::vector<double>{ 0.0, 35.5 });
				jpp.set("LIN_KD", std::vector<double>{ 1.0, 1.0 });
				jpp.popScope();
			}

			if (adEnabled)
				cadet::ad::setDirections(cadet::ad::getMaxDirections());

			cadet::IUnitOperation* const unit = createAndConfigureStirredTank(*mb, jpp, 1e-6);

			std::vector<double> y(unit->numDofs(), 0.0);
			cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, unit->numDofs());

			cadet::test::unitoperation::testConsistentInitialization(unit, adEnabled, y.data(), 1e-14, 1e-10);

			mb->destroyUnitOperation(unit);
			destroyModelBuilder(mb);
		}
	}
}

TEST_CASE("CSTR consistent sensitivity initialization with linear binding", "[CSTR],[ConsistentInit],[Sensitivity],[CI]")
{
	for (int bindingMode = 0; bindingMode < 2; ++bindingMode)
	{
		const bool isKinetic = (bindingMode == 0);

		SECTION(std::string(isKinetic ? "Kinetic binding" : "Quasi-stationary binding"))
		{
			cadet::IModelBuilder* const mb = cadet::createModelBuilder();
			REQUIRE(nullptr != mb);

			cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");
			cadet::test::setBindingMode(jpp, isKinetic);

			cadet::ad::setDirections(cadet::ad::getMaxDirections());

			cadet::IUnitOperation* const unit = createAndConfigureStirredTank(*mb, jpp, 1e-6);

			unit->setSensitiveParameter(cadet::makeParamId("INIT_C", 0, 0, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, cadet::SectionIndep), 0, 1.0);
			unit->setSensitiveParameter(cadet::makeParamId("LIQUID_VOLUME", 0, cadet::CompIndep, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, cadet::SectionIndep), 1, 1.0);
			REQUIRE(unit->numSensParams() == 2);

			std::vector<double> y(unit->numDofs(), 0.0);
			std::vector<double> yDot(unit->numDofs(), 0.0);
			cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, unit->numDofs());
			cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, unit->numDofs());

			cadet::test::unitoperation::testConsistentInitializationSensitivity(unit, true, y.data(), yDot.data(), 1e-10);

			mb->destroyUnitOperation(unit);
			destroyModelBuilder(mb);
		}
	}
}

TEST_CASE("CSTR batch uptake vs analytic solution", "[CSTR],[Simulation],[Analytic],[CI]")
{
	// A closed vessel without binding holds a single homogeneous particle type. The bulk and particle
	// liquid phase then form a linear 2x2 system,
	//     dc / dt   = -A (c - c_p),    A = (1 / eps_b - 1) * a_p * k_f
	//     dc_p / dt =  B (c - c_p),    B = a_p * k_f / eps_p
	// whose solution decays towards the common equilibrium c_inf = (B c_0 + A cp_0) / (A + B) as
	//     c(t)   = c_inf + A / (A + B) * (c_0 - cp_0) * exp(-(A + B) t)
	//     c_p(t) = c_inf - B / (A + B) * (c_0 - cp_0) * exp(-(A + B) t)
	const double parRadius = 4.5e-5;
	const double parPorosity = 0.75;
	const double bulkPorosity = 0.37;
	const double filmDiffusion = 6.9e-6;

	const double surfToVol = 3.0 / parRadius;
	const double A = (1.0 / bulkPorosity - 1.0) * surfToVol * filmDiffusion;
	const double B = surfToVol * filmDiffusion / parPorosity;

	const double c0 = 1.0;
	const double cp0 = 0.5;
	const double cInf = (B * c0 + A * cp0) / (A + B);

	cadet::JsonParameterProvider jpp = createFiniteBathBenchmark("HOMOGENEOUS_PARTICLE", "FV", false, 1, 1e-3, 1e-5);
	cadet::test::setSectionTimes(jpp, { 0.0, 1e-3 });
	setConnectionFlowRates(jpp, 0, 0.0, 0.0);

	cadet::Driver drv;
	drv.configure(jpp);
	drv.run();

	cadet::InternalStorageUnitOpRecorder const* const simData = drv.solution()->unitOperation(0);
	double const* outlet = simData->outlet();
	double const* time = drv.solution()->time();

	REQUIRE(simData->numDataPoints() > 1);

	for (unsigned int i = 0; i < simData->numDataPoints(); ++i, outlet += 2, ++time)
	{
		const double ref = cInf + A / (A + B) * (c0 - cp0) * std::exp(-(A + B) * (*time));

		CAPTURE(*time);
		CHECK(outlet[0] == cadet::test::makeApprox(ref, 1e-6, 1e-10));
	}
}

namespace
{
	const double closedTankLiquidVolume = 2.0;
	const double closedTankSolidVolume = 3.0;
	const double closedTankEndTime = 10.0;

	/**
	 * @brief Creates a closed stirred tank with two components and one bound state each
	 * @details The particle type is in rapid equilibrium with the bulk liquid and has no binding
	 *          model, so nothing moves unless a reaction moves it. The two volumes differ, which is
	 *          what makes a term that is weighted with the wrong one visible.
	 * @return Runnable model
	 */
	cadet::JsonParameterProvider createClosedTankWithBoundStates()
	{
		cadet::JsonParameterProvider jpp = createCSTRBenchmark(1, closedTankEndTime, closedTankEndTime / 100.0);

		for (unsigned int unit = 0; unit < 3; ++unit)
			cadet::test::setNumberOfComponents(jpp, unit, 2);

		cadet::test::setSectionTimes(jpp, { 0.0, closedTankEndTime });
		cadet::test::addBoundStates(jpp, { 1, 1 }, closedTankSolidVolume);
		cadet::test::setInitialConditions(jpp, { 2.0, 0.5 }, { 1.5, 0.25 }, closedTankLiquidVolume);
		cadet::test::setInletProfile(jpp, 0, 0, 0.0, 0.0, 0.0, 0.0);
		cadet::test::setInletProfile(jpp, 0, 1, 0.0, 0.0, 0.0, 0.0);
		cadet::test::setFlowRates(jpp, 0, 0.0, 0.0, 0.0);

		return jpp;
	}

	/**
	 * @brief Enters the reaction group of the only particle type of the only tank
	 */
	void pushParticleReactionScope(cadet::JsonParameterProvider& jpp, const std::string& phase, const std::string& type)
	{
		jpp.pushScope("model");
		jpp.pushScope("unit_000");
		jpp.pushScope("particle_type_000");

		jpp.set("NREAC_" + phase, 1);

		const std::string group = (phase == "SOLID") ? "solid_reaction_000" : "cross_phase_reaction_000";
		jpp.addScope(group);
		jpp.pushScope(group);

		jpp.set("TYPE", type);
	}

	void popParticleReactionScope(cadet::JsonParameterProvider& jpp)
	{
		jpp.popScope();
		jpp.popScope();
		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Compares the bulk liquid and the bound states against a reference
	 */
	void checkClosedTank(cadet::JsonParameterProvider& jpp, const std::function<double(double, unsigned int)>& solC, const std::function<double(double, unsigned int)>& solQ)
	{
		cadet::Driver drv;
		drv.configure(jpp);
		drv.run();

		cadet::InternalStorageUnitOpRecorder const* const simData = drv.solution()->unitOperation(0);
		double const* outlet = simData->outlet();
		double const* solid = simData->solid();
		double const* time = drv.solution()->time();

		REQUIRE(simData->numDataPoints() > 1);

		for (unsigned int i = 0; i < simData->numDataPoints(); ++i, outlet += 2, solid += 2, ++time)
		{
			CAPTURE(*time);

			for (unsigned int comp = 0; comp < 2; ++comp)
			{
				CAPTURE(comp);
				CHECK(outlet[comp] == cadet::test::makeApprox(solC(*time, comp), 1e-6, 1e-8));
				CHECK(solid[comp] == cadet::test::makeApprox(solQ(*time, comp), 1e-6, 1e-8));
			}
		}
	}
}

TEST_CASE("CSTR solid phase reaction vs analytic solution", "[CSTR],[Simulation],[Analytic],[ReactionModel],[CI]")
{
	// A reaction between bound states converts the first into the second,
	//     dq_0 / dt = -k q_0,    dq_1 / dt = k q_0.
	// It takes nothing from the bulk liquid, which the balance of the hosting tank has to reflect:
	// the uptake term V^p sum_m \dot{q}_m assumes that whatever the bound states gain comes out of
	// the liquid, so the reaction flux has to be subtracted there again. A tank that leaves that out
	// moves its liquid phase instead of holding it.
	const double k = 0.7;
	const double q0 = 1.5;
	const double q1 = 0.25;

	cadet::JsonParameterProvider jpp = createClosedTankWithBoundStates();

	pushParticleReactionScope(jpp, "SOLID", "MASS_ACTION_LAW");
	jpp.set("MAL_STOICHIOMETRY", std::vector<double>{ -1.0, 1.0 });
	jpp.set("MAL_KFWD", std::vector<double>{ k });
	jpp.set("MAL_KBWD", std::vector<double>{ 0.0 });
	popParticleReactionScope(jpp);

	checkClosedTank(jpp,
		[](double t, unsigned int comp) {
			return (comp == 0) ? 2.0 : 0.5;
		},
		[=](double t, unsigned int comp) {
			const double decay = q0 * std::exp(-k * t);
			return (comp == 0) ? decay : (q1 + q0 - decay);
	});
}

TEST_CASE("CSTR cross phase reaction vs analytic solution", "[CSTR],[Simulation],[Analytic],[ReactionModel],[CI]")
{
	// The reaction converts the first component into the second in both phases at once, at a rate of
	// its own per phase. The liquid phase part acts on the bulk liquid, because a particle in rapid
	// equilibrium has none of its own, so it is weighted with V^l rather than with V^p.
	const double kLiquid = 0.9;
	const double kSolid = 0.4;
	const double c0 = 2.0;
	const double c1 = 0.5;
	const double q0 = 1.5;
	const double q1 = 0.25;

	cadet::JsonParameterProvider jpp = createClosedTankWithBoundStates();

	pushParticleReactionScope(jpp, "CROSS_PHASE", "MASS_ACTION_LAW_CROSS_PHASE");
	jpp.set("MAL_STOICHIOMETRY_LIQUID", std::vector<double>{ -1.0, 1.0 });
	jpp.set("MAL_KFWD_LIQUID", std::vector<double>{ kLiquid });
	jpp.set("MAL_KBWD_LIQUID", std::vector<double>{ 0.0 });
	jpp.set("MAL_STOICHIOMETRY_SOLID", std::vector<double>{ -1.0, 1.0 });
	jpp.set("MAL_KFWD_SOLID", std::vector<double>{ kSolid });
	jpp.set("MAL_KBWD_SOLID", std::vector<double>{ 0.0 });
	popParticleReactionScope(jpp);

	checkClosedTank(jpp,
		[=](double t, unsigned int comp) {
			const double decay = c0 * std::exp(-kLiquid * t);
			return (comp == 0) ? decay : (c1 + c0 - decay);
		},
		[=](double t, unsigned int comp) {
			const double decay = q0 * std::exp(-kSolid * t);
			return (comp == 0) ? decay : (q1 + q0 - decay);
	});
}

TEST_CASE("CSTR batch uptake numerical Benchmark for rapid equilibrium Langmuir case", "[CSTR],[DG],[Simulation],[Reference],[CI]")
{
	// Uptake of m-xylene from solution by Y zeolite beads in a closed, well stirred vessel, taken from
	// S. Brandani, "Kinetics of liquid phase batch adsorption experiments", Adsorption 27 (2021) 353-368,
	// Fig. 2, the m-xylene curve at Sh = 2. The beads are general rate particles with a Langmuir isotherm
	// in local equilibrium with the pore liquid, and the vessel is closed, so both flow rates are zero.
	// The derivation of the CADET parameters from the parameters of the paper and the comparison against
	// the published curve are part of the CADET-Verification study of the same name.
	const std::string& modelFilePath = std::string("/data/model_finiteBath_GRP_reqLangmuir_1comp_Brandani2021.json");
	const std::string& refFilePath = std::string("/data/ref_finiteBath_GRP_reqLangmuir_1comp_Brandani2021_parDGP4Z8.h5");
	const std::vector<double> absTol = { 1e-10 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::DGParams disc(-1, 0, 0, 4, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("CSTR constant volume parameterization matches variable volume without particles", "[CSTR],[Simulation],[CI]")
{
	// A prescribed liquid volume has to give what a variable one whose flow rates cancel gives
	cadet::JsonParameterProvider jppBath = createCSTRBenchmark(1, 100.0, 1.0);
	{
		jppBath.pushScope("model");
		jppBath.pushScope("unit_000");

		jppBath.remove("INIT_LIQUID_VOLUME");
		jppBath.remove("CONST_SOLID_VOLUME");
		jppBath.remove("FLOWRATE_FILTER");
		jppBath.set("LIQUID_VOLUME", 1.0);
		jppBath.set("NPARTYPE", 0);
		jppBath.addScope("discretization");
		jppBath.pushScope("discretization");
		jppBath.set("USE_ANALYTIC_JACOBIAN", true);
		jppBath.popScope();

		jppBath.popScope();
		jppBath.popScope();
	}

	// The reference is the variable volume CSTR, whose volume stays constant because the flow rates cancel
	cadet::JsonParameterProvider jppCstr = createCSTRBenchmark(1, 100.0, 1.0);

	for (cadet::JsonParameterProvider* jpp : { &jppBath, &jppCstr })
	{
		cadet::test::setSectionTimes(*jpp, { 0.0, 100.0 });
		cadet::test::setInletProfile(*jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
		setConnectionFlowRates(*jpp, 0, 0.1, 0.1);
	}

	cadet::Driver drvBath;
	drvBath.configure(jppBath);
	drvBath.run();

	cadet::Driver drvCstr;
	drvCstr.configure(jppCstr);
	drvCstr.run();

	cadet::InternalStorageUnitOpRecorder const* const dataBath = drvBath.solution()->unitOperation(0);
	cadet::InternalStorageUnitOpRecorder const* const dataCstr = drvCstr.solution()->unitOperation(0);

	REQUIRE(dataBath->numDataPoints() == dataCstr->numDataPoints());
	REQUIRE(dataBath->numDataPoints() > 1);

	double const* outletBath = dataBath->outlet();
	double const* outletCstr = dataCstr->outlet();

	for (unsigned int i = 0; i < dataBath->numDataPoints(); ++i, ++outletBath, ++outletCstr)
	{
		CAPTURE(i);
		CHECK(*outletBath == cadet::test::makeApprox(*outletCstr, 1e-8, 1e-12));
	}
}

TEST_CASE("CSTR analytic Jacobian matches AD with a moving volume", "[CSTR],[Simulation],[AD],[CI]")
{
	// The flow rates do not cancel, so the liquid volume genuinely changes over the simulation and
	// the \dot{V} c term of the extensive bulk balance is exercised
	cadet::JsonParameterProvider jppBath = createCSTRBenchmark(1, 100.0, 1.0);
	{
		jppBath.pushScope("model");
		jppBath.pushScope("unit_000");

		jppBath.set("NPARTYPE", 0);
		jppBath.addScope("discretization");
		jppBath.pushScope("discretization");
		jppBath.set("USE_ANALYTIC_JACOBIAN", true);
		jppBath.popScope();

		jppBath.popScope();
		jppBath.popScope();
	}

	// The reference uses the AD Jacobian, which is what a CSTR input without a discretization group gets
	cadet::JsonParameterProvider jppCstr = createCSTRBenchmark(1, 100.0, 1.0);
	{
		jppCstr.pushScope("model");
		jppCstr.pushScope("unit_000");

		jppCstr.set("NPARTYPE", 0);
		jppCstr.addScope("discretization");
		jppCstr.pushScope("discretization");
		jppCstr.set("USE_ANALYTIC_JACOBIAN", false);
		jppCstr.popScope();

		jppCstr.popScope();
		jppCstr.popScope();
	}

	for (cadet::JsonParameterProvider* jpp : { &jppBath, &jppCstr })
	{
		cadet::test::setSectionTimes(*jpp, { 0.0, 100.0 });
		cadet::test::setInletProfile(*jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
		setConnectionFlowRates(*jpp, 0, 0.2, 0.1);

		jpp->pushScope("model");
		jpp->pushScope("unit_000");
		jpp->set("FLOWRATE_FILTER", std::vector<double>{ 0.05 });
		jpp->popScope();
		jpp->popScope();
	}

	cadet::Driver drvBath;
	drvBath.configure(jppBath);
	drvBath.run();

	cadet::Driver drvCstr;
	drvCstr.configure(jppCstr);
	drvCstr.run();

	cadet::InternalStorageUnitOpRecorder const* const dataBath = drvBath.solution()->unitOperation(0);
	cadet::InternalStorageUnitOpRecorder const* const dataCstr = drvCstr.solution()->unitOperation(0);

	REQUIRE(dataBath->numDataPoints() == dataCstr->numDataPoints());
	REQUIRE(dataBath->numDataPoints() > 1);

	double const* outletBath = dataBath->outlet();
	double const* outletCstr = dataCstr->outlet();
	double const* volumeBath = dataBath->volume();
	double const* volumeCstr = dataCstr->volume();

	REQUIRE(volumeBath != nullptr);
	REQUIRE(volumeCstr != nullptr);

	// The volume has to move, otherwise this test would not differ from the constant volume one
	CHECK(std::abs(volumeBath[dataBath->numDataPoints() - 1] - volumeBath[0]) > 1.0);

	for (unsigned int i = 0; i < dataBath->numDataPoints(); ++i, ++outletBath, ++outletCstr, ++volumeBath, ++volumeCstr)
	{
		CAPTURE(i);
		CHECK(*outletBath == cadet::test::makeApprox(*outletCstr, 1e-8, 1e-12));
		CHECK(*volumeBath == cadet::test::makeApprox(*volumeCstr, 1e-8, 1e-12));
	}
}

TEST_CASE("CSTR with particles in rapid equilibrium Jacobian vs AD", "[CSTR],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	// A particle without a film diffusion resistance is in rapid equilibrium with the bulk liquid and
	// holds bound states only, so the state vector of the tank has no particle liquid phase at all
	cadet::JsonParameterProvider jpp = createFiniteBath("HOMOGENEOUS_PARTICLE", "FV");

	jpp.pushScope("particle_type_000");
	jpp.set("HAS_FILM_DIFFUSION", false);
	jpp.popScope();

	jpp.remove("COL_POROSITY");
	jpp.set("TOTAL_POROSITY", 0.37);

	checkJacobianAD(jpp, 1e-6);
}

TEST_CASE("CSTR rejects mixed equilibrium and film diffusion particle types", "[CSTR],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("HOMOGENEOUS_PARTICLE", "FV");

	const double volFrac[] = { 0.3, 0.7 };
	const double parFactor[] = { 0.9 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 2, parFactor, volFrac);

	jpp.pushScope("particle_type_001");
	jpp.set("HAS_FILM_DIFFUSION", false);
	jpp.popScope();

	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	REQUIRE_THROWS_WITH(mb->createUnitOperation(jpp, 0), Catch::Contains("must not be mixed"));

	destroyModelBuilder(mb);
}

TEST_CASE("CSTR rejects unbalanced flow rates", "[CSTR],[CI]")
{
	// A model system already balances the flow rates of a unit operation that cannot accumulate, so
	// the guard of the unit operation is checked directly
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	cadet::JsonParameterProvider jpp = createFiniteBath("HOMOGENEOUS_PARTICLE", "FV");
	cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);

	const cadet::active in = 0.1;
	const cadet::active out = 0.05;
	unit->setFlowRates(&in, &out);

	std::vector<double> y(unit->numDofs(), 0.0);
	const cadet::AdJacobianParams noParams{ nullptr, nullptr, 0u };

	REQUIRE_THROWS_WITH(unit->notifyDiscontinuousSectionTransition(0.0, 0u, { y.data(), nullptr }, noParams), Catch::Contains("constant liquid volume"));

	mb->destroyUnitOperation(unit);
	destroyModelBuilder(mb);
}

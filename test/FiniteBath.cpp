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
	cadet::IUnitOperation* createAndConfigureFiniteBath(cadet::IModelBuilder& mb, cadet::JsonParameterProvider& jpp, double flowRate)
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

		cadet::IUnitOperation* const unitAna = createAndConfigureFiniteBath(*mb, jpp, flowRate);
		cadet::IUnitOperation* const unitAD = createAndConfigureFiniteBath(*mb, jpp, flowRate);
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

		cadet::IUnitOperation* const unit = createAndConfigureFiniteBath(*mb, jpp, flowRate);

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

TEST_CASE("FiniteBath with homogeneous particles Jacobian vs AD", "[FiniteBath],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("HOMOGENEOUS_PARTICLE", "FV");
	checkJacobianAD(jpp, 1e-6);
}

TEST_CASE("FiniteBath with general rate particles and FV discretization Jacobian vs AD", "[FiniteBath],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");
	checkJacobianAD(jpp, 1e-6);
}

TEST_CASE("FiniteBath with general rate particles and DG discretization Jacobian vs AD", "[FiniteBath],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "DG");
	checkJacobianAD(jpp, 1e-6);
}

TEST_CASE("FiniteBath with two particle types Jacobian vs AD", "[FiniteBath],[UnitOp],[Residual],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");

	const double volFrac[] = { 0.3, 0.7 };
	const double parFactor[] = { 0.9 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 2, parFactor, volFrac);

	checkJacobianAD(jpp, 1e-6);
}

TEST_CASE("FiniteBath without flow Jacobian vs AD", "[FiniteBath],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("HOMOGENEOUS_PARTICLE", "FV");
	checkJacobianAD(jpp, 0.0);
}

TEST_CASE("FiniteBath with homogeneous particles time derivative Jacobian vs FD", "[FiniteBath],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("HOMOGENEOUS_PARTICLE", "FV");
	checkTimeDerivativeJacobianFD(jpp, 1e-6, 1e-6, 0.0, 1e-4);
}

TEST_CASE("FiniteBath with general rate particles time derivative Jacobian vs FD", "[FiniteBath],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");
	checkTimeDerivativeJacobianFD(jpp, 1e-6, 1e-6, 0.0, 1e-4);
}

TEST_CASE("FiniteBath inlet DOF Jacobian", "[FiniteBath],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	cadet::JsonParameterProvider jpp = createFiniteBath("GENERAL_RATE_PARTICLE", "FV");
	cadet::IUnitOperation* const unit = createAndConfigureFiniteBath(*mb, jpp, 1e-6);

	cadet::test::unitoperation::testInletDofJacobian(unit, false);

	mb->destroyUnitOperation(unit);
	destroyModelBuilder(mb);
}

TEST_CASE("FiniteBath consistent initialization with linear binding", "[FiniteBath],[ConsistentInit],[CI]")
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

				cadet::IUnitOperation* const unit = createAndConfigureFiniteBath(*mb, jpp, 1e-6);

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

TEST_CASE("FiniteBath consistent sensitivity initialization with linear binding", "[FiniteBath],[ConsistentInit],[Sensitivity],[CI]")
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

			cadet::IUnitOperation* const unit = createAndConfigureFiniteBath(*mb, jpp, 1e-6);

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

TEST_CASE("FiniteBath batch uptake vs analytic solution", "[FiniteBath],[Simulation],[Analytic],[CI]")
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

TEST_CASE("FiniteBath batch uptake numerical Benchmark for rapid equilibrium Langmuir case", "[FiniteBath],[DG],[Simulation],[Reference],[CI]")
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

TEST_CASE("FiniteBath without particles matches constant volume CSTR", "[FiniteBath],[Simulation],[CI]")
{
	// Without particles, a finite bath is a constant volume CSTR
	cadet::JsonParameterProvider jppBath = createCSTRBenchmark(1, 100.0, 1.0);
	{
		jppBath.pushScope("model");
		jppBath.pushScope("unit_000");

		jppBath.set("UNIT_TYPE", "FINITE_BATH");
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

TEST_CASE("FiniteBath rejects particles in rapid equilibrium", "[FiniteBath],[CI]")
{
	cadet::JsonParameterProvider jpp = createFiniteBath("HOMOGENEOUS_PARTICLE", "FV");

	jpp.pushScope("particle_type_000");
	jpp.set("HAS_FILM_DIFFUSION", false);
	jpp.popScope();

	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	// The builder selects the CSTR for particles in rapid equilibrium, which rejects the
	// particle transport fields of a finite bath
	cadet::IModel* const unit = mb->createUnitOperation(jpp, 0);
	REQUIRE(nullptr != unit);
	REQUIRE(std::string(unit->unitOperationName()) == "CSTR");

	mb->destroyUnitOperation(reinterpret_cast<cadet::IUnitOperation*>(unit));
	destroyModelBuilder(mb);
}

TEST_CASE("FiniteBath rejects mixed equilibrium and film diffusion particle types", "[FiniteBath],[ParticleType],[CI]")
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

TEST_CASE("FiniteBath rejects unbalanced flow rates", "[FiniteBath],[CI]")
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

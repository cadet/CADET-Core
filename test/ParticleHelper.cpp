// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2019: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>
#include "Approx.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "ParticleHelper.hpp"
#include "SimHelper.hpp"
#include "ColumnTests.hpp"

#include "Utils.hpp"
#include "common/Driver.hpp"
#include "SimulationTypes.hpp"
#include "ParallelSupport.hpp"

#include "JsonTestModels.hpp"
#include "ModelBuilderImpl.hpp"
#include "UnitOperation.hpp"
#include "JacobianHelper.hpp"

namespace
{
	template <typename T>
	void replicateData(std::vector<T>& data, unsigned int nTimes)
	{
		data.reserve(data.size() * nTimes);
		const typename std::vector<T>::iterator itEnd = data.end();
		for (unsigned int i = 0; i < nTimes - 1; ++i)
		{
			data.insert(data.end(), data.begin(), itEnd);
		}
	}

	void replicateFieldDataDouble(cadet::JsonParameterProvider& jpp, const std::string& field, const std::vector<double>& factors)
	{
		if (!jpp.exists(field))
			return;

		std::vector<double> data = jpp.getDoubleArray(field);
		data.reserve(data.size() * (factors.size() + 1));
		const unsigned int numElements = data.size();

		for (unsigned int i = 0; i < factors.size(); ++i)
		{
			for (unsigned int j = 0; j < numElements; ++j)
				data.push_back(data[j] * factors[i]);
		}

		jpp.set(field, data);
	}

	void replicateFieldDataDouble(cadet::JsonParameterProvider& jpp, const std::string& field, unsigned int nTimes)
	{
		if (!jpp.exists(field))
			return;

		std::vector<double> data = jpp.getDoubleArray(field);
		replicateData(data, nTimes);
		jpp.set(field, data);
	}

	void replicateFieldDataInt(cadet::JsonParameterProvider& jpp, const std::string& field, unsigned int nTimes)
	{
		if (!jpp.exists(field))
			return;

		std::vector<int> data = jpp.getIntArray(field);
		replicateData(data, nTimes);
		jpp.set(field, data);
	}

	void replicateFieldDataString(cadet::JsonParameterProvider& jpp, const std::string& field, unsigned int nTimes)
	{
		if (!jpp.exists(field))
			return;

		std::vector<std::string> data = jpp.getStringArray(field);
		replicateData(data, nTimes);
		jpp.set(field, data);
	}

	/**
	 * @brief Creates a runnable column model
	 * @details Creates a column model and configures it using the given IParameterProvider @p jpp.
	 * @param [in] uoType Unit operation type
	 * @param [in] mb ModelBuilder
	 * @param [in] jpp Configuration of the model
	 * @return Runnable column model
	 */
	inline cadet::IUnitOperation* createAndConfigureUnit(const std::string& uoType, cadet::IModelBuilder& mb, cadet::JsonParameterProvider& jpp)
	{
		// Create a unit
		cadet::IModel* const iUnit = mb.createUnitOperation(uoType, 0);
		REQUIRE(nullptr != iUnit);

		cadet::IUnitOperation* const unit = reinterpret_cast<cadet::IUnitOperation*>(iUnit);

		// Configure
		cadet::ModelBuilder& temp = *reinterpret_cast<cadet::ModelBuilder*>(&mb);
		REQUIRE(unit->configureModelDiscretization(jpp, temp));
		REQUIRE(unit->configure(jpp));

		return unit;
	}

	template <typename factor_t>
	void extendModelToManyParticleTypesImpl(cadet::JsonParameterProvider& jpp, const factor_t& factors, unsigned int nTypes, double const* const volFrac)
	{
		{
			auto ds = cadet::test::util::makeOptionalGroupScope(jpp, "discretization");

			replicateFieldDataInt(jpp, "NBOUND", nTypes);
			replicateFieldDataInt(jpp, "NPAR", nTypes);
			replicateFieldDataString(jpp, "PAR_DISC_TYPE", nTypes);
			replicateFieldDataDouble(jpp, "PAR_DISC_VECTOR", nTypes);
		}

		replicateFieldDataDouble(jpp, "FILM_DIFFUSION", factors);
		replicateFieldDataDouble(jpp, "PAR_DIFFUSION", factors);
		replicateFieldDataDouble(jpp, "PAR_SURFDIFFUSION", factors);
		replicateFieldDataDouble(jpp, "PAR_RADIUS", factors);
		replicateFieldDataDouble(jpp, "PAR_CORERADIUS", factors);
		replicateFieldDataDouble(jpp, "PAR_POROSITY", factors);
		replicateFieldDataDouble(jpp, "PORE_ACCESSIBILITY", factors);

		replicateFieldDataDouble(jpp, "INIT_CP", nTypes);
		replicateFieldDataDouble(jpp, "INIT_Q", nTypes);

		replicateFieldDataString(jpp, "ADSORPTION_MODEL", nTypes);

		// Move group "adsorption" to "adsorption_000"
		if (jpp.exists("adsorption"))
		{
			jpp.copy("adsorption", "adsorption_000");
			jpp.remove("adsorption");
		}

		// Replicate "adsorption_000"
		std::ostringstream ss;
		for (unsigned int i = 1; i < nTypes; ++i)
		{
			ss.str("");
			ss << "adsorption_" << std::setfill('0') << std::setw(3) << i;
			jpp.copy("adsorption_000", ss.str());
		}

		if (volFrac)
			jpp.set("PAR_TYPE_VOLFRAC", std::vector<double>(volFrac, volFrac + nTypes));
	}

	void scrambleParticleTypeFractionsSpatially(cadet::JsonParameterProvider& jpp, unsigned int nParType)
	{
		auto ms = cadet::test::util::makeModelGroupScope(jpp);

		unsigned int nCol = 0;
		unsigned int nRad = 1;
		{
			auto ds = cadet::test::util::makeOptionalGroupScope(jpp, "discretization");
			nCol = jpp.getInt("NCOL");

			if (jpp.exists("NRAD"))
				nRad = jpp.getInt("NRAD");
		}
		
		const double baseFrac[] = {0.2, 0.45, 0.35};
		std::vector<double> volFrac(nCol * nRad * nParType, 0.0);
		for (unsigned int i = 0; i < nCol * nRad; ++i)
		{
			volFrac[i * nParType + 0] = baseFrac[(i+0) % nParType];
			volFrac[i * nParType + 1] = baseFrac[(i+1) % nParType];
			volFrac[i * nParType + 2] = baseFrac[(i+2) % nParType];
		}
		cadet::test::particle::setParticleTypeVolumeFractions(jpp, volFrac);		
	}
}

namespace cadet
{

namespace test
{

namespace particle
{
	void extendModelToManyParticleTypes(cadet::JsonParameterProvider& jpp, unsigned int nTypes, double const* const volFrac)
	{
		extendModelToManyParticleTypesImpl(jpp, nTypes, nTypes, volFrac);
	}

	void extendModelToManyParticleTypes(cadet::JsonParameterProvider& jpp, UnitOpIdx unit, unsigned int nTypes, double const* const volFrac)
	{
		auto ms = util::makeModelGroupScope(jpp, unit);
		extendModelToManyParticleTypes(jpp, nTypes, volFrac);
	}

	void extendModelToManyParticleTypes(cadet::JsonParameterProvider& jpp, const std::vector<double>& factors, double const* const volFrac)
	{
		extendModelToManyParticleTypesImpl(jpp, factors, factors.size() + 1, volFrac);
	}

	void setParticleTypeVolumeFractions(cadet::JsonParameterProvider& jpp, UnitOpIdx unit, const std::vector<double>& volFrac)
	{
		auto ms = util::makeModelGroupScope(jpp, unit);
		setParticleTypeVolumeFractions(jpp, volFrac);
	}

	void setParticleTypeVolumeFractions(cadet::JsonParameterProvider& jpp, const std::vector<double>& volFrac)
	{
		jpp.set("PAR_TYPE_VOLFRAC", volFrac);
	}

	void testOneVsTwoIdenticalParticleTypes(const char* uoType, double absTol, double relTol)
	{
		// Use Load-Wash-Elution test case
		cadet::JsonParameterProvider jpp = createLWE(uoType);
		testOneVsTwoIdenticalParticleTypes(jpp, absTol, relTol);
	}

	void testOneVsTwoIdenticalParticleTypes(cadet::JsonParameterProvider& jpp, double absTol, double relTol)
	{
		// Simulate without additional particle types
		cadet::Driver drv;
		drv.configure(jpp);
		drv.run();

		// Extend to two particle types
		const double volFrac[] = {1.0, 0.0};
		extendModelToManyParticleTypes(jpp, 0, 2, volFrac);

		// Simulate with first particle type only
		cadet::Driver drvP1;
		drvP1.configure(jpp);
		drvP1.run();

		cadet::InternalStorageUnitOpRecorder const* const data = drv.solution()->unitOperation(0);
		cadet::InternalStorageUnitOpRecorder const* const p1Data = drvP1.solution()->unitOperation(0);

		double const* outlet = data->outlet();
		double const* p1Outlet = p1Data->outlet();

		const unsigned int nComp = p1Data->numComponents();
		for (unsigned int i = 0; i < p1Data->numDataPoints() * p1Data->numInletPorts() * nComp; ++i, ++outlet, ++p1Outlet)
		{
			CAPTURE(i);
			CHECK((*outlet) == makeApprox(*p1Outlet, relTol, absTol));
		}
	}

	void testSeparateIdenticalParticleTypes(const char* uoType, double absTol, double relTol)
	{
		// Use Load-Wash-Elution test case
		cadet::JsonParameterProvider jpp = createLWE(uoType);
		testSeparateIdenticalParticleTypes(jpp, absTol, relTol);
	}

	void testSeparateIdenticalParticleTypes(cadet::JsonParameterProvider& jpp, double absTol, double relTol)
	{
		// Extend to two particle types
		const double volFrac[] = {1.0, 0.0};
		extendModelToManyParticleTypes(jpp, 0, 2, volFrac);

		// Simulate with first particle type only
		cadet::Driver drvP1;
		drvP1.configure(jpp);
		drvP1.run();

		// Simulate with second particle type only
		setParticleTypeVolumeFractions(jpp, 0, {0.0, 1.0});
		cadet::Driver drvP2;
		drvP2.configure(jpp);
		drvP2.run();

		cadet::InternalStorageUnitOpRecorder const* const p1Data = drvP1.solution()->unitOperation(0);
		cadet::InternalStorageUnitOpRecorder const* const p2Data = drvP2.solution()->unitOperation(0);

		double const* p1Outlet = p1Data->outlet();
		double const* p2Outlet = p2Data->outlet();

		// Check outlet
		const unsigned int nComp = p1Data->numComponents();
		for (unsigned int i = 0; i < p1Data->numDataPoints() * p1Data->numInletPorts() * nComp; ++i, ++p1Outlet, ++p2Outlet)
		{
			CAPTURE(i);
			CHECK((*p1Outlet) == makeApprox(*p2Outlet, relTol, absTol));
		}

		// Check volume if available
		if (p1Data->solutionConfig().storeVolume)
		{
			double const* p1Vol = p1Data->volume();
			double const* p2Vol = p2Data->volume();
			for (unsigned int i = 0; i < p1Data->numDataPoints(); ++i, ++p1Vol, ++p2Vol)
			{
				CAPTURE(i);
				CHECK((*p1Vol) == makeApprox(*p2Vol, relTol, absTol));
			}
		}
	}

	template <typename modifier_t>
	void testLinearMixedParticleTypesImpl(cadet::JsonParameterProvider& jpp, double absTol, double relTol, modifier_t modify)
	{
		for (int i = 0; i < 2; ++i)
		{
			const bool dynamicBinding = (i == 0);
			SECTION(std::string(dynamicBinding ? " Kinetic binding" : " Quasi-stationary binding"))
			{
				setBindingMode(jpp, dynamicBinding);

				// Simulate without additional particle types
				cadet::Driver drv;
				drv.configure(jpp);
				drv.run();

				// Extend to multiple particle types (such that we have a total of 3 types)
				const double volFrac[] = {0.3, 0.6, 0.1};
				extendModelToManyParticleTypes(jpp, 0, 3, volFrac);

				modify(jpp);

				// Simulate with mixed particle types
				cadet::Driver drvMP;
				drvMP.configure(jpp);
				drvMP.run();

				cadet::InternalStorageUnitOpRecorder const* const data = drv.solution()->unitOperation(0);
				cadet::InternalStorageUnitOpRecorder const* const mpData = drvMP.solution()->unitOperation(0);

				double const* outlet = data->outlet();
				double const* mpOutlet = mpData->outlet();

				// Check outlet
				const unsigned int nComp = mpData->numComponents();
				for (unsigned int i = 0; i < mpData->numDataPoints() * mpData->numInletPorts() * nComp; ++i, ++outlet, ++mpOutlet)
				{
					CAPTURE(i);
					CHECK((*outlet) == makeApprox(*mpOutlet, relTol, absTol));
				}

				// Check volume if available
				if (data->solutionConfig().storeVolume)
				{
					double const* vol = data->volume();
					double const* mpVol = mpData->volume();
					for (unsigned int i = 0; i < data->numDataPoints(); ++i, ++vol, ++mpVol)
					{
						CAPTURE(i);
						CHECK((*vol) == makeApprox(*mpVol, relTol, absTol));
					}
				}
			}
		}
	}

	void testLinearMixedParticleTypes(const char* uoType, double absTol, double relTol)
	{
		cadet::JsonParameterProvider jpp = createPulseInjectionColumn(uoType, true);
		testLinearMixedParticleTypesImpl(jpp, absTol, relTol, [](cadet::JsonParameterProvider& jpp) { });
	}

	void testLinearMixedParticleTypes(cadet::JsonParameterProvider& jpp, double absTol, double relTol)
	{
		testLinearMixedParticleTypesImpl(jpp, absTol, relTol, [](cadet::JsonParameterProvider& jpp) { });
	}

	void testLinearSpatiallyMixedParticleTypes(const char* uoType, double absTol, double relTol)
	{
		cadet::JsonParameterProvider jpp = createPulseInjectionColumn(uoType, true);
		testLinearMixedParticleTypesImpl(jpp, absTol, relTol, [](cadet::JsonParameterProvider& jpp) { scrambleParticleTypeFractionsSpatially(jpp, 3); });
	}

	void testJacobianMixedParticleTypesImpl(cadet::JsonParameterProvider& jpp)
	{
		cadet::IModelBuilder* const mb = cadet::createModelBuilder();
		REQUIRE(nullptr != mb);

		cadet::IUnitOperation* const unitAna = createAndConfigureUnit(jpp.getString("UNIT_TYPE"), *mb, jpp);
		cadet::IUnitOperation* const unitAD = createAndConfigureUnit(jpp.getString("UNIT_TYPE"), *mb, jpp);

		// Enable AD
		cadet::ad::setDirections(cadet::ad::getMaxDirections());
		unitAD->useAnalyticJacobian(false);

		cadet::active* adRes = new cadet::active[unitAD->numDofs()];
		cadet::active* adY = new cadet::active[unitAD->numDofs()];

		const AdJacobianParams noParams{nullptr, nullptr, 0u};
		const AdJacobianParams adParams{adRes, adY, 0u};

		unitAD->prepareADvectors(adParams);

		// Setup matrices
		unitAna->notifyDiscontinuousSectionTransition(0.0, 0u, noParams);
		unitAD->notifyDiscontinuousSectionTransition(0.0, 0u, adParams);

		// Obtain memory for state, Jacobian multiply direction, Jacobian column
		std::vector<double> y(unitAD->numDofs(), 0.0);
		std::vector<double> jacDir(unitAD->numDofs(), 0.0);
		std::vector<double> jacCol1(unitAD->numDofs(), 0.0);
		std::vector<double> jacCol2(unitAD->numDofs(), 0.0);
		cadet::util::ThreadLocalStorage<double> tls;
		tls.resize(unitAna->threadLocalMemorySize());

		// Fill state vector with some values
		util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, unitAna->numDofs());
//		util::populate(y.data(), [](unsigned int idx) { return 1.0; }, unitAna->numDofs());

		// Compute state Jacobian
		unitAna->residualWithJacobian(ActiveSimulationTime{0.0, 0u, 1.0}, ConstSimulationState{y.data(), nullptr}, jacDir.data(), noParams, tls);
		unitAD->residualWithJacobian(ActiveSimulationTime{0.0, 0u, 1.0}, ConstSimulationState{y.data(), nullptr}, jacDir.data(), adParams, tls);
		std::fill(jacDir.begin(), jacDir.end(), 0.0);

		// Compare Jacobians
		cadet::test::checkJacobianPatternFD(unitAna, unitAD, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data(), tls);
		cadet::test::checkJacobianPatternFD(unitAna, unitAna, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data(), tls);
		cadet::test::compareJacobian(unitAna, unitAD, nullptr, nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
//		cadet::test::compareJacobianFD(unitAna, unitAD, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());

		delete[] adRes;
		delete[] adY;
		mb->destroyUnitOperation(unitAna);
		mb->destroyUnitOperation(unitAD);
		destroyModelBuilder(mb);
	}

	void testJacobianMixedParticleTypes(const std::string& uoType)
	{
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding(uoType);
		testJacobianMixedParticleTypes(jpp);
	}

	void testJacobianMixedParticleTypes(cadet::JsonParameterProvider& jpp)
	{
		// Add more particle types (such that we have a total of 3 types)
		const double volFrac[] = {0.3, 0.6, 0.1};
		extendModelToManyParticleTypes(jpp, {0.9, 0.8}, volFrac);

		testJacobianMixedParticleTypesImpl(jpp);
	}

	void testJacobianSpatiallyMixedParticleTypes(const std::string& uoType)
	{
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding(uoType);

		// Add more particle types (such that we have a total of 3 types)
		extendModelToManyParticleTypes(jpp, {0.9, 0.8}, nullptr);

		// Spatially inhomogeneous
		scrambleParticleTypeFractionsSpatially(jpp, 3);

		testJacobianMixedParticleTypesImpl(jpp);
	}

	void testTimeDerivativeJacobianMixedParticleTypesFD(const std::string& uoType, double h, double absTol, double relTol)
	{
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding(uoType);
		testTimeDerivativeJacobianMixedParticleTypesFD(jpp, h, absTol, relTol);
	}

	void testTimeDerivativeJacobianMixedParticleTypesFD(cadet::JsonParameterProvider& jpp, double h, double absTol, double relTol)
	{
		cadet::IModelBuilder* const mb = cadet::createModelBuilder();
		REQUIRE(nullptr != mb);

		// Add more particle types (such that we have a total of 3 types)
		const double volFrac[] = {0.3, 0.6, 0.1};
		extendModelToManyParticleTypes(jpp, {0.9, 0.8}, volFrac);

		cadet::IUnitOperation* const unit = createAndConfigureUnit(jpp.getString("UNIT_TYPE"), *mb, jpp);

		// Setup matrices
		unit->notifyDiscontinuousSectionTransition(0.0, 0u, AdJacobianParams{nullptr, nullptr, 0u});

		// Obtain memory for state, Jacobian multiply direction, Jacobian column
		const unsigned int nDof = unit->numDofs();
		std::vector<double> y(nDof, 0.0);
		std::vector<double> yDot(nDof, 0.0);
		std::vector<double> jacDir(nDof, 0.0);
		std::vector<double> jacCol1(nDof, 0.0);
		std::vector<double> jacCol2(nDof, 0.0);
		cadet::util::ThreadLocalStorage<double> tls;
		tls.resize(unit->threadLocalMemorySize());

		// Fill state vectors with some values
		util::populate(y.data(), [=](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, nDof);
		util::populate(yDot.data(), [=](unsigned int idx) { return std::abs(std::sin((idx + nDof) * 0.13)) + 1e-4; }, nDof);

		// Compare Jacobians
		cadet::test::compareTimeDerivativeJacobianFD(unit, unit, y.data(), yDot.data(), jacDir.data(), jacCol1.data(), jacCol2.data(), tls, h, absTol, relTol);

		mb->destroyUnitOperation(unit);
		destroyModelBuilder(mb);
	}

	void testArrowHeadJacobianSpatiallyMixedParticleTypes(const std::string& uoType, double h, double absTol, double relTol)
	{
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding(uoType);

		// Add more particle types (such that we have a total of 3 types)
		extendModelToManyParticleTypes(jpp, {0.9, 0.8}, nullptr);

		// Spatially inhomogeneous
		scrambleParticleTypeFractionsSpatially(jpp, 3);

		cadet::IModelBuilder* const mb = cadet::createModelBuilder();
		REQUIRE(nullptr != mb);

		cadet::IUnitOperation* const unitAna = createAndConfigureUnit(uoType, *mb, jpp);
		cadet::IUnitOperation* const unitFD = createAndConfigureUnit(uoType, *mb, jpp);

		// Obtain offset to fluxes
		const unsigned int fluxOffset = column::fluxOffsetOfColumnUnitOp(unitFD);

		// Setup matrices
		const AdJacobianParams noParams{nullptr, nullptr, 0u};
		unitAna->notifyDiscontinuousSectionTransition(0.0, 0u, noParams);
		unitFD->notifyDiscontinuousSectionTransition(0.0, 0u, noParams);

		// Obtain memory for state, Jacobian multiply direction, Jacobian column
		std::vector<double> y(unitFD->numDofs(), 0.0);
		std::vector<double> jacDir(unitFD->numDofs(), 0.0);
		std::vector<double> jacCol1(unitFD->numDofs(), 0.0);
		std::vector<double> jacCol2(unitFD->numDofs(), 0.0);
		cadet::util::ThreadLocalStorage<double> tls;
		tls.resize(unitAna->threadLocalMemorySize());

		// Fill state vector with some values
		util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, unitAna->numDofs());
//		util::populate(y.data(), [](unsigned int idx) { return 1.0; }, unitAna->numDofs());

		// Compute state Jacobian
		unitAna->residualWithJacobian(ActiveSimulationTime{0.0, 0u, 1.0}, ConstSimulationState{y.data(), nullptr}, jacDir.data(), noParams, tls);
		unitFD->residualWithJacobian(ActiveSimulationTime{0.0, 0u, 1.0}, ConstSimulationState{y.data(), nullptr}, jacDir.data(), noParams, tls);
		std::fill(jacDir.begin(), jacDir.end(), 0.0);

		// Compare Jacobians
		cadet::test::compareJacobianArrowHeadFD(
			[=, &tls](double const* lDir, double* res) -> void { unitFD->residual(SimulationTime{0.0, 0u, 1.0}, ConstSimulationState{lDir, nullptr}, res, tls); }, 
			[&](double const* lDir, double* res) -> void { unitAna->multiplyWithJacobian(SimulationTime{0.0, 0u, 1.0}, ConstSimulationState{y.data(), nullptr}, lDir, 1.0, 0.0, res); }, 
			y.data(), jacDir.data(), jacCol1.data(), jacCol2.data(), unitFD->numDofs(), fluxOffset, h, absTol, relTol);

		mb->destroyUnitOperation(unitAna);
		mb->destroyUnitOperation(unitFD);
		destroyModelBuilder(mb);
	}

} // namespace particle
} // namespace test
} // namespace cadet

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

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "AutoDiff.hpp"
#include "linalg/Norms.hpp"
#include "UnitOperation.hpp"
#include "UnitOperationTests.hpp"
#include "SimulationTypes.hpp"
#include "ParallelSupport.hpp"

#include "Utils.hpp"

#include <vector>

namespace cadet
{

namespace test
{

namespace unitoperation
{

	void testConsistentInitialization(cadet::IUnitOperation* const unit, bool adEnabled, double* const y, double consTol, double absTol)
	{
		cadet::active* adRes = nullptr;
		cadet::active* adY = nullptr;

		// Enable AD
		if (adEnabled)
		{
			cadet::ad::setDirections(cadet::ad::getMaxDirections());
			unit->useAnalyticJacobian(false);

			adRes = new cadet::active[unit->numDofs()];
			adY = new cadet::active[unit->numDofs()];
			unit->prepareADvectors(AdJacobianParams{adRes, adY, 0});
		}
		else
		{
			unit->useAnalyticJacobian(true);
		}

		const AdJacobianParams adParams{adRes, adY, 0};

		// Setup matrices
		unit->notifyDiscontinuousSectionTransition(0.0, 0u, adParams);

		// Obtain memory
		std::vector<double> yDot(unit->numDofs(), 0.0);
		std::vector<double> res(unit->numDofs(), 0.0);
		cadet::util::ThreadLocalStorage tls;
		tls.resize(unit->threadLocalMemorySize());

		// Initialize algebraic variables in state vector
		unit->consistentInitialState(SimulationTime{0.0, 0u, 1.0}, y, adParams, consTol, tls);

		// Evaluate residual without yDot and update Jacobians
		// Store residual in yDot
		unit->residualWithJacobian(ActiveSimulationTime{0.0, 0u, 1.0}, ConstSimulationState{y, nullptr}, yDot.data(), adParams, tls);
		CAPTURE(yDot);

		// Initialize time derivative of state
		unit->consistentInitialTimeDerivative(SimulationTime{0.0, 0u, 1.0}, y, yDot.data(), tls);

		// Check norm of residual (but skip over connection DOFs at the beginning of the local state vector slice)
		unit->residual(SimulationTime{0.0, 0u, 1.0}, ConstSimulationState{y, yDot.data()}, res.data(), tls);
		INFO("Residual " << linalg::linfNorm(res.data() + unit->numComponents() * unit->numInletPorts(), unit->numDofs() - unit->numComponents() * unit->numInletPorts()));
		CAPTURE(res);
		CHECK(linalg::linfNorm(res.data() + unit->numComponents() * unit->numInletPorts(), unit->numDofs() - unit->numComponents() * unit->numInletPorts()) <= absTol);

		// Clean up
		delete[] adRes;
		delete[] adY;
	}

	void testConsistentInitializationSensitivity(cadet::IUnitOperation* const unit, bool adEnabled, double const* const y, double const* const yDot, double absTol)
	{
		const unsigned int nSens = unit->numSensParams();
		REQUIRE(nSens > 0);

		cadet::active* adRes = new cadet::active[unit->numDofs()];
		cadet::active* adY = nullptr;

		// Enable AD
		if (adEnabled)
		{
			cadet::ad::setDirections(cadet::ad::getMaxDirections());
			unit->useAnalyticJacobian(false);

			adY = new cadet::active[unit->numDofs()];
			unit->prepareADvectors(AdJacobianParams{adRes, adY, nSens});
		}
		else
		{
			unit->useAnalyticJacobian(true);
		}

		const AdJacobianParams adParams{adRes, adY, nSens};

		cadet::util::ThreadLocalStorage tls;
		tls.resize(unit->threadLocalMemorySize());

		// Setup matrices
		unit->notifyDiscontinuousSectionTransition(0.0, 0u, adParams);

		// Calculate dres / dp and Jacobians
		const ActiveSimulationTime simTime{0.0, 0u, 1.0};
		unit->residualSensFwdWithJacobian(simTime, ConstSimulationState{y, yDot}, adParams, tls);

		// Obtain memory
		std::vector<double> sensY(unit->numDofs() * nSens, 0.0);
		std::vector<double> sensYdot(unit->numDofs() * nSens, 0.0);
		std::vector<double> res(unit->numDofs() * nSens, 0.0);

		std::vector<double*> vecSensY(nSens, nullptr);
		std::vector<double*> vecSensYdot(nSens, nullptr);
		std::vector<double const*> cVecSensY(nSens, nullptr);
		std::vector<double const*> cVecSensYdot(nSens, nullptr);
		std::vector<double*> vecRes(nSens, nullptr);

		for (unsigned int i = 0; i < nSens; ++i)
		{
			vecSensY[i] = sensY.data() + i * unit->numDofs();
			vecSensYdot[i] = sensYdot.data() + i * unit->numDofs();
			cVecSensY[i] = sensY.data() + i * unit->numDofs();
			cVecSensYdot[i] = sensYdot.data() + i * unit->numDofs();
			vecRes[i] = res.data() + i * unit->numDofs();
		}

		std::vector<double> tmp1(unit->numDofs(), 0.0);
		std::vector<double> tmp2(unit->numDofs(), 0.0);
		std::vector<double> tmp3(unit->numDofs(), 0.0);

		// Fill sensY and sensYdot with some values
		util::populate(sensY.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, sensY.size());
		util::populate(sensYdot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, sensYdot.size());

		// Initialize sensitivity state vectors
		unit->consistentInitialSensitivity(simTime, ConstSimulationState{y, yDot}, vecSensY, vecSensYdot, adRes,tls);

		// Check norm of residual (but skip over connection DOFs at the beginning of the local state vector slice)
		unit->residualSensFwdAdOnly(simTime, ConstSimulationState{y, yDot}, adRes, tls);
		unit->residualSensFwdCombine(simTime, ConstSimulationState{y, yDot}, cVecSensY, cVecSensYdot, vecRes, adRes, tmp1.data(), tmp2.data(), tmp3.data());

		for (unsigned int i = 0; i < nSens; ++i)
		{
			INFO("Residual " << i << ": " << linalg::linfNorm(vecRes[i] + unit->numComponents() * unit->numInletPorts(), unit->numDofs() - unit->numComponents() * unit->numInletPorts()));
			CHECK(linalg::linfNorm(vecRes[i] + unit->numComponents() * unit->numInletPorts(), unit->numDofs() - unit->numComponents() * unit->numInletPorts()) <= absTol);
		}

		// Clean up
		delete[] adRes;
		delete[] adY;
	}

	void testInletDofJacobian(cadet::IUnitOperation* const unit, bool adEnabled)
	{
		cadet::active* adRes = nullptr;
		cadet::active* adY = nullptr;

		// Enable AD
		if (adEnabled)
		{
			cadet::ad::setDirections(cadet::ad::getMaxDirections());
			unit->useAnalyticJacobian(false);

			adRes = new cadet::active[unit->numDofs()];
			adY = new cadet::active[unit->numDofs()];
			unit->prepareADvectors(AdJacobianParams{adRes, adY, 0});
		}
		else
		{
			unit->useAnalyticJacobian(true);
		}

		const AdJacobianParams adParams{adRes, adY, 0};

		// Setup matrices
		unit->notifyDiscontinuousSectionTransition(0.0, 0u, adParams);

		// Obtain memory
		std::vector<double> y(unit->numDofs(), 0.0);
		std::vector<double> jac(unit->numDofs(), 0.0);
		cadet::util::ThreadLocalStorage tls;
		tls.resize(unit->threadLocalMemorySize());

		// Assemble Jacobian
		const ActiveSimulationTime simTime{0.0, 0u, 1.0};
		const ConstSimulationState simState{y.data(), nullptr};
		unit->residualWithJacobian(simTime, simState, jac.data(), adParams, tls);

		// Check ports and components
		const unsigned int inletBlockSize = unit->numInletPorts() * unit->numComponents();
		CHECK(inletBlockSize == unit->numDofs() - unit->numPureDofs());

		for (unsigned int port = 0; port < unit->numInletPorts(); ++port)
		{
			const unsigned int inletOffset = unit->localInletComponentIndex(port);
			const unsigned int inletStride = unit->localInletComponentStride(port);

			for (unsigned int comp = 0; comp < unit->numComponents(); ++comp)
			{
				// Get Jacobian column
				const unsigned int curItem = inletOffset + comp * inletStride;

				y[curItem] = 1.0;
				unit->multiplyWithJacobian(SimulationTime{0.0, 0u, 1.0}, simState, y.data(), 1.0, 0.0, jac.data());
				y[curItem] = 0.0;

				// Check inlet block of Jacobian column
				for (unsigned int i = 0; i < inletBlockSize; ++i)
				{
					CAPTURE(port);
					CAPTURE(comp);
					CAPTURE(curItem);
					if (i == curItem)
						CHECK(jac[i] == 1.0);
					else
						CHECK(jac[i] == 0.0);
				}
			}
		}

		// Clean up
		delete[] adRes;
		delete[] adY;
	}

} // namespace unitoperation
} // namespace test
} // namespace cadet

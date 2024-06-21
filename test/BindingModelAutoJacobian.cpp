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

#include <catch.hpp>
#include "Approx.hpp"
#include "cadet/cadet.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "model/binding/BindingModelBase.hpp"

#include "JacobianHelper.hpp"
#include "Dummies.hpp"

#include "linalg/DenseMatrix.hpp"
#include "AdUtils.hpp"
#include "AutoDiff.hpp"

namespace cadet
{
	class BindingWithoutJacobian : public model::BindingModelBase
	{
	public:
		BindingWithoutJacobian() { }

		virtual const char* name() const CADET_NOEXCEPT { return "BindingWithoutJacobian"; }
		virtual bool dependsOnTime() const CADET_NOEXCEPT { return false; }

		virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx) override { return true; }
		virtual bool hasSalt() const CADET_NOEXCEPT { return false; }
		virtual bool supportsMultistate() const CADET_NOEXCEPT { return true; }
		virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }

		CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE
	protected:
		virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return false; }

		template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
		int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
			CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
		{
			int bndIdx = 0;
			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				// Residual
				for (int j = 0; j < _nBoundStates[i]; ++j)
				{
					res[bndIdx] = (i+2+j) * y[bndIdx] - (i+1) * yCp[i];

					// Next bound state
					++bndIdx;
				}
			}

			return 0;
		}
	};

	class BindingWithJacobian : public model::BindingModelBase
	{
	public:
		BindingWithJacobian() { }

		virtual const char* name() const CADET_NOEXCEPT { return "BindingWithJacobian"; }
		virtual bool dependsOnTime() const CADET_NOEXCEPT { return false; }

		virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx) override { return true; }
		virtual bool hasSalt() const CADET_NOEXCEPT { return false; }
		virtual bool supportsMultistate() const CADET_NOEXCEPT { return true; }
		virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }

		CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE
	protected:
		virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

		template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
		int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
			CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
		{
			return 0;
		}
	};

} // namespace cadet

TEST_CASE("Automatic AD binding model Jacobian vs FD", "[BindingModel],[Jacobian]")
{
	cadet::BindingWithoutJacobian bm;
	const int nComp = 4;
	const int totalBoundStates = 5;
	const int numDofs = nComp + totalBoundStates;
	const unsigned int boundStates[] = {1, 0, 3, 1};
	const unsigned int boundOffset[] = {0, 1, 1, 4};

	const double yState[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, -8.0, -9.0};

	DummyParameterProvider pp;
	bm.configureModelDiscretization(pp, nComp, boundStates, nullptr);

	REQUIRE(bm.requiresWorkspace());
	REQUIRE(bm.requiredADdirs() == nComp + totalBoundStates);

	const unsigned int workspaceSize = bm.workspaceSize(4, totalBoundStates, boundOffset);
	std::vector<char> buffer(workspaceSize, 0);
	cadet::LinearBufferAllocator workSpace(buffer.data(), buffer.data() + workspaceSize);

	// Enable AD
	cadet::ad::setDirections(cadet::ad::getMaxDirections());

	cadet::linalg::DenseMatrix jacAna;
	jacAna.resize(numDofs, numDofs);
	bm.analyticJacobian(1.0, 0u, cadet::ColumnPosition{0.0, 0.0, 0.0}, yState, nComp, jacAna.row(nComp), workSpace);

	std::vector<double> y(numDofs, 0.0);
	std::vector<double> dir(numDofs, 0.0);
	std::vector<double> colA(totalBoundStates, 0.0);
	std::vector<double> colB(totalBoundStates, 0.0);

	// Compare against FD: Since our flux() is linear, the Jacobian should be exact
	cadet::test::compareJacobianFD(
		[=](double const* y, double* r) { std::fill_n(r, totalBoundStates, 0.0); bm.flux(0.0, 0u, cadet::ColumnPosition{0.0, 0.0, 0.0}, y + nComp, y, r, workSpace); },
		[&](double const* y, double* r) { jacAna.submatrixMultiplyVector(y, nComp, 0, totalBoundStates, numDofs, r); },
		y.data(), dir.data(), colA.data(), colB.data(), numDofs, totalBoundStates, 1e-7, 0.0, 1e-15
	);
}

TEST_CASE("Automatic AD disabled for binding model with Jacobian", "[BindingModel],[Jacobian]")
{
	cadet::BindingWithJacobian bm;
	const int nComp = 4;
	const int totalBoundStates = 5;
	const int numDofs = nComp + totalBoundStates;
	const unsigned int boundStates[] = {1, 0, 3, 1};
	const unsigned int boundOffset[] = {0, 1, 1, 4};

	const double yState[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, -8.0, -9.0};

	DummyParameterProvider pp;
	bm.configureModelDiscretization(pp, nComp, boundStates, nullptr);

	REQUIRE(!bm.requiresWorkspace());
	REQUIRE(bm.requiredADdirs() == 0);

	REQUIRE(bm.workspaceSize(4, totalBoundStates, boundOffset) == 0);
}

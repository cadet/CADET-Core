// =============================================================================
//  CADET
//  
//  Copyright © The CADET Authors
//            Please see the CONTRIBUTORS.md file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>
#include <cmath>
#include <vector>
#include <iostream>
#include "Approx.hpp"
#include "cadet/cadet.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"
#include "LoggingUtils.hpp"

//#include "common/JsonParameterProvider.hpp"
//#include "Dummies.hpp"
#include "linalg/DenseMatrix.hpp"
#include "AdUtils.hpp"
#include "AutoDiff.hpp"
#include "nonlin/PseudoTransientContinuation.hpp"

namespace
{
	struct PTCIterateOutputPolicy
	{
		inline static void iteration(unsigned int idxIter, double const* const x, double const* const fx, double residualNorm, double stepSize, unsigned int size)
        {
            std::cout << "Iter " << idxIter << " Step " << stepSize << " residual " << residualNorm << std::endl;
        }
	};

	template <typename T>
	bool residualProb1(T const* const in, T* const out)
	{
		using std::sin;
		out[0] = in[0] * in[1] - 4.0;
		out[1] = sin(in[1] + in[0]);
        return true;
	}
}

TEST_CASE("Pseudo-Transient Continuation Problem 1", "[NonlinearSolver],[PTC],[CI]")
{
	const int problemSize = 2;
	const int maxIter = 5000;
	const double resTol = 1e-8;
    const int numNonMonotone = 100;
    const double initStep = 1e-2;

	double x[] = {1.0, 3.0};

	// Initialize AD and allocate AD vectors
	cadet::ad::setDirections(problemSize);

	std::vector<cadet::active> adIn(problemSize, 0.0);
	std::vector<cadet::active> adOut(problemSize, 0.0);

	// Set seed vectors
	cadet::ad::prepareAdVectorSeedsForDenseMatrix(adIn.data(), 0, problemSize);
	cadet::ad::fillAd(adIn.data(), problemSize, 0.0);

	const auto jacobian = [&](double const* const point, cadet::linalg::detail::DenseMatrixBase& jac) -> bool
	{
		cadet::ad::copyToAd(point, adIn.data(), problemSize);
		residualProb1(adIn.data(), adOut.data());
		cadet::ad::extractDenseJacobianFromAd(adOut.data(), 0, jac);
        return true;
	};

	std::vector<double> mem(4 * problemSize, 0.0);
	cadet::linalg::DenseMatrix jacMat;
	jacMat.resize(problemSize, problemSize);

	const bool res = cadet::nonlin::pseudoTransientContinuation<PTCIterateOutputPolicy>(&residualProb1<double>, jacobian, maxIter, resTol, numNonMonotone, initStep, nullptr, false, x, jacMat, mem.data(), problemSize);
	CHECK(res);

	std::vector<double> residual(problemSize, 0.0);
	residualProb1(x, residual.data());

	CHECK(cadet::linalg::l2Norm(residual.data(), problemSize) <= resTol);
    std::cout << x[0] << " " << x[1] << std::endl;
    std::cout << residual[0] << " " << residual[1] << std::endl;
}

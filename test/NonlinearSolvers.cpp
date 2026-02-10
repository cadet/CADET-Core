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
#include <random>
#include <chrono>
#include <iostream>

#include "Approx.hpp"
#include "cadet/cadet.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"
#include "LoggingUtils.hpp"

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
            std::cout << "Iter " << idxIter << " Timestep " << stepSize << " residual " << residualNorm << std::endl;
			std::cout << "  Residual ";
			for (unsigned int i = 0; i < size; ++i)
			{
				std::cout << fx[i] << " ";
			}
			std::cout << std::endl;
        }
	};

	struct Problem1
	{
		static constexpr int problemSize = 2;
		static constexpr int maxIter = 100;
		static constexpr double resTol = 1e-8;
		static constexpr int numNonMonotone = 10;
		static constexpr double initStep = 20.0;

		template <typename T>
		static bool residual(T const* const in, T* const out)
		{
			using std::sin;
			out[0] = in[0] * in[1] - 4.0;
			out[1] = sin(in[1] + in[0]);
			return true;
		}

		static void checkInvariants(double const* start, double const* end) { }

		static void initialPoint(double* init)
		{
			init[0] = 1.0;
			init[1] = 3.0;
		}

		template <typename rand_eng_t>
		static void randomInitialPoint(double* init, rand_eng_t& generator)
		{
			std::normal_distribution<double> distribution(0.0, 10.0);
			init[0] = distribution(generator);
			init[1] = distribution(generator);
		}

		static void scaleFromInitialPoint(double* scale, double const* init)
		{
			scale[0] = 0.5;
			scale[1] = 2.0;
		}

		static bool shouldRunTest(bool applyScaling, bool useVariant) { return true; }
	};

	struct SingleThreeCompReaction
	{
		static constexpr int problemSize = 3;
		static constexpr int maxIter = 100;
		static constexpr double resTol = 1e-8;
		static constexpr int numNonMonotone = 10;
		static constexpr double initStep = 20.0;

		template <typename T>
		static bool residual(T const* const in, T* const out)
		{
			// Reaction A + B <-> C, rate function fwd * a * b - bwd * c
			const double fwd = 2.3;
			const double bwd = 1.1;

			const T rate = fwd * in[0] * in[1] - bwd * in[2];
			out[0] = -rate;
			out[1] = -rate;
			out[2] = rate;
			return true;
		}

		static void checkInvariants(double const* start, double const* end)
		{
			const double invariant1Start = start[0] + start[2];
			const double invariant2Start = start[1] + start[2];

			const double invariant1End = end[0] + end[2];
			const double invariant2End = end[1] + end[2];

			CHECK(std::abs(invariant1End - invariant1Start) <= resTol);
			CHECK(std::abs(invariant2End - invariant2Start) <= resTol);
		}

		static void initialPoint(double* init)
		{
			init[0] = 1.0;
			init[1] = 1.0;
			init[2] = 0.0;
		}

		template <typename rand_eng_t>
		static void randomInitialPoint(double* init, rand_eng_t& generator)
		{
			std::normal_distribution<double> distribution(0.0, 10.0);
			init[0] = std::abs(distribution(generator));
			init[1] = std::abs(distribution(generator));
			init[2] = std::abs(distribution(generator));
		}

		static void scaleFromInitialPoint(double* scale, double const* init)
		{
			scale[0] = 0.5;
			scale[1] = 1.0;
			scale[2] = 2.0;
		}

		static bool shouldRunTest(bool applyScaling, bool useVariant) { return true; }
	};

	struct ThreeCompMassActionLawReaction
	{
		static constexpr int problemSize = 3;
		static constexpr int maxIter = 100;
		static constexpr double resTol = 1e-8;
		static constexpr int numNonMonotone = 10;
		static constexpr double initStep = 20.0;

		template <typename T>
		static bool residual(T const* const in, T* const out)
		{
			// Reaction A + B <-> C, rate function fwd * a * b - bwd * c
			// Raction A <-> B, rate function fwd * a - bwd * b
			const double fwd1 = 2.3;
			const double bwd1 = 1.1;

			const double fwd2 = 1e-2;
			const double bwd2 = 2.0;

			const T rate1 = fwd1 * in[0] * in[1] - bwd1 * in[2];
			const T rate2 = fwd2 * in[0] - bwd2 * in[1];
			out[0] = -rate1 - rate2;
			out[1] = -rate1 + rate2;
			out[2] = rate1;
			return true;
		}

		static void checkInvariants(double const* start, double const* end)
		{
			const double invariant1Start = start[0] + start[1] + 2 * start[2];
			const double invariant1End = end[0] + end[1] + 2 * end[2];

			CHECK(std::abs(invariant1End - invariant1Start) <= resTol);
		}

		static void initialPoint(double* init)
		{
			init[0] = 1.0;
			init[1] = 1.0;
			init[2] = 0.0;
		}

		template <typename rand_eng_t>
		static void randomInitialPoint(double* init, rand_eng_t& generator)
		{
			std::normal_distribution<double> distribution(0.0, 10.0);
			init[0] = std::abs(distribution(generator));
			init[1] = std::abs(distribution(generator));
			init[2] = std::abs(distribution(generator));
		}

		static void scaleFromInitialPoint(double* scale, double const* init)
		{
			scale[0] = 0.5;
			scale[1] = 1.0;
			scale[2] = 2.0;
		}

		static bool shouldRunTest(bool applyScaling, bool useVariant) { return true; }
	};

	struct SMAProblem
	{
		static constexpr double yCp[4] = {5.8377002519964755e+01, 2.9352296732047269e-03, 1.5061023667222263e-02, 1.3523701213590386e-01};
		static constexpr double _kA[4] = {0.0, 35.5, 1.59, 7.7};
		static constexpr double _kD[4] = {0.0, 1000.0, 1000.0, 1000.0};
		static constexpr double _lambda = 1.2e3;
		static constexpr double _nu[4] = {0.0, 4.7, 5.29, 3.7};
		static constexpr double _sigma[4] = {0.0, 11.83, 10.6, 10.0};

		static constexpr int problemSize = 4;
		static constexpr int maxIter = 1000;
		static constexpr double resTol = 1e-8;
		static constexpr int numNonMonotone = 100;
		static constexpr double initStep = 1e-4;

		template <typename T>
		static bool residual(T const* const x, T* const res)
		{
			using std::pow;

			// Salt equation: q_0 - Lambda + Sum[nu_j * q_j, j] == 0 
			//           <=>  q_0 == Lambda - Sum[nu_j * q_j, j] 
			// Also compute \bar{q}_0 = q_0 - Sum[sigma_j * q_j, j]
			res[0] = x[0] - _lambda;
			T q0_bar = x[0];

			for (int j = 1; j < 4; ++j)
			{
				res[0] += _nu[j] * x[j];
				q0_bar -= _sigma[j] * x[j];
			}

			// Protein equations: dq_i / dt - ( k_{a,i} * c_{p,i} * \bar{q}_0^{nu_i} - k_{d,i} * q_i * c_{p,0}^{nu_i} ) == 0
			//               <=>  dq_i / dt == k_{a,i} * c_{p,i} * \bar{q}_0^{nu_i} - k_{d,i} * q_i * c_{p,0}^{nu_i}
			for (int i = 1; i < 4; ++i)
			{
				const T c0_pow_nu = pow(yCp[0], _nu[i]);
				const T q0_bar_pow_nu = pow(q0_bar, _nu[i]);

				// Residual
				res[i] = _kD[i] * x[i] * c0_pow_nu - _kA[i] * yCp[i] * q0_bar_pow_nu;
			}
			return true;
		}

		static void checkInvariants(double const* start, double const* end) { }

		static void initialPoint(double* init)
		{
			init[0] = 1.0485785488181000e+03;
			init[1] = 1.1604726694141368e+01;
			init[2] = 1.1469542586742687e+01;
			init[3] = 9.7852311988018670e+00;
		}

		template <typename rand_eng_t>
		static void randomInitialPoint(double* init, rand_eng_t& generator) { }

		static void scaleFromInitialPoint(double* scale, double const* init) { }

		static bool shouldRunTest(bool applyScaling, bool useVariant) { return !applyScaling && useVariant; }
	};


	template <typename Func_t>
	auto makeADJacobian(Func_t func, cadet::active* adIn, cadet::active* adOut, int problemSize)
	{
		const auto jacobian = [=](double const* const point, cadet::linalg::detail::DenseMatrixBase& jac) -> bool
		{
			cadet::ad::copyToAd(point, adIn, problemSize);
			func(adIn, adOut);
			cadet::ad::extractDenseJacobianFromAd(adOut, 0, jac);
			return true;
		};

		return jacobian;
	}
}

template <typename prob_t, typename init_func_t>
void testPTCFixedInitialPoint(init_func_t initFunc)
{
	const int problemSize = prob_t::problemSize;

	std::vector<double> initPoint(problemSize, 0.0);
	initFunc(initPoint.data());

	// Copy for checking invariants later
	std::vector<double> x = initPoint;

	std::vector<double> scaleVals(problemSize, 0.0);
	prob_t::scaleFromInitialPoint(scaleVals.data(), initPoint.data());

	// Initialize AD and allocate AD vectors
	cadet::ad::setDirections(problemSize);

	std::vector<cadet::active> adIn(problemSize, 0.0);
	std::vector<cadet::active> adOut(problemSize, 0.0);

	// Set seed vectors
	cadet::ad::prepareAdVectorSeedsForDenseMatrix(adIn.data(), 0, problemSize);
	cadet::ad::fillAd(adIn.data(), problemSize, 0.0);

	// Use AD to compute dense Jacobian
	const auto jacobian = makeADJacobian(prob_t::template residual<cadet::active>, adIn.data(), adOut.data(), problemSize);

	std::vector<double> mem(4 * problemSize, 0.0);
	cadet::linalg::DenseMatrix jacMat;
	jacMat.resize(problemSize, problemSize);

	for (int idxScale = 0; idxScale < 2; ++idxScale)
	{
		// Check with and without scaling
		double const* scale = (idxScale == 0) ? nullptr : scaleVals.data();
		CAPTURE(idxScale);
		SECTION("Scaling of rows")
		{
			for (int variant = 0; variant < 2; ++variant)
			{
				// Check both variants for solving the linear system
				const bool useVariant = variant;
				CAPTURE(useVariant);

				// Check if we want to run this particular test
				if (!prob_t::shouldRunTest(idxScale != 0, useVariant))
					continue;

				SECTION("Numerical stability")
				{
					const bool res = cadet::nonlin::pseudoTransientContinuation(
						&prob_t::template residual<double>,
						jacobian,
						prob_t::maxIter,
						prob_t::resTol,
						prob_t::numNonMonotone,
						prob_t::initStep,
						scale,
						useVariant,
						x.data(),
						jacMat,
						mem.data(),
						problemSize);
					CHECK(res);

					std::vector<double> residual(problemSize, 0.0);
					prob_t::residual(x.data(), residual.data());

					CHECK(cadet::linalg::l2Norm(residual.data(), problemSize) <= prob_t::resTol);

					prob_t::checkInvariants(initPoint.data(), x.data());
				}
			}
		}
	}
}

template <typename prob_t>
void testPTCFixedInitialPoint()
{
	testPTCFixedInitialPoint<prob_t>(prob_t::initialPoint);
}

template <typename prob_t>
void testPTCRandomInitialPoint(int numRounds)
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	CAPTURE(seed);

	for (int i = 0; i < numRounds; ++i)
		testPTCFixedInitialPoint<prob_t>([&](double* init) -> void { prob_t::randomInitialPoint(init, generator); });
}


TEST_CASE("Pseudo-Transient Continuation Problem 1 fixed init", "[NonlinearSolver],[PTC],[CI]")
{
	testPTCFixedInitialPoint<Problem1>();
}

TEST_CASE("Pseudo-Transient Continuation Problem 1 random inits", "[NonlinearSolver],[PTC],[Flaky]")
{
	testPTCRandomInitialPoint<Problem1>(5);
}

TEST_CASE("Pseudo-Transient Continuation Single Three-Component Reaction fixed init", "[NonlinearSolver],[PTC],[CI]")
{
	testPTCFixedInitialPoint<SingleThreeCompReaction>();
}

TEST_CASE("Pseudo-Transient Continuation Single Three-Component Reaction random inits", "[NonlinearSolver],[PTC],[Flaky]")
{
	testPTCRandomInitialPoint<SingleThreeCompReaction>(5);
}

TEST_CASE("Pseudo-Transient Continuation Three-Component Mass Action Law Reaction fixed init", "[NonlinearSolver],[PTC],[CI]")
{
	testPTCFixedInitialPoint<ThreeCompMassActionLawReaction>();
}

TEST_CASE("Pseudo-Transient Continuation Three-Component Mass Action Law Reaction random inits", "[NonlinearSolver],[PTC],[Flaky]")
{
	testPTCRandomInitialPoint<ThreeCompMassActionLawReaction>(5);
}

TEST_CASE("Pseudo-Transient Continuation SMA fixed init", "[NonlinearSolver],[PTC],[CI]")
{
	testPTCFixedInitialPoint<SMAProblem>();
}

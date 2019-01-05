// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2018: The CADET Authors
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

#include "JacobianHelper.hpp"
#include "BindingModelTests.hpp"

#include "common/JsonParameterProvider.hpp"

#include "BindingModelFactory.hpp"
#include "model/BindingModel.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "AdUtils.hpp"
#include "AutoDiff.hpp"

#include <cstring>
#include <algorithm>

namespace
{
	inline RelApprox makeApprox(double val, double relTol, double absTol)
	{
		return RelApprox(val).epsilon(relTol).margin(absTol);
	}

	inline cadet::model::IBindingModel* createBindingModel(const char* name)
	{
		cadet::BindingModelFactory bmf;
		cadet::model::IBindingModel* const bm = bmf.create(name);
		
		REQUIRE(nullptr != bm);
		return bm;
	}

	class ConstExternalFunction : public cadet::IExternalFunction
	{
	public:
		virtual bool configure(cadet::IParameterProvider* paramProvider) { return true; }
		virtual const char* name() const CADET_NOEXCEPT { return "CONSTFUN"; }
		virtual double externalProfile(double t, double z, double r, unsigned int sec) { return 1.0; }
		virtual double timeDerivative(double t, double z, double r, unsigned int sec) { return 0.0; }
		virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections) { }
	};

	class LinearExternalFunction : public cadet::IExternalFunction
	{
	public:
		virtual bool configure(cadet::IParameterProvider* paramProvider) { return true; }
		virtual const char* name() const CADET_NOEXCEPT { return "LINFUN"; }
		virtual double externalProfile(double t, double z, double r, unsigned int sec) { return t; }
		virtual double timeDerivative(double t, double z, double r, unsigned int sec) { return 1.0; }
		virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections) { }
	};

	class ConfiguredBindingModel
	{
	public:

		ConfiguredBindingModel(ConfiguredBindingModel&& cpy) CADET_NOEXCEPT 
			: _binding(cpy._binding), _nComp(cpy._nComp), _nBound(cpy._nBound), _boundOffset(cpy._boundOffset), _buffer(cpy._buffer), _extFuns(cpy._extFuns)
		{
			cpy._binding = nullptr;
			cpy._nBound = nullptr;
			cpy._boundOffset = nullptr;
			cpy._buffer = nullptr;
			cpy._extFuns = nullptr;
		}

		~ConfiguredBindingModel()
		{
			if (_buffer)
				delete[] _buffer;
			if (_extFuns)
				delete[] _extFuns;
			delete[] _boundOffset;
			delete _binding;
		}

		inline ConfiguredBindingModel& operator=(ConfiguredBindingModel&& cpy) CADET_NOEXCEPT
		{
			_binding = cpy._binding;
			_nComp = cpy._nComp;
			_nBound = cpy._nBound;
			_boundOffset = cpy._boundOffset;
			_buffer = cpy._buffer;
			_extFuns = cpy._extFuns;

			cpy._binding = nullptr;
			cpy._nBound = nullptr;
			cpy._boundOffset = nullptr;
			cpy._buffer = nullptr;
			cpy._extFuns = nullptr;

			return *this;			
		}

		static ConfiguredBindingModel create(const char* name, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config)
		{
			cadet::model::IBindingModel* const bm = createBindingModel(name);

			// Calculate offset of bound states
			unsigned int* boundOffset = new unsigned int[nComp];
			boundOffset[0] = 0;
			for (unsigned int i = 1; i < nComp; ++i)
			{
				boundOffset[i] = boundOffset[i-1] + nBound[i-1];
			}

			// Configure
			cadet::JsonParameterProvider jpp(config);
			bm->configureModelDiscretization(jpp, nComp, nBound, boundOffset);
			if (bm->requiresConfiguration())
			{
				jpp.set("IS_KINETIC", isKinetic);
				jpp.set("EXTFUN", std::vector<int>(1, 0));
				REQUIRE(bm->configure(jpp, 0, 0));
			}

			// Assign external functions
			cadet::IExternalFunction* extFuns = new LinearExternalFunction[50];
			bm->setExternalFunctions(&extFuns, 50);

			// Allocate memory buffer
			unsigned int requiredMem = 0;
			if (bm->requiresWorkspace())
				requiredMem = bm->workspaceSize();

			char* buffer = nullptr;
			if (requiredMem > 0)
			{
				buffer = new char[requiredMem];
				std::memset(buffer, 0, requiredMem);
			}

			return ConfiguredBindingModel(bm, nComp, nBound, boundOffset, buffer, extFuns);
		}

		inline cadet::model::IBindingModel& model() { return *_binding; }
		inline const cadet::model::IBindingModel& model() const { return *_binding; }

		inline void* buffer() { return _buffer; }
		inline unsigned int nComp() const { return _nComp; }
		inline unsigned int const* nBound() const { return _nBound; }
		inline unsigned int const* boundOffset() const { return _boundOffset; }

		inline unsigned int numBoundStates() const { return _boundOffset[_nComp - 1] + _nBound[_nComp - 1]; }

	private:

		ConfiguredBindingModel(cadet::model::IBindingModel* binding, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset, char* buffer, cadet::IExternalFunction* extFuns) 
			: _binding(binding), _nComp(nComp), _nBound(nBound), _boundOffset(boundOffset), _buffer(buffer), _extFuns(extFuns)
		{
		}

		cadet::model::IBindingModel* _binding;
		unsigned int _nComp;
		unsigned int const* _nBound;
		unsigned int const* _boundOffset;
		char* _buffer;
		cadet::IExternalFunction* _extFuns;
	};
}

namespace cadet
{

namespace test
{

namespace binding
{

void testJacobianAD(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, double const* point, double absTol, double relTol)
{
	ConfiguredBindingModel cbm = ConfiguredBindingModel::create(modelName, nComp, nBound, isKinetic, config);

	const unsigned int numDofs = cbm.nComp() + cbm.numBoundStates();
	const unsigned int numEq = cbm.numBoundStates();
	std::vector<double> yState(numDofs, 0.0);
	double* const yBound = yState.data() + cbm.nComp();
	std::copy_n(point, numDofs, yState.data());

	std::vector<double> yStateDot(numDofs, 0.0);
	double* const yBoundDot = yStateDot.data() + cbm.nComp();

	std::vector<double> dir(numDofs, 0.0);
	std::vector<double> colA(numEq, 0.0);
	std::vector<double> colB(numEq, 0.0);

	// Calculate analytic Jacobian
	cadet::linalg::DenseMatrix jacAna;
	jacAna.resize(numDofs, numDofs);
	cbm.model().analyticJacobian(1.0, 0.0, 0.0, 0u, yBound, jacAna.row(cbm.nComp()), cbm.buffer());

	// Enable AD
	cadet::ad::setDirections(cadet::ad::getMaxDirections());
	cadet::active* adRes = new cadet::active[numEq];
	cadet::active* adY = new cadet::active[numDofs];

	// Evaluate with AD
	ad::prepareAdVectorSeedsForDenseMatrix(adY, 0, numDofs);
	ad::copyToAd(yState.data(), adY, numDofs);
	cbm.model().residual(1.0, 0.0, 0.0, 0u, 1.0, adY + cbm.nComp(), yBoundDot, adRes, cbm.buffer());

	// Extract Jacobian
	cadet::linalg::DenseMatrix jacAD;
	jacAD.resize(numEq, numDofs);
	ad::extractDenseJacobianFromAd(adRes, 0, jacAD);

	delete[] adY;
	delete[] adRes;

	cadet::test::checkJacobianPatternFD(
		[&](double const* lDir, double* res) -> void { cbm.model().residual(1.0, 0.0, 0.0, 0u, 1.0, lDir + cbm.nComp(), yBoundDot, res, cbm.buffer()); },
		[&](double const* lDir, double* res) -> void { jacAna.submatrixMultiplyVector(lDir, cbm.nComp(), 0, numEq, numDofs, res); }, 
		yState.data(), dir.data(), colA.data(), colB.data(), numDofs, numEq);

	cadet::test::checkJacobianPatternFD(
		[&](double const* lDir, double* res) -> void { cbm.model().residual(1.0, 0.0, 0.0, 0u, 1.0, lDir + cbm.nComp(), yBoundDot, res, cbm.buffer()); }, 
		[&](double const* lDir, double* res) -> void { jacAD.multiplyVector(lDir, res); }, 
		yState.data(), dir.data(), colA.data(), colB.data(), numDofs, numEq);

	// Check Jacobians against each other
	for (unsigned int row = 0; row < numEq; ++row)
	{
		for (unsigned int col = 0; col < numDofs; ++col)
		{
			CAPTURE(row);
			CAPTURE(col);
			CHECK(jacAna.native(row + cbm.nComp(), col) == makeApprox(jacAD.native(row, col), absTol, relTol));
		}
	}
}

void testTimeDerivativeJacobianFD(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, double h, double absTol, double relTol)
{
	ConfiguredBindingModel cbm = ConfiguredBindingModel::create(modelName, nComp, nBound, isKinetic, config);

	const unsigned int numDofs = cbm.nComp() + cbm.numBoundStates();
	const unsigned int numEq = cbm.numBoundStates();
	std::vector<double> yState(numDofs, 0.0);
	double* const yBound = yState.data() + cbm.nComp();

	std::vector<double> dir(numDofs, 0.0);
	std::vector<double> colA(numEq, 0.0);
	std::vector<double> colB(numEq, 0.0);

	// Calculate time derivative Jacobian
	cadet::linalg::DenseMatrix jacMat;
	jacMat.resize(numDofs, numDofs);
	jacMat.setAll(0.0);
	cbm.model().jacobianAddDiscretized(1.0, jacMat.row(cbm.nComp()));

	cadet::test::compareJacobianFD(
		[&](double const* lDir, double* res) -> void { cbm.model().residual(1.0, 0.0, 0.0, 0u, 1.0, yBound, lDir + cbm.nComp(), res, cbm.buffer()); }, 
		[&](double const* lDir, double* res) -> void { jacMat.submatrixMultiplyVector(lDir, cbm.nComp(), 0, numEq, numDofs, res); }, 
		yState.data(), dir.data(), colA.data(), colB.data(), numDofs, numEq, h, absTol, relTol);
}

void testTimeDerivativeJacobianMultiplyFunction(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, double absTol, double relTol)
{
	ConfiguredBindingModel cbm = ConfiguredBindingModel::create(modelName, nComp, nBound, isKinetic, config);

	const unsigned int numDofs = cbm.nComp() + cbm.numBoundStates();
	const unsigned int numEq = cbm.numBoundStates();

	std::vector<double> dir(numDofs, 0.0);
	std::vector<double> colA(numEq, 0.0);
	std::vector<double> colB(numEq, 0.0);

	// Calculate time derivative Jacobian
	cadet::linalg::DenseMatrix jacMat;
	jacMat.resize(numDofs, numDofs);
	jacMat.setAll(0.0);
	cbm.model().jacobianAddDiscretized(1.0, jacMat.row(cbm.nComp()));

	cadet::test::compareJacobian(
		[&](double const* lDir, double* res) -> void { cbm.model().multiplyWithDerivativeJacobian(lDir + cbm.nComp(), res, 1.0); }, 
		[&](double const* lDir, double* res) -> void { jacMat.submatrixMultiplyVector(lDir, cbm.nComp(), 0, numEq, numDofs, res); }, 
		dir.data(), colA.data(), colB.data(), numDofs, numEq, absTol, relTol);
}

void testConsistentInitialization(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, bool useAD, double const* point, double consTol, double absTol)
{
	ConfiguredBindingModel cbm = ConfiguredBindingModel::create(modelName, nComp, nBound, isKinetic, config);

	if (!cbm.model().hasAlgebraicEquations())
		return;

	const unsigned int numDofs = cbm.nComp() + cbm.numBoundStates();
	const unsigned int numEq = cbm.numBoundStates();
	std::vector<double> yState(numDofs, 0.0);
	double* const yBound = yState.data() + cbm.nComp();
	std::copy_n(point, numDofs, yState.data());

	// Memory for matrix operations
	cadet::linalg::DenseMatrix jacMat;
	jacMat.resize(cbm.numBoundStates(), cbm.numBoundStates());

	// Perform consistent initialization
	if (useAD)
	{
		// Enable AD
		cadet::ad::setDirections(cadet::ad::getMaxDirections());
		cadet::active* adRes = new cadet::active[numEq];
		cadet::active* adY = new cadet::active[numDofs];
		ad::prepareAdVectorSeedsForDenseMatrix(adY, 0, numDofs);

		cbm.model().consistentInitialState(1.0, 0.0, 0.0, 0u, yBound, consTol, adRes, adY + cbm.nComp(), 0, nComp, ad::DenseJacobianExtractor(), reinterpret_cast<double*>(cbm.buffer()), jacMat);

		delete[] adY;
		delete[] adRes;
	}
	else
	{
		cbm.model().consistentInitialState(1.0, 0.0, 0.0, 0u, yBound, consTol, nullptr, nullptr, 0, 0, ad::DenseJacobianExtractor(), reinterpret_cast<double*>(cbm.buffer()), jacMat);
	}

	// Check
	std::vector<double> res(numEq, 0.0);
	cbm.model().residual(1.0, 0.0, 0.0, 0u, 1.0, yBound, nullptr, res.data(), cbm.buffer());

	unsigned int algStart = 0;
	unsigned int algLen = 0;
	cbm.model().getAlgebraicBlock(algStart, algLen);

	for (unsigned int i = algStart; i < algStart + algLen; ++i)
	{
		CAPTURE(i);
		CHECK(std::abs(res[i]) <= absTol);
	}
}

void testConsistencyOfAlgebraicEquationFunctions(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config)
{
	ConfiguredBindingModel cbm = ConfiguredBindingModel::create(modelName, nComp, nBound, isKinetic, config);

	const unsigned int numDofs = cbm.nComp() + cbm.numBoundStates();

	// Calculate time derivative Jacobian
	cadet::linalg::DenseMatrix jacMat;
	jacMat.resize(numDofs, numDofs);
	jacMat.setAll(0.0);
	cbm.model().jacobianAddDiscretized(1.0, jacMat.row(cbm.nComp()));

	if (cbm.model().hasAlgebraicEquations())
	{
		// There must be non-zero rows outside the algebraic block and all zero rows inside
		unsigned int algStart = 0;
		unsigned int algLen = 0;
		cbm.model().getAlgebraicBlock(algStart, algLen);

		// Before algebraic block
		for (unsigned int r = cbm.nComp(); r < cbm.nComp() + algStart; ++r)
		{
			bool timeJacobianRowIsNonZeroOutsideAlgBlock = false;
			for (unsigned int col = 0; col < numDofs; ++col)
			{
				if (jacMat.native(r, col) != 0.0)
				{
					timeJacobianRowIsNonZeroOutsideAlgBlock = true;
					break;
				}
			}
			const unsigned int row = r - cbm.nComp();
			CAPTURE(row);
			CHECK(timeJacobianRowIsNonZeroOutsideAlgBlock);
		}

		// In algebraic block
		for (unsigned int r = cbm.nComp() + algStart; r < cbm.nComp() + algStart + algLen; ++r)
		{
			bool timeJacobianRowIsAllZeroInAlgBlock = true;
			for (unsigned int col = 0; col < numDofs; ++col)
			{
				if (jacMat.native(r, col) != 0.0)
				{
					timeJacobianRowIsAllZeroInAlgBlock = false;
					break;
				}
			}
			const unsigned int row = r - cbm.nComp();
			CAPTURE(row);
			CHECK(timeJacobianRowIsAllZeroInAlgBlock);
		}

		// Behind algebraic block
		for (unsigned int r = cbm.nComp() + algStart + algLen; r < numDofs; ++r)
		{
			bool timeJacobianRowIsNonZeroOutsideAlgBlock = false;
			for (unsigned int col = 0; col < numDofs; ++col)
			{
				if (jacMat.native(r, col) != 0.0)
				{
					timeJacobianRowIsNonZeroOutsideAlgBlock = true;
					break;
				}
			}
			const unsigned int row = r - cbm.nComp();
			CAPTURE(row);
			CHECK(timeJacobianRowIsNonZeroOutsideAlgBlock);
		}
	}
	else
	{
		// There must not be a row with all zeros in jacMat
		for (unsigned int r = cbm.nComp(); r < numDofs; ++r)
		{
			bool timeJacobianRowIsNonZero = false;
			for (unsigned int col = 0; col < numDofs; ++col)
			{
				if (jacMat.native(r, col) != 0.0)
				{
					timeJacobianRowIsNonZero = true;
					break;
				}
			}
			const unsigned int row = r - cbm.nComp();
			CAPTURE(row);
			CHECK(timeJacobianRowIsNonZero);
		}
	}
}

void testNormalExternalConsistency(const char* modelName, const char* modelNameExt, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, double const* point)
{
	ConfiguredBindingModel cbm = ConfiguredBindingModel::create(modelName, nComp, nBound, isKinetic, config);
	ConfiguredBindingModel cbmExt = ConfiguredBindingModel::create(modelNameExt, nComp, nBound, isKinetic, config);

	const unsigned int numDofs = cbm.nComp() + cbm.numBoundStates();
	const unsigned int numEq = cbm.numBoundStates();

	std::vector<double> res(numEq, 0.0);
	std::vector<double> resExt(numEq, 0.0);
	std::vector<double> yState(numDofs, 0.0);
	double* const yBound = yState.data() + cbm.nComp();
	std::copy_n(point, numDofs, yState.data());

	// Evaluate residuals
	cbm.model().residual(1.0, 0.0, 0.0, 0u, 1.0, yBound, nullptr, res.data(), cbm.buffer());
	cbmExt.model().residual(1.0, 0.0, 0.0, 0u, 1.0, yBound, nullptr, resExt.data(), cbmExt.buffer());

	// Check residuals against each other
	for (unsigned int i = 0; i < numEq; ++i)
	{
		CAPTURE(i);
		CHECK(resExt[i] == RelApprox(res[i]));
	}

	// Calculate analytic Jacobians
	cadet::linalg::DenseMatrix jac;
	jac.resize(numDofs, numDofs);
	cbm.model().analyticJacobian(1.0, 0.0, 0.0, 0u, yBound, jac.row(cbm.nComp()), cbm.buffer());

	cadet::linalg::DenseMatrix jacExt;
	jacExt.resize(numDofs, numDofs);
	cbmExt.model().analyticJacobian(1.0, 0.0, 0.0, 0u, yBound, jacExt.row(cbmExt.nComp()), cbmExt.buffer());

	// Check Jacobians against each other
	for (unsigned int r = cbm.nComp(); r < numDofs; ++r)
	{
		const unsigned int row = r - cbm.nComp();
		for (unsigned int col = 0; col < numDofs; ++col)
		{
			CAPTURE(row);
			CAPTURE(col);
			CHECK(jac.native(r, col) == RelApprox(jacExt.native(r, col)));
		}
	}
}

void testNonBindingConsistency(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, bool useAD, double const* point)
{
	ConfiguredBindingModel cbm = ConfiguredBindingModel::create(modelName, nComp, nBound, isKinetic, config);

	const unsigned int numDofs = cbm.nComp() + cbm.numBoundStates();
	const unsigned int numEq = cbm.numBoundStates();

	std::vector<double> yState(numDofs, 0.0);
	double* const yBound = yState.data() + cbm.nComp();
	std::copy_n(point, numDofs, yState.data());

	std::vector<double> yStateDot(numDofs, 0.0);
	double* const yBoundDot = yStateDot.data() + cbm.nComp();

	// Calculate Jacobian
	cadet::linalg::DenseMatrix jac;
	if (useAD)
	{
		jac.resize(numEq, numDofs);

		// Enable AD
		cadet::ad::setDirections(cadet::ad::getMaxDirections());
		cadet::active* adRes = new cadet::active[numEq];
		cadet::active* adY = new cadet::active[numDofs];

		// Evaluate with AD
		ad::prepareAdVectorSeedsForDenseMatrix(adY, 0, numDofs);
		ad::copyToAd(yState.data(), adY, numDofs);
		cbm.model().residual(1.0, 0.0, 0.0, 0u, 1.0, adY + cbm.nComp(), yBoundDot, adRes, cbm.buffer());

		// Extract Jacobian
		ad::extractDenseJacobianFromAd(adRes, 0, jac);

		delete[] adY;
		delete[] adRes;
	}
	else
	{
		jac.resize(numDofs, numDofs);
		cbm.model().analyticJacobian(1.0, 0.0, 0.0, 0u, yBound, jac.row(cbm.nComp()), cbm.buffer());
	}

	// Check that columns of non-binding liquid phase components are all zero
	for (unsigned int col = 0; col < nComp; ++col)
	{
		if (nBound[col] > 0)
			continue;

		for (unsigned int row = 0; row < numEq; ++row)
		{
			if (useAD)
			{
				CAPTURE(row);
				CAPTURE(col);
				CHECK(jac.native(row, col) == 0.0);
			}
			else
			{
				CAPTURE(row);
				CAPTURE(col);
				CHECK(jac.native(row + cbm.nComp(), col) == 0.0);
			}
		}
	}
}

void testNonbindingBindingConsistency(const char* modelName, unsigned int nCompBnd, unsigned int nCompNonBnd, unsigned int const* nBound, unsigned int const* nBoundNonBnd, bool isKinetic, const char* configBnd, const char* configNonBnd, bool useAD, double const* pointBnd, double const* pointNonBnd)
{
	ConfiguredBindingModel cbmBnd = ConfiguredBindingModel::create(modelName, nCompBnd, nBound, isKinetic, configBnd);
	ConfiguredBindingModel cbmNonBnd = ConfiguredBindingModel::create(modelName, nCompNonBnd, nBoundNonBnd, isKinetic, configNonBnd);

	// Setup all binding
	const unsigned int numDofsBnd = cbmBnd.nComp() + cbmBnd.numBoundStates();
	const unsigned int numEqBnd = cbmBnd.numBoundStates();

	std::vector<double> yStateBnd(numDofsBnd, 0.0);
	double* const yBoundBnd = yStateBnd.data() + cbmBnd.nComp();
	std::copy_n(pointBnd, numDofsBnd, yStateBnd.data());

	std::vector<double> yStateBndDot(numDofsBnd, 0.0);
	double* const yBoundBndDot = yStateBndDot.data() + cbmBnd.nComp();

	// Setup with nonbinding components
	const unsigned int numDofsNonBnd = cbmNonBnd.nComp() + cbmNonBnd.numBoundStates();
	const unsigned int numEqNonBnd = cbmNonBnd.numBoundStates();

	std::vector<double> yStateNonBnd(numDofsNonBnd, 0.0);
	double* const yBoundNonBnd = yStateNonBnd.data() + cbmNonBnd.nComp();
	std::copy_n(pointNonBnd, numDofsNonBnd, yStateNonBnd.data());

	std::vector<double> yStateNonBndDot(numDofsNonBnd, 0.0);
	double* const yBoundNonBndDot = yStateNonBndDot.data() + cbmNonBnd.nComp();

	// Number of equations has to be the same
	REQUIRE(numEqBnd == numEqNonBnd);

	// Check residual
	std::vector<double> res1(numEqBnd, 0.0);
	std::vector<double> res2(numEqBnd, 0.0);

	cbmBnd.model().residual(1.0, 0.0, 0.0, 0u, 1.0, yBoundBnd, yBoundBndDot, res1.data(), cbmBnd.buffer());
	cbmNonBnd.model().residual(1.0, 0.0, 0.0, 0u, 1.0, yBoundNonBnd, yBoundNonBndDot, res2.data(), cbmNonBnd.buffer());

	for (unsigned int i = 0; i < numEqBnd; ++i)
	{
		CAPTURE(i);
		CHECK(res1[i] == res2[i]);
	}

	// Calculate Jacobians
	cadet::linalg::DenseMatrix jacBnd;
	cadet::linalg::DenseMatrix jacNonBnd;
	if (useAD)
	{
		jacBnd.resize(numEqBnd, numDofsBnd);
		jacNonBnd.resize(numEqBnd, numDofsNonBnd);

		// Enable AD
		cadet::ad::setDirections(cadet::ad::getMaxDirections());
		cadet::active* adRes = new cadet::active[numEqBnd];
		cadet::active* adY = new cadet::active[numDofsNonBnd];

		// Evaluate with AD, all binding
		ad::prepareAdVectorSeedsForDenseMatrix(adY, 0, numDofsBnd);
		ad::copyToAd(yStateBnd.data(), adY, numDofsBnd);
		cbmBnd.model().residual(1.0, 0.0, 0.0, 0u, 1.0, adY + cbmBnd.nComp(), yBoundBndDot, adRes, cbmBnd.buffer());

		// Extract Jacobian, all binding
		ad::extractDenseJacobianFromAd(adRes, 0, jacBnd);

		// Evaluate with AD, with nonbinding
		ad::resetAd(adRes, numEqBnd);
		ad::resetAd(adY, numDofsNonBnd);
		ad::prepareAdVectorSeedsForDenseMatrix(adY, 0, numDofsNonBnd);
		ad::copyToAd(yStateNonBnd.data(), adY, numDofsNonBnd);
		cbmNonBnd.model().residual(1.0, 0.0, 0.0, 0u, 1.0, adY + cbmNonBnd.nComp(), yBoundNonBndDot, adRes, cbmNonBnd.buffer());

		// Extract Jacobian, with nonbinding
		ad::extractDenseJacobianFromAd(adRes, 0, jacNonBnd);

		delete[] adY;
		delete[] adRes;
	}
	else
	{
		jacBnd.resize(numDofsBnd, numDofsBnd);
		jacNonBnd.resize(numDofsNonBnd, numDofsNonBnd);
		cbmBnd.model().analyticJacobian(1.0, 0.0, 0.0, 0u, yBoundBnd, jacBnd.row(cbmBnd.nComp()), cbmBnd.buffer());
		cbmNonBnd.model().analyticJacobian(1.0, 0.0, 0.0, 0u, yBoundNonBnd, jacNonBnd.row(cbmNonBnd.nComp()), cbmNonBnd.buffer());
	}

	// Compare Jacobians
	const unsigned int rowOffsetBnd = useAD ? 0 : cbmBnd.nComp();
	const unsigned int rowOffsetNonBnd = useAD ? 0 : cbmNonBnd.nComp();
	unsigned int colBnd = 0;
	for (unsigned int col = 0; col < numDofsNonBnd; ++col)
	{
		// Skip non-binding liquid phase columns
		if ((col < nCompNonBnd) && (nBoundNonBnd[col] == 0))
			continue;

		for (unsigned int row = 0; row < numEqBnd; ++row)
		{
			CAPTURE(col);
			CAPTURE(colBnd);
			CAPTURE(row);
			CHECK(jacBnd.native(row + rowOffsetBnd, colBnd) == jacNonBnd.native(row + rowOffsetNonBnd, col));
		}
		++colBnd;
	}
}

} // namespace binding
} // namespace test
} // namespace cadet

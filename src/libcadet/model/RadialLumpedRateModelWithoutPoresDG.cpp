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

#include "model/RadialLumpedRateModelWithoutPoresDG.hpp"
#include "BindingModelFactory.hpp"
#include "ParamReaderHelper.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/ExternalFunction.hpp"
#include "cadet/SolutionRecorder.hpp"
#include "ConfigurationHelper.hpp"
#include "model/BindingModel.hpp"
#include "SimulationTypes.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/Norms.hpp"
#include "linalg/Subset.hpp"
#include "model/parts/BindingCellKernel.hpp"

#include "AdUtils.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <functional>

#include "ParallelSupport.hpp"
#ifdef CADET_PARALLELIZE
#include <tbb/parallel_for.h>
#endif

#define EIGEN_USE_MKL_ALL

using namespace Eigen;

namespace cadet
{

	namespace model
	{

		RadialLumpedRateModelWithoutPoresDG::RadialLumpedRateModelWithoutPoresDG(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
			_analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr), _initC(0),
			_initCs(0), _initState(0), _initStateDot(0)
		{
			// Multiple particle types are not supported
			_singleBinding = true;
			_singleDynReaction = true;
		}

		RadialLumpedRateModelWithoutPoresDG::~RadialLumpedRateModelWithoutPoresDG() CADET_NOEXCEPT
		{
			delete[] _tempState;

			delete _linearSolver;

			for (IDynamicReactionModel* p : _dynReaction) delete p;
		}

		unsigned int RadialLumpedRateModelWithoutPoresDG::numDofs() const CADET_NOEXCEPT
		{
			// Column bulk DOFs: nElem * nNodes * nComp mobile phase and nElem * nNodes * (sum boundStates) solid phase
			// Inlet DOFs: nComp
			return _disc.nPoints * (_disc.nComp + _disc.strideBound) + _disc.nComp;
		}

		unsigned int RadialLumpedRateModelWithoutPoresDG::numPureDofs() const CADET_NOEXCEPT
		{
			// Column bulk DOFs: nElem * nNodes * nComp mobile phase and nElem * nNodes * (sum boundStates) solid phase
			return _disc.nPoints * (_disc.nComp + _disc.strideBound);
		}


		bool RadialLumpedRateModelWithoutPoresDG::usesAD() const CADET_NOEXCEPT
		{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
			// We always need AD if we want to check the analytical Jacobian
			return true;
#else
			// We only need AD if we are not computing the Jacobian analytically
			return !_analyticJac;
#endif
		}

		bool RadialLumpedRateModelWithoutPoresDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper)
		{
			const bool firstConfigCall = _tempState == nullptr;

			// Read discretization
			_disc.nComp = paramProvider.getInt("NCOMP");

			if (paramProvider.exists("NPARTYPE") ? paramProvider.getInt("NPARTYPE") != 1 : false)
				throw InvalidParameterException("Number of particle types must be 1 for EQUILIBRIUM particles, i.e. RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES");

			paramProvider.pushScope("particle_type_000");

			if (paramProvider.getBool("HAS_FILM_DIFFUSION"))
				throw InvalidParameterException("HAS_FILM_DIFFUSION must be false for RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES");
			if (paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false)
				throw InvalidParameterException("HAS_PORE_DIFFUSION must be false for RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES");
			if (paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false)
				throw InvalidParameterException("HAS_SURFACE_DIFFUSION must be false for RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES");

			std::vector<int> nBound;
			nBound = paramProvider.getIntArray("NBOUND");
			if (nBound.size() < _disc.nComp)
				throw InvalidParameterException("Field NBOUND contains too few elements (NCOMP = " + std::to_string(_disc.nComp) + " required)");

			if (firstConfigCall)
				_disc.nBound = new unsigned int[_disc.nComp];
			std::copy_n(nBound.begin(), _disc.nComp, _disc.nBound);

			paramProvider.popScope();

			paramProvider.pushScope("discretization");

			if (firstConfigCall)
				_linearSolver = cadet::linalg::setLinearSolver(paramProvider.exists("LINEAR_SOLVER") ? paramProvider.getString("LINEAR_SOLVER") : "SparseLU");

			if (paramProvider.exists("POLYDEG"))
				_disc.polyDeg = paramProvider.getInt("POLYDEG");
			else
				_disc.polyDeg = 4u; // default value
			if (paramProvider.getInt("POLYDEG") < 1)
				throw InvalidParameterException("Polynomial degree must be at least 1!");
			else if (_disc.polyDeg < 3)
				LOG(Warning) << "Polynomial degree > 2 in bulk discretization (cf. POLYDEG) is always recommended for performance reasons.";

			_disc.nNodes = _disc.polyDeg + 1;

			if (paramProvider.exists("NELEM"))
				_disc.nElem = paramProvider.getInt("NELEM");
			else if (paramProvider.exists("NCOL"))
				_disc.nElem = std::max(1u, paramProvider.getInt("NCOL") / _disc.nNodes);
			else
				throw InvalidParameterException("Specify field NELEM (or NCOL)");

			if (_disc.nElem < 1)
				throw InvalidParameterException("Number of radial elements must be at least 1!");

			_disc.nPoints = _disc.nNodes * _disc.nElem;

			// Radial DG always uses exact integration
			int polynomial_integration_mode = 1;

			// Precompute offsets and total number of bound states
			if (firstConfigCall)
				_disc.boundOffset = new unsigned int[_disc.nComp];
			_disc.boundOffset[0] = 0;
			for (unsigned int i = 1; i < _disc.nComp; ++i)
			{
				_disc.boundOffset[i] = _disc.boundOffset[i - 1] + _disc.nBound[i - 1];
			}
			_disc.strideBound = _disc.boundOffset[_disc.nComp - 1] + _disc.nBound[_disc.nComp - 1];

			// Determine whether analytic Jacobian should be used
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
			const bool analyticJac = paramProvider.getBool("USE_ANALYTIC_JACOBIAN");
#else
			const bool analyticJac = false;
#endif

			// Allocate space for initial conditions
			_initC.resize(_disc.nElem * _disc.nNodes * _disc.nComp);
			_initCs.resize(_disc.nElem * _disc.nNodes * _disc.strideBound);

			// Create nonlinear solver for consistent initialization
			configureNonlinearSolver(paramProvider);

			paramProvider.popScope();

			const unsigned int strideNode = _disc.nComp + _disc.strideBound;
			const bool transportSuccess = _convDispOp.configureModelDiscretization(paramProvider, helper, _disc.nComp, polynomial_integration_mode, _disc.nElem, _disc.polyDeg, strideNode);

			_disc.curSection = -1;

			// Allocate memory
			Indexer idxr(_disc);

			// Radial DG always uses exact integration
			_jacInlet.resize(_disc.nNodes, 1);
			_jac.resize((_disc.nComp + _disc.strideBound) * _disc.nPoints, (_disc.nComp + _disc.strideBound) * _disc.nPoints);
			_jacDisc.resize((_disc.nComp + _disc.strideBound) * _disc.nPoints, (_disc.nComp + _disc.strideBound) * _disc.nPoints);

			// Set whether analytic Jacobian is used
			useAnalyticJacobian(analyticJac);

			// ==== Construct and configure binding model

			clearBindingModels();
			_binding.push_back(nullptr);

			paramProvider.pushScope("particle_type_000");

			if (paramProvider.exists("ADSORPTION_MODEL"))
				_binding[0] = helper.createBindingModel(paramProvider.getString("ADSORPTION_MODEL"));
			else
				_binding[0] = helper.createBindingModel("NONE");

			if (!_binding[0])
				throw InvalidParameterException("Unknown binding model " + paramProvider.getString("ADSORPTION_MODEL"));

			bool bindingConfSuccess = true;
			if (_binding[0]->usesParamProviderInDiscretizationConfig())
				paramProvider.pushScope("adsorption");

			bindingConfSuccess = _binding[0]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound, _disc.boundOffset);

			if (_binding[0]->usesParamProviderInDiscretizationConfig())
				paramProvider.popScope();

			paramProvider.popScope();

			// ==== Construct and configure dynamic reaction model
			bool reactionConfSuccess = true;
			clearDynamicReactionModels();
			_reaction.getDynReactionVector("liquid").clear();
			for (IDynamicReactionModel* p : _dynReaction) delete p;
			_dynReaction.clear();
			_dynReaction.push_back(nullptr);

			if (paramProvider.exists("REACTION_MODEL"))
			{
				_dynReaction[0] = helper.createDynamicReactionModel(paramProvider.getString("REACTION_MODEL"));
				if (!_dynReaction[0])
					throw InvalidParameterException("Unknown dynamic reaction model " + paramProvider.getString("REACTION_MODEL"));

				if (_dynReaction[0]->usesParamProviderInDiscretizationConfig())
					paramProvider.pushScope("reaction");

				reactionConfSuccess = _dynReaction[0]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound, _disc.boundOffset);

				if (_dynReaction[0]->usesParamProviderInDiscretizationConfig())
					paramProvider.popScope();
			}

			// Sync _reaction for CellParameters (non-owning refs)
			if (_dynReaction[0]) _reaction.getDynReactionVector("liquid").push_back(_dynReaction[0]);

			// Setup the memory for tempState based on state vector
			if (firstConfigCall)
				_tempState = new double[numDofs()];

			return bindingConfSuccess && reactionConfSuccess && transportSuccess;
		}

		bool RadialLumpedRateModelWithoutPoresDG::configure(IParameterProvider& paramProvider)
		{
			_parameters.clear();

			_convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

			// Read geometry parameters
			_totalPorosity = paramProvider.getDouble("TOTAL_POROSITY");

			// Add parameters to map
			_parameters[makeParamId(hashString("TOTAL_POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_totalPorosity;

			// Register initial conditions parameters
			for (unsigned int i = 0; i < _disc.nComp; ++i)
				_parameters[makeParamId(hashString("INIT_C"), _unitOpIdx, i, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = _initC.data() + i;

			if (_binding[0])
			{
				std::vector<ParameterId> initParams(_disc.strideBound);
				_binding[0]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, cadet::ParTypeIndep);

				for (unsigned int i = 0; i < _disc.strideBound; ++i)
					_parameters[initParams[i]] = _initCs.data() + i;
			}

			// Reconfigure binding model
			paramProvider.pushScope("particle_type_000");

			bool bindingConfSuccess = true;
			if (_binding[0] && paramProvider.exists("adsorption") && _binding[0]->requiresConfiguration())
			{
				paramProvider.pushScope("adsorption");
				bindingConfSuccess = _binding[0]->configure(paramProvider, _unitOpIdx, cadet::ParTypeIndep);
				paramProvider.popScope();
			}

			paramProvider.popScope();

			// Reconfigure dynamic reaction model
			bool reactionConfSuccess = true;
			if (_dynReaction[0] && paramProvider.exists("reaction") && _dynReaction[0]->requiresConfiguration())
			{
				paramProvider.pushScope("reaction");
				reactionConfSuccess = _dynReaction[0]->configure(paramProvider, _unitOpIdx, cadet::ParTypeIndep);
				paramProvider.popScope();
			}

			setPattern(_jac, true, _dynReaction[0] && (_dynReaction[0]->numReactionsCombined() > 0));
			setPattern(_jacDisc, true, _dynReaction[0] && (_dynReaction[0]->numReactionsCombined() > 0));

			_linearSolver->analyzePattern(_jacDisc);

			return bindingConfSuccess && reactionConfSuccess;
		}

		unsigned int RadialLumpedRateModelWithoutPoresDG::threadLocalMemorySize() const CADET_NOEXCEPT
		{
			LinearMemorySizer lms;

			if (_binding[0] && _binding[0]->requiresWorkspace())
				lms.addBlock(_binding[0]->workspaceSize(_disc.nComp, _disc.strideBound, _disc.nBound));

			if (_dynReaction[0])
			{
				lms.addBlock(_dynReaction[0]->workspaceSize(_disc.nComp, _disc.strideBound, _disc.nBound));
				lms.add<active>(_disc.strideBound);
				lms.add<double>(_disc.strideBound * (_disc.strideBound + _disc.nComp));
			}

			return lms.bufferSize();
		}

		unsigned int RadialLumpedRateModelWithoutPoresDG::requiredADdirs() const CADET_NOEXCEPT
		{
			return _jacobianAdDirs;
		}

		void RadialLumpedRateModelWithoutPoresDG::useAnalyticJacobian(const bool analyticJac)
		{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
			_analyticJac = analyticJac;
			if (!_analyticJac)
				_jacobianAdDirs = _convDispOp.requiredADdirs();
			else
				_jacobianAdDirs = 0;
#else
			_analyticJac = false;
			_jacobianAdDirs = _convDispOp.requiredADdirs();
#endif
		}

		void RadialLumpedRateModelWithoutPoresDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac)
		{
			_convDispOp.notifyDiscontinuousSectionTransition(t, secIdx, _jacInlet);

			updateSection(secIdx);
		}

		void RadialLumpedRateModelWithoutPoresDG::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
		{
			_convDispOp.setFlowRates(in[0], out[0], _totalPorosity);
		}

		void RadialLumpedRateModelWithoutPoresDG::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
		{
			Exporter expr(_disc, *this, solution);
			recorder.beginUnitOperation(_unitOpIdx, *this, expr);
			recorder.endUnitOperation();
		}

		void RadialLumpedRateModelWithoutPoresDG::reportSolutionStructure(ISolutionRecorder& recorder) const
		{
			Exporter expr(_disc, *this, nullptr);
			recorder.unitOperationStructure(_unitOpIdx, *this, expr);
		}

		int RadialLumpedRateModelWithoutPoresDG::Exporter::writeMobilePhase(double* buffer) const
		{
			const int blockSize = numMobilePhaseDofs();
			std::copy_n(&_data[_idx.offsetC()], blockSize, buffer);
			return blockSize;
		}

		int RadialLumpedRateModelWithoutPoresDG::Exporter::writeSolidPhase(double* buffer) const
		{
			int numWritten = 0;
			for (unsigned int i = 0; i < _disc.nPoints; ++i)
			{
				const int localOffset = _idx.offsetC() + i * _idx.strideColNode() + _idx.strideColLiquid();
				std::copy_n(&_data[localOffset], _disc.strideBound, buffer);
				buffer += _disc.strideBound;
				numWritten += _disc.strideBound;
			}
			return numWritten;
		}

		int RadialLumpedRateModelWithoutPoresDG::Exporter::writeSolidPhase(unsigned int parType, double* buffer) const
		{
			cadet_assert(parType == 0);
			return writeSolidPhase(buffer);
		}

		int RadialLumpedRateModelWithoutPoresDG::Exporter::writeInlet(unsigned int port, double* buffer) const
		{
			cadet_assert(port == 0);
			std::copy_n(_data, _disc.nComp, buffer);
			return _disc.nComp;
		}

		int RadialLumpedRateModelWithoutPoresDG::Exporter::writeInlet(double* buffer) const
		{
			std::copy_n(_data, _disc.nComp, buffer);
			return _disc.nComp;
		}

		int RadialLumpedRateModelWithoutPoresDG::Exporter::writeOutlet(unsigned int port, double* buffer) const
		{
			cadet_assert(port == 0);

			if (_model._convDispOp.forwardFlow())
			{
				// Forward flow: Outlet is at the end (last node of last element)
				const unsigned int offset = _idx.offsetC() + (_disc.nPoints - 1) * _idx.strideColNode();
				for (unsigned int i = 0; i < _disc.nComp; ++i)
					buffer[i] = _data[offset + i * _idx.strideColComp()];
			}
			else
			{
				// Backward flow: Outlet is at the beginning (first node)
				for (unsigned int i = 0; i < _disc.nComp; ++i)
					buffer[i] = _data[_idx.offsetC() + i * _idx.strideColComp()];
			}

			return _disc.nComp;
		}

		int RadialLumpedRateModelWithoutPoresDG::Exporter::writeOutlet(double* buffer) const
		{
			return writeOutlet(0, buffer);
		}

		void RadialLumpedRateModelWithoutPoresDG::prepareADvectors(const AdJacobianParams& adJac) const
		{
			Indexer idxr(_disc);
			// Radial DG always uses exact integration, so bandwidth is 2 * nNodes * stride
			int lowerBandwidth = 2 * _disc.nNodes * idxr.strideColNode();
			int upperBandwidth = lowerBandwidth;

			ad::prepareAdVectorSeedsForBandMatrix(adJac.adY + _disc.nComp, adJac.adDirOffset, _jac.rows(), lowerBandwidth, upperBandwidth, lowerBandwidth);
		}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

		void RadialLumpedRateModelWithoutPoresDG::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
		{
			Indexer idxr(_disc);

			LOG(Debug) << "AD dir offset: " << adDirOffset << " DiagBlockSize: " << numPureDofs() << " numPureDofs: " << numPureDofs();

			const double* const adVec = adRes + idxr.offsetC();

			for (int row = 0; row < _jac.rows(); row++)
			{
				for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(_jac, row); it; ++it)
				{
					const int col = it.col();

					const double relErr = std::abs(it.value() - adVec[row].getADValue(adDirOffset + col)) / std::max(1.0, std::max(std::abs(it.value()), std::abs(adVec[row].getADValue(adDirOffset + col))));

					if (relErr > 1e-6)
					{
						LOG(Debug) << "MISMATCH (" << row << ", " << col << "): " << it.value() << " vs AD " << adVec[row].getADValue(adDirOffset + col) << " relErr: " << relErr;
					}
				}
			}
		}

#endif

		int RadialLumpedRateModelWithoutPoresDG::jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerResidual);

			_factorizeJacobian = true;

			if (_analyticJac)
				return residualImpl<double, double, double, true, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
			else
				return residualWithJacobian(simTime, ConstSimulationState{ simState.vecStateY, nullptr }, nullptr, adJac, threadLocalMem);
		}

		int RadialLumpedRateModelWithoutPoresDG::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerResidual);

			return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
		}

		int RadialLumpedRateModelWithoutPoresDG::residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerResidual);

			return residual(simTime, simState, res, adJac, threadLocalMem, true, false);
		}

		int RadialLumpedRateModelWithoutPoresDG::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem, bool updateJacobian, bool paramSensitivity)
		{
			if (updateJacobian)
			{
				_factorizeJacobian = true;

#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
				if (_analyticJac)
				{
					if (paramSensitivity)
					{
						const int retCode = residualImpl<double, active, active, true>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

						if (res)
							ad::copyFromAd(adJac.adRes, res, numDofs());

						return retCode;
					}
					else
						return residualImpl<double, double, double, true>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
				}
				else
				{
					ad::copyToAd(simState.vecStateY, adJac.adY, numDofs());
					ad::resetAd(adJac.adRes, numDofs());

					int retCode = 0;
					if (paramSensitivity)
						retCode = residualImpl<active, active, active, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);
					else
						retCode = residualImpl<active, active, double, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

					if (res)
						ad::copyFromAd(adJac.adRes, res, numDofs());

					extractJacobianFromAD(adJac.adRes, adJac.adDirOffset);

					return retCode;
				}
#else

				ad::copyToAd(simState.vecStateY, adJac.adY, numDofs());
				ad::resetAd(adJac.adRes, numDofs());

				int retCode = 0;
				if (paramSensitivity)
					retCode = residualImpl<active, active, active, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);
				else
					retCode = residualImpl<active, active, double, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

				if (res)
				{
					retCode = residualImpl<double, double, double, true>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);

					checkAnalyticJacobianAgainstAd(adJac.adRes, adJac.adDirOffset);
				}

				extractJacobianFromAD(adJac.adRes, adJac.adDirOffset);

				return retCode;
#endif
			}
			else
			{
				if (paramSensitivity)
				{
					ad::resetAd(adJac.adRes, numDofs());

					const int retCode = residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

					if (res)
						ad::copyFromAd(adJac.adRes, res, numDofs());

					return retCode;
				}
				else
					return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
			}
		}

		template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
		int RadialLumpedRateModelWithoutPoresDG::residualImpl(double t, unsigned int secIdx, StateType const* const y_, double const* const yDot_, ResidualType* const res_, util::ThreadLocalStorage& threadLocalMem)
		{
			Indexer idxr(_disc);

			bool success = 1;

			if (wantJac) {

				if (!wantRes || _disc.newStaticJac) {

					success = _convDispOp.calcTransportJacobian(_jac, _jacInlet);

					_disc.newStaticJac = false;
				}

				if (cadet_unlikely(!success))
					LOG(Error) << "Jacobian pattern did not fit the Jacobian estimation";

			}

			/*	Compute bulk convection dispersion residual	*/
			if (wantRes)
				_convDispOp.residual(*this, t, secIdx, y_, yDot_, res_, typename cadet::ParamSens<ParamType>::enabled());

			/* Compute binding, reaction residual */

#ifdef CADET_PARALLELIZE
			tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nPoints), [&](std::size_t blk)
#else
			for (unsigned int blk = 0; blk < _disc.nPoints; ++blk)
#endif
			{
				linalg::BandedEigenSparseRowIterator jacIt(_jac, blk * idxr.strideColNode());
				StateType const* const localY = y_ + idxr.offsetC() + idxr.strideColNode() * blk;

				const parts::cell::CellParameters cellResParams
				{
					_disc.nComp,
					_disc.nBound,
					_disc.boundOffset,
					_disc.strideBound,
					_binding[0]->reactionQuasiStationarity(),
					_totalPorosity,
					nullptr,
					_binding[0],
					_reaction.hasReactions() ? &_reaction : nullptr
				};

				// Relative position of current node - needed in externally dependent adsorption kinetic
				double z = _convDispOp.relativeCoordinate(blk);
				if (wantRes)
				{
					ResidualType* const localRes = res_ + idxr.offsetC() + idxr.strideColNode() * blk;
					double const* const localYdot = yDot_ ? yDot_ + idxr.offsetC() + idxr.strideColNode() * blk : nullptr;
					parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantJac, false>(
						t, secIdx, ColumnPosition{ z, 0.0, 0.0 }, localY, localYdot, localRes, jacIt, cellResParams, threadLocalMem.get()
					);
				}
				else
				{
					parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantJac, false, false>(
						t, secIdx, ColumnPosition{ z, 0.0, 0.0 }, localY, nullptr, nullptr, jacIt, cellResParams, threadLocalMem.get()
					);
				}
			} CADET_PARFOR_END;

			// Handle inlet DOFs, which are simply copied to res
			for (unsigned int i = 0; i < _disc.nComp; ++i)
				res_[i] = y_[i];

			return 0;
		}

		int RadialLumpedRateModelWithoutPoresDG::residualSensFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerResidualSens);

			return residual(simTime, simState, nullptr, adJac, threadLocalMem, true, true);
		}

		int RadialLumpedRateModelWithoutPoresDG::residualSensFwdAdOnly(const SimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerResidualSens);

			return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adRes, threadLocalMem);
		}

		int RadialLumpedRateModelWithoutPoresDG::residualSensFwdCombine(const SimulationTime& simTime, const ConstSimulationState& simState,
			const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes,
			double* const tmp1, double* const tmp2, double* const tmp3)
		{
			BENCH_SCOPE(_timerResidualSens);

			for (std::size_t param = 0; param < yS.size(); ++param)
			{
				multiplyWithJacobian(SimulationTime{ 0.0, 0u }, ConstSimulationState{ nullptr, nullptr }, yS[param], 1.0, 0.0, tmp1);
				multiplyWithDerivativeJacobian(SimulationTime{ 0.0, 0u }, ConstSimulationState{ nullptr, nullptr }, ySdot[param], tmp2);

				double* const ptrResS = resS[param];

				BENCH_START(_timerResidualSensPar);

#ifdef CADET_PARALLELIZE
				tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(numDofs()), [&](std::size_t i)
#else
				for (unsigned int i = 0; i < numDofs(); ++i)
#endif
				{
					ptrResS[i] = tmp1[i] + tmp2[i] + adRes[i].getADValue(param);
				} CADET_PARFOR_END;

				BENCH_STOP(_timerResidualSensPar);
			}

			return 0;
		}

		void RadialLumpedRateModelWithoutPoresDG::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
		{
			Indexer idxr(_disc);

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				ret[comp] = alpha * yS[comp] + beta * ret[comp];
			}

			Eigen::Map<Eigen::VectorXd> ret_vec(ret + idxr.offsetC(), numPureDofs());
			Eigen::Map<const Eigen::VectorXd> yS_vec(yS + idxr.offsetC(), numPureDofs());
			ret_vec = alpha * _jac * yS_vec + beta * ret_vec;

			// Map inlet DOFs to the column inlet (first bulk cells for forward flow)
			unsigned int offInlet = _convDispOp.forwardFlow() ? 0 : (_disc.nElem - 1u) * idxr.strideColCell();

			for (unsigned int comp = 0; comp < _disc.nComp; comp++) {
				for (unsigned int node = 0; node < _disc.nNodes; node++) {
					ret[idxr.offsetC() + offInlet + comp * idxr.strideColComp() + node * idxr.strideColNode()] += alpha * _jacInlet(node, 0) * yS[comp];
				}
			}
		}

		void RadialLumpedRateModelWithoutPoresDG::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
		{
			Indexer idxr(_disc);
			const double invBeta = (1.0 / static_cast<double>(_totalPorosity) - 1.0);

			_convDispOp.multiplyWithDerivativeJacobian(simTime, sDot, ret);

			for (unsigned int node = 0; node < _disc.nPoints; ++node)
			{
				const unsigned int localOffset = idxr.offsetC() + node * idxr.strideColNode();
				double const* const localSdot = sDot + localOffset;
				double* const localRet = ret + localOffset;

				parts::cell::multiplyWithDerivativeJacobianKernel<false>(localSdot, localRet, _disc.nComp, _disc.nBound, _disc.boundOffset, _disc.strideBound, _binding[0]->reactionQuasiStationarity(), 1.0, invBeta);
			}

			std::fill_n(ret, _disc.nComp, 0.0);
		}

		void RadialLumpedRateModelWithoutPoresDG::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
		{
			if (_binding[0])
				_binding[0]->setExternalFunctions(extFuns, size);
		}

		unsigned int RadialLumpedRateModelWithoutPoresDG::localOutletComponentIndex(unsigned int port) const CADET_NOEXCEPT
		{
			if (static_cast<double>(_convDispOp.currentVelocity()) >= 0.0)
				return _disc.nComp + (_disc.nPoints - 1) * (_disc.nComp + _disc.strideBound);
			else
				return _disc.nComp;
		}

		unsigned int RadialLumpedRateModelWithoutPoresDG::localInletComponentIndex(unsigned int port) const CADET_NOEXCEPT
		{
			return 0;
		}

		unsigned int RadialLumpedRateModelWithoutPoresDG::localOutletComponentStride(unsigned int port) const CADET_NOEXCEPT
		{
			return 1;
		}

		unsigned int RadialLumpedRateModelWithoutPoresDG::localInletComponentStride(unsigned int port) const CADET_NOEXCEPT
		{
			return 1;
		}

		void RadialLumpedRateModelWithoutPoresDG::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
		{
			// @todo Write this function
		}

		int RadialLumpedRateModelWithoutPoresDG::linearSolve(double t, double alpha, double outerTol, double* const rhs, double const* const weight,
			const ConstSimulationState& simState)
		{
			BENCH_SCOPE(_timerLinearSolve);

			Indexer idxr(_disc);

			bool success = true;
			bool result = true;

			assembleDiscretizedJacobian(alpha, idxr);

			Eigen::Map<VectorXd> r(rhs, numDofs());

			if (_factorizeJacobian)
			{
				_linearSolver->factorize(_jacDisc);

				if (_linearSolver->info() != Success) {
					LOG(Error) << "factorization failed";
					success = false;
				}
			}

			r.segment(idxr.offsetC(), numPureDofs()) = _linearSolver->solve(r.segment(idxr.offsetC(), numPureDofs()));

			if (_linearSolver->info() != Success) {
				LOG(Error) << "solve() failed";
				result = false;
			}

			unsigned int offInlet = _convDispOp.forwardFlow() ? 0 : (_disc.nElem - 1u) * idxr.strideColCell();

			for (unsigned int comp = 0; comp < _disc.nComp; comp++) {
				for (unsigned int node = 0; node < _disc.nNodes; node++) {
					r[idxr.offsetC() + offInlet + comp * idxr.strideColComp() + node * idxr.strideColNode()] -= _jacInlet(node, 0) * r[comp];
				}
			}

			return (success && result) ? 0 : 1;
		}

		void RadialLumpedRateModelWithoutPoresDG::assembleDiscretizedJacobian(double alpha, const Indexer& idxr)
		{
			_jacDisc = _jac;

			_convDispOp.addTimeDerivativeToJacobian(alpha, _jacDisc);

			const double invBeta = 1.0 / static_cast<double>(_totalPorosity) - 1.0;
			linalg::BandedEigenSparseRowIterator jac(_jacDisc, 0);
			for (unsigned int j = 0; j < _disc.nPoints; ++j)
			{
				addTimeDerivativeToJacobianNode(jac, idxr, alpha, invBeta);
			}
		}

		void RadialLumpedRateModelWithoutPoresDG::addTimeDerivativeToJacobianNode(linalg::BandedEigenSparseRowIterator& jac, const Indexer& idxr, double alpha, double invBeta) const
		{
			for (int comp = 0; comp < static_cast<int>(_disc.nComp); ++comp, ++jac)
			{
				for (int i = 0; i < static_cast<int>(_disc.nBound[comp]); ++i)
				{
					jac[idxr.strideColLiquid() - comp + idxr.offsetBoundComp(comp) + i] += alpha * invBeta;
				}
			}

			int const* const qsReaction = _binding[0]->reactionQuasiStationarity();
			for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd, ++jac)
			{
				if (qsReaction[bnd])
					continue;

				jac[0] += alpha;
			}
		}

		void RadialLumpedRateModelWithoutPoresDG::applyInitialCondition(const SimulationState& simState) const
		{
			Indexer idxr(_disc);

			if (!_initState.empty())
			{
				std::fill(simState.vecStateY, simState.vecStateY + idxr.offsetC(), 0.0);
				std::copy(_initState.data(), _initState.data() + numPureDofs(), simState.vecStateY + idxr.offsetC());

				if (!_initStateDot.empty())
				{
					std::fill(simState.vecStateYdot, simState.vecStateYdot + idxr.offsetC(), 0.0);
					std::copy(_initStateDot.data(), _initStateDot.data() + numPureDofs(), simState.vecStateYdot + idxr.offsetC());
				}

				return;
			}

			double* const stateYbulk = simState.vecStateY + idxr.offsetC();

			for (unsigned int point = 0; point < _disc.nPoints; ++point)
			{
				double* const stateYbulkNode = stateYbulk + point * idxr.strideColNode();
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					stateYbulkNode[comp * idxr.strideColComp()] = static_cast<double>(_initC[comp]);
				}

				double* const stateYbulkSolid = stateYbulkNode + idxr.strideColLiquid();
				for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
					stateYbulkSolid[bnd] = static_cast<double>(_initCs[bnd]);
			}
		}

		void RadialLumpedRateModelWithoutPoresDG::readInitialCondition(IParameterProvider& paramProvider)
		{
			_initState.clear();
			_initStateDot.clear();

			if (paramProvider.exists("INIT_STATE"))
			{
				const std::vector<double> initState = paramProvider.getDoubleArray("INIT_STATE");
				_initState = std::vector<double>(initState.begin(), initState.begin() + numPureDofs());

				if (paramProvider.exists("INIT_STATE_DOT"))
				{
					const std::vector<double> initStateDot = paramProvider.getDoubleArray("INIT_STATE_DOT");
					_initStateDot = std::vector<double>(initStateDot.begin(), initStateDot.begin() + numPureDofs());
				}
				return;
			}

			const std::vector<double> initC = paramProvider.getDoubleArray("INIT_C");

			if (initC.size() < _disc.nComp)
				throw InvalidParameterException("INIT_C does not contain enough values for all components");

			ad::copyToAd(initC.data(), _initC.data(), _disc.nComp);

			if (!_binding[0] || (_disc.strideBound == 0))
				return;

			paramProvider.pushScope("particle_type_000");

			if (paramProvider.exists("INIT_Q"))
			{
				const std::vector<double> initQ = paramProvider.getDoubleArray("INIT_Q");

				if (initQ.size() < _disc.strideBound)
					throw InvalidParameterException("INIT_Q does not contain enough values");

				ad::copyToAd(initQ.data(), _initCs.data(), _disc.strideBound);
			}

			paramProvider.popScope();
		}

		void RadialLumpedRateModelWithoutPoresDG::consistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerConsistentInit);

			Indexer idxr(_disc);

			// Step 1: Solve algebraic equations
			if (!_binding[0]->hasQuasiStationaryReactions())
				return;

			double* const c = vecStateY + idxr.offsetC();

			for (unsigned int point = 0; point < _disc.nPoints; ++point)
			{
				LinearBufferAllocator tlmAlloc = threadLocalMem.get();

				const unsigned int localOffset = point * idxr.strideColNode();
				double* const localC = c + localOffset;
				double* const localQ = localC + idxr.strideColLiquid();

				const double rho = _convDispOp.relativeCoordinate(point);
				const ColumnPosition colPos{rho, 0.0, 0.0};

				// Use pre/post consistent initial state from binding model
				_binding[0]->preConsistentInitialState(simTime.t, simTime.secIdx, colPos, localQ, localC, tlmAlloc);
				_binding[0]->postConsistentInitialState(simTime.t, simTime.secIdx, colPos, localQ, localC, tlmAlloc);
			}
		}

		void RadialLumpedRateModelWithoutPoresDG::consistentInitialTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerConsistentInit);

			Indexer idxr(_disc);

			double* const cDot = vecStateYdot + idxr.offsetC();
			double* const qDot = cDot + idxr.strideColLiquid();

			double const* const c = vecStateY + idxr.offsetC();
			double const* const q = c + idxr.strideColLiquid();

			LinearBufferAllocator tlmAlloc = threadLocalMem.get();

			for (unsigned int point = 0; point < _disc.nPoints; ++point)
			{
				const unsigned int localOffset = point * idxr.strideColNode();
				double const* const localC = c + localOffset;
				double const* const localQ = q + localOffset;
				double* const localCdot = cDot + localOffset;
				double* const localQdot = qDot + localOffset;

				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					localCdot[comp] = 0.0;

				if (!_binding[0]->hasQuasiStationaryReactions())
				{
					std::fill(localQdot, localQdot + _disc.strideBound, 0.0);
					continue;
				}

				const ColumnPosition colPos{_convDispOp.relativeCoordinate(point), 0.0, 0.0};

				// timeDerivativeQuasiStationaryFluxes(t, secIdx, colPos, yCp, y, dResDt, workSpace)
				_binding[0]->timeDerivativeQuasiStationaryFluxes(simTime.t, simTime.secIdx, colPos, localC, localQ, localQdot, tlmAlloc);
			}
		}

		void RadialLumpedRateModelWithoutPoresDG::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
		{
			Indexer idxr(_disc);

			for (std::size_t param = 0; param < vecSensY.size(); ++param)
			{
				double* const stateYbulk = vecSensY[param] + idxr.offsetC();

				std::fill(stateYbulk, stateYbulk + numPureDofs(), 0.0);
			}
		}

		void RadialLumpedRateModelWithoutPoresDG::consistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
			std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
		{
			// For now, initialize sensitivities to zero
			// TODO: Implement full sensitivity calculation for quasi-stationary binding
			Indexer idxr(_disc);

			for (std::size_t param = 0; param < vecSensY.size(); ++param)
			{
				double* const sensYdot = vecSensYdot[param];
				double* const sensCdot = sensYdot + idxr.offsetC();
				double* const sensQdot = sensCdot + idxr.strideColLiquid();

				for (unsigned int point = 0; point < _disc.nPoints; ++point)
				{
					const unsigned int localOffset = point * idxr.strideColNode();
					double* const localSensCdot = sensCdot + localOffset;
					double* const localSensQdot = sensQdot + localOffset;

					for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
						localSensCdot[comp] = 0.0;

					std::fill(localSensQdot, localSensQdot + _disc.strideBound, 0.0);
				}
			}
		}

		void RadialLumpedRateModelWithoutPoresDG::leanConsistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
		{
			consistentInitialState(simTime, vecStateY, adJac, errorTol, threadLocalMem);
		}

		void RadialLumpedRateModelWithoutPoresDG::leanConsistentInitialTimeDerivative(double t, double const* const vecStateY, double* const vecStateYdot, double* const res, util::ThreadLocalStorage& threadLocalMem)
		{
			Indexer idxr(_disc);

			double* const cDot = vecStateYdot + idxr.offsetC();
			double* const qDot = cDot + idxr.strideColLiquid();

			for (unsigned int point = 0; point < _disc.nPoints; ++point)
			{
				const unsigned int localOffset = point * idxr.strideColNode();
				double* const localCdot = cDot + localOffset;
				double* const localQdot = qDot + localOffset;

				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					localCdot[comp] = 0.0;

				std::fill(localQdot, localQdot + _disc.strideBound, 0.0);
			}
		}

		void RadialLumpedRateModelWithoutPoresDG::leanConsistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
			std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
		{
			consistentInitialSensitivity(simTime, simState, vecSensY, vecSensYdot, adRes, threadLocalMem);
		}

		bool RadialLumpedRateModelWithoutPoresDG::setParameter(const ParameterId& pId, double value)
		{
			if (pId.unitOperation == _unitOpIdx)
			{
				if (_convDispOp.setParameter(pId, value))
					return true;
			}

			const bool found = UnitOperationBase::setParameter(pId, value);
			if (!found && (pId.unitOperation == _unitOpIdx))
			{
				if (_binding[0] && _binding[0]->setParameter(pId, value))
					return true;
				if (_dynReaction[0] && _dynReaction[0]->setParameter(pId, value))
					return true;
			}

			return found;
		}

		void RadialLumpedRateModelWithoutPoresDG::setSensitiveParameterValue(const ParameterId& pId, double value)
		{
			if (pId.unitOperation == _unitOpIdx)
			{
				if (_convDispOp.setSensitiveParameterValue(_sensParams, pId, value))
					return;
			}

			UnitOperationBase::setSensitiveParameterValue(pId, value);
		}

		bool RadialLumpedRateModelWithoutPoresDG::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
		{
			if (pId.unitOperation == _unitOpIdx)
			{
				if (_convDispOp.setSensitiveParameter(_sensParams, pId, adDirection, adValue))
					return true;
			}

			return UnitOperationBase::setSensitiveParameter(pId, adDirection, adValue);
		}

		void RadialLumpedRateModelWithoutPoresDG::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
		{
			Indexer idxr(_disc);
			const int nDOFs = numPureDofs();
			const double* const adVec = reinterpret_cast<const double*>(adRes) + idxr.offsetC();

			for (int row = 0; row < _jac.rows(); row++)
			{
				for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(_jac, row); it; ++it)
				{
					const int col = it.col();
					it.valueRef() = adVec[row * (adDirOffset + nDOFs + 1) + adDirOffset + col];
				}
			}
		}

	}  // namespace model

}  // namespace cadet

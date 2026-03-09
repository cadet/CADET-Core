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

#include "model/RadialLumpedRateModelWithPoresDG.hpp"
#include "BindingModelFactory.hpp"
#include "ParamReaderHelper.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/ExternalFunction.hpp"
#include "cadet/SolutionRecorder.hpp"
#include "ConfigurationHelper.hpp"
#include "model/BindingModel.hpp"
#include "model/ParameterDependence.hpp"
#include "SimulationTypes.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/Norms.hpp"
#include "linalg/Subset.hpp"
#include "model/parts/BindingCellKernel.hpp"
#include "model/parts/DGToolbox.hpp"

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

		RadialLumpedRateModelWithPoresDG::RadialLumpedRateModelWithPoresDG(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
			_filmDiffDep(nullptr), _variableFilmDiff(false),
			_analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr),
			_initC(0), _initCp(0), _initCs(0), _initState(0), _initStateDot(0)
		{
			// Multiple particle types are not supported
			_singleBinding = true;
			_singleDynReaction = true;
		}

		RadialLumpedRateModelWithPoresDG::~RadialLumpedRateModelWithPoresDG() CADET_NOEXCEPT
		{
			delete[] _tempState;
			delete _linearSolver;
			delete _filmDiffDep;
			for (IDynamicReactionModel* p : _dynReaction) delete p;
		}

		unsigned int RadialLumpedRateModelWithPoresDG::numDofs() const CADET_NOEXCEPT
		{
			// Inlet DOFs: nComp
			// Bulk DOFs: nPoints * nComp (mobile phase c)
			// Pore DOFs: nPoints * (nComp + strideBound) (pore liquid cp + solid q)
			return _disc.nComp + _disc.nPoints * _disc.nComp + _disc.nPoints * (_disc.nComp + _disc.strideBound);
		}

		unsigned int RadialLumpedRateModelWithPoresDG::numPureDofs() const CADET_NOEXCEPT
		{
			// Bulk DOFs: nPoints * nComp
			// Pore DOFs: nPoints * (nComp + strideBound)
			return _disc.nPoints * _disc.nComp + _disc.nPoints * (_disc.nComp + _disc.strideBound);
		}


		bool RadialLumpedRateModelWithPoresDG::usesAD() const CADET_NOEXCEPT
		{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
			// We always need AD if we want to check the analytical Jacobian
			return true;
#else
			// We only need AD if we are not computing the Jacobian analytically
			return !_analyticJac;
#endif
		}

		bool RadialLumpedRateModelWithPoresDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper)
		{
			const bool firstConfigCall = _tempState == nullptr;

			// Read discretization
			_disc.nComp = paramProvider.getInt("NCOMP");

			if (paramProvider.exists("NPARTYPE") ? paramProvider.getInt("NPARTYPE") != 1 : false)
				throw InvalidParameterException("Number of particle types must be 1 for RADIAL_LUMPED_RATE_MODEL_WITH_PORES_DG");

			paramProvider.pushScope("particle_type_000");

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
			_initCp.resize(_disc.nElem * _disc.nNodes * _disc.nComp);
			_initCs.resize(_disc.nElem * _disc.nNodes * _disc.strideBound);

			// Create nonlinear solver for consistent initialization
			configureNonlinearSolver(paramProvider);

			paramProvider.popScope();

			// Configure convection-dispersion operator (bulk phase stride is just nComp)
			const unsigned int strideNode = _disc.nComp;
			const bool transportSuccess = _convDispOp.configureModelDiscretization(paramProvider, helper, _disc.nComp, polynomial_integration_mode, _disc.nElem, _disc.polyDeg, strideNode);

			// Configure film diffusion parameter dependency
			if (paramProvider.exists("film_diffusion"))
			{
				paramProvider.pushScope("film_diffusion");
				if (paramProvider.exists("FILM_DIFFUSION_DEP"))
				{
					const std::string paramDepName = paramProvider.getString("FILM_DIFFUSION_DEP");
					_filmDiffDep = helper.createParameterParameterDependence(paramDepName);
					if (!_filmDiffDep)
						throw InvalidParameterException("Unknown parameter dependence " + paramDepName + " in FILM_DIFFUSION_DEP");

					_filmDiffDep->configureModelDiscretization(paramProvider);
					_variableFilmDiff = true;
				}
				else
				{
					_filmDiffDep = helper.createParameterParameterDependence("CONSTANT_ONE");
					_variableFilmDiff = false;
				}
				paramProvider.popScope();
			}
			else
			{
				_filmDiffDep = helper.createParameterParameterDependence("CONSTANT_ONE");
				_variableFilmDiff = false;
			}

			_disc.curSection = -1;

			// Allocate memory
			Indexer idxr(_disc);

			// Allocate film diffusion matrices
			_M_K.resize(_disc.nElem);
			_filmDiffCoupling.resize(_disc.nElem);
			_filmDiffAtNodes.resize(_disc.nElem);
			for (unsigned int cell = 0; cell < _disc.nElem; ++cell)
			{
				_M_K[cell].resize(_disc.nNodes, _disc.nNodes);
				_filmDiffCoupling[cell].resize(_disc.nNodes, _disc.nNodes);
				_filmDiffAtNodes[cell].resize(_disc.nNodes);
			}

			// Jacobian allocation
			_jacInlet.resize(_disc.nNodes, 1);

			// Total Jacobian size: bulk + pore DOFs
			const unsigned int jacSize = _disc.nPoints * _disc.nComp + _disc.nPoints * (_disc.nComp + _disc.strideBound);
			_jac.resize(jacSize, jacSize);
			_jacDisc.resize(jacSize, jacSize);

			// Set whether analytic Jacobian is used
			useAnalyticJacobian(analyticJac);

			// ==== Construct and configure binding model

			clearBindingModels();
			_binding.push_back(nullptr);

			paramProvider.pushScope("particle_type_000");

			if (paramProvider.exists("adsorption_model"))
			{
				bool bindingConfSuccess = true;
				paramProvider.pushScope("adsorption");

				const std::string bindModelName = paramProvider.getString("ADSORPTION_MODEL");
				_binding[0] = helper.createBindingModel(bindModelName);
				if (!_binding[0])
					throw InvalidParameterException("Unknown binding model " + bindModelName);

				bindingConfSuccess = _binding[0]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound, _disc.boundOffset);

				paramProvider.popScope();

				if (!bindingConfSuccess)
					throw InvalidParameterException("Failed to configure binding model for particle type 000");
			}
			else
			{
				_binding[0] = helper.createBindingModel("NONE");
				_binding[0]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound, _disc.boundOffset);
			}

			// ==== Configure dynamic reaction model
			clearDynamicReactionModels();
			_reaction.getDynReactionVector("liquid").clear();
			for (IDynamicReactionModel* p : _dynReaction) delete p;
			_dynReaction.clear();
			_dynReaction.push_back(nullptr);

			if (paramProvider.exists("reaction_model"))
			{
				paramProvider.pushScope("reaction");

				const std::string dynReactionName = paramProvider.getString("REACTION_MODEL");
				_dynReaction[0] = helper.createDynamicReactionModel(dynReactionName);
				if (!_dynReaction[0])
					throw InvalidParameterException("Unknown dynamic reaction model " + dynReactionName);

				if (!_dynReaction[0]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound, _disc.boundOffset))
					throw InvalidParameterException("Failed to configure dynamic reaction model");

				paramProvider.popScope();
			}

			// Sync _reaction for CellParameters (non-owning refs)
			if (_dynReaction[0]) _reaction.getDynReactionVector("liquid").push_back(_dynReaction[0]);

			paramProvider.popScope(); // particle_type_000

			// Allocate temporaries
			const unsigned int maxStride = std::max({
				_disc.nComp + _disc.strideBound,
				_binding[0] ? _binding[0]->workspaceSize(_disc.nComp, _disc.strideBound, _disc.nBound) : 0u,
				_dynReaction[0] ? _dynReaction[0]->workspaceSize(_disc.nComp, _disc.strideBound, _disc.nBound) : 0u
				});

			const unsigned int nDof = numDofs();
			if (firstConfigCall || (_tempState == nullptr))
			{
				delete[] _tempState;
				_tempState = new double[nDof + _disc.nPoints * maxStride];
			}

			return transportSuccess;
		}

		bool RadialLumpedRateModelWithPoresDG::configure(IParameterProvider& paramProvider)
		{
			_parameters.clear();

			// Read geometry parameters
			_colPorosity = paramProvider.getDouble("COL_POROSITY");

			// Register column porosity
			_parameters[makeParamId(hashString("COL_POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colPorosity;

			// Read particle parameters
			paramProvider.pushScope("particle_type_000");

			_parPorosity = paramProvider.getDouble("PAR_POROSITY");
			_parRadius = paramProvider.getDouble("PAR_RADIUS");
			_parGeomSurfToVol = paramProvider.exists("PAR_GEOM_SURF_TO_VOL") ? paramProvider.getDouble("PAR_GEOM_SURF_TO_VOL") : 3.0;

			_parameters[makeParamId(hashString("PAR_POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parPorosity;
			_parameters[makeParamId(hashString("PAR_RADIUS"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parRadius;

			// Read film diffusion
			readScalarParameterOrArray(_filmDiffusion, paramProvider, "FILM_DIFFUSION", _disc.nComp);
			_filmDiffusionMode = readAndRegisterMultiplexCompSecParam(paramProvider, _parameters, _filmDiffusion, "FILM_DIFFUSION", _disc.nComp, 0, false, _unitOpIdx);

			// Read pore accessibility factor
			if (paramProvider.exists("PORE_ACCESSIBILITY"))
			{
				readScalarParameterOrArray(_poreAccessFactor, paramProvider, "PORE_ACCESSIBILITY", _disc.nComp);
				_poreAccessFactorMode = readAndRegisterMultiplexCompSecParam(paramProvider, _parameters, _poreAccessFactor, "PORE_ACCESSIBILITY", _disc.nComp, 0, false, _unitOpIdx);
			}
			else
			{
				_poreAccessFactor.resize(_disc.nComp, 1.0);
				_poreAccessFactorMode = MultiplexMode::Independent;
			}

			paramProvider.popScope();

			// Configure binding model
			if (_binding[0] && _binding[0]->requiresConfiguration())
			{
				paramProvider.pushScope("particle_type_000");
				paramProvider.pushScope("adsorption");

				const bool bindSuccess = _binding[0]->configure(paramProvider, _unitOpIdx, ParTypeIndep);

				paramProvider.popScope();
				paramProvider.popScope();

				if (!bindSuccess)
					throw InvalidParameterException("Failed to configure binding model");
			}

			// Configure dynamic reaction model
			if (_dynReaction[0])
			{
				paramProvider.pushScope("particle_type_000");
				paramProvider.pushScope("reaction");

				const bool dynReactSuccess = _dynReaction[0]->configure(paramProvider, _unitOpIdx, ParTypeIndep);

				paramProvider.popScope();
				paramProvider.popScope();

				if (!dynReactSuccess)
					throw InvalidParameterException("Failed to configure dynamic reaction model");
			}

			// Configure film diffusion dependency
			if (_filmDiffDep && _variableFilmDiff && paramProvider.exists("film_diffusion"))
			{
				paramProvider.pushScope("film_diffusion");
				if (!_filmDiffDep->configure(paramProvider, _unitOpIdx, ParTypeIndep, BoundStateIndep, "FILM_DIFFUSION_DEP"))
					throw InvalidParameterException("Failed to configure film diffusion dependency");
				paramProvider.popScope();
			}

			// Configure convection-dispersion operator
			const bool transportSuccess = _convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

			// Compute initial film diffusion matrices (constant k_f case)
			computeFilmDiffusionMatrices();

			// Set Jacobian pattern
			setPattern(_jac, false, _dynReaction[0] && (_dynReaction[0]->numReactionsCombined() > 0));
			setPattern(_jacDisc, true, _dynReaction[0] && (_dynReaction[0]->numReactionsCombined() > 0));

			return transportSuccess;
		}

		void RadialLumpedRateModelWithPoresDG::computeFilmDiffusionMatrices()
		{
			const unsigned int nNodes = _disc.nNodes;
			const double deltaRho = static_cast<double>(_convDispOp.columnLength()) / _disc.nElem;

			// Get film diffusion coefficient (use first component for now)
			const active kf = _filmDiffusion.size() > 0 ? _filmDiffusion[0] : active(1.0);

			for (unsigned int cell = 0; cell < _disc.nElem; ++cell)
			{
				const double rho_left = _convDispOp.elemLeftBound(cell);

				if (!_variableFilmDiff)
				{
					// Constant k_f: M_K = k_f * M_rho
					_M_K[cell] = static_cast<double>(kf) * _convDispOp.MRho(cell);
				}
				else
				{
					// Variable k_f: use radialFilmDiffusionMatrix from DGToolbox
					Eigen::VectorXd kfAtNodes(_disc.nNodes);
					for (unsigned int node = 0; node < nNodes; ++node)
					{
						kfAtNodes[node] = static_cast<double>(kf); // Will be updated per-component
					}
					// Convert LGLnodes to Eigen::VectorXd
					Eigen::VectorXd LGLnodes_eigen(_disc.nNodes);
					for (unsigned int i = 0; i < _disc.nNodes; ++i)
						LGLnodes_eigen[i] = _convDispOp.LGLnodes()[i];
					_M_K[cell] = parts::dgtoolbox::radialFilmDiffusionMatrix(
						_disc.polyDeg, LGLnodes_eigen, rho_left, deltaRho, kfAtNodes);
				}

				// Precompute M_ρ^{-1} * M_K for each cell
				_filmDiffCoupling[cell] = _convDispOp.invMRho(cell) * _M_K[cell];
			}
		}

		void RadialLumpedRateModelWithPoresDG::updateFilmDiffusionValues(unsigned int secIdx, unsigned int comp)
		{
			if (!_variableFilmDiff)
				return;

			const unsigned int nNodes = _disc.nNodes;
			const double deltaRho = static_cast<double>(_convDispOp.columnLength()) / _disc.nElem;
			const active* const kf = getSectionDependentSlice(_filmDiffusion, _disc.nComp, secIdx);

			// Convert LGLnodes to Eigen::VectorXd once
			Eigen::VectorXd LGLnodes_eigen(_disc.nNodes);
			for (unsigned int i = 0; i < _disc.nNodes; ++i)
				LGLnodes_eigen[i] = _convDispOp.LGLnodes()[i];

			for (unsigned int cell = 0; cell < _disc.nElem; ++cell)
			{
				const double rho_left = _convDispOp.elemLeftBound(cell);

				for (unsigned int node = 0; node < nNodes; ++node)
				{
					const double relPos = (rho_left + 0.5 * deltaRho * (1.0 + _convDispOp.LGLnodes()[node])) /
						static_cast<double>(_convDispOp.columnLength());
					const active modifier = _filmDiffDep->getValue(*this, ColumnPosition{ relPos, 0.0, 0.0 }, comp, ParTypeIndep, BoundStateIndep, 0.0);
					_filmDiffAtNodes[cell][node] = static_cast<double>(kf[comp]) * static_cast<double>(modifier);
				}

				// Recompute M_K for this cell using variable k_f
				_M_K[cell] = parts::dgtoolbox::radialFilmDiffusionMatrix(
					_disc.polyDeg, LGLnodes_eigen, rho_left, deltaRho, _filmDiffAtNodes[cell]);

				// Recompute coupling matrix
				_filmDiffCoupling[cell] = _convDispOp.invMRho(cell) * _M_K[cell];
			}
		}

		unsigned int RadialLumpedRateModelWithPoresDG::threadLocalMemorySize() const CADET_NOEXCEPT
		{
			return 0;
		}

		unsigned int RadialLumpedRateModelWithPoresDG::requiredADdirs() const CADET_NOEXCEPT
		{
			return _jacobianAdDirs;
		}

		void RadialLumpedRateModelWithPoresDG::useAnalyticJacobian(const bool analyticJac)
		{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
			_analyticJac = analyticJac;
			if (!_analyticJac)
			{
				// AD Jacobian required - compute from conv-disp bandwidth
				// The bulk phase uses conv-disp structure, pore phase has local coupling
				// Total bandwidth: conv-disp bandwidth for bulk + pore local coupling
				const unsigned int bulkBandwidth = _convDispOp.requiredADdirs();
				const unsigned int parStride = _disc.nComp + _disc.strideBound;
				// Film diffusion couples bulk and pore within each cell, doubling effective bandwidth
				_jacobianAdDirs = std::max(bulkBandwidth, 2 * _disc.nNodes * parStride + 1);
			}
			else
				_jacobianAdDirs = 0;
#else
			// Use AD Jacobian if analytic Jacobian is to be checked
			_analyticJac = false;
			const unsigned int bulkBandwidth = _convDispOp.requiredADdirs();
			const unsigned int parStride = _disc.nComp + _disc.strideBound;
			_jacobianAdDirs = std::max(bulkBandwidth, 2 * _disc.nNodes * parStride + 1);
#endif
		}

		void RadialLumpedRateModelWithPoresDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac)
		{
			// Update section index
			updateSection(secIdx);

			// Update convection-dispersion operator for new section
			_convDispOp.notifyDiscontinuousSectionTransition(t, secIdx, _jacInlet);

			// Recompute film diffusion matrices if section-dependent
			if (isSectionDependent(_filmDiffusionMode))
			{
				computeFilmDiffusionMatrices();
			}
		}

		void RadialLumpedRateModelWithPoresDG::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
		{
			_convDispOp.setFlowRates(in[0], out[0], _colPorosity);
		}

		void RadialLumpedRateModelWithPoresDG::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
		{
			Exporter expr(_disc, *this, solution);
			recorder.beginUnitOperation(_unitOpIdx, *this, expr);
			recorder.endUnitOperation();
		}

		void RadialLumpedRateModelWithPoresDG::reportSolutionStructure(ISolutionRecorder& recorder) const
		{
			Exporter expr(_disc, *this, nullptr);
			recorder.unitOperationStructure(_unitOpIdx, *this, expr);
		}

		int RadialLumpedRateModelWithPoresDG::jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerResidual);

			if (_analyticJac)
				return residualImpl<double, double, double, true, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
			else
				return residualWithJacobian(simTime, ConstSimulationState{ simState.vecStateY, nullptr }, nullptr, adJac, threadLocalMem);
		}

		int RadialLumpedRateModelWithPoresDG::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerResidual);

			return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
		}

		int RadialLumpedRateModelWithPoresDG::residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerResidual);

			return residual(simTime, simState, res, adJac, threadLocalMem, true, false);
		}

		int RadialLumpedRateModelWithPoresDG::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem, bool updateJacobian, bool paramSensitivity)
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
		int RadialLumpedRateModelWithPoresDG::residualImpl(double t, unsigned int secIdx, StateType const* const y_, double const* const yDot_, ResidualType* const res_, util::ThreadLocalStorage& threadLocalMem)
		{
			Indexer idxr(_disc);

			// Compute Jacobian if requested
			if (wantJac) {
				if (!wantRes || _disc.newStaticJac) {
					// TODO: Implement analytic Jacobian computation
					// For now, just rely on AD
					_disc.newStaticJac = false;
				}
			}

			// Initialize residual to zero
			if (wantRes)
			{
				std::fill(res_, res_ + numDofs(), ResidualType(0.0));
			}

			// Compute bulk convection dispersion residual
			// Note: The residual function operates on bulk phase c only
			if (wantRes)
				_convDispOp.residual(*this, t, secIdx, y_, yDot_, res_, typename cadet::ParamSens<ParamType>::enabled());

			// Compute film transfer and binding residuals
			const ParamType Fc = static_cast<ParamType>((1.0 - static_cast<double>(_colPorosity)) / static_cast<double>(_colPorosity));
			const ParamType Q = static_cast<ParamType>(_parGeomSurfToVol / static_cast<double>(_parRadius));
			const ParamType epsP = static_cast<ParamType>(_parPorosity);

			const active* const poreAccFactor = _poreAccessFactor.data();

			// Process each component
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				// Update film diffusion if variable
				if (_variableFilmDiff)
					updateFilmDiffusionValues(secIdx, comp);

				// Process each cell
				for (unsigned int cell = 0; cell < _disc.nElem; ++cell)
				{
					const unsigned int cellOffset = cell * _disc.nNodes;

					// Extract c and cp values into local vectors for this cell
					Eigen::Matrix<StateType, Eigen::Dynamic, 1> c_cell(_disc.nNodes);
					Eigen::Matrix<StateType, Eigen::Dynamic, 1> cp_cell(_disc.nNodes);

					for (unsigned int node = 0; node < _disc.nNodes; ++node)
					{
						c_cell[node] = y_[idxr.offsetC() + comp + (cellOffset + node) * idxr.strideColNode()];
						cp_cell[node] = y_[idxr.offsetCp() + comp + (cellOffset + node) * idxr.strideParNode()];
					}

					// Compute film diffusion term: filmDiff = Q * M_ρ^{-1} * M_K * (c - cp)
					Eigen::Matrix<StateType, Eigen::Dynamic, 1> diff_cp = c_cell - cp_cell;
					Eigen::Matrix<StateType, Eigen::Dynamic, 1> filmTermState = _filmDiffCoupling[cell].template cast<StateType>() * diff_cp;

					// Scale by Q and convert to ResidualType (element-wise to handle AD types)
					Eigen::Matrix<ResidualType, Eigen::Dynamic, 1> filmTerm(_disc.nNodes);
					for (unsigned int node = 0; node < _disc.nNodes; ++node)
						filmTerm[node] = filmTermState[node] * Q;

					if (wantRes)
					{
						// Add film transfer to bulk phase residual
						// dc/dt = convDisp - Fc * filmTerm
						for (unsigned int node = 0; node < _disc.nNodes; ++node)
						{
							const unsigned int bulkIdx = idxr.offsetC() + comp + (cellOffset + node) * idxr.strideColNode();
							res_[bulkIdx] -= Fc * filmTerm[node];
						}

						// Add film transfer to pore phase residual
						// dcp/dt = filmTerm / eps_p - Fp * dq/dt
						for (unsigned int node = 0; node < _disc.nNodes; ++node)
						{
							const unsigned int poreIdx = idxr.offsetCp() + comp + (cellOffset + node) * idxr.strideParNode();
							res_[poreIdx] += filmTerm[node] / (epsP * static_cast<ParamType>(poreAccFactor[comp]));
						}
					}
				}
			}

			// Compute binding kinetics at each node
#ifdef CADET_PARALLELIZE
			tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nPoints), [&](std::size_t blk)
#else
			for (unsigned int blk = 0; blk < _disc.nPoints; ++blk)
#endif
			{
				linalg::BandedEigenSparseRowIterator jacIt(_jac, _disc.nPoints * _disc.nComp + blk * idxr.strideParNode());

				// Get pointer to particle state at this node (cp and q)
				StateType const* const localY = y_ + idxr.offsetCp() + idxr.strideParNode() * blk;

				const parts::cell::CellParameters cellResParams
				{
					_disc.nComp,
					_disc.nBound,
					_disc.boundOffset,
					_disc.strideBound,
					_binding[0]->reactionQuasiStationarity(),
					_parPorosity,
					poreAccFactor,
					_binding[0],
					(_dynReaction[0] && (_dynReaction[0]->numReactionsCombined() > 0)) ? &_reaction : nullptr
				};

				// Relative position of current node
				double z = _convDispOp.relativeCoordinate(blk);

				if (wantRes)
				{
					ResidualType* const localRes = res_ + idxr.offsetCp() + idxr.strideParNode() * blk;
					double const* const localYdot = yDot_ ? yDot_ + idxr.offsetCp() + idxr.strideParNode() * blk : nullptr;

					parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantJac, true>(
						t, secIdx, ColumnPosition{ z, 0.0, 0.0 }, localY, localYdot, localRes, jacIt, cellResParams, threadLocalMem.get()
					);
				}
				else if (wantJac)
				{
					parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantJac, true, false>(
						t, secIdx, ColumnPosition{ z, 0.0, 0.0 }, localY, nullptr, nullptr, jacIt, cellResParams, threadLocalMem.get()
					);
				}

			} CADET_PARFOR_END;

			// Add time derivatives
			if (wantRes && yDot_)
			{
				// Bulk phase: dc/dt
				for (unsigned int node = 0; node < _disc.nPoints; ++node)
				{
					for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					{
						const unsigned int bulkIdx = idxr.offsetC() + comp + node * idxr.strideColNode();
						res_[bulkIdx] += yDot_[bulkIdx];
					}
				}

				// Pore phase: dcp/dt + Fp * dq/dt is already handled in bindingKernel
				for (unsigned int node = 0; node < _disc.nPoints; ++node)
				{
					for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					{
						const unsigned int poreIdx = idxr.offsetCp() + comp + node * idxr.strideParNode();
						res_[poreIdx] += yDot_[poreIdx];
					}
				}
			}

			return 0;
		}

		int RadialLumpedRateModelWithPoresDG::residualSensFwdAdOnly(const SimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerResidualSens);

			return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adRes, threadLocalMem);
		}

		int RadialLumpedRateModelWithPoresDG::residualSensFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerResidualSens);

			return residual(simTime, simState, nullptr, adJac, threadLocalMem, true, true);
		}

		int RadialLumpedRateModelWithPoresDG::residualSensFwdCombine(const SimulationTime& simTime, const ConstSimulationState& simState,
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

		void RadialLumpedRateModelWithPoresDG::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
		{
			Indexer idxr(_disc);

			// Inlet DOFs
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				ret[comp] = alpha * yS[comp] + beta * ret[comp];
			}

			// Pure DOFs (bulk + pore)
			Eigen::Map<Eigen::VectorXd> ret_vec(ret + idxr.offsetC(), numPureDofs());
			Eigen::Map<const Eigen::VectorXd> yS_vec(yS + idxr.offsetC(), numPureDofs());
			ret_vec = alpha * _jac * yS_vec + beta * ret_vec;

			// Inlet contribution
			unsigned int offInlet = _convDispOp.forwardFlow() ? 0 : (_disc.nElem - 1u) * idxr.strideColCell();

			for (unsigned int comp = 0; comp < _disc.nComp; comp++) {
				for (unsigned int node = 0; node < _disc.nNodes; node++) {
					ret[idxr.offsetC() + offInlet + node * idxr.strideColNode() + comp] += alpha * _jacInlet.coeff(node, 0) * yS[comp];
				}
			}
		}

		void RadialLumpedRateModelWithPoresDG::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
		{
			Indexer idxr(_disc);
			const double invBetaP = 1.0 / static_cast<double>(_parPorosity) - 1.0;

			std::fill_n(ret, numDofs(), 0.0);

			// Bulk phase: dc/dt
			for (unsigned int node = 0; node < _disc.nPoints; ++node)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					const unsigned int idx = idxr.offsetC() + comp + node * idxr.strideColNode();
					ret[idx] = sDot[idx];
				}
			}

			// Pore phase: dcp/dt + Fp * sum_j dq_j/dt
			for (unsigned int node = 0; node < _disc.nPoints; ++node)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					const unsigned int cpIdx = idxr.offsetCp() + comp + node * idxr.strideParNode();
					ret[cpIdx] = sDot[cpIdx];

					// Add bound state contributions
					if (_binding[0] && !_binding[0]->reactionQuasiStationarity()[comp])
					{
						for (unsigned int bnd = 0; bnd < _disc.nBound[comp]; ++bnd)
						{
							const unsigned int qIdx = idxr.offsetCp() + idxr.strideParLiquid() + _disc.boundOffset[comp] + bnd + node * idxr.strideParNode();
							ret[cpIdx] += invBetaP * sDot[qIdx];
						}
					}
				}

				// Bound states: dq/dt
				for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
				{
					const unsigned int qIdx = idxr.offsetCp() + idxr.strideParLiquid() + bnd + node * idxr.strideParNode();
					ret[qIdx] = sDot[qIdx];
				}
			}
		}

		void RadialLumpedRateModelWithPoresDG::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
		{
			if (_binding[0])
				_binding[0]->setExternalFunctions(extFuns, size);
			if (_dynReaction[0])
				_dynReaction[0]->setExternalFunctions(extFuns, size);
		}

		unsigned int RadialLumpedRateModelWithPoresDG::localOutletComponentIndex(unsigned int port) const CADET_NOEXCEPT
		{
			Indexer idxr(_disc);
			// Outlet is at end of column for forward flow, beginning for backward
			if (_convDispOp.forwardFlow())
				return idxr.offsetC() + (_disc.nPoints - 1) * idxr.strideColNode();
			else
				return idxr.offsetC();
		}

		unsigned int RadialLumpedRateModelWithPoresDG::localOutletComponentStride(unsigned int port) const CADET_NOEXCEPT
		{
			return 1;
		}

		unsigned int RadialLumpedRateModelWithPoresDG::localInletComponentIndex(unsigned int port) const CADET_NOEXCEPT
		{
			return 0;
		}

		unsigned int RadialLumpedRateModelWithPoresDG::localInletComponentStride(unsigned int port) const CADET_NOEXCEPT
		{
			return 1;
		}

		void RadialLumpedRateModelWithPoresDG::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
		{
			// TODO: Implement proper error tolerance expansion
			std::fill_n(expandOut, numDofs(), errorSpec[0]);
		}

		int RadialLumpedRateModelWithPoresDG::linearSolve(double t, double alpha, double tol, double* const rhs, double const* const weight, const ConstSimulationState& simState)
		{
			BENCH_SCOPE(_timerLinearSolve);

			Indexer idxr(_disc);

			if (_factorizeJacobian)
			{
				// Assemble discretized Jacobian
				assembleDiscretizedJacobian(alpha, idxr);

				// Factorize
				_linearSolver->analyzePattern(_jacDisc);
				_linearSolver->factorize(_jacDisc);

				_factorizeJacobian = false;
			}

			// Solve
			Eigen::Map<Eigen::VectorXd> r(rhs + idxr.offsetC(), numPureDofs());
			r = _linearSolver->solve(r);

			return 0;
		}

		void RadialLumpedRateModelWithPoresDG::assembleDiscretizedJacobian(double alpha, const Indexer& idxr)
		{
			_jacDisc = _jac;

			const double invBetaP = 1.0 / static_cast<double>(_parPorosity) - 1.0;

			// Add time derivatives to bulk phase
			for (unsigned int node = 0; node < _disc.nPoints; ++node)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					const unsigned int idx = comp + node * idxr.strideColNode();
					_jacDisc.coeffRef(idx, idx) += alpha;
				}
			}

			// Add time derivatives to pore phase
			const unsigned int poreOffset = _disc.nPoints * _disc.nComp;
			for (unsigned int node = 0; node < _disc.nPoints; ++node)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					const unsigned int cpIdx = poreOffset + comp + node * idxr.strideParNode();
					_jacDisc.coeffRef(cpIdx, cpIdx) += alpha;

					// Coupling with bound states
					for (unsigned int bnd = 0; bnd < _disc.nBound[comp]; ++bnd)
					{
						const unsigned int qIdx = poreOffset + idxr.strideParLiquid() + _disc.boundOffset[comp] + bnd + node * idxr.strideParNode();
						_jacDisc.coeffRef(cpIdx, qIdx) += alpha * invBetaP;
					}
				}

				// Bound states
				for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
				{
					const unsigned int qIdx = poreOffset + idxr.strideParLiquid() + bnd + node * idxr.strideParNode();
					_jacDisc.coeffRef(qIdx, qIdx) += alpha;
				}
			}
		}

		void RadialLumpedRateModelWithPoresDG::prepareADvectors(const AdJacobianParams& adJac) const
		{
			ad::prepareAdVectorSeedsForBandMatrix(adJac.adY + Indexer(_disc).offsetC(), adJac.adDirOffset, numPureDofs(), _jac.outerSize(), 0, numPureDofs());
		}

		void RadialLumpedRateModelWithPoresDG::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
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

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
		void RadialLumpedRateModelWithPoresDG::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
		{
			// TODO: Implement Jacobian check
		}
#endif

		void RadialLumpedRateModelWithPoresDG::applyInitialCondition(const SimulationState& simState) const
		{
			Indexer idxr(_disc);

			// Inlet DOFs
			std::fill_n(simState.vecStateY, _disc.nComp, 0.0);

			// Bulk phase
			for (unsigned int i = 0; i < _disc.nPoints * _disc.nComp; ++i)
				simState.vecStateY[idxr.offsetC() + i] = static_cast<double>(_initC[i]);

			// Pore phase
			for (unsigned int i = 0; i < _disc.nPoints * _disc.nComp; ++i)
				simState.vecStateY[idxr.offsetCp() + i] = static_cast<double>(_initCp[i]);

			// Solid phase
			for (unsigned int node = 0; node < _disc.nPoints; ++node)
			{
				for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
				{
					simState.vecStateY[idxr.offsetCp() + idxr.strideParLiquid() + bnd + node * idxr.strideParNode()] =
						static_cast<double>(_initCs[node * _disc.strideBound + bnd]);
				}
			}

			// Time derivatives
			if (simState.vecStateYdot)
				std::fill_n(simState.vecStateYdot, numDofs(), 0.0);
		}

		void RadialLumpedRateModelWithPoresDG::readInitialCondition(IParameterProvider& paramProvider)
		{
			Indexer idxr(_disc);

			// Bulk phase initial conditions
			std::vector<double> initC = paramProvider.getDoubleArray("INIT_C");
			if (initC.size() < _disc.nComp)
				throw InvalidParameterException("INIT_C has fewer than NCOMP (" + std::to_string(_disc.nComp) + ") entries");

			// Fill initial conditions for all nodes
			if (initC.size() >= _disc.nPoints * _disc.nComp)
			{
				// Per-node initialization
				for (unsigned int i = 0; i < _disc.nPoints * _disc.nComp; ++i)
					_initC[i] = initC[i];
			}
			else
			{
				// Per-component initialization
				for (unsigned int node = 0; node < _disc.nPoints; ++node)
				{
					for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
						_initC[node * _disc.nComp + comp] = initC[comp];
				}
			}

			// Pore phase initial conditions
			if (paramProvider.exists("INIT_CP"))
			{
				std::vector<double> initCp = paramProvider.getDoubleArray("INIT_CP");
				if (initCp.size() >= _disc.nPoints * _disc.nComp)
				{
					for (unsigned int i = 0; i < _disc.nPoints * _disc.nComp; ++i)
						_initCp[i] = initCp[i];
				}
				else
				{
					for (unsigned int node = 0; node < _disc.nPoints; ++node)
					{
						for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
							_initCp[node * _disc.nComp + comp] = initCp[comp % initCp.size()];
					}
				}
			}
			else
			{
				// Default: pore phase = bulk phase
				_initCp = _initC;
			}

			// Solid phase initial conditions
			if (_disc.strideBound > 0)
			{
				std::vector<double> initCs = paramProvider.getDoubleArray("INIT_Q");
				if (initCs.size() >= _disc.nPoints * _disc.strideBound)
				{
					for (unsigned int i = 0; i < _disc.nPoints * _disc.strideBound; ++i)
						_initCs[i] = initCs[i];
				}
				else
				{
					for (unsigned int node = 0; node < _disc.nPoints; ++node)
					{
						for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
							_initCs[node * _disc.strideBound + bnd] = initCs[bnd % initCs.size()];
					}
				}
			}
		}

		void RadialLumpedRateModelWithPoresDG::consistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerConsistentInit);

			// Apply initial conditions
			applyInitialCondition(SimulationState{ vecStateY, nullptr });

			// Solve for consistent initial state if binding is quasi-stationary
			if (_binding[0] && _binding[0]->hasQuasiStationaryReactions())
			{
				// TODO: Implement quasi-stationary binding initialization
			}
		}

		void RadialLumpedRateModelWithPoresDG::consistentInitialTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerConsistentInit);

			// Set time derivatives to zero initially
			std::fill_n(vecStateYdot, numDofs(), 0.0);
		}

		void RadialLumpedRateModelWithPoresDG::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
		{
			for (double* sensY : vecSensY)
				std::fill_n(sensY, numDofs(), 0.0);
		}

		void RadialLumpedRateModelWithPoresDG::consistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
			std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
		{
			for (unsigned int i = 0; i < vecSensY.size(); ++i)
			{
				std::fill_n(vecSensY[i], numDofs(), 0.0);
				std::fill_n(vecSensYdot[i], numDofs(), 0.0);
			}
		}

		void RadialLumpedRateModelWithPoresDG::leanConsistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
		{
			consistentInitialState(simTime, vecStateY, adJac, errorTol, threadLocalMem);
		}

		void RadialLumpedRateModelWithPoresDG::leanConsistentInitialTimeDerivative(double t, double const* const vecStateY, double* const vecStateYdot, double* const res, util::ThreadLocalStorage& threadLocalMem)
		{
			std::fill_n(vecStateYdot, numDofs(), 0.0);
		}

		void RadialLumpedRateModelWithPoresDG::leanConsistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
			std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
		{
			consistentInitialSensitivity(simTime, simState, vecSensY, vecSensYdot, adRes, threadLocalMem);
		}

		bool RadialLumpedRateModelWithPoresDG::setParameter(const ParameterId& pId, double value)
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

		bool RadialLumpedRateModelWithPoresDG::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
		{
			if (pId.unitOperation == _unitOpIdx)
			{
				if (_convDispOp.setSensitiveParameter(_sensParams, pId, adDirection, adValue))
					return true;
			}

			return UnitOperationBase::setSensitiveParameter(pId, adDirection, adValue);
		}

		void RadialLumpedRateModelWithPoresDG::setSensitiveParameterValue(const ParameterId& pId, double value)
		{
			if (pId.unitOperation == _unitOpIdx)
			{
				if (_convDispOp.setSensitiveParameterValue(_sensParams, pId, value))
					return;
			}

			UnitOperationBase::setSensitiveParameterValue(pId, value);
		}

		// Jacobian pattern functions
		void RadialLumpedRateModelWithPoresDG::setPattern(Eigen::SparseMatrix<double, RowMajor>& mat, bool stateDer, bool has_reaction)
		{
			std::vector<T> tripletList;

			Indexer idxr(_disc);

			// Estimate number of entries
			unsigned int conv_disp_entries = _convDispOp.nJacEntries(false);
			unsigned int film_diff_entries = 2 * _disc.nPoints * _disc.nComp; // c-cp coupling
			unsigned int binding_entries = _disc.nPoints * _disc.strideBound * (_disc.strideBound + _disc.nComp);
			unsigned int reaction_entries = has_reaction ? _disc.nPoints * _disc.nComp * (_disc.strideBound + _disc.nComp) : 0;

			tripletList.reserve(conv_disp_entries + film_diff_entries + binding_entries + reaction_entries);

			// Add convection-dispersion pattern (bulk phase)
			_convDispOp.convDispJacPattern(tripletList);

			// Add film diffusion coupling pattern
			filmDiffusionPattern(tripletList);

			// Add binding pattern
			bindingAndReactionPattern(tripletList, has_reaction);

			if (stateDer)
				stateDerPattern(tripletList);

			mat.setFromTriplets(tripletList.begin(), tripletList.end());
		}

		void RadialLumpedRateModelWithPoresDG::filmDiffusionPattern(std::vector<T>& tripletList)
		{
			Indexer idxr(_disc);
			const unsigned int poreOffset = _disc.nPoints * _disc.nComp;

			for (unsigned int cell = 0; cell < _disc.nElem; ++cell)
			{
				for (unsigned int node = 0; node < _disc.nNodes; ++node)
				{
					const unsigned int point = cell * _disc.nNodes + node;

					for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					{
						const unsigned int bulkIdx = comp + point * idxr.strideColNode();
						const unsigned int poreIdx = poreOffset + comp + point * idxr.strideParNode();

						// Film diffusion couples c and cp within each cell (due to M_K matrix)
						for (unsigned int nodeInner = 0; nodeInner < _disc.nNodes; ++nodeInner)
						{
							const unsigned int pointInner = cell * _disc.nNodes + nodeInner;
							const unsigned int bulkIdxInner = comp + pointInner * idxr.strideColNode();
							const unsigned int poreIdxInner = poreOffset + comp + pointInner * idxr.strideParNode();

							// dc/dt depends on c and cp in the cell
							tripletList.push_back(T(bulkIdx, bulkIdxInner, 0.0));
							tripletList.push_back(T(bulkIdx, poreIdxInner, 0.0));

							// dcp/dt depends on c and cp in the cell
							tripletList.push_back(T(poreIdx, bulkIdxInner, 0.0));
							tripletList.push_back(T(poreIdx, poreIdxInner, 0.0));
						}
					}
				}
			}
		}

		void RadialLumpedRateModelWithPoresDG::stateDerPattern(std::vector<T>& tripletList)
		{
			Indexer idxr(_disc);
			const unsigned int poreOffset = _disc.nPoints * _disc.nComp;

			// Bulk phase time derivative
			for (unsigned int point = 0; point < _disc.nPoints; point++)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; comp++)
				{
					const unsigned int idx = comp + point * idxr.strideColNode();
					tripletList.push_back(T(idx, idx, 0.0));
				}
			}

			// Pore phase time derivative and coupling to bound states
			for (unsigned int point = 0; point < _disc.nPoints; point++)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; comp++)
				{
					const unsigned int cpIdx = poreOffset + comp + point * idxr.strideParNode();
					tripletList.push_back(T(cpIdx, cpIdx, 0.0));

					// Coupling with bound states
					for (unsigned int bnd = 0; bnd < _disc.nBound[comp]; ++bnd)
					{
						const unsigned int qIdx = poreOffset + idxr.strideParLiquid() + _disc.boundOffset[comp] + bnd + point * idxr.strideParNode();
						tripletList.push_back(T(cpIdx, qIdx, 0.0));
					}
				}

				// Bound state time derivatives
				for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
				{
					const unsigned int qIdx = poreOffset + idxr.strideParLiquid() + bnd + point * idxr.strideParNode();
					tripletList.push_back(T(qIdx, qIdx, 0.0));
				}
			}
		}

		void RadialLumpedRateModelWithPoresDG::bindingAndReactionPattern(std::vector<T>& tripletList, bool has_reaction)
		{
			Indexer idxr(_disc);
			const unsigned int poreOffset = _disc.nPoints * _disc.nComp;

			for (unsigned int point = 0; point < _disc.nPoints; point++)
			{
				// Binding: q depends on cp and q
				for (unsigned int solid = 0; solid < _disc.strideBound; solid++)
				{
					const unsigned int qIdx = poreOffset + idxr.strideParLiquid() + solid + point * idxr.strideParNode();

					for (unsigned int conc = 0; conc < _disc.nComp + _disc.strideBound; conc++)
					{
						const unsigned int depIdx = poreOffset + conc + point * idxr.strideParNode();
						tripletList.push_back(T(qIdx, depIdx, 0.0));
					}
				}

				// Reaction in pore phase
				if (has_reaction)
				{
					for (unsigned int liquid = 0; liquid < _disc.nComp; liquid++)
					{
						const unsigned int cpIdx = poreOffset + liquid + point * idxr.strideParNode();

						for (unsigned int conc = 0; conc < _disc.nComp + _disc.strideBound; conc++)
						{
							const unsigned int depIdx = poreOffset + conc + point * idxr.strideParNode();
							tripletList.push_back(T(cpIdx, depIdx, 0.0));
						}
					}
				}
			}
		}

		// Exporter implementations
		int RadialLumpedRateModelWithPoresDG::Exporter::writeMobilePhase(double* buffer) const
		{
			const double* data = _data + _idx.offsetC();
			for (unsigned int i = 0; i < _disc.nPoints; ++i)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					buffer[comp * _disc.nPoints + i] = data[i * _idx.strideColNode() + comp];
			}
			return _disc.nComp * _disc.nPoints;
		}

		int RadialLumpedRateModelWithPoresDG::Exporter::writeSolidPhase(double* buffer) const
		{
			const double* data = _data + _idx.offsetCp() + _idx.strideParLiquid();
			for (unsigned int i = 0; i < _disc.nPoints; ++i)
			{
				for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
					buffer[bnd * _disc.nPoints + i] = data[i * _idx.strideParNode() + bnd];
			}
			return _disc.strideBound * _disc.nPoints;
		}

		int RadialLumpedRateModelWithPoresDG::Exporter::writeParticleMobilePhase(double* buffer) const
		{
			const double* data = _data + _idx.offsetCp();
			for (unsigned int i = 0; i < _disc.nPoints; ++i)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					buffer[comp * _disc.nPoints + i] = data[i * _idx.strideParNode() + comp];
			}
			return _disc.nComp * _disc.nPoints;
		}

		int RadialLumpedRateModelWithPoresDG::Exporter::writeSolidPhase(unsigned int parType, double* buffer) const
		{
			return writeSolidPhase(buffer);
		}

		int RadialLumpedRateModelWithPoresDG::Exporter::writeParticleMobilePhase(unsigned int parType, double* buffer) const
		{
			return writeParticleMobilePhase(buffer);
		}

		int RadialLumpedRateModelWithPoresDG::Exporter::writeInlet(unsigned int port, double* buffer) const
		{
			std::copy_n(_data, _disc.nComp, buffer);
			return _disc.nComp;
		}

		int RadialLumpedRateModelWithPoresDG::Exporter::writeInlet(double* buffer) const
		{
			return writeInlet(0, buffer);
		}

		int RadialLumpedRateModelWithPoresDG::Exporter::writeOutlet(unsigned int port, double* buffer) const
		{
			if (_model._convDispOp.forwardFlow())
			{
				const double* data = _data + _idx.offsetC() + (_disc.nPoints - 1) * _idx.strideColNode();
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					buffer[comp] = data[comp];
			}
			else
			{
				const double* data = _data + _idx.offsetC();
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					buffer[comp] = data[comp];
			}
			return _disc.nComp;
		}

		int RadialLumpedRateModelWithPoresDG::Exporter::writeOutlet(double* buffer) const
		{
			return writeOutlet(0, buffer);
		}

	} // namespace model
} // namespace cadet

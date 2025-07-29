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


#include "model/MultiChannelTransportModel.hpp"
#include "AdUtils.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "ParamReaderHelper.hpp"
#include "model/parts/BindingCellKernel.hpp"
#include "SimulationTypes.hpp"
#include "SensParamUtil.hpp"
#include "linalg/Subset.hpp"

#include <algorithm>
#include <functional>
#include <set>

#include "LoggingUtils.hpp"
#include "Logging.hpp"
#include "ParallelSupport.hpp"
#include "iostream"

namespace cadet
{

namespace model
{

int MultiChannelTransportModel::multiplexInitialConditions(const cadet::ParameterId& pId, unsigned int adDirection, double adValue)
{
	if (pId.name == hashString("INIT_C") && (pId.section == SectionIndep) && (pId.boundState == BoundStateIndep) && (pId.particleType == ParTypeIndep) && (pId.component != CompIndep))
	{
		if ((pId.reaction == ReactionIndep) && _singleRadiusInitC)
		{
			_sensParams.insert(&_initC[pId.component]);
			for (unsigned int r = 0; r < _disc.nChannel; ++r)
				_initC[r * _disc.nComp + pId.component].setADValue(adDirection, adValue);

			return 1;
		}
		else if ((pId.reaction != ReactionIndep) && !_singleRadiusInitC)
		{
			_sensParams.insert(&_initC[pId.reaction * _disc.nComp + pId.component]);
			_initC[pId.reaction * _disc.nComp + pId.component].setADValue(adDirection, adValue);
			return 1;
		}

		return -1;
	}
	else if (pId.name == hashString("INIT_C"))
		return -1;

	return 0;
}

int MultiChannelTransportModel::multiplexInitialConditions(const cadet::ParameterId& pId, double val, bool checkSens)
{
	if (pId.name == hashString("INIT_C") && (pId.section == SectionIndep) && (pId.boundState == BoundStateIndep) && (pId.particleType == ParTypeIndep) && (pId.component != CompIndep))
	{
		if ((pId.reaction == ReactionIndep) && _singleRadiusInitC)
		{
			if (checkSens && !contains(_sensParams, &_initC[pId.component]))
				return -1;

			for (unsigned int r = 0; r < _disc.nChannel; ++r)
				_initC[r * _disc.nComp + pId.component].setValue(val);

			return 1;
		}
		else if ((pId.reaction != ReactionIndep) && !_singleRadiusInitC)
		{
			if (checkSens && !contains(_sensParams, &_initC[pId.reaction * _disc.nComp + pId.component]))
				return -1;

			_initC[pId.reaction * _disc.nComp + pId.component].setValue(val);
			return 1;
		}
		else
			return -1;
	}

	return 0;
}

void MultiChannelTransportModel::applyInitialCondition(const SimulationState& simState) const
{
	Indexer idxr(_disc);

	// Check whether full state vector is available as initial condition
	if (!_initState.empty())
	{
		std::fill(simState.vecStateY, simState.vecStateY + idxr.offsetC(), 0.0);
		std::copy(_initState.data(), _initState.data() + numPureDofs(), simState.vecStateY + idxr.offsetC());

		if (!_initStateDot.empty())
		{
			std::fill(simState.vecStateYdot, simState.vecStateYdot + idxr.offsetC(), 0.0);
			std::copy(_initStateDot.data(), _initStateDot.data() + numPureDofs(), simState.vecStateYdot + idxr.offsetC());
		}
		else
			std::fill(simState.vecStateYdot, simState.vecStateYdot + numDofs(), 0.0);

		return;
	}

	double* const stateYbulk = simState.vecStateY + idxr.offsetC();

	// Loop over axial cells
	for (unsigned int col = 0; col < _disc.nCol; ++col)
	{
		// Loop over radial cells
		for (unsigned int rad = 0; rad < _disc.nChannel; ++rad)
		{
			// Loop over components in cell
			for (unsigned comp = 0; comp < _disc.nComp; ++comp)
				stateYbulk[col * idxr.strideColAxialCell() + rad * idxr.strideChannelCell() + comp * idxr.strideColComp()] = static_cast<double>(_initC[comp + rad * _disc.nComp]);
		}
	}
}

void MultiChannelTransportModel::readInitialCondition(IParameterProvider& paramProvider)
{
	_initState.clear();
	_initStateDot.clear();

	// Check if INIT_STATE is present
	if (paramProvider.exists("INIT_STATE"))
	{
		const std::vector<double> initState = paramProvider.getDoubleArray("INIT_STATE");
		_initState = std::vector<double>(initState.begin(), initState.begin() + numPureDofs());

		// Check if INIT_STATE contains the full state and its time derivative
		if (initState.size() >= 2 * numPureDofs())
			_initStateDot = std::vector<double>(initState.begin() + numPureDofs(), initState.begin() + 2 * numPureDofs());
		return;
	}

	const std::vector<double> initC = paramProvider.getDoubleArray("INIT_C");
	_singleRadiusInitC = (initC.size() < _disc.nComp * _disc.nChannel);

	if (((initC.size() < _disc.nComp) && _singleRadiusInitC) || ((initC.size() < _disc.nComp * _disc.nChannel) && !_singleRadiusInitC))
		throw InvalidParameterException("INIT_C does not contain enough values for all components (and radial zones)");

	if (!_singleRadiusInitC)
		ad::copyToAd(initC.data(), _initC.data(), _disc.nComp * _disc.nChannel);
	else
	{
		for (unsigned int r = 0; r < _disc.nChannel; ++r)
			ad::copyToAd(initC.data(), _initC.data() + r * _disc.nComp, _disc.nComp);
	}
}

/**
 * @brief Creates quasi-stationary masks and conservation groups for a specific component
 * @param [in] comp Component index
 * @param [in] exchange Exchange model to query for quasi-stationary pairs
 * @param [in] nChannel Number of channels
 * @param [in] nComp Number of components
 * @param [out] qsMask Output mask vector (will be resized and populated)
 * @param [out] conservationGroups Groups of channels that form separate conservation constraints
 * @return Number of active entries in the mask
 */
static int createQuasiStationaryMaskWithGroups(unsigned int comp, IExchangeModel* exchange,
	unsigned int nChannel, unsigned int nComp,
	std::vector<int>& qsMask,
	std::vector<std::vector<unsigned int>>& conservationGroups)
{
	// Get quasi-stationary channel pairs for this component
	std::vector<std::pair<unsigned int, unsigned int>> qsChannelPairs;
	exchange->quasiStationarityMap(comp, qsChannelPairs);

	// Initialize mask with zeros
	qsMask.assign(nChannel * nComp, 0);
	conservationGroups.clear();

	if (qsChannelPairs.empty()) {
		return 0;
	}

	// Build connected components (groups of channels that are connected via quasi stationary exchange)
	std::vector<std::set<unsigned int>> groups;

	for (const auto& channelPair : qsChannelPairs)
	{
		const auto source = channelPair.first;
		const auto destination = channelPair.second;

		if (source >= nChannel || destination >= nChannel) continue;

		// Find if either channel is already in a group
		int sourceGroup = -1;
		int destGroup = -1;

		for (size_t i = 0; i < groups.size(); ++i)
		{
			if (groups[i].count(source)) sourceGroup = i;
			if (groups[i].count(destination)) destGroup = i;
		}

		//either sourse or destionation are in the group
		if (sourceGroup == -1 && destGroup == -1)
		{
			// Create new group
			groups.emplace_back();
			groups.back().insert(source);
			groups.back().insert(destination);
		}
		else if (sourceGroup != -1 && destGroup == -1)
		{
			// Add destination to source's group
			groups[sourceGroup].insert(destination);
		}
		else if (sourceGroup == -1 && destGroup != -1)
		{
			// Add source to destination's group
			groups[destGroup].insert(source);
		}
		else if (sourceGroup != destGroup)
		{
			// Merge two groups
			groups[sourceGroup].insert(groups[destGroup].begin(), groups[destGroup].end());
			groups.erase(groups.begin() + destGroup);
		}
		// If sourceGroup == destGroup, both are already in the same group
	}

	// Convert to output format and set mask
	for (const auto& group : groups)
	{
		std::vector<unsigned int> channelGroup(group.begin(), group.end());
		conservationGroups.push_back(channelGroup);

		// Set mask entries for this group
		for (unsigned int channel : group)
		{
			qsMask[channel * nComp + comp] = 1;
		}
	}

	// Count active entries
	return std::count(qsMask.begin(), qsMask.end(), 1);
}

/**
 * @brief Computes consistent initial values (state variables without their time derivatives)
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *
 *          The process works in two steps:
 *          <ol>
 *              <li>Solve all algebraic equations in the model (e.g., quasi-stationary isotherms, reaction equilibria).
 *                 Once all @f$ c_i @f$, @f$ c_{p,i} @f$, and @f$ q_i^{(j)} @f$ have been computed, solve for the
 *                 fluxes @f$ j_{f,i} @f$ (only linear equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0.
 *                 However, because of the algebraic equations, we need additional conditions to fully determine
 *                 @f$ \dot{y}@f$. By differentiating the algebraic equations with respect to time, we get the
 *                 missing linear equations (recall that the state vector @f$ y @f$ is fixed). The resulting system
 *                 has a similar structure as the system Jacobian.
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                                & \dot{J}_1     &        &           &   \\
 *                                &         & \ddots &           &   \\
 *                                &         &        & \dot{J}_{N_z}   &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_i @f$ denotes the Jacobian with respect to @f$ \dot{y}@f$. Note that the
 *                 @f$ J_{i,f} @f$ matrices in the right column are missing.
 *
 *     The right hand side of the linear system is given by the negative residual without contribution
 *     of @f$ \dot{y} @f$ for differential equations and 0 for algebraic equations
 *     (@f$ -\frac{\partial F}{\partial t}@f$, to be more precise).
 *
 *     The linear system is solved by backsubstitution. First, the diagonal blocks are solved in parallel.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     diagonal blocks.</li>
 *          </ol>
 *     This function performs step 1. See consistentInitialTimeDerivative() for step 2.
 *
 * 	   This function is to be used with consistentInitialTimeDerivative(). Do not mix normal and lean
 *     consistent initialization!
 *
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in,out] vecStateY State vector with initial values that are to be updated for consistency
 * @param [in,out] adJac Jacobian information for AD (AD vectors for residual and state, direction offset)
 * @param [in] errorTol Error tolerance for algebraic equations
 * @todo Decrease amount of allocated memory by partially using temporary vectors (state and Schur complement)
 */
void MultiChannelTransportModel::consistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	// Step 1: Solve algebraic equations
	// Step 1a: Compute quasi-stationary exchange model state
	
	// state vector structure: 
	// y = [channel1, channel2,..., channelN] with channeli = [c1, c2, ..., cNComp]

	// initialize evey compoent in every channel consistent
	for (auto comp = 0; comp < _disc.nComp; comp++)
	{
		if (!_exchange[0]->hasQuasiStationary(comp))
			return;

		// Create quasi-stationary mask for this component with conservation groups
		std::vector<int> qsMask;
		std::vector<std::vector<unsigned int>> conservationGroups;
		const int numActiveMask = createQuasiStationaryMaskWithGroups(comp, _exchange[0], _disc.nChannel, _disc.nComp, qsMask, conservationGroups);

		const linalg::ConstMaskArray mask{ qsMask.data(), static_cast<int>(_disc.nChannel * _disc.nComp) };
		const int probSize = linalg::numMaskActive(mask);
		
		//Problem capturing variables here
	#ifdef CADET_PARALLELIZE
		BENCH_SCOPE(_timerConsistentInitPar);
		tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nCol), [&](std::size_t pblk)
	#else
		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	#endif
		{
			LinearBufferAllocator tlmAlloc = threadLocalMem.get();

			// Reuse memory of band matrix for dense matrix
			linalg::DenseMatrixView fullJacobianMatrix(_convDispOp.jacobian().data(), nullptr, mask.len, mask.len);
			linalg::CompressedSparseMatrix& testJac = _convDispOp.jacobian();

			// Midpoint of current column cell (z coordinate) - needed in externally dependent adsorption kinetic
			const double z = (0.5 + static_cast<double>(pblk)) / static_cast<double>(_disc.nCol);

			// Get workspace memory
			BufferedArray<double> nonlinMemBuffer = tlmAlloc.array<double>(_nonlinearSolver->workspaceSize(probSize));
			double* const nonlinMem = static_cast<double*>(nonlinMemBuffer);

			BufferedArray<double> solutionBuffer = tlmAlloc.array<double>(probSize);
			double* const solution = static_cast<double*>(solutionBuffer);

			BufferedArray<double> fullResidualBuffer = tlmAlloc.array<double>(numDofs());
			double* const fullResidual = static_cast<double*>(fullResidualBuffer);

			BufferedArray<double> fullXBuffer = tlmAlloc.array<double>(numDofs());
			double* const fullX = static_cast<double*>(fullXBuffer);

			BufferedArray<double> jacobianMemBuffer = tlmAlloc.array<double>(probSize * probSize);
			linalg::DenseMatrixView jacobianMatrix(static_cast<double*>(jacobianMemBuffer), new lapackInt_t[probSize], probSize, probSize);

			// Get pointer to c variables in the channels
			auto cShellOfset = _disc.nComp * _disc.nChannel * (1 + pblk);
			double* const cShell = vecStateY + cShellOfset;
			active* const localAdRes = adJac.adRes ? adJac.adRes : nullptr;
			active* const localAdY = adJac.adY ? adJac.adY : nullptr;

			const ColumnPosition colPos{ z, 0.0, 0.0 };

			// Determine whether nonlinear solver is required
			//todo 

			// Extract initial values from current state
			linalg::selectVectorSubset(cShell, mask, solution);	// Save values of conserved moieties for this component
			
			// The total amount of every component in quasi-stationary channels must be constant
			std::vector<std::pair<unsigned int, unsigned int>> qsChannelPairs;
			_exchange[0]->quasiStationarityMap(comp, qsChannelPairs);
			
			// Calculate conserved amounts for each group
			std::vector<double> conservedQuants;
			for (const auto& group : conservationGroups)
			{
				double totalConservedAmount = 0.0;
				for (unsigned int channel : group)
				{
					totalConservedAmount += cShell[channel * _disc.nComp + comp];
				}
				conservedQuants.push_back(totalConservedAmount);
			}

			std::function<bool(double const* const, linalg::detail::DenseMatrixBase&)> jacFunc;
			if (localAdY && localAdRes)
			{
				jacFunc = [&](double const* const x, linalg::detail::DenseMatrixBase& mat)
				{
					// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
					// and initialize residuals with zero (also resetting directional values)
					ad::copyToAd(cShell, localAdY + cShellOfset, mask.len);
					ad::resetAd(localAdRes, numDofs());

					// Prepare input vector by overwriting masked items
					linalg::applyVectorSubset(x, mask, localAdY + cShellOfset);

					// Call residual function with AD types
					residualImpl<active, active, double, false>(simTime.t, simTime.secIdx, localAdY, nullptr, localAdRes, threadLocalMem);

					// Clear the matrix
					mat.setAll(0.0);

					// Extract Jacobian submatrix corresponding to masked variables
					const int numExchangeEq = probSize - conservedQuants.size();
					
					// Part 1: Extract exchange equation Jacobians from AD derivatives
					// Map global indices to local indices for masked variables
					std::vector<int> globalToLocal(mask.len, -1);
					int localIdx = 0;
					for (int i = 0; i < mask.len; ++i)
					{
						if (mask.mask[i])
						{
							globalToLocal[i] = localIdx++;
						}
					}
					
					// Extract relevant entries from AD Jacobian
					for (int localRowIdx = 0; localRowIdx < numExchangeEq; ++localRowIdx)
					{
						int globalRowIdx = -1;
						int localCount = 0;
						for (int i = 0; i < mask.len; ++i)
						{
							if (mask.mask[i])
							{
								if (localCount == localRowIdx)
								{
									globalRowIdx = i; 
									break;
								}
								localCount++;
							}
						}
						
						if (globalRowIdx >= 0 && globalRowIdx < mask.len)
						{
							// Get AD residual at global row index (offset by cShellOfset)
							const active& adResVal = localAdRes[cShellOfset + globalRowIdx];
							
							// Extract derivatives for masked columns
							for (int j = 0; j < mask.len; ++j)
							{
								if (mask.mask[j])
								{
									const int localColIdx = globalToLocal[j];
									if (localColIdx >= 0)
									{
										// Extract derivative with respect to variable at position (cShellOfset + j)
										mat.native(localRowIdx, localColIdx) = adResVal.getADValue(adJac.adDirOffset + j);
									}
								}
							}
						}
					}
					
					// Part 2: Set up conservation constraint equations in lower rows
					if (conservedQuants.size() > 0)
					{
						const int conservationStartRow = numExchangeEq;
						
						// For each conservation group, set up the constraint equation
						for (unsigned int conIdx = 0; conIdx < conservedQuants.size(); ++conIdx)
						{
							const int rowIdx = conservationStartRow + conIdx;
							const auto& currentGroup = conservationGroups[conIdx];
							
							// Set coefficients for channels in this conservation group
							for (unsigned int channel : currentGroup)
							{
								// Find the column index in the reduced system (mask space)
								const unsigned int maskIdx = channel * _disc.nComp + comp;
								
								if (maskIdx < mask.len && mask.mask[maskIdx])
								{
									const int localColIdx = globalToLocal[maskIdx];
									if (localColIdx >= 0)
									{
										mat.native(rowIdx, localColIdx) = -1.0;
									}
								}
							}
						}
					}

					return true;
				};
			}
			else
			{
				jacFunc = [&](double const* const x, linalg::detail::DenseMatrixBase& mat)
				{
					// Prepare input vector by overwriting masked items
					std::copy_n(cShell, mask.len, fullX + cShellOfset);
					linalg::applyVectorSubset(x, mask, fullX + cShellOfset);

					// Call residual function to compute Jacobian
					residualImpl<double, double, double, true>(simTime.t, simTime.secIdx, fullX, nullptr, fullResidual, threadLocalMem);

					// Clear the matrix
					mat.setAll(0.0);

					// Extract Jacobian submatrix corresponding to masked variables
					const int numExchangeEq = probSize - conservedQuants.size();
					
					// Part 1: Manually extract exchange equation Jacobians from sparse matrix
					// Map global indices to local indices for masked variables
					std::vector<int> globalToLocal(mask.len, -1);
					int localIdx = 0;
					for (int i = 0; i < mask.len; ++i)
					{
						if (mask.mask[i])
						{
							globalToLocal[i] = localIdx++;
						}
					}
					
					// Extract relevant entries from sparse Jacobian
					int extractedEntries = 0;
					for (int localRowIdx = 0; localRowIdx < numExchangeEq; ++localRowIdx)
					{
						int globalRowIdx = -1;
						int localCount = 0;
						for (int i = 0; i < mask.len; ++i)
						{
							if (mask.mask[i])
							{
								if (localCount == localRowIdx)
								{
									globalRowIdx = i; 
									break;
								}
								localCount++;
							}
						}
						
						//std::cout << "Local row " << localRowIdx << " -> Global row " << globalRowIdx << "\n";
						
						if (globalRowIdx >= 0 && globalRowIdx < testJac.rows())
						{
							// Get sparse row data from the compressed sparse matrix
							double const* const vals = testJac.valuesOfRow(globalRowIdx);
							const int nnz = testJac.numNonZerosInRow(globalRowIdx);
							int const* colIdx = testJac.columnIndicesOfRow(globalRowIdx);
							

							//std::cout << "  Row " << globalRowIdx << " has " << nnz << " non-zeros\n";
							
							// Copy relevant entries to dense matrix
							for (int j = 0; j < nnz; ++j)
							{
								const int globalColIdx = colIdx[j];
								//std::cout << "    Col[" << j << "]: global=" << colIdx[j] << " local=" << globalColIdx << "\n";
								
								if (globalColIdx >= 0 && globalColIdx < mask.len && mask.mask[globalColIdx])
								{
									const int localColIdx = globalToLocal[globalColIdx];
									if (localColIdx >= 0)
									{
										mat.native(localRowIdx, localColIdx) = vals[j];
										extractedEntries++;
										//std::cout << "      Setting mat(" << localRowIdx << "," << localColIdx << ") = " << vals[j] << "\n";
									}
								}
							}
						}
					}
					
					//std::cout << "Total extracted entries: " << extractedEntries << "\n";
					
					// Part 2: Set up conservation constraint equations in lower rows
					if (conservedQuants.size() > 0)
					{
						//std::cout << "Setting up conservation constraints...\n";
						const int conservationStartRow = numExchangeEq;
						
						// For each conservation group, set up the constraint equation
						for (unsigned int conIdx = 0; conIdx < conservedQuants.size(); ++conIdx)
						{
							const int rowIdx = conservationStartRow + conIdx;
							const auto& currentGroup = conservationGroups[conIdx];
							
							//std::cout << "  Conservation group " << conIdx << " has channels: ";
							//for (unsigned int ch : currentGroup) std::cout << ch << " ";
							//std::cout << "\n";
							
							// Set coefficients for channels in this conservation group
							for (unsigned int channel : currentGroup)
							{
								// Find the column index in the reduced system (mask space)
								const unsigned int maskIdx = channel * _disc.nComp + comp;
								
								if (maskIdx < mask.len && mask.mask[maskIdx])
								{
									const int localColIdx = globalToLocal[maskIdx];
									if (localColIdx >= 0)
									{
										mat.native(rowIdx, localColIdx) = -1.0;
										//std::cout << "    Setting conservation mat(" << rowIdx << "," << localColIdx << ") = 1.0\n";
									}
								}
							}
						}
					}

					//std::cout<< "Final Jacobian matrix for component " << comp << " in column " << pblk << ":\n";
					//for (int i = 0; i < mat.rows(); ++i)
					//{
					//	for (int j = 0; j < mat.columns(); ++j)
					//	{
					//		std::cout << mat.native(i, j) << " ";
					//	}
					//	std::cout << "\n";
					//}
					return true;
				};
			}
			std::copy_n(vecStateY, numDofs(), fullX);
			// Apply nonlinear solver
			_nonlinearSolver->solve(
				[&](double const* const x, double* const r)
				{
					// Prepare input vector by overwriting masked items
					std::copy_n(cShell, mask.len, fullX + cShellOfset);
					linalg::applyVectorSubset(x, mask, fullX + cShellOfset);

					// Call residual function
					residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, fullX, nullptr, fullResidual, threadLocalMem);

					// Extract values from residual 
					// r = [ci_channel1, ci_channel2,..., ci_channelNeq, cm1, ..., cmNEQ]
					// where ci_channeli is the value of the i-th component in the channel in rapid equilibrium
					// and cm1, ..., cmNEQ are the values of the conservation relations
					linalg::selectVectorSubset(fullResidual + cShellOfset, mask, r);					
					
					// and save them in the lower part of r
					unsigned int rIdx = probSize - conservedQuants.size();
					std::fill_n(r + rIdx, conservedQuants.size(), 0.0);
					
					// For this component, the conservation constraint is:
					// sum of concentrations in each quasi-stationary group = constant
					if (conservedQuants.size() > 0)
					{
						for (unsigned int conIdx = 0; conIdx < conservedQuants.size(); ++conIdx)
						{
							r[rIdx] = conservedQuants[conIdx]; // target conserved amount for this group
						
							// Subtract current concentrations in quasi-stationary channels for this specific group
							const auto& currentGroup = conservationGroups[conIdx];
							for (unsigned int channel : currentGroup)
							{
								// Find the corresponding index in the mask
								unsigned int maskIdx = channel * _disc.nComp + comp;
								if (maskIdx < mask.len && mask.mask[maskIdx])
								{
									// Find corresponding index in solution vector x
									unsigned int solIdx = 0;
									for (unsigned int j = 0; j <= maskIdx; ++j)
									{
										if (mask.mask[j]) solIdx++;
									}
									r[rIdx] -= x[solIdx - 1];
								}
							}
							++rIdx;
						}
					}

					return true;
				},
				jacFunc, errorTol, solution, nonlinMem, jacobianMatrix, probSize);

			// Apply solution
			linalg::applyVectorSubset(solution, mask, cShell);

		}
	}
	// Step 1b: Compute fluxes j_f
	// Reset j_f to 0.0
	/*double* const jf = vecStateY + idxr.offsetJf();
	std::fill(jf, jf + _disc.nComp * _disc.nCol * _disc.nParType, 0.0);

	solveForFluxes(vecStateY, idxr);*/
}

/**
 * @brief Computes consistent initial time derivatives
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *
 *          The process works in two steps:
 *          <ol>
 *              <li>Solve all algebraic equations in the model (e.g., quasi-stationary isotherms, reaction equilibria).
 *                 Once all @f$ c_i @f$, @f$ c_{p,i} @f$, and @f$ q_i^{(j)} @f$ have been computed, solve for the
 *                 fluxes @f$ j_{f,i} @f$ (only linear equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0.
 *                 However, because of the algebraic equations, we need additional conditions to fully determine
 *                 @f$ \dot{y}@f$. By differentiating the algebraic equations with respect to time, we get the
 *                 missing linear equations (recall that the state vector @f$ y @f$ is fixed). The resulting system
 *                 has a similar structure as the system Jacobian.
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                                & \dot{J}_1     &        &           &   \\
 *                                &         & \ddots &           &   \\
 *                                &         &        & \dot{J}_{N_z}   &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_i @f$ denotes the Jacobian with respect to @f$ \dot{y}@f$. Note that the
 *                 @f$ J_{i,f} @f$ matrices in the right column are missing.
 *
 *     The right hand side of the linear system is given by the negative residual without contribution
 *     of @f$ \dot{y} @f$ for differential equations and 0 for algebraic equations
 *     (@f$ -\frac{\partial F}{\partial t}@f$, to be more precise).
 *
 *     The linear system is solved by backsubstitution. First, the diagonal blocks are solved in parallel.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     diagonal blocks.</li>
 *          </ol>
 *     This function performs step 2. See consistentInitialState() for step 1.
 *
 * 	   This function is to be used with consistentInitialState(). Do not mix normal and lean
 *     consistent initialization!
 *
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in] vecStateY Consistently initialized state vector
 * @param [in,out] vecStateYdot On entry, residual without taking time derivatives into account. On exit, consistent state time derivatives.
 */
void MultiChannelTransportModel::consistentInitialTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	// Step 2: Compute the correct time derivative of the state vector

	// Step 2a: Assemble, factorize, and solve diagonal blocks of linear system

	// Note that the residual has not been negated, yet. We will do that now.
	for (unsigned int i = 0; i < numDofs(); ++i)
		vecStateYdot[i] = -vecStateYdot[i];

	// Handle bulk column block
	_convDispOp.solveTimeDerivativeSystem(simTime, vecStateYdot + idxr.offsetC());

	// Process quasi-stationary exchange constraints for each component
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
	{
		if (!_exchange[0]->hasQuasiStationary(comp))
			continue;

		// Create quasi-stationary mask for this component with conservation groups
		std::vector<int> qsMask;
		std::vector<std::vector<unsigned int>> conservationGroups;
		const int numActiveMask = createQuasiStationaryMaskWithGroups(comp, _exchange[0], _disc.nChannel, _disc.nComp, qsMask, conservationGroups);

		if (conservationGroups.empty())
			continue;

#ifdef CADET_PARALLELIZE
		BENCH_START(_timerConsistentInitPar);
		tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nCol),
			[comp, qsMask, conservationGroups, mask, probSize, dFluxDt, &simTime, &vecStateY, &vecStateYdot, &threadLocalMem, this](std::size_t pblk)
#else
		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
#endif
		{
			LinearBufferAllocator tlmAlloc = threadLocalMem.get();
			const double z = (0.5 + static_cast<double>(pblk)) / static_cast<double>(_disc.nCol);


			// Get offset for this column block
			auto cShellOffset = _disc.nComp * _disc.nChannel * (1 + pblk);
			double* const cShell = const_cast<double*>(vecStateY) + cShellOffset;
			double* const cShellDot = vecStateYdot + cShellOffset;

			const linalg::ConstMaskArray mask{ qsMask.data(), static_cast<int>(_disc.nChannel * _disc.nComp) };
			const int probSize = linalg::numMaskActive(mask);

			// Initialize Jacobian matrix for this component/column block
			if (probSize <= 0 || conservationGroups.empty())
				continue;


			// Compute time derivative of quasi-stationary fluxes (if exchange depends on time)
			BufferedArray<double> dFluxDtBuffer = tlmAlloc.array<double>(mask.len);
			double* const dFluxDt = static_cast<double*>(dFluxDtBuffer);
			std::fill_n(dFluxDt, mask.len, 0.0);

			if (_exchange[0]->dependsOnTime())
			{
				// TODO: Call exchange model's time derivative function if it exists
				// This is analogous to _binding[type]->timeDerivativeQuasiStationaryFluxes in LRM
				_exchange[0]->timeDerivativeQuasiStationaryExchange(simTime.t, simTime.secIdx,
					ColumnPosition{ z, 0.0, 0.0 }, cShell, dFluxDt, tlmAlloc);
			}

			// Allocate memory for factorizable matrix (similar to LRM's _jacPdisc)
			BufferedArray<double> jacMatrixBuffer = tlmAlloc.array<double>(probSize * probSize);
			linalg::DenseMatrixView jacMatrix(static_cast<double*>(jacMatrixBuffer),
				new lapackInt_t[probSize], probSize, probSize);

			// Set up matrix similar to LRM's approach
			jacMatrix.setAll(0.0);

			// Midpoint of current column cell (z coordinate) - needed in externally dependent exchange

			// Map global indices to local indices for masked variables
			std::vector<int> globalToLocal(mask.len, -1);
			int localIdx = 0;
			for (int i = 0; i < mask.len; ++i)
			{
				if (mask.mask[i])
				{
					globalToLocal[i] = localIdx++;
				}
			}

			// Add time derivative contributions to Jacobian matrix
			// This is the REAL Jacobian ∂F/∂ẏ, NOT identity matrix!
			const int numExchangeEq = probSize - conservationGroups.size();

			// Part 1: Calculate time derivative jacobian
			BufferedArray<double> tempJacBuffer = tlmAlloc.array<double>(mask.len * mask.len);
			linalg::DenseMatrixView tempJac(static_cast<double*>(tempJacBuffer),
				new lapackInt_t[mask.len], mask.len, mask.len);
			tempJac.setAll(0.0);


			ColumnPosition colPos{ z, 0.0, 0.0 };

			//jacobianTimeDerivative(simTime.t, simTime.secIdx, colPos, comp, 
			//     cShell, cShellDot, tempJac.data(), tlmAlloc);

			// for now, calculate explizit
			for (int i = 0; i < numExchangeEq; ++i)
			{
				// Finde entsprechende globale Indizes
				int globalRow = -1, globalCol = -1;
				int count = 0;
				for (int j = 0; j < mask.len; ++j)
				{
					if (mask.mask[j])
					{
						if (count == i) { globalRow = j; break; }
						count++;
					}
				}

				if (globalRow >= 0)
				{
					jacMatrix.native(i, i) = 1.0;
				}
			}

			if (conservationGroups.size() > 0)
			{
				// Set up constraint rows in Jacobian (ersetze quasi-stationäre Zeilen)
				for (unsigned int conIdx = 0; conIdx < conservationGroups.size(); ++conIdx)
				{
					const int rowIdx = numExchangeEq + conIdx;
					const auto& currentGroup = conservationGroups[conIdx];

					// Clear the row first
					for (int col = 0; col < probSize; ++col)
						jacMatrix.native(rowIdx, col) = 0.0;

					// Set conservation constraint: sum of time derivatives in group = 0
					for (unsigned int channel : currentGroup)
					{
						const unsigned int maskIdx = channel * _disc.nComp + comp;
						if (maskIdx < mask.len && mask.mask[maskIdx])
						{
							const int localColIdx = globalToLocal[maskIdx];
							if (localColIdx >= 0)
							{
								jacMatrix.native(rowIdx, localColIdx) = 1.0;
							}
						}
					}
				}

				// Set right hand side for conservation constraints
				for (unsigned int conIdx = 0; conIdx < conservationGroups.size(); ++conIdx)
				{
					const auto& currentGroup = conservationGroups[conIdx];

					// Calculate the constraint RHS: time derivative of exchange flux for this group
					double rhsValue = 0.0;
					for (unsigned int channel : currentGroup)
					{
						const unsigned int maskIdx = channel * _disc.nComp + comp;
						if (maskIdx < mask.len)
						{
							rhsValue -= dFluxDt[maskIdx]; // Add time derivative of exchange flux
						}
					}

					// Apply constraint RHS to representative channel
					if (!currentGroup.empty())
					{
						unsigned int reprChannel = currentGroup[0];
						const unsigned int reprMaskIdx = reprChannel * _disc.nComp + comp;
						if (reprMaskIdx < mask.len && mask.mask[reprMaskIdx])
						{
							cShellDot[reprMaskIdx] = rhsValue;
						}
					}
				}
			}

			// Precondition, factorize, and solve
			BufferedArray<double> scaleFactorsBuffer = tlmAlloc.array<double>(probSize);
			double* const scaleFactors = static_cast<double*>(scaleFactorsBuffer);

			jacMatrix.rowScaleFactors(scaleFactors);
			jacMatrix.scaleRows(scaleFactors);

			// Factorize
			const bool result = jacMatrix.factorize();
			if (!result)
			{
				LOG(Error) << "Factorize() failed for component " << comp << " column block " << pblk;
			}

			// Extract RHS vector from current time derivatives
			BufferedArray<double> rhsBuffer = tlmAlloc.array<double>(probSize);
			double* const rhs = static_cast<double*>(rhsBuffer);
			linalg::selectVectorSubset(cShellDot, mask, rhs);

			// Solve
			const bool result2 = jacMatrix.solve(scaleFactors, rhs);
			if (!result2)
			{
				LOG(Error) << "Solve() failed for component " << comp << " column block " << pblk;
			}

			// Apply solution back to time derivatives
			linalg::applyVectorSubset(rhs, mask, cShellDot);

		} CADET_PARFOR_END;

#ifdef CADET_PARALLELIZE
		BENCH_STOP(_timerConsistentInitPar);
#endif
	}
}

/**
 * @brief Computes approximately / partially consistent initial values (state variables without their time derivatives)
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *
 *          This function performs a relaxed consistent initialization: Only parts of the vectors are updated
 *          and, hence, consistency is not guaranteed. Since there is less work to do, it is probably faster than
 *          the standard process represented by consistentInitialState().
 *
 *          The process works in two steps:
 *          <ol>
 *              <li>Keep state and time derivative vectors as they are (i.e., do not solve algebraic equations).
 *                 Only solve for the fluxes @f$ j_{f,i} @f$ (only linear equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0 in the column
 *                 bulk and flux blocks. The resulting equations are stated below:
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_0 @f$ denotes the bulk block Jacobian with respect to @f$ \dot{y}@f$.
 *
 *     The right hand side of the linear system is given by the negative residual without contribution
 *     of @f$ \dot{y} @f$ for the bulk block and 0 for the flux block.
 *
 *     The linear system is solved by backsubstitution. First, the bulk block is solved.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     bulk block and the unchanged particle block time derivative vectors.</li>
 *          </ol>
 *     This function performs step 1. See leanConsistentInitialTimeDerivative() for step 2.
 *
 * 	   This function is to be used with leanConsistentInitialTimeDerivative(). Do not mix normal and lean
 *     consistent initialization!
 *
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in,out] vecStateY State vector with initial values that are to be updated for consistency
 * @param [in,out] adJac Jacobian information for AD (AD vectors for residual and state, direction offset)
 * @param [in] errorTol Error tolerance for algebraic equations
 */
void MultiChannelTransportModel::leanConsistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	consistentInitialState(simTime, vecStateY, adJac, errorTol, threadLocalMem);
}

/**
 * @brief Computes approximately / partially consistent initial time derivatives
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *
 *          This function performs a relaxed consistent initialization: Only parts of the vectors are updated
 *          and, hence, consistency is not guaranteed. Since there is less work to do, it is probably faster than
 *          the standard process represented by consistentInitialTimeDerivative().
 *
 *          The process works in two steps:
 *          <ol>
 *              <li>Keep state and time derivative vectors as they are (i.e., do not solve algebraic equations).
 *                 Only solve for the fluxes @f$ j_{f,i} @f$ (only linear equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0 in the column
 *                 bulk and flux blocks. The resulting equations are stated below:
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_0 @f$ denotes the bulk block Jacobian with respect to @f$ \dot{y}@f$.
 *
 *     The right hand side of the linear system is given by the negative residual without contribution
 *     of @f$ \dot{y} @f$ for the bulk block and 0 for the flux block.
 *
 *     The linear system is solved by backsubstitution. First, the bulk block is solved.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     bulk block and the unchanged particle block time derivative vectors.</li>
 *          </ol>
 *     This function performs step 2. See leanConsistentInitialState() for step 1.
 *
 * 	   This function is to be used with leanConsistentInitialState(). Do not mix normal and lean
 *     consistent initialization!
 *
 * @param [in] t Current time point
 * @param [in] vecStateY (Lean) consistently initialized state vector
 * @param [in,out] vecStateYdot On entry, inconsistent state time derivatives. On exit, partially consistent state time derivatives.
 * @param [in] res On entry, residual without taking time derivatives into account. The data is overwritten during execution of the function.
 */
void MultiChannelTransportModel::leanConsistentInitialTimeDerivative(double t, double const* const vecStateY, double* const vecStateYdot, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	// Step 2: Compute the correct time derivative of the state vector

	// Step 2a: Assemble, factorize, and solve column bulk block of linear system

	// Note that the residual is not negated as required at this point. We will fix that later.

	double* const resSlice = res + idxr.offsetC();

	// Handle bulk block
	_convDispOp.solveTimeDerivativeSystem(SimulationTime{t, 0u}, resSlice);

	// Note that we have solved with the *positive* residual as right hand side
	// instead of the *negative* one. Fortunately, we are dealing with linear systems,
	// which means that we can just negate the solution.
	double* const yDotSlice = vecStateYdot + idxr.offsetC();
	for (unsigned int i = 0; i < _disc.nCol * _disc.nChannel * _disc.nComp; ++i)
		yDotSlice[i] = -resSlice[i];
}

void MultiChannelTransportModel::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
{
	Indexer idxr(_disc);
	for (std::size_t param = 0; param < vecSensY.size(); ++param)
	{
		double* const stateYbulk = vecSensY[param] + idxr.offsetC();

		// Loop over axial cells
		for (unsigned int col = 0; col < _disc.nCol; ++col)
		{
			// Loop over radial cells
			for (unsigned int rad = 0; rad < _disc.nChannel; ++rad)
			{
				// Loop over components in cell
				for (unsigned comp = 0; comp < _disc.nComp; ++comp)
					stateYbulk[col * idxr.strideColAxialCell() + rad * idxr.strideChannelCell() + comp * idxr.strideColComp()] = _initC[comp + rad * _disc.nComp].getADValue(param);
			}
		}
	}
}

/**
 * @brief Computes consistent initial values and time derivatives of sensitivity subsystems
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] and initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$,
 *          the sensitivity system for a parameter @f$ p @f$ reads
 *          \f[ \frac{\partial F}{\partial y}(t, y, \dot{y}) s + \frac{\partial F}{\partial \dot{y}}(t, y, \dot{y}) \dot{s} + \frac{\partial F}{\partial p}(t, y, \dot{y}) = 0. \f]
 *          The initial values of this linear DAE, @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p} @f$
 *          have to be consistent with the sensitivity DAE. This functions updates the initial sensitivity\f$ s_0 \f$ and \f$ \dot{s}_0 \f$ such that they are consistent.
 *
 *          The process follows closely the one of consistentInitialConditions() and, in fact, is a linearized version of it.
 *          This is necessary because the initial conditions of the sensitivity system \f$ s_0 \f$ and \f$ \dot{s}_0 \f$ are
 *          related to the initial conditions \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ of the original DAE by differentiating them
 *          with respect to @f$ p @f$: @f$ s_0 =
	for (std::size_t param = 0; param < vecSensY.size(); ++param)
	{
		double* const sensY = vecSensY[param];
		double* const sensYdot = vecSensYdot[param];

		// Copy parameter derivative dF / dp from AD and negate it
		for (unsigned int i = _disc.nComp * _disc.nChannel; i < numDofs(); ++i)
			sensYdot[i] = -adRes[i].getADValue(param);

		// Step 1: Solve algebraic equations

		// Step 2: Compute the correct time derivative of the state vector

		// Step 2a: Assemble, factorize, and solve diagonal blocks of linear system

		// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in sensYdot

		// Note that we have correctly negated the right hand side

		// Handle bulk block
		_convDispOp.solveTimeDerivativeSystem(simTime, sensYdot + idxr.offsetC());
	}
}

/**
 * @brief Computes approximately / partially consistent initial values and time derivatives of sensitivity subsystems
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] and initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$,
 *          the sensitivity system for a parameter @f$ p @f$ reads
 *          \f[ \frac{\partial F}{\partial y}(t, y, \dot{y}) s + \frac{\partial F}{\partial \dot{y}}(t, y, \dot{y}) \dot{s} + \frac{\partial F}{\partial p}(t, y, \dot{y}) = 0. \f]
 *          The initial values of this linear DAE, @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p} @f$
 *          have to be consistent with the sensitivity DAE. This functions updates the initial sensitivity\f$ s_0 \f$ and \f$ \dot{s}_0 \f$ such that they are consistent.
 *
 *          The process follows closely the one of leanConsistentInitialConditions() and, in fact, is a linearized version of it.
 *          This is necessary because the initial conditions of the sensitivity system \f$ s_0 \f$ and \f$ \dot{s}_0 \f$ are
 *          related to the initial conditions \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ of the original DAE by differentiating them
 *          with respect to @f$ p @f$: @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p}. @f$
 *          <ol>
 *              <li>Keep state and time derivative vectors as they are (i.e., do not solve algebraic equations).
 *                 Only solve for the fluxes @f$ j_{f,i} @f$ (only linear equations).</li>
 *              <li>Compute the time derivatives of the sensitivity @f$ \dot{s} @f$ such that the differential equations hold.
 *                 However, because of the algebraic equations, we need additional conditions to fully determine
 *                 @f$ \dot{s}@f$. By differentiating the algebraic equations with respect to time, we get the
 *                 missing linear equations (recall that the sensitivity vector @f$ s @f$ is fixed). The resulting
 *                 equations are stated below:
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_0 @f$ denotes the bulk block Jacobian with respect to @f$ \dot{y}@f$.
 *
 *     Let @f$ \mathcal{I}_d @f$ denote the index set of differential equations.
 *     The right hand side of the linear system is given by @f[ -\frac{\partial F}{\partial y}(t, y, \dot{y}) s - \frac{\partial F}{\partial p}(t, y, \dot{y}), @f]
 *     which is 0 for algebraic equations (@f$ -\frac{\partial^2 F}{\partial t \partial p}@f$, to be more precise).
 *
 *     The linear system is solved by backsubstitution. First, the bulk block is solved.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     bulk block and the unchanged particle block time derivative vectors.</li>
 *          </ol>
 *     This function requires the parameter sensitivities to be computed beforehand and up-to-date Jacobians.
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in] simState Consistent state of the simulation (state vector and its time derivative)
 * @param [in,out] vecSensY Sensitivity subsystem state vectors
 * @param [in,out] vecSensYdot Time derivative state vectors of the sensitivity subsystems to be initialized
 * @param [in] adRes Pointer to residual vector of AD datatypes with parameter sensitivities
 * @todo Decrease amount of allocated memory by partially using temporary vectors (state and Schur complement)
 */
void MultiChannelTransportModel::leanConsistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	for (std::size_t param = 0; param < vecSensY.size(); ++param)
	{
		double* const sensY = vecSensY[param];
		double* const sensYdot = vecSensYdot[param];

		// Copy parameter derivative from AD to tempState and negate it
		// We need to use _tempState in order to keep sensYdot unchanged at this point
		for (unsigned int i = 0; i < numDofs(); ++i)
			_tempState[i] = -adRes[i].getADValue(param);

		// Step 2: Compute the correct time derivative of the state vector

		// Step 2a: Assemble, factorize, and solve bulk block of linear system

		// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in _tempState
		multiplyWithJacobian(simTime, simState, sensY, -1.0, 1.0, _tempState);

		// Copy relevant parts to sensYdot for use as right hand sides
		std::copy(_tempState + idxr.offsetC(), _tempState + numDofs(), sensYdot + idxr.offsetC());

		// Handle bulk block
		_convDispOp.solveTimeDerivativeSystem(simTime, sensYdot + idxr.offsetC());
	}
}

void MultiChannelTransportModel::consistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
{}

}  // namespace model

}  // namespace cadet

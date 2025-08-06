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

#include "model/parts/ParticleDiffusionOperatorFV.hpp"
#include "cadet/Exceptions.hpp"
#include "ParamReaderHelper.hpp"
#include "ParamReaderScopes.hpp"
#include "AdUtils.hpp"
#include "SimulationTypes.hpp"
#include "model/ParameterDependence.hpp"
#include "SensParamUtil.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

namespace cadet
{

namespace model
{

namespace parts
{
	/**
	 * @brief Creates a ParticleDiffusionOperatorFV
	 */
	ParticleDiffusionOperatorFV::ParticleDiffusionOperatorFV()
	{
	}

	ParticleDiffusionOperatorFV::~ParticleDiffusionOperatorFV() CADET_NOEXCEPT
	{
	}

	bool ParticleDiffusionOperatorFV::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
	{
		const bool baseConfigSuccess = ParticleDiffusionOperatorBase::configureModelDiscretization(paramProvider, helper, nComp, parTypeIdx, nParType, strideBulkComp);

		paramProvider.pushScope("discretization");

		_nParPoints = paramProvider.getInt("NCELLS");
		if (_nParPoints < 1)
			throw InvalidParameterException("Field NCELLS in discretization of particle type " + std::to_string(_parTypeIdx) + " must be > 0 but is " + std::to_string(_nParPoints));
		// Default boundary order handling to second order
		_boundaryOrderFV = 2;
		if (paramProvider.exists("FV_BOUNDARY_ORDER"))
		{
			_boundaryOrderFV = paramProvider.getInt("FV_BOUNDARY_ORDER");
			if ((_boundaryOrderFV < 1) || (_boundaryOrderFV > 2))
				throw InvalidParameterException("Field FV_BOUNDARY_ORDER is out of valid range (1 or 2)");
		}

		// Configure particle discretization
		_parCenterRadius.resize(_nParPoints);
		_parOuterSurfAreaPerVolume.resize(_nParPoints);
		_parInnerSurfAreaPerVolume.resize(_nParPoints);

		// Read particle discretization mode and default to "EQUIDISTANT_PAR"
		_parDiscMode = ParticleDiscretizationMode::Equidistant;
		std::string pdt = paramProvider.getString("PAR_DISC_TYPE");

		if (pdt == "EQUIVOLUME_PAR")
			_parDiscMode = ParticleDiscretizationMode::Equivolume;
		else if (pdt == "USER_DEFINED_PAR")
			_parDiscMode = ParticleDiscretizationMode::UserDefined;

		if (paramProvider.exists("PAR_DISC_VECTOR"))
		{
			_parDiscVector = paramProvider.getDoubleArray("PAR_DISC_VECTOR");
			if (_parDiscVector.size() < _nParPoints)
				throw InvalidParameterException("Field PAR_DISC_VECTOR contains too few elements (" + std::to_string(_nParPoints) + " required)");
		}

		paramProvider.popScope(); // discretization

		_discParFlux.resize(sizeof(active) * _nComp);

		return baseConfigSuccess;
	}

	bool ParticleDiffusionOperatorFV::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound, const int* reqBinding)
	{
		const bool baseConfigSuccess = ParticleDiffusionOperatorBase::configure(unitOpIdx, paramProvider, parameters, nParType, nBoundBeforeType, nTotalBound, reqBinding);

		// Compute particle metrics
		_deltaR.resize(_nParPoints);
		updateRadialDisc();

		return baseConfigSuccess;
	}

	/**
	 * @brief Computes equidistant radial nodes in the beads
	 */
	void ParticleDiffusionOperatorFV::setEquidistantRadialDisc()
	{
		active* const ptrCenterRadius = _parCenterRadius.data();
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data();
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data();

		const active radius = _parRadius - _parCoreRadius;
		const active dr = radius / static_cast<double>(_nParPoints);
		std::fill(_deltaR.data(), _deltaR.data() + _nParPoints, dr);

		if (_parGeomSurfToVol == _SurfVolRatioSphere)
		{
			for (int cell = 0; cell < _nParPoints; ++cell)
			{
				const active r_out = _parCoreRadius + static_cast<double>(cell + 1) * dr;
				const active r_in = _parCoreRadius + static_cast<double>(cell) * dr;

				ptrCenterRadius[cell] = _parCoreRadius + (0.5 + static_cast<double>(cell)) * dr;

				// Compute denominator -> corresponding to cell volume
				const active vol = pow(r_out, 3.0) - pow(r_in, 3.0);

				ptrOuterSurfAreaPerVolume[cell] = 3.0 * sqr(r_out) / vol;
				ptrInnerSurfAreaPerVolume[cell] = 3.0 * sqr(r_in) / vol;
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioCylinder)
		{
			for (int cell = 0; cell < _nParPoints; ++cell)
			{
				const active r_out = _parCoreRadius + static_cast<double>(cell + 1) * dr;
				const active r_in = _parCoreRadius + static_cast<double>(cell) * dr;

				ptrCenterRadius[cell] = _parCoreRadius + (0.5 + static_cast<double>(cell)) * dr;

				// Compute denominator -> corresponding to cell volume
				const active vol = sqr(r_out) - sqr(r_in);

				ptrOuterSurfAreaPerVolume[cell] = 3.0 * sqr(r_out) / vol;
				ptrInnerSurfAreaPerVolume[cell] = 3.0 * sqr(r_in) / vol;
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioSlab)
		{
			for (int cell = 0; cell < _nParPoints; ++cell)
			{
				const active r_out = _parCoreRadius + static_cast<double>(cell + 1) * dr;
				const active r_in = _parCoreRadius + static_cast<double>(cell) * dr;

				ptrCenterRadius[cell] = _parCoreRadius + (0.5 + static_cast<double>(cell)) * dr;

				// Compute denominator -> corresponding to cell volume
				const active vol = r_out - r_in;

				ptrOuterSurfAreaPerVolume[cell] = 1.0 / vol;
				ptrInnerSurfAreaPerVolume[cell] = 1.0 / vol;
			}
		}
	}

	/**
	 * @brief Computes the radial nodes in the beads in such a way that all shells have the same volume
	 */
	void ParticleDiffusionOperatorFV::setEquivolumeRadialDisc()
	{
		active* const ptrCellSize = _deltaR.data();
		active* const ptrCenterRadius = _parCenterRadius.data();
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data();
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data();

		if (_parGeomSurfToVol == _SurfVolRatioSphere)
		{
			active r_out = _parRadius;
			active r_in = _parCoreRadius;
			const active volumePerShell = (pow(_parRadius, 3.0) - pow(_parCoreRadius, 3.0)) / static_cast<double>(_nParPoints);

			for (int cell = 0; cell < _nParPoints; ++cell)
			{
				if (cell != _nParPoints)
					r_out = pow(volumePerShell + pow(r_in, 3.0), (1.0 / 3.0));
				else
					r_out = _parRadius;

				ptrCellSize[cell] = r_out - r_in;
				ptrCenterRadius[cell] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[cell] = 3.0 * sqr(r_out) / volumePerShell;
				ptrInnerSurfAreaPerVolume[cell] = 3.0 * sqr(r_in) / volumePerShell;

				// For the next cell: r_in == r_out of the current cell
				r_in = r_out;
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioCylinder)
		{
			active r_out = _parRadius;
			active r_in = _parCoreRadius;
			const active volumePerShell = (sqr(_parRadius) - sqr(_parCoreRadius)) / static_cast<double>(_nParPoints);

			for (int cell = 0; cell < _nParPoints; ++cell)
			{
				if (cell != (_nParPoints - 1))
					r_out = sqrt(volumePerShell + sqr(r_in));
				else
					r_out = _parRadius;

				ptrCellSize[cell] = r_out - r_in;
				ptrCenterRadius[cell] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[cell] = 2.0 * r_out / volumePerShell;
				ptrInnerSurfAreaPerVolume[cell] = 2.0 * r_in / volumePerShell;

				// For the next cell: r_in == r_out of the current cell
				r_in = r_out;
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioSlab)
		{
			active r_out = _parRadius;
			active r_in = _parCoreRadius;
			const active volumePerShell = (_parRadius - _parCoreRadius) / static_cast<double>(_nParPoints);

			for (int cell = 0; cell < _nParPoints; ++cell)
			{
				if (cell != (_nParPoints - 1))
					r_out = volumePerShell + r_in;
				else
					r_out = _parRadius;

				ptrCellSize[cell] = r_out - r_in;
				ptrCenterRadius[cell] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[cell] = 1.0 / volumePerShell;
				ptrInnerSurfAreaPerVolume[cell] = 1.0 / volumePerShell;

				// For the next cell: r_in == r_out of the current cell
				r_in = r_out;
			}
		}
	}

	/**
	 * @brief Computes all helper quantities for radial bead discretization from given radial cell boundaries
	 * @details Calculates surface areas per volume for every shell and the radial shell centers.
	 */
	void ParticleDiffusionOperatorFV::setUserdefinedRadialDisc()
	{
		active* const ptrCellSize = _deltaR.data();
		active* const ptrCenterRadius = _parCenterRadius.data();
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data();
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data();

		// Care for the right ordering and include 0.0 / 1.0 if not already in the vector.
		std::vector<active> orderedInterfaces = std::vector<active>(_parDiscVector.begin(), _parDiscVector.begin() + _nParPoints + 1);

		// Sort in ascending order
		std::sort(orderedInterfaces.begin(), orderedInterfaces.end(), std::less<active>());

		// Force first and last element to be 1.0 and 0.0, respectively
		orderedInterfaces[0] = 0.0;
		orderedInterfaces.back() = 1.0;

		// Map [0, 1] -> [core radius, particle radius] via linear interpolation
		for (int cell = 0; cell < _nParPoints; ++cell)
			orderedInterfaces[cell] = static_cast<double>(orderedInterfaces[cell]) * (_parRadius - _parCoreRadius) + _parCoreRadius;

		if (_parGeomSurfToVol == _SurfVolRatioSphere)
		{
			for (int cell = 0; cell < _nParPoints; ++cell)
			{
				ptrCellSize[cell] = orderedInterfaces[cell + 1] - orderedInterfaces[cell];
				ptrCenterRadius[cell] = (orderedInterfaces[cell + 1] + orderedInterfaces[cell]) * 0.5;

				// Compute denominator -> corresponding to cell volume
				const active vol = pow(orderedInterfaces[cell + 1], 3.0) - pow(orderedInterfaces[cell], 3.0);

				ptrOuterSurfAreaPerVolume[cell] = 3.0 * sqr(orderedInterfaces[cell + 1]) / vol;
				ptrInnerSurfAreaPerVolume[cell] = 3.0 * sqr(orderedInterfaces[cell]) / vol;
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioCylinder)
		{
			for (int cell = 0; cell < _nParPoints; ++cell)
			{
				ptrCellSize[cell] = orderedInterfaces[cell + 1] - orderedInterfaces[cell];
				ptrCenterRadius[cell] = (orderedInterfaces[cell + 1] + orderedInterfaces[cell]) * 0.5;

				// Compute denominator -> corresponding to cell volume
				const active vol = sqr(orderedInterfaces[cell + 1]) - sqr(orderedInterfaces[cell]);

				ptrOuterSurfAreaPerVolume[cell] = 2.0 * orderedInterfaces[cell + 1] / vol;
				ptrInnerSurfAreaPerVolume[cell] = 2.0 * orderedInterfaces[cell] / vol;
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioSlab)
		{
			for (int cell = 0; cell < _nParPoints; ++cell)
			{
				ptrCellSize[cell] = orderedInterfaces[cell + 1] - orderedInterfaces[cell];
				ptrCenterRadius[cell] = (orderedInterfaces[cell + 1] + orderedInterfaces[cell]) * 0.5;

				// Compute denominator -> corresponding to cell volume
				const active vol = orderedInterfaces[cell + 1] - orderedInterfaces[cell];

				ptrOuterSurfAreaPerVolume[cell] = 1.0 / vol;
				ptrInnerSurfAreaPerVolume[cell] = 1.0 / vol;
			}
		}
	}

	void ParticleDiffusionOperatorFV::updateRadialDisc()
	{
		if (_parDiscMode == ParticleDiscretizationMode::Equidistant)
			setEquidistantRadialDisc();
		else if (_parDiscMode == ParticleDiscretizationMode::Equivolume)
			setEquivolumeRadialDisc();
		else if (_parDiscMode == ParticleDiscretizationMode::UserDefined)
			setUserdefinedRadialDisc();
	}

	bool ParticleDiffusionOperatorFV::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx)
	{
		return ParticleDiffusionOperatorBase::notifyDiscontinuousSectionTransition(t, secIdx);
	}
	/**
	 * @brief calculates the physical radial/particle coordinates of the DG discretization with double! interface nodes
	 */
	int ParticleDiffusionOperatorFV::writeParticleCoordinates(double* coords) const
	{
		active const* const pcr = _parCenterRadius.data();
		for (int i = 0; i < _nParPoints; ++i)
			coords[i] = static_cast<double>(pcr[i]);
		return _nParPoints;
	}

	int ParticleDiffusionOperatorFV::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithoutParamSensitivity)
	{
		if (resPar)
		{
			if (jacIt.data())
				return residualImpl<double, double, double, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, jacIt);
			else
				return residualImpl<double, double, double, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, jacIt);
		}
		else if (jacIt.data())
			return residualImpl<double, double, double, true, false>(t, secIdx, yPar, yBulk, yDotPar, resPar, jacIt);
		else
			return -1;
	}
	int ParticleDiffusionOperatorFV::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithParamSensitivity)
	{
		 if (jacIt.data())
			return residualImpl<double, active, active, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, jacIt);
		else
			 return residualImpl<double, active, active, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, jacIt);
	}
	int ParticleDiffusionOperatorFV::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithoutParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<active, active, double, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, jacIt);
		else
			return residualImpl<active, active, double, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, jacIt);
	}
	int ParticleDiffusionOperatorFV::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<active, active, active, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, jacIt);
		else
			return residualImpl<active, active, active, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, jacIt);
	}

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
	int ParticleDiffusionOperatorFV::residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, linalg::BandedEigenSparseRowIterator& jacBase)
	{
		const active* const filmDiff = getSectionDependentSlice(_filmDiffusion, _nComp, secIdx);
		active const* const parDiff = getSectionDependentSlice(_parDiffusion, _nComp, secIdx);
		active const* const parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound, secIdx);

		// get temporary pointer to move through for particle diffusion
		StateType const* y = yPar;
		ResidualType* res = resPar;

		active const* const outerSurfPerVol = _parOuterSurfAreaPerVolume.data();
		active const* const innerSurfPerVol = _parInnerSurfAreaPerVolume.data();
		active const* const parCenterRadius = _parCenterRadius.data();

		// Loop over particle cells
		for (int par = 0; par < _nParPoints; ++par)
		{
			// Geometry
			const ParamType outerAreaPerVolume = static_cast<ParamType>(outerSurfPerVol[par]);
			const ParamType innerAreaPerVolume = static_cast<ParamType>(innerSurfPerVol[par]);

			// Mobile phase
			for (int comp = 0; comp < _nComp; ++comp, ++y, ++jacBase)
			{
				const unsigned int nBound = _nBound[comp];
				const ParamType invBetaP = (1.0 - static_cast<ParamType>(_parPorosity)) / (static_cast<ParamType>(_poreAccessFactor[comp]) * static_cast<ParamType>(_parPorosity));

				const ParamType dp = static_cast<ParamType>(parDiff[comp]);

				// Add flow through outer surface
				// Note that inflow boundary conditions are handled in residualFlux().
				if (cadet_likely(par != _nParPoints - 1))
				{
					// Difference between two cell-centers
					const ParamType dr = static_cast<ParamType>(parCenterRadius[par + 1]) - static_cast<ParamType>(parCenterRadius[par]);

					// Molecular diffusion contribution
					const ResidualType gradCp = (y[strideParPoint()] - y[0]) / dr;
					if (wantRes)
						*res -= outerAreaPerVolume * dp * gradCp;

					// Surface diffusion contribution for quasi-stationary bound states
					if (cadet_unlikely(_hasSurfaceDiffusion))
					{
						for (int i = 0; i < nBound; ++i)
						{
							// Index explanation:
							//   - comp go back to beginning of liquid phase
							//   + strideParLiquid skip over liquid phase to solid phase
							//   + offsetBoundComp jump to component comp (skips all bound states of previous components)
							//   + i go to current bound state
							const int curIdx = strideParLiquid() - comp + offsetBoundComp()[comp] + i;
							const ResidualType gradQ = (y[strideParPoint() + curIdx] - y[curIdx]) / dr;
							if (wantRes)
								*res -= outerAreaPerVolume * static_cast<ParamType>(parSurfDiff[offsetBoundComp()[comp] + i]) * invBetaP * gradQ;
						}

						if (wantJac)
						{
							const double localInvBetaP = static_cast<double>(invBetaP);
							const double ouApV = static_cast<double>(outerAreaPerVolume);
							const double ldr = static_cast<double>(dr);

							// Liquid phase
							jacBase[0] += ouApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j)
							jacBase[strideParPoint()] += -ouApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j-1)

							// Solid phase
							for (int i = 0; i < nBound; ++i)
							{
								// See above for explanation of curIdx value
								const int curIdx = strideParLiquid() - comp + offsetBoundComp()[comp] + i;
								jacBase[curIdx] += ouApV * localInvBetaP * static_cast<double>(parSurfDiff[offsetBoundComp()[comp] + i]) / ldr; // dres / dq_i^(p,j)
								jacBase[strideParPoint() + curIdx] += -ouApV * localInvBetaP * static_cast<double>(parSurfDiff[offsetBoundComp()[comp] + i]) / ldr; // dres / dq_i^(p,j-1)
							}
						}
					}
					else if (wantJac)
					{
						// No surface diffusion
						// Liquid phase
						const double ouApV = static_cast<double>(outerAreaPerVolume);
						const double ldr = static_cast<double>(dr);

						jacBase[0] += ouApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j)
						jacBase[strideParPoint()] += -ouApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j+1)
					}
				}

				// Add flow through inner surface
				// Note that this term vanishes for the most inner shell due to boundary conditions
				if (cadet_likely(par != 0))
				{
					// Difference between two cell-centers
					const ParamType dr = static_cast<ParamType>(parCenterRadius[par]) - static_cast<ParamType>(parCenterRadius[par - 1]);

					// Molecular diffusion contribution
					const ResidualType gradCp = (y[0] - y[-strideParPoint()]) / dr;
					if (wantRes)
						*res += innerAreaPerVolume * dp * gradCp;

					// Surface diffusion contribution
					if (cadet_unlikely(_hasSurfaceDiffusion))
					{
						for (int i = 0; i < nBound; ++i)
						{
							// See above for explanation of curIdx value
							const int curIdx = strideParLiquid() - comp + offsetBoundComp()[comp] + i;
							const ResidualType gradQ = (y[curIdx] - y[-strideParPoint() + curIdx]) / dr;
							if (wantRes)
								*res += innerAreaPerVolume * static_cast<ParamType>(parSurfDiff[offsetBoundComp()[comp] + i]) * invBetaP * gradQ;
						}

						if (wantJac)
						{
							const double localInvBetaP = static_cast<double>(invBetaP);
							const double inApV = static_cast<double>(innerAreaPerVolume);
							const double ldr = static_cast<double>(dr);

							// Liquid phase
							jacBase[0] += inApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j)
							jacBase[-strideParPoint()] += -inApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j+1)

							// Solid phase
							for (unsigned int i = 0; i < nBound; ++i)
							{
								// See above for explanation of curIdx value
								const int curIdx = strideParLiquid() - comp + offsetBoundComp()[comp] + i;
								jacBase[curIdx] += inApV * localInvBetaP * static_cast<double>(parSurfDiff[offsetBoundComp()[comp] + i]) / ldr; // dres / dq_i^(p,j)
								jacBase[-strideParPoint() + curIdx] += -inApV * localInvBetaP * static_cast<double>(parSurfDiff[offsetBoundComp()[comp] + i]) / ldr; // dres / dq_i^(p,j-1)
							}
						}
					}
					else if (wantJac)
					{
						// No surface diffusion
						// Liquid phase
						const double inApV = static_cast<double>(innerAreaPerVolume);
						const double ldr = static_cast<double>(dr);

						jacBase[0] += inApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j)
						jacBase[-strideParPoint()] += -inApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j-1)
					}
				}

				if (wantRes)
					res++;
			}

			// Solid phase
			if (cadet_unlikely(_hasSurfaceDiffusion && _hasDynamicReactions))
			{
				for (int bnd = 0; bnd < _strideBound; ++bnd, ++res, ++y, ++jacBase)
				{
					// Skip quasi-stationary bound states
					if (_reqBinding[bnd])
						continue;

					// Add flow through outer surface
					// Note that this term vanishes for the most outer shell due to boundary conditions
					if (cadet_likely(par != _nParPoints - 1))
					{
						// Difference between two cell-centers
						const ParamType dr = static_cast<ParamType>(parCenterRadius[par + 1]) - static_cast<ParamType>(parCenterRadius[par]);

						const ResidualType gradQ = (y[strideParPoint()] - y[0]) / dr;

						if (wantRes)
							*res -= outerAreaPerVolume * static_cast<ParamType>(parSurfDiff[bnd]) * gradQ;

						if (wantJac)
						{
							const double ouApV = static_cast<double>(outerAreaPerVolume);
							const double ldr = static_cast<double>(dr);

							jacBase[0] += ouApV * static_cast<double>(parSurfDiff[bnd]) / ldr; // dres / dq_i^(p,j)
							jacBase[strideParPoint()] += -ouApV * static_cast<double>(parSurfDiff[bnd]) / ldr; // dres / dq_i^(p,j-1)
						}
					}

					// Add flow through inner surface
					// Note that this term vanishes for the most inner shell due to boundary conditions
					if (cadet_likely(par != 0))
					{
						// Difference between two cell-centers
						const ParamType dr = static_cast<ParamType>(parCenterRadius[par]) - static_cast<ParamType>(parCenterRadius[par - 1]);

						const ResidualType gradQ = (y[0] - y[-strideParPoint()]) / dr;

						if (wantRes)
							*res += innerAreaPerVolume * static_cast<ParamType>(parSurfDiff[bnd]) * gradQ;

						if (wantJac)
						{
							const double inApV = static_cast<double>(innerAreaPerVolume);
							const double ldr = static_cast<double>(dr);

							jacBase[0] += inApV * static_cast<double>(parSurfDiff[bnd]) / ldr; // dres / dq_i^(p,j)
							jacBase[-strideParPoint()] += -inApV * static_cast<double>(parSurfDiff[bnd]) / ldr; // dres / dq_i^(p,j-1)
						}
					}
				}
			}
			else
			{
				// Advance pointers over solid phase
				if (wantRes)
					res += _strideBound;
				y += _strideBound;
				jacBase += _strideBound;
			}
		}

		/* Film diffusion */
		// note that bulk equation part is treated outside this operator; here we only handle the particle equation film diffusion

		// Discretized film diffusion kf for finite volumes
		ParamType* const kf_FV = _discParFlux.create<ParamType>(_nComp);
		const ParamType epsP = static_cast<ParamType>(_parPorosity);
		const ParamType surfaceToVolumeRatio = _parGeomSurfToVol / static_cast<ParamType>(_parRadius);
		const ParamType outerAreaPerVolume = static_cast<ParamType>(_parOuterSurfAreaPerVolume[_nParPoints - 1]);

		const ParamType jacPF_val = -outerAreaPerVolume / epsP;

		// Discretized film diffusion kf for finite volumes
		if (cadet_likely(_boundaryOrderFV == 2))
		{
			const ParamType absOuterShellHalfRadius = 0.5 * static_cast<ParamType>(_deltaR[_nParPoints - 1]);
			for (int comp = 0; comp < _nComp; ++comp)
				kf_FV[comp] = 1.0 / (absOuterShellHalfRadius / epsP / static_cast<ParamType>(_poreAccessFactor[comp]) / static_cast<ParamType>(parDiff[comp]) + 1.0 / static_cast<ParamType>(filmDiff[comp]));
		}
		else
		{
			for (int comp = 0; comp < _nComp; ++comp)
				kf_FV[comp] = static_cast<ParamType>(filmDiff[comp]);
		}

		// bead boundary condition in outer bead shell equation
		for (int comp = 0; comp < _nComp; ++comp)
		{
			ResidualType flux = kf_FV[comp] * (yBulk[comp * strideBulkComp()] - yPar[(_nParPoints - 1) * strideParPoint() + comp]);

			if (cadet_unlikely(_hasSurfaceDiffusion && _hasReqReactions && (_nParPoints > 1)))
			{
				// Ordering of particle surface diffusion:
				// bnd0comp0, bnd1comp0, bnd0comp1, bnd1comp1, bnd0comp2, bnd1comp2
				active const* const parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound, secIdx);
				active const* const parCenterRadius = _parCenterRadius.data();
				const ParamType absOuterShellHalfRadius = 0.5 * static_cast<ParamType>(_deltaR[_nParPoints - 1]);

				for (int comp = 0; comp < _nComp; ++comp)
					kf_FV[comp] = (1.0 - static_cast<ParamType>(_parPorosity)) / (1.0 + epsP * static_cast<ParamType>(_poreAccessFactor[comp]) * static_cast<ParamType>(parDiff[comp]) / (absOuterShellHalfRadius * static_cast<ParamType>(filmDiff[comp])));

				const ParamType dr = static_cast<ParamType>(parCenterRadius[_nParPoints - 1]) - static_cast<ParamType>(parCenterRadius[_nParPoints - 2]);

				const unsigned int nBound = _nBound[comp];

				for (int i = 0; i < nBound; ++i)
				{
					const int idxBnd = offsetBoundComp()[comp] + i;

					// Skip quasi-stationary bound states
					if (!_reqBinding[idxBnd])
						continue;

					// Evaluate surface diffusion coefficient and apply weighted arithmetic mean
					const int curIdx = strideParLiquid() + idxBnd;

					const ResidualType gradQ = (yPar[(_nParPoints - 1) * strideParPoint() + curIdx] - yPar[(_nParPoints - 2) * strideParPoint() + curIdx]) / dr;
					flux += kf_FV[comp] * static_cast<ParamType>(parSurfDiff[comp]) * gradQ;
				}
			}

			resPar[(_nParPoints - 1) * strideParPoint() + comp] += jacPF_val / static_cast<ParamType>(_poreAccessFactor[comp]) * flux;
		}

		_discParFlux.destroy<ParamType>();

		return true;
	}

	// ==========================================================================================================================================================  //
	// ========================================						FV particle Jacobian							=============================================  //
	// ==========================================================================================================================================================  //
	
	/**
	 * @brief calculates the bandwidth of the banded particle Jaocbian structure and returns an approximate (>=) number of Jaocbian entries
	 */
	int ParticleDiffusionOperatorFV::particleJacobianBandwidth(unsigned int& lowerBandwidth, unsigned int& upperBandwidth) const
	{
		// Base case: No surface diffusion -> Need to reach same element in previous and next cell
		const unsigned int cellSize = _nComp + _strideBound;
		lowerBandwidth = cellSize;
		upperBandwidth = cellSize;

		return _nParPoints * strideParPoint() * (lowerBandwidth + 1 + upperBandwidth);
	}

	/**
	 * @brief calculates the particle dispersion jacobian Pattern, including entries for the dependence of particle entries on bulk entries through film diffusion boundary condition
	 * @detail Does NOT add film diffusion entries for the dependence of bulk conc. on particle conc.
	*/
	void ParticleDiffusionOperatorFV::setParticleJacobianPattern(std::vector<ParticleDiffusionOperatorFV::T>& tripletList, unsigned int offsetPar, unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx)
	{
		/* particle binding and diffusion entries */
		
		unsigned int lowerBandwidth = 0;
		unsigned int upperBandwidth = 0;
		particleJacobianBandwidth(lowerBandwidth, upperBandwidth);

		for (int parEntry = 0; parEntry < strideParBlock(); parEntry++)
		{
			const int localOffset = offsetPar + parEntry;

			for (int band = 0; band < lowerBandwidth + strideParPoint() + upperBandwidth; band++)
			{
				if (offsetPar > localOffset - lowerBandwidth + band || localOffset - lowerBandwidth + band >= offsetPar + strideParBlock()) // stay within particle block
					continue;

				tripletList.push_back(T(localOffset, localOffset - lowerBandwidth + band, 0.0));
			}
		}

		/* film diffusion entries */
		
		for (int comp = 0; comp < _nComp; comp++)
		{
			// cb on cp dependence through source term
			tripletList.push_back(Eigen::Triplet<double>(offsetBulk + comp, offsetPar + strideParBlock() - strideParPoint() + comp, 0.0));

			// cp on cb dependence through BC
			// row: add particle offset to current parType and particle, go to last cell and current cell and add component offset
			// col: add flux offset to current component, jump over previous cells and components
			tripletList.push_back(T(offsetPar + (_nParPoints - 1) * strideParPoint() + comp * strideParComp(),
				offsetBulk + comp,
				0.0));
		}
	}
	
	unsigned int ParticleDiffusionOperatorFV::jacobianNNZperParticle() const
	{
		unsigned int lowerBandwidth = 0;
		unsigned int upperBandwidth = 0;
		int nBndParDiffEntries = particleJacobianBandwidth(lowerBandwidth, upperBandwidth);

		const int fdNNZ = 2 * _nComp * (1 + 2);

		return fdNNZ + upperBandwidth;
	}
	/**
	 * @brief does nothing, the particle diffusion Jacobian is computed in the residual
	 */
	int ParticleDiffusionOperatorFV::calcParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac)
	{
		return 1;
	}

	int ParticleDiffusionOperatorFV::calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool outliersOnly)
	{
		// Ordering of diffusion:
		// sec0type0comp0, sec0type0comp1, sec0type0comp2, sec0type1comp0, sec0type1comp1, sec0type1comp2,
		// sec1type0comp0, sec1type0comp1, sec1type0comp2, sec1type1comp0, sec1type1comp1, sec1type1comp2, ...
		active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _nComp, secIdx);
		active const* const parDiff = getSectionDependentSlice(_parDiffusion, _nComp, secIdx);
		double* const kf_FV = _discParFlux.create<double>(_nComp);
		const double epsP = static_cast<double>(_parPorosity);
		const double outerAreaPerVolume = static_cast<double>(_parOuterSurfAreaPerVolume[_nParPoints - 1]);
		const double jacPF_val = -outerAreaPerVolume / epsP;

		linalg::BandedEigenSparseRowIterator jacCl(globalJac, offsetC);
		linalg::BandedEigenSparseRowIterator jacCp(globalJac, offsetCp + (_nParPoints - 1) * strideParPoint()); // iterator at the outer particle boundary

		for (int blk = 0; blk < nBulkPoints; blk++)
		{
			for (int comp = 0; comp < _nComp; comp++, ++jacCp, ++jacCl) {
				// add Cb on Cb entries (added since these entries are also touched by bulk jacobian)
				// row: already at bulk phase. already at current node and component.
				// col: already at bulk phase. already at current node and component.
				if (!outliersOnly)
					jacCl[0] += static_cast<double>(filmDiff[comp]) * (1.0 - colPorosity) / colPorosity
					* _parGeomSurfToVol / static_cast<double>(_parRadius)
					* static_cast<double>(parTypeVolFrac[_parTypeIdx + blk * nParType]);
				// add Cb on Cp entries
				// row: already at bulk phase. already at current node and component.
				// col: go to current particle phase entry.
				jacCl[jacCp.row() - jacCl.row()] = -static_cast<double>(filmDiff[comp]) * (1.0 - colPorosity) / colPorosity
					* _parGeomSurfToVol / static_cast<double>(_parRadius)
					* static_cast<double>(parTypeVolFrac[_parTypeIdx + blk * nParType]);


				// Discretized film diffusion kf for finite volumes
				if (cadet_likely(_boundaryOrderFV == 2))
				{
					const double absOuterShellHalfRadius = 0.5 * static_cast<double>(_deltaR[_nParPoints - 1]);
					for (int comp = 0; comp < _nComp; ++comp)
						kf_FV[comp] = 1.0 / (absOuterShellHalfRadius / epsP / static_cast<double>(_poreAccessFactor[comp]) / static_cast<double>(parDiff[comp]) + 1.0 / static_cast<double>(filmDiff[comp]));
				}
				else
				{
					for (int comp = 0; comp < _nComp; ++comp)
						kf_FV[comp] = static_cast<double>(filmDiff[comp]);
				}


				// Cp on Cb entry
				jacCp[jacCl.row() - jacCp.row()] = kf_FV[comp] * jacPF_val / static_cast<double>(_poreAccessFactor[comp]);

				// Cp on Cp entry
				if (!outliersOnly)
					jacCp[0] -= kf_FV[comp] * jacPF_val / static_cast<double>(_poreAccessFactor[comp]);

				// Cp on Cs entries
				if (cadet_unlikely(_hasSurfaceDiffusion && _hasReqReactions && (_nParPoints > 1)))
				{
					// Ordering of particle surface diffusion:
					// bnd0comp0, bnd1comp0, bnd0comp1, bnd1comp1, bnd0comp2, bnd1comp2
					active const* const parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound, secIdx);
					active const* const parCenterRadius = _parCenterRadius.data();
					const double absOuterShellHalfRadius = 0.5 * static_cast<double>(_deltaR[_nParPoints - 1]);

					kf_FV[comp] = (1.0 - static_cast<double>(_parPorosity)) / (1.0 + epsP * static_cast<double>(_poreAccessFactor[comp]) * static_cast<double>(parDiff[comp]) / (absOuterShellHalfRadius * static_cast<double>(filmDiff[comp])));

					const double dr = static_cast<double>(parCenterRadius[0]) - static_cast<double>(parCenterRadius[1]);

					for (int i = 0; i < _nBound[comp]; ++i)
					{
						const int idxBnd = offsetBoundComp()[comp] + i;

						// Skip quasi-stationary bound states
						if (!_reqBinding[idxBnd])
							continue;

						const double v = kf_FV[comp] * static_cast<double>(parSurfDiff[idxBnd]) / dr;
						const int curIdx = strideParLiquid() - comp + idxBnd;

						// Cp on Cs entries
						jacCp[curIdx] -= v;
						jacCp[curIdx - strideParPoint()] += v;
					}
				}
			}

			if (blk < nBulkPoints - 1) // go to last cell of next particle block
				jacCp += _strideBound + (_nParPoints - 1) * strideParPoint();
		}

		_discParFlux.destroy<double>();

		return 1;
	}

	bool ParticleDiffusionOperatorFV::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
	{
		return ParticleDiffusionOperatorBase::setSensitiveParameter(sensParams, pId, adDirection, adValue);
	}

}  // namespace parts

}  // namespace model

}  // namespace cadet

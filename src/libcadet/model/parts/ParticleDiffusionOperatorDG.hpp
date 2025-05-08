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

/**
 * @file
 * Defines the particle dispersion transport operator according to the discontinuous Galerkin discretization.
 */

#ifndef LIBCADET_PARTCICLEDIFFUSIONOPERATORDG_HPP_
#define LIBCADET_PARTCICLEDIFFUSIONOPERATORDG_HPP_

#include "cadet/StrongTypes.hpp"
#include "ParamIdUtil.hpp"
#include "AutoDiff.hpp"
#include "Memory.hpp"
#include "SimulationTypes.hpp"
#include "ParamReaderHelper.hpp"
#include "linalg/BandedEigenSparseRowIterator.hpp"
#include "model/parts/DGToolbox.hpp"
#include "model/ParameterMultiplexing.hpp"

#include <unordered_map>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <vector>

using namespace Eigen;

namespace cadet
{

class IParameterProvider;
class IConfigHelper;
struct AdJacobianParams;
struct SimulationTime;
class IModel;

namespace model
{

class IParameterStateDependence;

namespace parts
{

	constexpr double _SurfVolRatioSphere = 3.0; //!< Surface to volume ratio for a spherical particle
	constexpr double _SurfVolRatioCylinder = 2.0; //!< Surface to volume ratio for a cylindrical particle
	constexpr double _SurfVolRatioSlab = 1.0; //!< Surface to volume ratio for a slab-shaped particle

	/**
	 * @brief Particle dispersion transport operator
	 * @details Implements the equation
	 *
	 * @f[\begin{align}
	 \frac{\partial c^{\mathrm{p}}_{i}}{\partial t} + \frac{1 - \varepsilon^{\mathrm{p}}}{\varepsilon^{\mathrm{p}}} \frac{\partial c^{\mathrm{s}}_{i}}{\partial t} &=
	 \frac{1}{r^2} \frac{\partial }{\partial r} \left( r^2 D_{i}^{\mathrm{p}} \frac{\partial c^{\mathrm{p}}_{i}}{\partial r} \right)
	 - \frac{1 - \varepsilon^{\mathrm{p}}}{\varepsilon^{\mathrm{p}}} \frac{1}{r^2} \frac{\partial }{\partial r} \left( r^2 D_{i}^{\mathrm{s}} \frac{\partial c^{\mathrm{s}}_{i}}{\partial r} \right) \\
	 \end{align} @f]
	 * with Danckwerts boundary conditions (see @cite Danckwerts1953)
	 @f[ \begin{align}
	 - \left. \left( \varepsilon^{\mathrm{p}} D^{\mathrm{p}}_{i} \frac{\partial c^{\mathrm{p}}_{i}}{\partial r} + (1 - \varepsilon^{\mathrm{p}}) D^{\mathrm{s}}_{i} \frac{\partial c^{\mathrm{s}}_{i}}{\partial r} \right) \right|_{r=0}
	 &= 0, \\
	 \varepsilon^{\mathrm{p}} \left. \left( \varepsilon^{\mathrm{p}}  D^{\mathrm{p}}_{i} \frac{\partial c^{\mathrm{p}}_{i}}{\partial r} + (1 - \varepsilon^{\mathrm{p}} ) D^{\mathrm{s}}_{i} \frac{\partial c^{\mathrm{s}}_{i}}{\partial r} \right)\right|_{r = R^{\mathrm{p}}_{}}
	 &= k^{\mathrm{f}}_{i} \left. \left( c^{\mathrm{b}}_i - c^{\mathrm{p}}_{i} \right|_{r = R^{\mathrm{p}}_{}} \right)
	 \end{align} @f]
	 * Additionally implements the variants for cylindrical and slab-shaped particles
	 * Methods are described in @cite Breuer2023
	 *
	 * This class does not store the Jacobian. It only fills existing matrices given to its residual() functions.
	 * It assumes that there is no offset to the particle entries in the local state vector
	 */
	class ParticleDiffusionOperatorDG
	{
	public:

		ParticleDiffusionOperatorDG();
		~ParticleDiffusionOperatorDG() CADET_NOEXCEPT;

		bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int nParType, const int strideBulkComp);
		bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters);

		void setEquidistantRadialDisc(unsigned int parType);
		void setEquivolumeRadialDisc(unsigned int parType);
		void setUserdefinedRadialDisc(unsigned int parType);
		void updateRadialDisc();

		void clearParDepSurfDiffusion();

		bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active const* const filmDiff, active const* const poreAccessFactor);

		int residual(double t, unsigned int secIdx, unsigned int parType, unsigned int colNode, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, int const* const qsBinding, WithoutParamSensitivity);
		int residual(double t, unsigned int secIdx, unsigned int parType, unsigned int colNode, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithParamSensitivity);
		int residual(double t, unsigned int secIdx, unsigned int parType, unsigned int colNode, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithParamSensitivity);
		int residual(double t, unsigned int secIdx, unsigned int parType, unsigned int colNode, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithoutParamSensitivity);


		/* Physical model parameters */

		std::vector<active> _parRadius; //!< Particle radius \f$ r_p \f$
		bool _singleParRadius;
		std::vector<active> _parCoreRadius; //!< Particle core radius \f$ r_c \f$
		bool _singleParCoreRadius;
		std::vector<active> _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$
		bool _singleParPorosity;
		std::vector<double> _parGeomSurfToVol; //!< Particle surface to volume ratio factor (i.e., 3.0 for spherical, 2.0 for cylindrical, 1.0 for hexahedral)
		std::vector<active> _parOuterSurfAreaPerVolume; //!< Particle element outer sphere surface to volume ratio
		std::vector<active> _parInnerSurfAreaPerVolume; //!< Particle element inner sphere surface to volume ratio

		std::vector<active> _filmDiffusion; //!< Particle diffusion coefficient \f$ D_p \f$
		//MultiplexMode _filmDiffusionMode;
		std::vector<active> _poreAccessFactor; //!< Pore accessibility factor \f$ F_{\text{acc}} \f$
		//MultiplexMode _poreAccessFactorMode;
		std::vector<active> _parDiffusion; //!< Particle diffusion coefficient \f$ D_p \f$
		MultiplexMode _parDiffusionMode;
		std::vector<active> _parSurfDiffusion; //!< Particle surface diffusion coefficient \f$ D_s \f$
		MultiplexMode _parSurfDiffusionMode;
		std::vector<IParameterStateDependence*> _parDepSurfDiffusion; //!< Parameter dependencies for particle surface diffusion
		bool _singleParDepSurfDiffusion; //!< Determines whether a single parameter dependence for particle surface diffusion is used
		bool _hasParDepSurfDiffusion; //!< Determines whether particle surface diffusion parameter dependencies are present
		std::vector<bool> _hasSurfaceDiffusion; //!< Determines whether surface diffusion is present in each particle type

		unsigned int _nComp; //!< Number of components
		unsigned int _nParType; //!< Number of particle types

		/* Model discretization */

		enum class ParticleDiscretizationMode : int
		{
			/**
			 * Equidistant distribution of element edges
			 */
			Equidistant,

			/**
			 * Volumes of elements are uniform
			 */
			Equivolume,

			/**
			 * element edges specified by user
			 */
			UserDefined
		};

		// todo when separating discretization and particle model
		//class Indexer
		//{
		//public:

		//	Indexer(const Discretization& disc) : _disc(disc) { }

		//	inline int strideColNode() const CADET_NOEXCEPT { return static_cast<int>(_nComp); }

		//	protected:
		//		_disc;
		//};

		int _strideBulkComp;
		inline int strideBulkComp() const CADET_NOEXCEPT { return _strideBulkComp; }
		inline int strideParComp() const CADET_NOEXCEPT { return 1; }
		inline int strideParLiquid() const CADET_NOEXCEPT { return static_cast<int>(_nComp); }
		inline int strideParBound(int parType) const CADET_NOEXCEPT { return static_cast<int>(_strideBound[parType]); }
		inline int strideParNode(int parType) const CADET_NOEXCEPT { return strideParLiquid() + strideParBound(parType); }
		inline int strideParElem(int parType) const CADET_NOEXCEPT { return strideParNode(parType) * _nParNode[parType]; }
		inline int strideParBlock(int parType) const CADET_NOEXCEPT { return static_cast<int>(_nParPoints[parType]) * strideParNode(parType); }
		inline int offsetBoundComp(ParticleTypeIndex pti, ComponentIndex comp) const CADET_NOEXCEPT { return _boundOffset[pti.value * _nComp + comp.value]; }
		/**
		 * @brief returns the offset between state order and parameter storage (bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2, ...) for one components bound state for a certain particle type
		 * @todo find a different more elegant/easy way?
		 */
		unsigned int getOffsetSurfDiff(unsigned int parType, unsigned int comp, unsigned int bnd) {

			unsigned int offNextBound = 0;

			// we need to estimate the offset to the next parameter of current components next bound state
			// Ordering of particle surface diffusion: bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2
			for (unsigned int _comp = 0; _comp < _nComp; _comp++) {
				if (_comp < comp) // if its a component that occurs before comp, add all bound states of that component up to bnd + 1 (note that bound index starts at 0 -> +1).
					offNextBound += std::min(bnd + 1u, _nBound[parType * _nComp + _comp]);
				else // Otherwise, only add all previous (i.e. up to bnd) bound states of that component. This includes the current component itself (comp == _comp).
					offNextBound += std::min(bnd, _nBound[parType * _nComp + _comp]);
			}

			return offNextBound;
		}
		/**
		 * @brief calculate offsets between surface diffusion parameter storage and state ordering
		 */
		void orderSurfDiff() {

			for (unsigned int type = 0; type < _nParType; type++) {
				for (unsigned int comp = 0; comp < _nComp; comp++) {
					for (unsigned int bnd = 0; bnd < _nBound[type * _nComp + comp]; bnd++) {
						_offsetSurfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd] = getOffsetSurfDiff(type, comp, bnd);
					}
				}
			}
		}

		std::vector<ParticleDiscretizationMode> _parDiscMode; //!< Particle discretization mode

		std::vector<double> _parDiscVector; //!< Particle discretization element edges

		unsigned int* _nParElem; //!< Array with number of radial elements in each particle type
		unsigned int* _nParPointsBeforeType; //!< Array with total number of radial points before a particle type (cumulative sum of nParPoints), additional last element contains total number of particle elements
		unsigned int* _parPolyDeg; //!< polynomial degree of particle elements
		unsigned int* _nParNode; //!< Array with number of radial nodes per element in each particle type
		unsigned int* _nParPoints; //!< Array with number of radial nodes per element in each particle type
		bool* _parExactInt; //!< 1 for exact integration, 0 for inexact LGL quadrature for each particle type
		bool* _parGSM; //!< specifies whether (single element) Galerkin spectral method should be used in particles
		unsigned int* _offsetSurfDiff; //!< particle surface diffusion (may be section and component dependent)
		//unsigned int* _parTypeOffset; //!< Array with offsets (in particle block) to particle type, additional last element contains total number of particle DOFs
		unsigned int* _nBound; //!< Array with number of bound states for each component and particle type (particle type major ordering)
		unsigned int* _boundOffset; //!< Array with offset to the first bound state of each component in the solid phase (particle type major ordering)
		unsigned int* _strideBound; //!< Total number of bound states for each particle type, additional last element contains total number of bound states for all types
		unsigned int* _nBoundBeforeType; //!< Array with number of bound states before a particle type (cumulative sum of strideBound)

		std::vector<active> _parElementSize; //!< Particle element size
		std::vector<active> _parCenterRadius; //!< Particle node-centered position for each particle node

		/* DG specific operators */

		active* _deltaR; //!< equidistant particle element spacing for each particle type
		Eigen::VectorXd* _parNodes; //!< Array with positions of nodes in radial reference element for each particle
		Eigen::MatrixXd* _parPolyDerM; //!< Array with polynomial derivative Matrix for each particle
		Eigen::MatrixXd* _minus_InvMM_ST; //!< equals minus inverse mass matrix times transposed stiffness matrix. Required solely for exact integration DG discretization of particle equation
		Eigen::VectorXd* _parInvWeights; //!< Array with weights for LGL quadrature of size nNodes for each particle
		Eigen::MatrixXd* _parInvMM; //!< dense inverse mass matrix for exact integration of integrals with metrics, for each particle
		Eigen::MatrixXd* _parInvMM_Leg; //!< dense inverse mass matrix (Legendre) for exact integration of integral without metric, for each particle
		Eigen::MatrixXd* _secondOrderStiffnessM; //!< specific second order stiffness matrix
		Eigen::MatrixXd* _minus_parInvMM_Ar; //!< inverse mass matrix times specific second order stiffness matrix
		Eigen::Vector<active, Dynamic>* _Ir; //!< metric part for each particle type and element, particle type major ordering
		Eigen::MatrixXd* _Dr; //!< derivative matrices including metrics for each particle type and element, particle type major ordering
		Eigen::VectorXi _offsetMetric; //!< offset required to access metric dependent DG operator storage of Ir, Dr -> summed up nCells of all previous parTypes

		Eigen::MatrixXd* _DGjacParDispBlocks; //!< particle dispersion blocks of DG jacobian

		Eigen::Vector<active, Dynamic>* _g_p; //!< auxiliary variable g = dc_p / dr
		Eigen::Vector<active, Dynamic>* _g_pSum; //!< auxiliary variable g = sum_{k \in p, s_i} dc_k / dr
		Eigen::Vector<active, Dynamic>* _surfaceFluxParticle; //!< stores the surface flux values for each particle
		active* _localFlux; //!< stores the local (at respective particle) film diffusion flux


		void initializeDG();

		void initializeDGjac(std::vector<double> parGeomSurfToVol);
		Eigen::MatrixXd DGjacobianParDispBlock(unsigned int elemIdx, unsigned int parType, double parGeomSurfToVol);
		Eigen::MatrixXd GSMjacobianParDispBlock(unsigned int parType, double parGeomSurfToVol);
		Eigen::MatrixXd getParBMatrix(int parType, int element, double parGeomSurfToVol);
		Eigen::MatrixXd parAuxBlockGstar(unsigned int elemIdx, unsigned int parType, MatrixXd leftG, MatrixXd middleG, MatrixXd rightG);
		Eigen::MatrixXd getParGBlock(unsigned int elemIdx, unsigned int parType);


	protected:

		template <typename StateType, typename ResidualType, typename ParamType>
		int residualImpl(double t, unsigned int secIdx, unsigned int parType, unsigned int colNode, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, int const* const qsBinding);

		/**
		 * @brief promotes doubles to actives
		 * @detail promotes consecutive doubles to consecutive actives (with zero gradients) based on input double pointer
		 */
		void vectorPromoter(double* state, const unsigned int nVals) {

			const int nDirs = ad::getDirections();
			const int stride = (1 + nDirs);
			const int ADsize = stride * nVals;
			double buff = 0.0;

			for (int val = 1; val <= nVals; val++) // start with last entry to avoid overwriting
			{
				buff = state[nVals - val];
				std::fill(state + ADsize - val * stride, state + ADsize - (val - 1) * stride, 0.0);
				state[ADsize - val * stride] = buff;
			}
		}

		template<typename ResidualType, typename ParamType>
		void applyParInvMap(Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>>& state, unsigned int parType) {
			for (int cell = 0; cell < _nParElem[parType]; cell++) {
				state.segment(cell * _nParNode[parType], _nParNode[parType]) *= 2.0 / static_cast<ParamType>(_deltaR[_offsetMetric[parType] + cell]) * 2.0 / static_cast<ParamType>(_deltaR[_offsetMetric[parType] + cell]);
			}
		}

		template<typename StateType, typename ResidualType>
		void parGSMVolumeIntegral(const int parType, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer) {

			int nNodes = _nParNode[parType];

			stateDer.segment(0, nNodes)
				-= (_minus_parInvMM_Ar[parType].template cast<StateType>() * state.segment(0, nNodes)).template cast<ResidualType>();
		}

		template<typename StateType, typename ResidualType>
		void parVolumeIntegral(const int parType, const bool aux, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer) {

			int nNodes = _nParNode[parType];

			/* no additional metric term for auxiliary equation or particle equation with exact integration scheme
			   -> res = - D * (d_p * c^p + invBeta_p sum_mi d_s c^s) */
			if (aux || (_parExactInt[parType] && _parGeomSurfToVol[parType] == _SurfVolRatioSlab)) {
				// comp-cell-node state vector: use of Eigen lib performance
				for (unsigned int Cell = 0; Cell < _nParElem[parType]; Cell++) {
					stateDer.segment(Cell * nNodes, nNodes)
						-= (_parPolyDerM[parType].template cast<StateType>() * state.segment(Cell * nNodes, nNodes)).template cast<ResidualType>();
				}
			}
			else if (_parExactInt[parType] && _parGeomSurfToVol[parType] != _SurfVolRatioSlab) {
				// comp-cell-node state vector: use of Eigen lib performance
				for (unsigned int Cell = 0; Cell < _nParElem[parType]; Cell++) {
					stateDer.segment(Cell * nNodes, nNodes)
						-= (_minus_InvMM_ST[_offsetMetric[parType] + Cell].template cast<StateType>() * state.segment(Cell * nNodes, nNodes)).template cast<ResidualType>();
				}
			}
			/* include metrics for main particle equation -> res = - D * (d_p * c^p + invBeta_p sum_mi d_s c^s) */
			else { // inexact integration, main equation

				int Cell0 = 0; // auxiliary variable to distinguish special case

				// special case for non slab-shaped particles without core => r(xi_0) = 0
				if (_parGeomSurfToVol[parType] != _SurfVolRatioSlab && _parCoreRadius[parType] == 0.0) {
					Cell0 = 1;

					// compute volume integral except for boundary node
					stateDer.segment(1, nNodes - 1) -= (_Dr[_offsetMetric[parType]].block(1, 1, nNodes - 1, nNodes - 1).template cast<StateType>() * state.segment(1, nNodes - 1)).template cast<ResidualType>();
					// estimate volume integral for boundary node: sum_{j=1}^N state_j * w_j * D_{j,0} * r_j
					stateDer[0] += static_cast<ResidualType>(
						(state.segment(1, nNodes - 1).array()
							* _parInvWeights[parType].segment(1, nNodes - 1).array().cwiseInverse().template cast<StateType>()
							* _parPolyDerM[parType].block(1, 0, nNodes - 1, 1).array().template cast<StateType>()
							* _Ir[_offsetMetric[parType]].segment(1, nNodes - 1).array().template cast<StateType>()
							).sum()
						);
				}

				// "standard" computation for remaining cells
				for (int cell = Cell0; cell < _nParElem[parType]; cell++) {
					stateDer.segment(cell * nNodes, nNodes) -= (_Dr[_offsetMetric[parType] + cell].template cast<StateType>() * state.segment(cell * nNodes, nNodes)).template cast<ResidualType>();
				}
			}
		}
		/*
		 * @brief calculates the interface fluxes g* of particle mass balance equation and implements the respective boundary conditions
		 * @param [in] aux bool if interface flux for auxiliary equation
		 * @param [in] addParDisc bool if interface flux for additional particle DG-discretized equation
		*/
		template<typename StateType>
		void InterfaceFluxParticle(int parType, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state,
			const unsigned int strideCell, const unsigned int strideNode, const bool aux, const int comp, const bool addParDisc = false) {

			Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> _surfFluxPar(reinterpret_cast<StateType*>(&_surfaceFluxParticle[parType][0]), _nParElem[parType] + 1, InnerStride<>(1));

			// reset surface flux storage as it is used multiple times
			_surfFluxPar.setZero();

			// numerical flux: state* = 0.5 (state^+ + state^-)

			// calculate inner interface fluxes
			for (unsigned int Cell = 1u; Cell < _nParElem[parType]; Cell++) {
				_surfFluxPar[Cell] // left interfaces
					= 0.5 * (state[Cell * strideCell - strideNode] + // outer/left node
						state[Cell * strideCell]); // inner/right node
			}

			// calculate boundary interface fluxes.
			if (aux) { // ghost nodes given by state^- := state^+ for auxiliary equation
				_surfFluxPar[0] = state[0];

				_surfFluxPar[_nParElem[parType]] = state[_nParElem[parType] * strideCell - strideNode];
			}
			else if (addParDisc) {
				_surfFluxPar[0] = 0.0;

				_surfFluxPar[_nParElem[parType]] = 0.0;
			}
			else {

				// film diffusion BC
				_surfFluxPar[_nParElem[parType]] = static_cast<StateType>(_localFlux[comp])
					/ (static_cast<double>(_parPorosity[parType]) * static_cast<double>(_poreAccessFactor[parType * _nComp + comp]))
					* (2.0 / static_cast<double>(_deltaR[_offsetMetric[parType]])); // inverse squared mapping was also applied, so we apply Map * invMap^2 = invMap

				// inner particle BC
				_surfFluxPar[0] = 0.0;

			}
		}
		/**
		 * @brief calculates the particle surface Integral (type- and component-wise)
		 * @param [in] parType current particle type
		 * @param [in] state relevant state vector
		 * @param [in] stateDer state derivative vector the solution is added to
		 * @param [in] aux true for auxiliary equation, false for main equation
		 * @param [in] strideCell component-wise cell stride
		 * @param [in] strideNodecomponent-wise node stride
		 * @param [in] comp current component
		*/
		template<typename StateType, typename ResidualType>
		void parSurfaceIntegral(int parType, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state,
			Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer, unsigned const int strideCell, unsigned const int strideNode,
			const bool aux, const int comp = 0, const bool addParDisc = false) {

			Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> _surfFluxPar(reinterpret_cast<StateType*>(&_surfaceFluxParticle[parType][0]), _nParElem[parType] + 1, InnerStride<>(1));

			// calc numerical flux values
			InterfaceFluxParticle<StateType>(parType, state, strideCell, strideNode, aux, comp, addParDisc);

			// strong surface integral -> M^-1 B [state - state*]
			if (!_parExactInt[parType]) { // inexact integration approach -> diagonal mass matrix
				int Cell0 = 0; // auxiliary variable to distinguish special case
				// special case for sphere and cylinder if particle core = 0.0 -> leave out inner particle boundary flux
				if (_parGeomSurfToVol[parType] != _SurfVolRatioSlab && _parCoreRadius[parType] == 0.0) {

					Cell0 = 1;

					stateDer[_parPolyDeg[parType] * strideNode] // last cell node
						+= _parInvWeights[parType][_parPolyDeg[parType]] * (state[_parPolyDeg[parType] * strideNode] - _surfFluxPar[1]);
				}

				for (unsigned int Cell = Cell0; Cell < _nParElem[parType]; Cell++) {

					stateDer[Cell * strideCell] // first cell node
						-= _parInvWeights[parType][0] * (state[Cell * strideCell] - _surfFluxPar[Cell]);

					stateDer[Cell * strideCell + _parPolyDeg[parType] * strideNode] // last cell node
						+= _parInvWeights[parType][_parPolyDeg[parType]] * (state[Cell * strideCell + _parPolyDeg[parType] * strideNode] - _surfFluxPar[Cell + 1u]);
				}
			}
			else { // exact integration approach -> dense mass matrix
				for (unsigned int Cell = 0; Cell < _nParElem[parType]; Cell++) {

					for (unsigned int Node = 0; Node < _nParNode[parType]; Node++) {
						if (aux) { // strong surface integral -> M^-1 B [state - state*] and B has two non-zero entries, -1 and 1
							stateDer[Cell * strideCell + Node * strideNode]
								-= _parInvMM_Leg[parType](Node, 0) * (state[Cell * strideCell] - _surfFluxPar[Cell])
								- _parInvMM_Leg[parType](Node, _parPolyDeg[parType]) * (state[Cell * strideCell + _parPolyDeg[parType] * strideNode] - _surfFluxPar[Cell + 1u]);
						}
						else {
							if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab) { // strong surface integral -> M^-1 B [state - state*] and B has two non-zero entries, -1 and 1
								stateDer[Cell * strideCell + Node * strideNode]
									-= static_cast<ResidualType>(
										_parInvMM[parType](Node, 0) * (state[Cell * strideCell] - _surfFluxPar[Cell])
										- _parInvMM[parType](Node, _parPolyDeg[parType]) * (state[Cell * strideCell + _parPolyDeg[parType] * strideNode] - _surfFluxPar[Cell + 1u])
										);
							}
							else if (_parGeomSurfToVol[parType] == _SurfVolRatioCylinder) { // weak surface integral -> M^-1 B [- state*] and B has two non-zero entries, which depend on metrics
								stateDer[Cell * strideCell + Node * strideNode]
									-= static_cast<ResidualType>(
										_Ir[_offsetMetric[parType] + Cell][0] * _parInvMM[_offsetMetric[parType] + Cell](Node, 0) * (-_surfFluxPar[Cell])
										+ _Ir[_offsetMetric[parType] + Cell][_nParNode[parType] - 1] * _parInvMM[_offsetMetric[parType] + Cell](Node, _parPolyDeg[parType]) * _surfFluxPar[Cell + 1u]
										);
							}
							else if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere) { // weak surface integral -> M^-1 B [- state*] and B has two non-zero entries, which depend on metrics
								stateDer[Cell * strideCell + Node * strideNode]
									-= static_cast<ResidualType>(
										_Ir[_offsetMetric[parType] + Cell][0] * _parInvMM[_offsetMetric[parType] + Cell](Node, 0) * (-_surfFluxPar[Cell])
										+ _Ir[_offsetMetric[parType] + Cell][_nParNode[parType] - 1] * _parInvMM[_offsetMetric[parType] + Cell](Node, _parPolyDeg[parType]) * _surfFluxPar[Cell + 1u]
										);
							}
						}
					}
				}
			}
		}
		/**
		 * @brief solves the auxiliary system g = d c / d xi
		 * @detail computes g = Dc - M^-1 B [c - c^*] and stores this in _g_p
		*/
		template<typename StateType>
		void solve_auxiliary_DG(int parType, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<>>& conc, unsigned int strideCell, unsigned int strideNode, int comp) {

			Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> g_p(reinterpret_cast<StateType*>(&_g_p[parType][0]), _nParPoints[parType], InnerStride<>(1));
			Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> _surfFluxPar(reinterpret_cast<StateType*>(&_surfaceFluxParticle[parType][0]), _nParElem[parType] + 1, InnerStride<>(1));
			_surfFluxPar.setZero(); // reset surface flux storage as it is used multiple times
			g_p.setZero(); // reset auxiliary variable g

			// ========================================================================================//
			// solve auxiliary systems g = d c / d xi	 =>		g_p = Dc - M^-1 B [c - c^*]			   //
			// ========================================================================================//

			parVolumeIntegral<StateType, StateType>(parType, true, conc, g_p); // volumne integral in strong DG form: - D c

			parSurfaceIntegral<StateType>(parType, conc, g_p, strideCell, strideNode, true, comp); // surface integral in strong DG form: M^-1 B [c - c^*]

			g_p *= -1.0; // auxiliary factor -1
		}

		// ==========================================================================================================================================================  //
		// ========================================						DG particle Jacobian							=============================================  //
		// ==========================================================================================================================================================  //

		//typedef Eigen::Triplet<double> T;
		///**
		// * @brief calculates the particle dispersion jacobian Pattern of the exact/inexact integration DG scheme for the given particle type and bead
		//*/
		//void calcParticleJacobianPattern(std::vector<T>& tripletList, unsigned int parType, unsigned int colNode, unsigned int secIdx) {

		//	// Ordering of particle surface diffusion:
		//	// bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2
		//	active const* const _parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound[_nParType], secIdx) + _nBoundBeforeType[parType];

		//	Indexer idxr(_disc);

		//	// (global) strides
		//	unsigned int sCell = _nParNode[parType] * strideParNode(parType);
		//	unsigned int sNode = strideParNode(parType);
		//	unsigned int sComp = 1u;
		//	unsigned int offset = offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode });
		//	unsigned int nNodes = _nParNode[parType];

		//	// case: one cell  -> diffBlock \in R^(nParNodes x nParNodes), GBlock = parPolyDerM
		//	if (_nParElem[parType] == 1) {

		//		// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
		//		for (unsigned int comp = 0; comp < _nComp; comp++) {
		//			for (unsigned int i = 0; i < nNodes; i++) {
		//				for (unsigned int j = 0; j < nNodes; j++) {
		//					// handle liquid state
		//					// row: add component offset and go node strides from there for each dispersion block entry
		//					// col: add component offset and go node strides from there for each dispersion block entry
		//					tripletList.push_back(T(offset + comp * sComp + i * sNode,
		//						offset + comp * sComp + j * sNode, 0.0));

		//					// handle surface diffusion of bound states.
		//					if (_hasSurfaceDiffusion[parType]) {

		//						int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//						for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//							if (_parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) {
		//								// row: add current component offset and go node strides from there for each dispersion block entry
		//								// col: jump oover liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//								tripletList.push_back(T(offset + comp * sComp + i * sNode,
		//									offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode, 0.0));

		//								/* add surface diffusion dispersion block to solid */
		//								if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//									// row: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//									// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//									tripletList.push_back(T(offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//										offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode, 0.0));

		//								}
		//							}
		//						}
		//					}
		//				}
		//			}
		//		}
		//	}
		//	else {

		//		if (!_parExactInt[parType]) {

		//			/*			 left boundary cell				*/

		//			// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
		//			for (unsigned int comp = 0; comp < _nComp; comp++) {
		//				for (unsigned int i = 0; i < nNodes; i++) {
		//					for (unsigned int j = nNodes; j < 3 * nNodes; j++) {
		//						// pattern is more sparse than a nNodes x 2*nNodes block.
		//						if ((j >= nNodes - 1 && j <= 2 * nNodes) ||
		//							(i == 0 && j <= 2 * nNodes) ||
		//							(i == nNodes - 1 && j >= nNodes - 1)) {
		//							// handle liquid state
		//							// row: add component offset and go node strides from there for each dispersion block entry
		//							// col: add component offset and go node strides from there for each dispersion block entry
		//							tripletList.push_back(T(offset + comp * sComp + i * sNode,
		//								offset + comp * sComp + (j - nNodes) * sNode,
		//								0.0));

		//							// handle surface diffusion of bound states. binding is handled in residualKernel().
		//							if (_hasSurfaceDiffusion[parType]) {

		//								int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//								for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//									if (_parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) {
		//										// row: add current component offset and go node strides from there for each dispersion block entry
		//										// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//										tripletList.push_back(T(offset + comp * sComp + i * sNode,
		//											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + (j - nNodes) * sNode,
		//											0.0));

		//										/* add surface diffusion dispersion block to solid */
		//										if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//											// row: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//											// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//											tripletList.push_back(T(offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//												offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + (j - nNodes) * sNode,
		//												0.0));

		//										}
		//									}
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}

		//			/*			 right boundary cell				*/

		//			// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
		//			for (unsigned int comp = 0; comp < _nComp; comp++) {
		//				for (unsigned int i = 0; i < nNodes; i++) {
		//					for (unsigned int j = 0; j < 2 * nNodes; j++) {
		//						// pattern is more sparse than a nNodes x 2*nNodes block.
		//						if ((j >= nNodes - 1 && j <= 2 * nNodes) ||
		//							(i == 0 && j <= 2 * nNodes) ||
		//							(i == nNodes - 1 && j >= nNodes - 1)) {
		//							// handle liquid state
		//							// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//							// col: add component offset and jump over previous cells. Go back one cell and go node strides from there for each dispersion block entry
		//							tripletList.push_back(T(offset + comp * sComp + (_nParElem[parType] - 1) * sCell + i * sNode,
		//								offset + comp * sComp + (_nParElem[parType] - 2) * sCell + j * sNode,
		//								0.0));
		//							// handle surface diffusion of bound states. binding is handled in residualKernel().
		//							if (_hasSurfaceDiffusion[parType]) {

		//								int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//								for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//									if (_parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) {
		//										// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//										// col: jump over liquid states, add current bound state offset and jump over previous cells. Go back one cell and go node strides from there for each dispersion block entry
		//										tripletList.push_back(T(offset + comp * sComp + (_nParElem[parType] - 1) * sCell + i * sNode,
		//											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + (_nParElem[parType] - 2) * sCell + j * sNode,
		//											0.0));

		//										/* add surface diffusion dispersion block to solid */
		//										if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//											// row: jump over previous cells and liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//											// col: jump over previous cells and liquid states, go back one cell, add current bound state offset and go node strides from there for each dispersion block entry
		//											tripletList.push_back(T(offset + (_nParElem[parType] - 1) * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//												offset + (_nParElem[parType] - 2) * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode,
		//												0.0));

		//										}
		//									}
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}

		//			/*				inner cells				*/

		//			for (int cell = 1; cell < _nParElem[parType] - 1; cell++) {

		//				// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
		//				for (unsigned int comp = 0; comp < _nComp; comp++) {
		//					for (unsigned int i = 0; i < nNodes; i++) {
		//						for (unsigned int j = 0; j < 3 * nNodes; j++) {
		//							// pattern is more sparse than a nNodes x 3*nNodes block.
		//							if ((j >= nNodes - 1 && j <= 2 * nNodes) ||
		//								(i == 0 && j <= 2 * nNodes) ||
		//								(i == nNodes - 1 && j >= nNodes - 1)) {
		//								// handle liquid state
		//								// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//								// col: add component offset and jump over previous cells. Go back one cell and go node strides from there for each dispersion block entry
		//								tripletList.push_back(T(offset + comp * sComp + cell * sCell + i * sNode,
		//									offset + comp * sComp + (cell - 1) * sCell + j * sNode, 0.0));
		//								// handle surface diffusion of bound states. binding is handled in residualKernel().
		//								if (_hasSurfaceDiffusion[parType]) {

		//									int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//									for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//										if (_parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) {
		//											// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//											// col: jump over liquid states, add current bound state offset and jump over previous cells. Go back one cell and go node strides from there for each dispersion block entry
		//											tripletList.push_back(T(offset + comp * sComp + cell * sCell + i * sNode,
		//												offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + (cell - 1) * sCell + j * sNode,
		//												0.0));

		//											/* add surface diffusion dispersion block to solid */
		//											if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//												// row: jump over previous cells and liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//												// col: jump over previous cells and liquid states, go back one cell, add current bound state offset and go node strides from there for each dispersion block entry
		//												tripletList.push_back(T(offset + cell * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//													offset + (cell - 1) * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode,
		//													0.0));

		//											}
		//										}
		//									}
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}
		//		}
		//		else { //exact integration

		//			/*			boundary cells			*/

		//			/*			 left boundary cell				*/

		//			unsigned int special = 0u; if (_nParElem[parType] < 3u) special = 1u; // limits the iterator for special case nCells = 3 (dependence on additional entry)
		//			// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
		//			for (unsigned int comp = 0; comp < _nComp; comp++) {
		//				for (unsigned int i = 0; i < nNodes; i++) {
		//					for (unsigned int j = nNodes + 1; j < 3 * nNodes + 2 - special; j++) {
		//						// handle liquid state
		//						// row: add component offset and go node strides from there for each dispersion block entry
		//						// col: add component offset and go node strides from there for each dispersion block entry. adjust for j start
		//						tripletList.push_back(T(offset + comp * sComp + i * sNode,
		//							offset + comp * sComp + j * sNode - (nNodes + 1) * sNode,
		//							0.0));

		//						// handle surface diffusion of bound states. binding is handled in residualKernel().
		//						if (_hasSurfaceDiffusion[parType]) {

		//							int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//							for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//								if (_parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) {
		//									// row: add current component offset and go node strides from there for each dispersion block entry
		//									// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry. adjust for j start
		//									tripletList.push_back(T(offset + comp * sComp + i * sNode,
		//										offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - (nNodes + 1) * sNode,
		//										0.0));

		//									/* add surface diffusion dispersion block to solid */
		//									if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//										// row: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//										// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry. adjust for j start
		//										tripletList.push_back(T(offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - (nNodes + 1) * sNode,
		//											0.0));

		//									}
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}

		//			/*			 right boundary cell				*/

		//			// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
		//			for (unsigned int comp = 0; comp < _nComp; comp++) {
		//				for (unsigned int i = 0; i < nNodes; i++) {

		//					for (unsigned int j = special; j < 2 * nNodes + 1; j++) {
		//						// handle liquid state
		//						// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//						// col: add component offset and jump over previous cells. Go back one cell (and node or adjust for start) and go node strides from there for each dispersion block entry.
		//						tripletList.push_back(T(offset + comp * sComp + (_nParElem[parType] - 1) * sCell + i * sNode,
		//							offset + comp * sComp + (_nParElem[parType] - 1) * sCell - sCell - sNode + j * sNode,
		//							0.0));

		//						// handle surface diffusion of bound states. binding is handled in residualKernel().
		//						if (_hasSurfaceDiffusion[parType]) {

		//							int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//							for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//								if (_parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) {
		//									// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//									// col: jump over liquid states, add current bound state offset and jump over previous cells. Go back one cell (and node or adjust for start) and go node strides from there for each dispersion block entry.
		//									tripletList.push_back(T(offset + comp * sComp + (_nParElem[parType] - 1) * sCell + i * sNode,
		//										offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + (_nParElem[parType] - 2) * sCell - sNode + j * sNode,
		//										0.0));

		//									/* add surface diffusion dispersion block to solid */
		//									if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//										// row: jump over previous cells and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
		//										// col: jump over previous cells and over liquid states, add current bound state offset. go back one cell (and node or adjust for start) and go node strides from there for each dispersion block entry.
		//										tripletList.push_back(T(offset + (_nParElem[parType] - 1) * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//											offset + (_nParElem[parType] - 2) * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) - sNode + bnd + j * sNode,
		//											0.0));
		//									}
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}
		//			if (_nParElem[parType] == 3) {
		//				for (unsigned int comp = 0; comp < _nComp; comp++) {
		//					for (unsigned int i = 0; i < nNodes; i++) {
		//						for (unsigned int j = 1; j < 3 * nNodes + 2 - 1; j++) {
		//							// handle liquid state
		//							// row: add component offset and jump over previous cell. Go node strides from there for each dispersion block entry
		//							// col: add component offset. Go node strides from there for each dispersion block entry. adjust for j start
		//							tripletList.push_back(T(offset + comp * sComp + sCell + i * sNode,
		//								offset + comp * sComp + j * sNode - sNode,
		//								0.0));

		//							// handle surface diffusion of bound states. binding is handled in residualKernel().
		//							if (_hasSurfaceDiffusion[parType]) {

		//								int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//								for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//									if (_parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) {
		//										// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//										// col: jump over liquid states, add current bound state offset and jump over previous cell. go back one cell and go node strides from there for each dispersion block entry. adjust for j start
		//										tripletList.push_back(T(offset + comp * sComp + sCell + i * sNode,
		//											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - sNode,
		//											0.0));

		//										/* add surface diffusion dispersion block to solid */
		//										if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//											// row: jump over previous cells and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
		//											// col: jump over liquid states, add current bound state offset and jump over previous cell. go node strides from there for each dispersion block entry. adjust for j start
		//											tripletList.push_back(T(offset + sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//												offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - sNode,
		//												0.0));

		//										}
		//									}
		//								}
		//							}

		//						}
		//					}
		//				}
		//			}// special case nCells == 3
		//			/*	boundary cell neighbours (exist only if nCells >= 4)	*/
		//			if (_nParElem[parType] >= 4) {

		//				for (unsigned int comp = 0; comp < _nComp; comp++) {
		//					for (unsigned int i = 0; i < nNodes; i++) {
		//						for (unsigned int j = 1; j < 3 * nNodes + 2; j++) {
		//							// handle liquid state
		//							// row: add component offset and jump over previous cell. Go node strides from there for each dispersion block entry
		//							// col: add component offset. Go node strides from there for each dispersion block entry. adjust for j start
		//							tripletList.push_back(T(offset + comp * sComp + sCell + i * sNode,
		//								offset + comp * sComp + j * sNode - sNode,
		//								0.0));

		//							// handle surface diffusion of bound states. binding is handled in residualKernel().
		//							if (_hasSurfaceDiffusion[parType]) {

		//								int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//								for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//									if (_parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) {
		//										// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//										// col: jump over liquid states, add current bound state offset and jump over previous cells. Go back one cell and go node strides from there for each dispersion block entry. adjust for j start
		//										tripletList.push_back(T(offset + comp * sComp + sCell + i * sNode,
		//											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - sNode,
		//											0.0));

		//										/* add surface diffusion dispersion block to solid */
		//										if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//											// row: jump over previous cells and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
		//											// col: jump over previous cells and over liquid states, add current bound state offset. go back one cell and go node strides from there for each dispersion block entry. adjust for j start
		//											tripletList.push_back(T(offset + sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//												offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - sNode,
		//												0.0));

		//										}
		//									}
		//								}
		//							}

		//						}
		//					}
		//				}

		//				for (unsigned int comp = 0; comp < _nComp; comp++) {
		//					for (unsigned int i = 0; i < nNodes; i++) {
		//						for (unsigned int j = 0; j < 3 * nNodes + 2 - 1; j++) {
		//							// handle liquid state
		//							// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//							// col: add component offset and jump over previous cells. Go back one cell and node. Go node strides from there for each dispersion block entry.
		//							tripletList.push_back(T(offset + comp * sComp + (_nParElem[parType] - 2) * sCell + i * sNode,
		//								offset + comp * sComp + (_nParElem[parType] - 2) * sCell - sCell - sNode + j * sNode,
		//								0.0));

		//							// handle surface diffusion of bound states. binding is handled in residualKernel().
		//							if (_hasSurfaceDiffusion[parType]) {

		//								int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//								for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//									if (_parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) {
		//										// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//										// col: jump over liquid states, add current bound state offset and jump over previous cells. Go back one cell and node and go node strides from there for each dispersion block entry
		//										tripletList.push_back(T(offset + comp * sComp + (_nParElem[parType] - 2) * sCell + i * sNode,
		//											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + (_nParElem[parType] - 2) * sCell - sCell - sNode + j * sNode,
		//											0.0));

		//										/* add surface diffusion dispersion block to solid */
		//										if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//											// row: jump over previous cells and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
		//											// col: jump over previous cells and over liquid states, add current bound state offset. go back one cell and node and go node strides from there for each dispersion block entry
		//											tripletList.push_back(T(offset + (_nParElem[parType] - 2) * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//												offset + (_nParElem[parType] - 2) * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) - sCell - sNode + bnd + j * sNode,
		//												0.0));

		//										}
		//									}
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}

		//			/* Inner cells (exist only if nCells >= 5) */

		//			if (_nParElem[parType] >= 5) {

		//				for (unsigned int cell = 2; cell < _nParElem[parType] - 2; cell++) {

		//					for (unsigned int comp = 0; comp < _nComp; comp++) {
		//						for (unsigned int i = 0; i < nNodes; i++) {
		//							for (unsigned int j = 0; j < 3 * nNodes + 2; j++) {
		//								// handle liquid state
		//								// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//								// col: add component offset and jump over previous cells. Go back one cell and node. Go node strides from there for each dispersion block entry.
		//								tripletList.push_back(T(offset + comp * sComp + cell * sCell + i * sNode,
		//									offset + comp * sComp + cell * sCell - sCell - sNode + j * sNode,
		//									0.0));

		//								// handle surface diffusion of bound states. binding is handled in residualKernel().
		//								if (_hasSurfaceDiffusion[parType]) {

		//									int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//									for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//										if (_parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) {
		//											// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//											// col: jump over liquid states, add current bound state offset and jump over previous cells. Go back one cell and node and go node strides from there for each dispersion block entry
		//											tripletList.push_back(T(offset + comp * sComp + cell * sCell + i * sNode,
		//												offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + cell * sCell - sCell - sNode + j * sNode,
		//												0.0));

		//											/* add surface diffusion dispersion block to solid */
		//											if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//												// row: jump over previous cells and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
		//												// col: jump over previous cells and over liquid states, add current bound state offset. go back one cell and node and go node strides from there for each dispersion block entry
		//												tripletList.push_back(T(offset + cell * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//													offset + cell * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd - sCell - sNode + j * sNode,
		//													0.0));

		//											}
		//										}
		//									}
		//								}
		//							}
		//						}
		//					}
		//				}

		//			}

		//		} // parExactInt
		//	} // if nCells > 1
		//}
		//unsigned int calcParDispNNZ(int parType) {

		//	if (_parExactInt[parType]) {
		//		return _nComp * ((3u * _nParElem[parType] - 2u) * _nParNode[parType] * _nParNode[parType] + (2u * _nParElem[parType] - 3u) * _nParNode[parType]);
		//	}
		//	else {
		//		return _nComp * (_nParElem[parType] * _nParNode[parType] * _nParNode[parType] + 8u * _nParNode[parType]);
		//	}
		//}
		///**
		// * @brief sets the sparsity pattern of the binding Jacobian
		// */
		//void parBindingPattern_GRM(std::vector<T>& tripletList, unsigned int parType, unsigned int colNode) {

		//	int offset = offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode });

		//	// every bound state might depend on every bound and liquid state
		//	for (int parNode = 0; parNode < _nParPoints[parType]; parNode++) {
		//		for (int bnd = 0; bnd < _strideBound[parType]; bnd++) {
		//			for (int conc = 0; conc < strideParNode(parType); conc++) {
		//				// row: jump over previous nodes and liquid states and add current bound state offset
		//				// col: jump over previous nodes and add current concentration offset (liquid and bound)
		//				tripletList.push_back(T(offset + parNode * strideParNode(parType) + strideParLiquid() + bnd,
		//					offset + parNode * strideParNode(parType) + conc, 0.0));
		//			}
		//		}
		//	}
		//}
		///**
		// *@brief adds the time derivative entries from particle equations
		// *@detail since the main diagonal entries are already set, we actually only set the solid phase time derivative entries for the discretized particle mass balance equations
		// */
		//void parTimeDerJacPattern_GRM(std::vector<T>& tripletList, unsigned int parType, unsigned int colNode, unsigned int secIdx) {

		//	active const* const _parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound[_nParType], secIdx) + _nBoundBeforeType[parType];
		//	unsigned int offset = offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode });

		//	for (unsigned int parNode = 0; parNode < _nParPoints[parType]; parNode++) {

		//		// discretization special case: we get an algebraic equation at inner particle boundary
		//		if (!_parExactInt[parType] && parNode == 0u && _parGeomSurfToVol[parType] != _SurfVolRatioSlab && _parCoreRadius[parType] == 0.0)
		//			continue;

		//		for (unsigned int comp = 0; comp < _nComp; comp++) {

		//			for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//				// row: jump over previous nodes add current component offset
		//				// col: jump over previous nodes, liquid phase and previous bound states
		//				tripletList.push_back(T(offset + parNode * strideParNode(parType) + comp,
		//					offset + parNode * strideParNode(parType) + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd,
		//					0.0));
		//			}
		//		}
		//	}
		//}
		///**
		// * @brief sets the sparsity pattern of the global Jacobian
		// */
		//void setParJacPattern(std::vector<T>& tripletList, unsigned int parType, unsigned int colNode, unsigned int secIdx) {

		//	calcParticleJacobianPattern(tripletList, parType, colNode, secIdx);

		//	parTimeDerJacPattern_GRM(tripletList, parType, colNode, secIdx);

		//	parBindingPattern_GRM(tripletList, parType, colNode);
		//}
		///**
		// * @brief returns the offset between state order and parameter storage (bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2, ...) for one components bound state for a certain particle type
		// * @todo find a different more elegant/easy way?
		// */
		//unsigned int getOffsetSurfDiff(unsigned int parType, unsigned int comp, unsigned int bnd) {

		//	unsigned int offNextBound = 0;

		//	// we need to estimate the offset to the next parameter of current components next bound state
		//	// Ordering of particle surface diffusion: bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2
		//	for (unsigned int _comp = 0; _comp < _nComp; _comp++) {
		//		if (_comp < comp) // if its a component that occurs before comp, add all bound states of that component up to bnd + 1 (note that bound index starts at 0 -> +1).
		//			offNextBound += std::min(bnd + 1u, _nBound[parType * _nComp + _comp]);
		//		else // Otherwise, only add all previous (i.e. up to bnd) bound states of that component. This includes the current component itself (comp == _comp).
		//			offNextBound += std::min(bnd, _nBound[parType * _nComp + _comp]);
		//	}

		//	return offNextBound;
		//}
		///**
		// * @brief calculate offsets between surface diffusion parameter storage and state ordering
		// */
		//void orderSurfDiff() {

		//	for (unsigned int type = 0; type < _nParType; type++) {
		//		for (unsigned int comp = 0; comp < _nComp; comp++) {
		//			for (unsigned int bnd = 0; bnd < _nBound[type * _nComp + comp]; bnd++) {
		//				_offsetSurfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd] = getOffsetSurfDiff(type, comp, bnd);
		//			}
		//		}
		//	}
		//}
		///**
		// * @brief analytically calculates the particle dispersion jacobian of the DGSEM (exact integration) for a single particle type and bead
		// */
		//int calcParticleDGSEMJacobian(unsigned int parType, unsigned int colNode, const active* const parDiff, const active* const parSurfDiff, const double* const invBetaP, linalg::BandedEigenSparseRowIterator jac) {

		//	// (global) strides
		//	unsigned int sCell = _nParNode[parType] * strideParNode(parType);
		//	unsigned int sNode = strideParNode(parType);
		//	unsigned int sComp = 1u;
		//	unsigned int offset = offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode });
		//	unsigned int nNodes = _nParNode[parType];

		//	/* Special case */
		//	if (_nParElem[parType] == 1) {

		//		linalg::BandedEigenSparseRowIterator jacIt(_globalJac, offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode })); // row iterator starting at first cell, first component

		//		insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 0].block(0, nNodes + 1, nNodes, nNodes), jacIt, idxr, parDiff, parSurfDiff, invBetaP, _binding[parType]->reactionQuasiStationarity(), parType, 1u, 0);
		//		return 1;
		//	}

		//	/* Special case */
		//	if (_nParElem[parType] == 2) {

		//		linalg::BandedEigenSparseRowIterator jacIt(_globalJac, offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode })); // row iterator starting at first cell, first component

		//		insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 0].block(0, nNodes + 1, nNodes, 2 * nNodes), jacIt, idxr, parDiff, parSurfDiff, invBetaP, _binding[parType]->reactionQuasiStationarity(), parType, 1u, 0);
		//		// right Bacobian block, iterator is already moved to second cell
		//		insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 1].block(0, 1, nNodes, 2 * nNodes), jacIt, idxr, parDiff, parSurfDiff, invBetaP, _binding[parType]->reactionQuasiStationarity(), parType, 1u, -strideParShell(parType));
		//		return 1;
		//	}

		//	/* Special case */
		//	if (_nParElem[parType] == 3) {

		//		linalg::BandedEigenSparseRowIterator jacIt(_globalJac, offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }) + strideParShell(parType)); // row iterator starting at first cell, first component

		//		insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 1].block(0, 1, nNodes, 3 * nNodes), jacIt, idxr, parDiff, parSurfDiff, invBetaP, _binding[parType]->reactionQuasiStationarity(), parType, 1u, -strideParShell(parType));
		//	}

		//	/* Inner cells (exist only if nCells >= 5) */
		//	if (_nParElem[parType] >= 5) {

		//		linalg::BandedEigenSparseRowIterator jacIt(_globalJac, offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }) + strideParShell(parType) * 2); // row iterator starting at third cell, first component

		//		// insert all (nElem - 4) inner cell blocks
		//		for (unsigned int cell = 2; cell < _nParElem[parType] - 2; cell++)
		//			insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + cell], jacIt, idxr, parDiff, parSurfDiff, invBetaP, _binding[parType]->reactionQuasiStationarity(), parType, 1u, -(strideParShell(parType) + strideParNode(parType)));
		//	}

		//	/*	boundary cell neighbours (exist only if nCells >= 4)	*/
		//	if (_nParElem[parType] >= 4) {

		//		linalg::BandedEigenSparseRowIterator jacIt(_globalJac, offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }) + strideParShell(parType)); // row iterator starting at second cell, first component

		//		insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 1].block(0, 1, nNodes, 3 * nNodes + 1), jacIt, idxr, parDiff, parSurfDiff, invBetaP, _binding[parType]->reactionQuasiStationarity(), parType, 1u, -strideParShell(parType));

		//		jacIt += (_nParElem[parType] - 4) * strideParShell(parType); // move iterator to preultimate cell (already at third cell)
		//		insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + _nParElem[parType] - 2u].block(0, 0, nNodes, 3 * nNodes + 1), jacIt, idxr, parDiff, parSurfDiff, invBetaP, _binding[parType]->reactionQuasiStationarity(), parType, 1u, -(strideParShell(parType) + strideParNode(parType)));
		//	}

		//	/*			boundary cells (exist only if nCells >= 3)			*/
		//	if (_nParElem[parType] >= 3) {

		//		linalg::BandedEigenSparseRowIterator jacIt(_globalJac, offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode })); // row iterator starting at first cell, first component

		//		insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 0].block(0, nNodes + 1, nNodes, 2 * nNodes + 1), jacIt, idxr, parDiff, parSurfDiff, invBetaP, _binding[parType]->reactionQuasiStationarity(), parType, 1u, 0);

		//		jacIt += (_nParElem[parType] - 2) * strideParShell(parType); // move iterator to last cell (already at second cell)
		//		insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + _nParElem[parType] - 1u].block(0, 0, nNodes, 2 * nNodes + 1), jacIt, idxr, parDiff, parSurfDiff, invBetaP, _binding[parType]->reactionQuasiStationarity(), parType, 1u, -(strideParShell(parType) + strideParNode(parType)));
		//	}

		//	return 1;
		//}
		///**
		// * @brief analytically calculates the particle dispersion jacobian of the collocation DGSEM (inexact integration) for one particle type and bead
		// * @note deprecated, not further development
		// */
		//int calcParticleCollocationDGSEMJacobian(unsigned int parType, unsigned int colNode, const active* const parDiff, const active* const parSurfDiff, const double* const invBetaP, linalg::BandedEigenSparseRowIterator jac) {

		//	Indexer idxr(_disc);

		//	// (global) strides
		//	unsigned int sCell = _nParNode[parType] * strideParNode(parType);
		//	unsigned int sNode = strideParNode(parType);
		//	unsigned int sComp = 1u;
		//	unsigned int offset = offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode });
		//	unsigned int nNodes = _nParNode[parType];

		//	// blocks to compute jacobian
		//	Eigen::MatrixXd dispBlock;
		//	Eigen::MatrixXd B = MatrixXd::Zero(nNodes, nNodes);
		//	B(0, 0) = -1.0; B(nNodes - 1, nNodes - 1) = 1.0;

		//	// special case: one cell -> diffBlock \in R^(nParNodes x nParNodes), GBlock = parPolyDerM
		//	if (_nParElem[parType] == 1) {

		//		double invMap = (2.0 / static_cast<double>(_deltaR[_offsetMetric[parType]]));

		//		if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab || _parCoreRadius[parType] != 0.0)
		//			dispBlock = invMap * invMap * (_Dr[parType] - _parInvWeights[parType].asDiagonal() * B) * _parPolyDerM[parType];

		//		else { // special treatment of inner boundary node for spherical and cylindrical particles without particle core

		//			dispBlock = MatrixXd::Zero(nNodes, nNodes);

		//			// reduced system
		//			dispBlock.block(1, 0, nNodes - 1, nNodes)
		//				= (_Dr[parType].block(1, 1, nNodes - 1, nNodes - 1)
		//					- _parInvWeights[parType].segment(1, nNodes - 1).asDiagonal() * B.block(1, 1, nNodes - 1, nNodes - 1))
		//				* _parPolyDerM[parType].block(1, 0, nNodes - 1, nNodes);

		//			// inner boundary node
		//			dispBlock.block(0, 0, 1, nNodes)
		//				= -(_Ir[parType].segment(1, nNodes - 1).template cast<double>().cwiseProduct(
		//					_parInvWeights[parType].segment(1, nNodes - 1).cwiseInverse()).cwiseProduct(
		//						_parPolyDerM[parType].block(1, 0, nNodes - 1, 1))).transpose()
		//				* _parPolyDerM[parType].block(1, 0, nNodes - 1, nNodes);

		//			dispBlock *= invMap * invMap;
		//		}

		//		// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
		//		for (unsigned int comp = 0; comp < _nComp; comp++) {
		//			for (unsigned int i = 0; i < dispBlock.rows(); i++) {
		//				for (unsigned int j = 0; j < dispBlock.cols(); j++) {
		//					// handle liquid state
		//					// row: add component offset and go node strides from there for each dispersion block entry
		//					// col: add component offset and go node strides from there for each dispersion block entry
		//					_globalJac.coeffRef(offset + comp * sComp + i * sNode,
		//						offset + comp * sComp + j * sNode)
		//						= -(static_cast<double>(parDiff[comp])) * dispBlock(i, j); // - D_p * (Delta r / 2)^2 * (D_r D - M^-1 B D)

		//					// handle surface diffusion of bound states. binding is handled in residualKernel().
		//					if (_hasSurfaceDiffusion[parType]) {

		//						int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//						for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//							if (static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) != 0.0) {
		//								/* add surface diffusion dispersion block to liquid */
		//								// row: add current component offset and go node strides from there for each dispersion block entry
		//								// col: jump oover liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//								_globalJac.coeffRef(offset + comp * sComp + i * sNode,
		//									offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode)
		//									= -(static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) * invBetaP[comp]) * dispBlock(i, j); // -  D_s * (1 / Beta_p) * (Delta r / 2)^2 * (D_r D - M^-1 B D)

		//								/* add surface diffusion dispersion block to solid */
		//								if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//									// row: jump oover liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//									// col: jump oover liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//									_globalJac.coeffRef(offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//										offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode)
		//										= -(static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)])) * dispBlock(i, j); // -  D_s * (Delta r / 2)^2 * (D_r D - M^-1 B D)

		//								}
		//							}
		//						}
		//					}
		//				}
		//			}
		//		}
		//	}
		//	else {

		//		/*			boundary cells			*/

		//		// initialize dispersion and metric block matrices
		//		MatrixXd bnd_dispBlock = MatrixXd::Zero(nNodes, 2 * nNodes); // boundary cell specific
		//		dispBlock = MatrixXd::Zero(nNodes, 3 * nNodes);

		//		// compute blocks used for inexact integration scheme
		//		// auxiliary block [ d g(c) / d c ] for left boundary cell
		//		MatrixXd GBlock_l = MatrixXd::Zero(nNodes, nNodes + 1);
		//		GBlock_l.block(0, 0, nNodes, nNodes) = _parPolyDerM[parType];
		//		GBlock_l(nNodes - 1, nNodes - 1) -= 0.5 * _parInvWeights[parType][nNodes - 1];
		//		GBlock_l(nNodes - 1, nNodes) += 0.5 * _parInvWeights[parType][nNodes - 1];
		//		// auxiliary block [ d g(c) / d c ] for right boundary cell
		//		MatrixXd GBlock_r = MatrixXd::Zero(nNodes, nNodes + 1);
		//		GBlock_r.block(0, 1, nNodes, nNodes) = _parPolyDerM[parType];
		//		GBlock_r(0, 0) -= 0.5 * _parInvWeights[parType][0];
		//		GBlock_r(0, 1) += 0.5 * _parInvWeights[parType][0];
		//		// numerical flux contribution for right interface of left boundary cell -> d f^*_N / d cp
		//		MatrixXd bnd_gStarDC = MatrixXd::Zero(nNodes, 2 * nNodes);
		//		bnd_gStarDC.block(nNodes - 1, 0, 1, nNodes + 1) = GBlock_l.block(nNodes - 1, 0, 1, nNodes + 1);
		//		bnd_gStarDC.block(nNodes - 1, nNodes - 1, 1, nNodes + 1) += GBlock_r.block(0, 0, 1, nNodes + 1);
		//		bnd_gStarDC *= 0.5;

		//		/*			 left boundary cell				*/
		//		double invMap = (2.0 / static_cast<double>(_deltaR[_offsetMetric[parType]]));

		//		// "standard" computation for slab-shaped particles and spherical, cylindrical particles without core
		//		if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab || _parCoreRadius[parType] != 0.0) {
		//			// dispBlock <- invMap^2 * ( D * G_l - M^-1 * B * [G_l - g^*] )
		//			bnd_dispBlock.block(0, 0, nNodes, nNodes + 1) = (_Dr[_offsetMetric[parType]] - _parInvWeights[parType].asDiagonal() * B) * GBlock_l;
		//			bnd_dispBlock.block(0, 0, nNodes, 2 * nNodes) += _parInvWeights[parType].asDiagonal() * B * bnd_gStarDC;
		//			bnd_dispBlock *= invMap * invMap;
		//		}
		//		else { // special treatment of inner boundary node for spherical and cylindrical particles without particle core

		//			// inner boundary node
		//			bnd_dispBlock.block(0, 0, 1, nNodes + 1)
		//				= -(_Ir[_offsetMetric[parType]].template cast<double>().segment(1, nNodes - 1).cwiseProduct(
		//					_parInvWeights[parType].segment(1, nNodes - 1).cwiseInverse()).cwiseProduct(
		//						_parPolyDerM[parType].block(1, 0, nNodes - 1, 1))).transpose()
		//				* GBlock_l.block(1, 0, nNodes - 1, nNodes + 1);

		//			// reduced system for remaining nodes
		//			bnd_dispBlock.block(1, 0, nNodes - 1, nNodes + 1)
		//				= (_Dr[_offsetMetric[parType]].block(1, 1, nNodes - 1, nNodes - 1)
		//					- _parInvWeights[parType].segment(1, nNodes - 1).asDiagonal() * B.block(1, 1, nNodes - 1, nNodes - 1)
		//					) * GBlock_l.block(1, 0, nNodes - 1, nNodes + 1);

		//			bnd_dispBlock.block(1, 0, nNodes - 1, 2 * nNodes)
		//				+= _parInvWeights[parType].segment(1, nNodes - 1).asDiagonal() * B.block(1, 1, nNodes - 1, nNodes - 1) * bnd_gStarDC.block(1, 0, nNodes - 1, 2 * nNodes);

		//			// mapping
		//			bnd_dispBlock *= invMap * invMap;
		//		}

		//		// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
		//		for (unsigned int comp = 0; comp < _nComp; comp++) {
		//			for (unsigned int i = 0; i < bnd_dispBlock.rows(); i++) {
		//				for (unsigned int j = 0; j < bnd_dispBlock.cols(); j++) {
		//					// handle liquid state
		//					// inexact integration pattern is more sparse than a nNodes x 2*nNodes block.
		//					if ((j <= nNodes) || (i == nNodes - 1)) {
		//						// row: add component offset and go node strides from there for each dispersion block entry
		//						// col: add component offset and go node strides from there for each dispersion block entry
		//						_globalJac.coeffRef(offset + comp * sComp + i * sNode,
		//							offset + comp * sComp + j * sNode)
		//							= -static_cast<double>(parDiff[comp]) * bnd_dispBlock(i, j); // dispBlock <- D_p * [ M^-1 * M_r * G_l +  invMap * ( D * G_l - M^-1 * B * [G_l - g^*] ) ]

		//						// handle surface diffusion of bound states. binding is handled in residualKernel().
		//						if (_hasSurfaceDiffusion[parType]) {

		//							int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//							for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//								if (static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) != 0.0) {
		//									// row: add current component offset and go node strides from there for each dispersion block entry
		//									// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//									_globalJac.coeffRef(offset + comp * sComp + i * sNode,
		//										offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode)
		//										= -static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) * invBetaP[comp] * bnd_dispBlock(i, j); // dispBlock <- D_s * invBeta * [ M^-1 * M_r * G_l +  invMap * ( D * G_l - M^-1 * B * [G_l - g^*] ) ]

		//									/* add surface diffusion dispersion block to solid */
		//									if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//										// row: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//										// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//										_globalJac.coeffRef(offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode)
		//											= -(static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)])) * bnd_dispBlock(i, j); // -  D_s * (Delta r / 2)^2 * (D_r D - M^-1 B D)

		//									}
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}
		//		}

		//		/*			 right boundary cell				*/
		//		invMap = (2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + _nParElem[parType] - 1]));

		//		// numerical flux contribution for left interface of right boundary cell -> d f^*_0 / d cp
		//		bnd_gStarDC.setZero();
		//		bnd_gStarDC.block(0, nNodes - 1, 1, nNodes + 1) = GBlock_r.block(0, 0, 1, nNodes + 1);
		//		bnd_gStarDC.block(0, 0, 1, nNodes + 1) += GBlock_l.block(nNodes - 1, 0, 1, nNodes + 1);
		//		bnd_gStarDC *= 0.5;
		//		// dispBlock <- invMap * ( D_r * G_r - M^-1 * B * [G_r - g^*] )
		//		bnd_dispBlock.setZero();
		//		bnd_dispBlock.block(0, nNodes - 1, nNodes, nNodes + 1) = (_Dr[_offsetMetric[parType] + _nParElem[parType] - 1] - _parInvWeights[parType].asDiagonal() * B) * GBlock_r;
		//		bnd_dispBlock.block(0, 0, nNodes, 2 * nNodes) += _parInvWeights[parType].asDiagonal() * B * bnd_gStarDC;
		//		bnd_dispBlock *= invMap * invMap;

		//		// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
		//		for (unsigned int comp = 0; comp < _nComp; comp++) {
		//			for (unsigned int i = 0; i < bnd_dispBlock.rows(); i++) {
		//				for (unsigned int j = 0; j < bnd_dispBlock.cols(); j++) {
		//					// handle liquid state
		//					// inexact integration pattern is more sparse than a nNodes x 2*nNodes block.
		//					if ((j <= nNodes) || (i == nNodes - 1)) {
		//						// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//						// col: add component offset and jump over previous cells. Go back one cell and go node strides from there for each dispersion block entry
		//						_globalJac.coeffRef(offset + comp * sComp + (_nParElem[parType] - 1) * sCell + i * sNode,
		//							offset + comp * sComp + (_nParElem[parType] - 2) * sCell + j * sNode)
		//							= -static_cast<double>(parDiff[comp]) * bnd_dispBlock(i, j); // dispBlock <- D_p * [ M^-1 * M_r * G_l +  invMap * ( D * G_l - M^-1 * B * [G_l - g^*] ) ]

		//						// handle surface diffusion of bound states. binding is handled in residualKernel().
		//						if (_hasSurfaceDiffusion[parType]) {

		//							int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//							for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//								if (static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) != 0.0) {
		//									// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//									// col: jump over liquid states, add current bound state offset and jump over previous cells. Go back one cell and go node strides from there for each dispersion block entry
		//									_globalJac.coeffRef(offset + comp * sComp + (_nParElem[parType] - 1) * sCell + i * sNode,
		//										offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + (_nParElem[parType] - 2) * sCell + j * sNode)
		//										= -static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) * invBetaP[comp] * bnd_dispBlock(i, j); // dispBlock <- D_s * invBeta * [ M^-1 * M_r * G_l +  invMap * ( D * G_l - M^-1 * B * [G_l - g^*] ) ]

		//									/* add surface diffusion dispersion block to solid */
		//									if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//										// row: jump over previous cells and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
		//										// col: jump over previous cells and over liquid states, add current bound state offset. go back one cell and go node strides from there for each dispersion block entry
		//										_globalJac.coeffRef(offset + (_nParElem[parType] - 1) * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//											offset + (_nParElem[parType] - 2) * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode)
		//											= -(static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)])) * bnd_dispBlock(i, j); // -  D_s * (Delta r / 2)^2 * (D_r D - M^-1 B D)

		//									}
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}
		//		}

		//		/*				inner cells				*/

		//		// auxiliary block [ d g(c) / d c ] for inner cells
		//		MatrixXd GBlock = MatrixXd::Zero(nNodes, nNodes + 2);
		//		GBlock.block(0, 1, nNodes, nNodes) = _parPolyDerM[parType];
		//		GBlock(0, 0) -= 0.5 * _parInvWeights[parType][0];
		//		GBlock(0, 1) += 0.5 * _parInvWeights[parType][0];
		//		GBlock(nNodes - 1, nNodes) -= 0.5 * _parInvWeights[parType][nNodes - 1];
		//		GBlock(nNodes - 1, nNodes + 1) += 0.5 * _parInvWeights[parType][nNodes - 1];

		//		// numerical flux contribution
		//		MatrixXd gStarDC = MatrixXd::Zero(nNodes, 3 * nNodes);
		//		gStarDC.block(0, nNodes - 1, 1, nNodes + 2) = GBlock.block(0, 0, 1, nNodes + 2);
		//		gStarDC.block(0, 0, 1, nNodes + 1) += GBlock.block(nNodes - 1, 1, 1, nNodes + 1);
		//		gStarDC.block(nNodes - 1, nNodes - 1, 1, nNodes + 2) += GBlock.block(nNodes - 1, 0, 1, nNodes + 2);
		//		gStarDC.block(nNodes - 1, 2 * nNodes - 1, 1, nNodes + 1) += GBlock.block(0, 0, 1, nNodes + 1);
		//		gStarDC *= 0.5;

		//		dispBlock.setZero();
		//		// dispersion block part without metrics
		//		dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) = -1.0 * _parInvWeights[parType].asDiagonal() * B * GBlock;
		//		dispBlock.block(0, 0, nNodes, 3 * nNodes) += _parInvWeights[parType].asDiagonal() * B * gStarDC;
		//		dispBlock *= invMap * invMap;

		//		for (int cell = 1; cell < _nParElem[parType] - 1; cell++) {
		//			double invMap = (2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + cell]));

		//			// add metric part, dependent on current cell
		//			dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) += _Dr[_offsetMetric[parType] + cell] * GBlock * invMap * invMap;

		//			// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
		//			for (unsigned int comp = 0; comp < _nComp; comp++) {
		//				for (unsigned int i = 0; i < dispBlock.rows(); i++) {
		//					for (unsigned int j = 0; j < dispBlock.cols(); j++) {
		//						// handle liquid state
		//						// pattern is more sparse than a nNodes x 3*nNodes block.
		//						if ((j >= nNodes - 1 && j <= 2 * nNodes) ||
		//							(i == 0 && j <= 2 * nNodes) ||
		//							(i == nNodes - 1 && j >= nNodes - 1)) {
		//							// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//							// col: add component offset and jump over previous cells. Go back one cell and go node strides from there for each dispersion block entry
		//							_globalJac.coeffRef(offset + comp * sComp + cell * sCell + i * sNode,
		//								offset + comp * sComp + (cell - 1) * sCell + j * sNode)
		//								= -static_cast<double>(parDiff[comp]) * dispBlock(i, j); // dispBlock <- D_p * [ M^-1 * M_r * G_l +  invMap * ( D * G_l - M^-1 * B * [G_l - g^*] ) ]

		//							// handle surface diffusion of bound states. binding is handled in residualKernel().
		//							if (_hasSurfaceDiffusion[parType]) {

		//								int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();

		//								for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
		//									if (static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) != 0.0) {
		//										// row: add component offset and jump over previous cells. Go node strides from there for each dispersion block entry
		//										// col: jump over liquid states, add current bound state offset and jump over previous cells. Go back one cell and go node strides from there for each dispersion block entry
		//										_globalJac.coeffRef(offset + comp * sComp + cell * sCell + i * sNode,
		//											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + (cell - 1) * sCell + j * sNode)
		//											= -static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) * invBetaP[comp] * dispBlock(i, j); // dispBlock <- D_s * invBeta * [ M^-1 * M_r * G_l +  invMap * ( D * G_l - M^-1 * B * [G_l - g^*] ) ]

		//										/* add surface diffusion dispersion block to solid */
		//										if (!qsReaction[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
		//											// row: jump over previous cells and liquid states, add current bound state offset and go node strides from there for each dispersion block entry
		//											// col: jump over previous cells and liquid states, go back one cell and add current bound state offset and go node strides from there for each dispersion block entry
		//											_globalJac.coeffRef(offset + cell * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
		//												offset + (cell - 1) * sCell + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode)
		//												= -(static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)])) * dispBlock(i, j); // -  D_s * (Delta r / 2)^2 * (D_r D - M^-1 B D)

		//										}
		//									}
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}
		//			// substract metric part in preparation of next iteration
		//			dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) -= _Dr[_offsetMetric[parType] + cell] * GBlock * invMap * invMap;
		//		}

		//	} // if nCells > 1

		//	return 1;
		//}
		///**
		// * @brief adds jacobian entries which have been overwritten by the binding kernel (only use for surface diffusion combined with kinetic binding)
		// * @detail only adds the entries d RHS_i / d c^s_i, which lie on the diagonal
		// * @parType[in] current particle type
		// * @parSurfDiff[in] pointer to particle surface diffusion at current section and particle type
		// */
		//int addSolidDGentries(unsigned int parType, const active* const parSurfDiff) {

		//	if (!_parExactInt[parType])
		//		return addSolidDGentries_inexInt(parType, parSurfDiff);

		//	Indexer idxr(_disc);

		//	for (unsigned int col = 0; col < _nPoints; col++) {
		//		// Get jacobian iterator at first solid entry of first particle of current type
		//		linalg::BandedEigenSparseRowIterator jac(_globalJac, offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ col }) + strideParLiquid());

		//		for (unsigned int cell = 0; cell < _nParElem[parType]; cell++)
		//			addDiagonalSolidJacobianEntries(_DGjacParDispBlocks[_offsetMetric[parType] + cell].block(0, _nParNode[parType] + 1, _nParNode[parType], _nParNode[parType]),
		//				jac, idxr, parSurfDiff, _binding[parType]->reactionQuasiStationarity(), parType);
		//	}

		//	return 1;
		//}
		///**
		// * @brief adds jacobian entries which have been overwritten by the binding kernel (only use for surface diffusion combined with kinetic binding)
		// * @detail only adds the entries d RHS_i / d c^s_i, which lie on the diagonal
		// * @parType[in] current particle type
		// * @parSurfDiff[in] pointer to particle surface diffusion at current section and particle type
		// */
		//int addSolidDGentries_inexInt(unsigned int parType, const active* const parSurfDiff) {

		//	Indexer idxr(_disc);

		//	// (global) strides
		//	unsigned int sCell = _nParNode[parType] * strideParNode(parType);
		//	unsigned int sNode = strideParNode(parType);
		//	unsigned int sComp = 1u;
		//	unsigned int nNodes = _nParNode[parType];

		//	// blocks to compute jacobian
		//	Eigen::MatrixXd dispBlock;
		//	Eigen::MatrixXd B = MatrixXd::Zero(nNodes, nNodes);
		//	B(0, 0) = -1.0; B(nNodes - 1, nNodes - 1) = 1.0;

		//	// special case: one cell -> diffBlock \in R^(nParNodes x nParNodes), GBlock = parPolyDerM
		//	if (_nParElem[parType] == 1) {

		//		double invMap = (2.0 / static_cast<double>(_deltaR[_offsetMetric[parType]]));

		//		if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab || _parCoreRadius[parType] != 0.0)
		//			dispBlock = invMap * invMap * (_Dr[parType] - _parInvWeights[parType].asDiagonal() * B) * _parPolyDerM[parType];

		//		else { // special treatment of inner boundary node for spherical and cylindrical particles without particle core

		//			dispBlock = MatrixXd::Zero(nNodes, nNodes);

		//			// reduced system
		//			dispBlock.block(1, 0, nNodes - 1, nNodes)
		//				= (_Dr[parType].block(1, 1, nNodes - 1, nNodes - 1)
		//					- _parInvWeights[parType].segment(1, nNodes - 1).asDiagonal() * B.block(1, 1, nNodes - 1, nNodes - 1))
		//				* _parPolyDerM[parType].block(1, 0, nNodes - 1, nNodes);

		//			// inner boundary node
		//			dispBlock.block(0, 0, 1, nNodes)
		//				= -(_Ir[parType].segment(1, nNodes - 1).template cast<double>().cwiseProduct(
		//					_parInvWeights[parType].segment(1, nNodes - 1).cwiseInverse()).cwiseProduct(
		//						_parPolyDerM[parType].block(1, 0, nNodes - 1, 1))).transpose()
		//				* _parPolyDerM[parType].block(1, 0, nNodes - 1, nNodes);

		//			dispBlock *= invMap * invMap;
		//		}

		//		for (unsigned int colNode = 0; colNode < _nPoints; colNode++)
		//		{
		//			unsigned int offset = offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode });
		//			// start at first solid entry
		//			linalg::BandedEigenSparseRowIterator jac(_globalJac, offset + strideParLiquid());

		//			for (unsigned int node = 0; node < _nParNode[parType]; node++, jac += strideParLiquid())
		//			{
		//				for (unsigned int comp = 0; comp < _nComp; comp++)
		//				{
		//					for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++, ++jac)
		//					{
		//						if (static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) != 0.0)
		//						{
		//							jac[0] += -(static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)])) * dispBlock(node, node);
		//						}
		//					}
		//				}
		//			}
		//		}
		//	}
		//	else {

		//		/*			boundary cells			*/
		//		// initialize dispersion and metric block matrices
		//		MatrixXd bnd_dispBlock = MatrixXd::Zero(nNodes, 2 * nNodes); // boundary cell specific
		//		dispBlock = MatrixXd::Zero(nNodes, 3 * nNodes);

		//		// auxiliary block [ d g(c) / d c ] for left boundary cell
		//		MatrixXd GBlock_l = MatrixXd::Zero(nNodes, nNodes + 1);
		//		GBlock_l.block(0, 0, nNodes, nNodes) = _parPolyDerM[parType];
		//		GBlock_l(nNodes - 1, nNodes - 1) -= 0.5 * _parInvWeights[parType][nNodes - 1];
		//		GBlock_l(nNodes - 1, nNodes) += 0.5 * _parInvWeights[parType][nNodes - 1];
		//		// auxiliary block [ d g(c) / d c ] for right boundary cell
		//		MatrixXd GBlock_r = MatrixXd::Zero(nNodes, nNodes + 1);
		//		GBlock_r.block(0, 1, nNodes, nNodes) = _parPolyDerM[parType];
		//		GBlock_r(0, 0) -= 0.5 * _parInvWeights[parType][0];
		//		GBlock_r(0, 1) += 0.5 * _parInvWeights[parType][0];

		//		/*			 left boundary cell				*/
		//		int _cell = 0;
		//		double invMap = (2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + _cell]));

		//		// numerical flux contribution for right interface of left boundary cell -> d f^*_N / d cp
		//		MatrixXd bnd_gStarDC = MatrixXd::Zero(nNodes, 2 * nNodes);
		//		bnd_gStarDC.block(nNodes - 1, 0, 1, nNodes + 1) = GBlock_l.block(nNodes - 1, 0, 1, nNodes + 1);
		//		bnd_gStarDC.block(nNodes - 1, nNodes - 1, 1, nNodes + 1) += GBlock_r.block(0, 0, 1, nNodes + 1);
		//		bnd_gStarDC *= 0.5;

		//		// "standard" computation for slab-shaped particles and spherical, cylindrical particles without core
		//		if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab || _parCoreRadius[parType] != 0.0) {
		//			// dispBlock <- invMap^2 * ( D * G_l - M^-1 * B * [G_l - g^*] )
		//			bnd_dispBlock.block(0, 0, nNodes, nNodes + 1) = (_Dr[_offsetMetric[parType]] - _parInvWeights[parType].asDiagonal() * B) * GBlock_l;
		//			bnd_dispBlock.block(0, 0, nNodes, 2 * nNodes) += _parInvWeights[parType].asDiagonal() * B * bnd_gStarDC;
		//			bnd_dispBlock *= invMap * invMap;
		//		}
		//		else { // special treatment of inner boundary node for spherical and cylindrical particles without particle core

		//			// inner boundary node
		//			bnd_dispBlock.block(0, 0, 1, nNodes + 1)
		//				= -(_Ir[_offsetMetric[parType]].template cast<double>().segment(1, nNodes - 1).cwiseProduct(
		//					_parInvWeights[parType].segment(1, nNodes - 1).cwiseInverse()).cwiseProduct(
		//						_parPolyDerM[parType].block(1, 0, nNodes - 1, 1))).transpose()
		//				* GBlock_l.block(1, 0, nNodes - 1, nNodes + 1);

		//			// reduced system for remaining nodes
		//			bnd_dispBlock.block(1, 0, nNodes - 1, nNodes + 1)
		//				= (_Dr[_offsetMetric[parType]].block(1, 1, nNodes - 1, nNodes - 1)
		//					- _parInvWeights[parType].segment(1, nNodes - 1).asDiagonal() * B.block(1, 1, nNodes - 1, nNodes - 1)
		//					) * GBlock_l.block(1, 0, nNodes - 1, nNodes + 1);

		//			bnd_dispBlock.block(1, 0, nNodes - 1, 2 * nNodes)
		//				+= _parInvWeights[parType].segment(1, nNodes - 1).asDiagonal() * B.block(1, 1, nNodes - 1, nNodes - 1) * bnd_gStarDC.block(1, 0, nNodes - 1, 2 * nNodes);

		//			// mapping
		//			bnd_dispBlock *= invMap * invMap;
		//		}

		//		for (unsigned int colNode = 0; colNode < _nPoints; colNode++) {

		//			unsigned int offset = offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode });
		//			// start at first solid entry of first cell
		//			linalg::BandedEigenSparseRowIterator jac_left(_globalJac, offset + strideParLiquid());

		//			for (unsigned int node = 0; node < _nParNode[parType]; node++, jac_left += strideParLiquid()) {
		//				for (unsigned int comp = 0; comp < _nComp; comp++) {
		//					for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++, ++jac_left) {
		//						if (static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) != 0.0) {
		//							jac_left[0] += -(static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)])) * bnd_dispBlock(node, node);
		//						}
		//					}
		//				}
		//			}
		//		}

		//		/*			 right boundary cell				*/
		//		_cell = _nParElem[parType] - 1;
		//		invMap = (2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + _cell]));

		//		bnd_gStarDC = MatrixXd::Zero(nNodes, 2 * nNodes);
		//		// numerical flux contribution for left interface of right boundary cell -> d f^*_0 / d cp
		//		bnd_gStarDC.setZero();
		//		bnd_gStarDC.block(0, nNodes - 1, 1, nNodes + 1) = GBlock_r.block(0, 0, 1, nNodes + 1);
		//		bnd_gStarDC.block(0, 0, 1, nNodes + 1) += GBlock_l.block(nNodes - 1, 0, 1, nNodes + 1);
		//		bnd_gStarDC *= 0.5;
		//		// dispBlock <- invMap * ( D_r * G_r - M^-1 * B * [G_r - g^*] )
		//		bnd_dispBlock.setZero();
		//		bnd_dispBlock.block(0, nNodes - 1, nNodes, nNodes + 1) = (_Dr[_offsetMetric[parType] + _cell] - _parInvWeights[parType].asDiagonal() * B) * GBlock_r;
		//		bnd_dispBlock.block(0, 0, nNodes, 2 * nNodes) += _parInvWeights[parType].asDiagonal() * B * bnd_gStarDC;
		//		bnd_dispBlock *= invMap * invMap;

		//		for (unsigned int colNode = 0; colNode < _nPoints; colNode++) {

		//			unsigned int offset = offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode });
		//			// start at first solid entry of last cell
		//			linalg::BandedEigenSparseRowIterator jac_right(_globalJac, offset + (_nParElem[parType] - 1) * sCell + strideParLiquid());

		//			for (unsigned int node = 0; node < _nParNode[parType]; node++, jac_right += strideParLiquid()) {
		//				for (unsigned int comp = 0; comp < _nComp; comp++) {
		//					for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++, ++jac_right) {
		//						if (static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) != 0.0) {
		//							jac_right[0] += -(static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)])) * bnd_dispBlock(node, _nParNode[parType] + node);
		//						}
		//					}
		//				}
		//			}
		//		}

		//		/*				inner cells				*/

		//			// auxiliary block [ d g(c) / d c ] for inner cells
		//		MatrixXd GBlock = MatrixXd::Zero(nNodes, nNodes + 2);
		//		GBlock.block(0, 1, nNodes, nNodes) = _parPolyDerM[parType];
		//		GBlock(0, 0) -= 0.5 * _parInvWeights[parType][0];
		//		GBlock(0, 1) += 0.5 * _parInvWeights[parType][0];
		//		GBlock(nNodes - 1, nNodes) -= 0.5 * _parInvWeights[parType][nNodes - 1];
		//		GBlock(nNodes - 1, nNodes + 1) += 0.5 * _parInvWeights[parType][nNodes - 1];

		//		// numerical flux contribution
		//		MatrixXd gStarDC = MatrixXd::Zero(nNodes, 3 * nNodes);
		//		gStarDC.block(0, nNodes - 1, 1, nNodes + 2) = GBlock.block(0, 0, 1, nNodes + 2);
		//		gStarDC.block(0, 0, 1, nNodes + 1) += GBlock.block(nNodes - 1, 1, 1, nNodes + 1);
		//		gStarDC.block(nNodes - 1, nNodes - 1, 1, nNodes + 2) += GBlock.block(nNodes - 1, 0, 1, nNodes + 2);
		//		gStarDC.block(nNodes - 1, 2 * nNodes - 1, 1, nNodes + 1) += GBlock.block(0, 0, 1, nNodes + 1);
		//		gStarDC *= 0.5;

		//		dispBlock.setZero();
		//		// dispersion block part without metrics
		//		dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) = -1.0 * _parInvWeights[parType].asDiagonal() * B * GBlock;
		//		dispBlock.block(0, 0, nNodes, 3 * nNodes) += _parInvWeights[parType].asDiagonal() * B * gStarDC;

		//		for (int cell = 1; cell < _nParElem[parType] - 1; cell++) {

		//			// add metric part, dependent on current cell
		//			dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) += _Dr[_offsetMetric[parType] + cell] * GBlock;
		//			invMap = (2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + cell]));
		//			dispBlock *= invMap * invMap;

		//			for (unsigned int colNode = 0; colNode < _nPoints; colNode++) {

		//				unsigned int offset = offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode });
		//				// start at first solid entry of current inner cell
		//				linalg::BandedEigenSparseRowIterator jac_inner(_globalJac, offset + cell * sCell + strideParLiquid());

		//				for (unsigned int node = 0; node < _nParNode[parType]; node++, jac_inner += strideParLiquid())
		//				{
		//					for (unsigned int comp = 0; comp < _nComp; comp++)
		//					{
		//						for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++, ++jac_inner)
		//						{
		//							if (static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) != 0.0)
		//							{
		//								jac_inner[0] += -(static_cast<double>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)])) * dispBlock(node, _nParNode[parType] + node);
		//							}
		//						}
		//					}
		//				}
		//			}

		//			// substract metric part in preparation of next iteration
		//			dispBlock /= invMap * invMap;
		//			dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) -= _Dr[_offsetMetric[parType] + cell] * GBlock;
		//		}
		//	} // if nCells > 1

		//	return 1;
		//}
		///**
		// * @brief adds a state block into the system jacobian.
		// * @param [in] block (sub)block to be added
		// * @param [in] jac row iterator at first (i.e. upper) entry
		// * @param [in] offRowToCol column to row offset (i.e. start at upper left corner of block)
		// * @param [in] idxr Indexer
		// * @param [in] nCells determines how often the block is added (diagonally)
		// * @param [in] stateFactor state dependend factors
		// * @param [in] strideDead how many (dead) states to be jumped over after each state block
		// * @param [in] nStates how many states are concerned, defaults to nComp
		// */
		//void addJacBlock(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, int offRowToCol, Indexer& idxr, unsigned int nCells, unsigned int strideNode, unsigned int nStates, unsigned int strideDead = 0) {

		//	for (unsigned int cell = 0; cell < nCells; cell++) {
		//		for (unsigned int i = 0; i < block.rows(); i++, jac += strideDead) {
		//			for (unsigned int state = 0; state < nStates; state++, ++jac) {
		//				for (unsigned int j = 0; j < block.cols(); j++) {
		//					// row: at current node component
		//					// col: jump to node j
		//					jac[(j - i) * strideNode + offRowToCol] += block(i, j);
		//				}
		//			}
		//		}
		//	}
		//}
		///**
		// * @brief adds a state block into the system jacobian.
		// * @param [in] block (sub)block whose diagonal entries are to be added
		// * @param [in] jac row iterator at first (i.e. upper) entry
		// * @param [in] idxr Indexer
		// * @param [in] surfDiff pointer to surfaceDiffusion storage
		// * @param [in] nonKinetic pointer to binding kinetics
		// * @param [in] type particle type
		// */
		//template<typename ParamType>
		//void addDiagonalSolidJacobianEntries(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, Indexer& idxr, ParamType* surfDiff, const int* nonKinetic, unsigned int type) {

		//	for (unsigned int i = 0; i < block.rows(); i++, jac += strideParLiquid()) {
		//		for (unsigned int comp = 0; comp < _nComp; comp++) {
		//			for (unsigned int bnd = 0; bnd < _nBound[type * _nComp + comp]; bnd++, ++jac) {
		//				if (static_cast<double>(surfDiff[_offsetSurfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]]) != 0.0
		//					&& !nonKinetic[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]) {
		//					// row, col: at current node and bound state
		//					jac[0] += block(i, i)
		//						* static_cast<double>(surfDiff[_offsetSurfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]]);
		//				}
		//			}
		//		}
		//	}
		//}
		///**
		// * @brief adds a state block into the particle jacobian.
		// * @param [in] block (sub)block to be added
		// * @param [in] jac row iterator at first (i.e. upper) entry
		// * @param [in] idxr Indexer
		// * @param [in] idxr parDiff pointer to particle diffusion parameters
		// * @param [in] idxr surfDiff pointer to particle surface diffusion parameters
		// * @param [in] idxr beta_p pointer to particle porosity parameters
		// * @param [in] idxr nonKinetic pointer to binding kinetics parameters
		// * @param [in] type particle type
		// * @param [in] nBlocks number of blocks, i.e. cells/elements, to be inserted
		// * @param [in] offRowToCol column to row offset (i.e. start at upper left corner of block)
		// */
		//void insertParJacBlock(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, Indexer& idxr, const active* const parDiff, const active* const surfDiff, const double* const beta_p, const int* nonKinetic, unsigned int type, unsigned int nBlocks, int offRowToCol) {

		//	for (unsigned int cell = 0; cell < nBlocks; cell++) {
		//		for (unsigned int i = 0; i < block.rows(); i++) {
		//			for (unsigned int comp = 0; comp < _nComp; comp++, ++jac) {
		//				for (unsigned int j = 0; j < block.cols(); j++) {
		//					/* liquid on liquid blocks */
		//					// row: at current node and component; col: jump to node j
		//					jac[(j - i) * strideParNode(type) + offRowToCol] = block(i, j) * static_cast<double>(parDiff[comp]);
		//				}
		//				/* liquid on solid blocks */
		//				for (unsigned int bnd = 0; bnd < _nBound[type * _nComp + comp]; bnd++) {
		//					if (static_cast<double>(surfDiff[_offsetSurfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]]) != 0.0) {
		//						for (unsigned int j = 0; j < block.cols(); j++) {
		//							// row: at current node and component; col: jump to node j and to current bound state
		//							jac[(j - i) * strideParNode(type) + offRowToCol + strideParLiquid() - comp
		//								+ offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd
		//							]
		//								= block(i, j) * static_cast<double>(beta_p[comp])
		//								* static_cast<double>(surfDiff[_offsetSurfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]]);
		//						}
		//					}
		//				}
		//			}
		//			/* solid on solid blocks */
		//			for (unsigned int comp = 0; comp < _nComp; comp++) {
		//				for (unsigned int bnd = 0; bnd < _nBound[type * _nComp + comp]; bnd++, ++jac) {
		//					if (static_cast<double>(surfDiff[_offsetSurfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]]) != 0.0
		//						&& !nonKinetic[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]) {
		//						for (unsigned int j = 0; j < block.cols(); j++) {
		//							// row: at current node and bound state; col: jump to node j
		//							jac[(j - i) * strideParNode(type) + offRowToCol + bnd]
		//								= block(i, j)
		//								* static_cast<double>(surfDiff[_offsetSurfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]]);
		//						}
		//					}
		//				}
		//			}
		//		}
		//	}
		//}
		///**
		// * @brief analytically calculates the static (per section) particle jacobian
		// * @return 1 if jacobain calculation fits the predefined pattern of the jacobian, 0 if not.
		// */
		//int calcStaticAnaParticleDiffJacobian(unsigned int parType, unsigned int colNode, const active* const parDiff, const active* const parSurfDiff, const double* const invBetaP, linalg::BandedEigenSparseRowIterator jac) {

		//	// DG particle dispersion Jacobian
		//	if (_parExactInt[parType])
		//		calcParticleDGSEMJacobian(parType, colNode, parDiff, parSurfDiff, invBetaP, jac);
		//	else // deprecated
		//		calcParticleCollocationDGSEMJacobian(parType, colNode, parDiff, parSurfDiff, invBetaP, jac);

		//	return true;
		//}


		//void setJacobianPattern_GRM(SparseMatrix<double, RowMajor>& globalJ, unsigned int secIdx, bool hasBulkReaction) {

		//	Indexer idxr(_disc);

		//	std::vector<T> tripletList;
		//	// reserve space for all entries
		//	int bulkEntries = _convDispOp.nConvDispEntries(false);
		//	if (hasBulkReaction)
		//		bulkEntries += _nPoints * _nComp * _nComp; // add nComp entries for every component at each discrete bulk point

		//	// particle
		//	int addTimeDer = 0; // additional time derivative entries: bound states in particle dispersion equation
		//	int isothermNNZ = 0;
		//	int particleEntries = 0;
		//	for (int type = 0; type < _nParType; type++) {
		//		isothermNNZ = (strideParNode(type)) * _nParPoints[type] * _strideBound[type]; // every bound satte might depend on every bound and liquid state
		//		addTimeDer = _nParPoints[type] * _strideBound[type];
		//		particleEntries += calcParDispNNZ(type) + addTimeDer + isothermNNZ;
		//	}

		//	int fluxEntries = 4 * _nParType * _nPoints * _nComp;

		//	tripletList.reserve(fluxEntries + bulkEntries + particleEntries);

		//	// NOTE: inlet and jacF flux jacobian are set in calc jacobian function (identity matrices)
		//	// Note: flux jacobian (identity matrix) is handled in calc jacobian function
		//	unsigned int bulkOffset = offsetC();
		//	_convDispOp.convDispJacPattern(tripletList, bulkOffset);

		//	// bulk reaction jacobian
		//	if (hasBulkReaction) {
		//		for (unsigned int colNode = 0; colNode < _nPoints; colNode++) {
		//			for (unsigned int comp = 0; comp < _nComp; comp++) {
		//				for (unsigned int toComp = 0; toComp < _nComp; toComp++) {
		//					tripletList.push_back(T(offsetC() + colNode * strideColNode() + comp * strideColComp(),
		//						offsetC() + colNode * strideColNode() + toComp * strideColComp(),
		//						0.0));
		//				}
		//			}
		//		}
		//	}

		//	// particle jacobian (including isotherm and time derivative)
		//	for (int colNode = 0; colNode < _nPoints; colNode++) {
		//		for (int type = 0; type < _nParType; type++) {
		//			setParJacPattern(tripletList, type, colNode, secIdx);
		//		}
		//	}

		//	// flux jacobians
		//	for (unsigned int type = 0; type < _nParType; type++) {
		//		for (unsigned int colNode = 0; colNode < _nPoints; colNode++) {
		//			for (unsigned int comp = 0; comp < _nComp; comp++) {
		//				// add Cl on Cp entries
		//				// row: add bulk offset, jump over previous nodes and components
		//				// col: add flux offset to current parType, jump over previous nodes and components
		//				tripletList.push_back(T(offsetC() + colNode * strideColNode() + comp * strideColComp(),
		//					offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ colNode }) + (_nParPoints[type] - 1) * strideParNode(type) + comp * strideParComp(), 0.0));

		//				// add Cp on Cl entries
		//				if (!_parExactInt[type])
		//					// row: add particle offset to current parType and particle, go to last node and add component offset
		//					// col: add flux offset to current component, jump over previous nodes and components
		//					tripletList.push_back(T(offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ colNode })
		//						+ (_nParPoints[type] - 1) * strideParNode(type) + comp * strideParComp(),
		//						offsetC() + colNode * strideColNode() + comp, 0.0));
		//				else {
		//					for (unsigned int node = 0; node < _nParNode[type]; node++) {
		//						// row: add particle offset to current parType and particle, go to last cell and current node and add component offset
		//						// col: add flux offset to current component, jump over previous nodes and components
		//						tripletList.push_back(T(offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ colNode }) + (_nParElem[type] - 1) * _nParNode[type] * strideParNode(type)
		//							+ node * strideParNode(type) + comp * strideParComp(),
		//							offsetC() + colNode * strideColNode() + comp, 0.0));
		//					}
		//				}
		//			}
		//		}
		//	}

		//	globalJ.setFromTriplets(tripletList.begin(), tripletList.end());
		//}

		//int calcFluxJacobians(unsigned int secIdx, bool outliersOnly = false) {

		//	Indexer idxr(_disc);

		//	for (unsigned int type = 0; type < _nParType; type++) {

		//		// lifting matrix entry for exact integration scheme depends on metrics for sphere and cylinder
		//		double exIntLiftContribution = static_cast<double>(_Ir[_offsetMetric[type] + _nParElem[type] - 1][_nParNode[type] - 1]);
		//		if (_parGeomSurfToVol[type] == _SurfVolRatioSlab)
		//			exIntLiftContribution = 1.0;

		//		// Ordering of diffusion:
		//		// sec0type0comp0, sec0type0comp1, sec0type0comp2, sec0type1comp0, sec0type1comp1, sec0type1comp2,
		//		// sec1type0comp0, sec1type0comp1, sec1type0comp2, sec1type1comp0, sec1type1comp1, sec1type1comp2, ...
		//		active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _nComp * _nParType, secIdx) + type * _nComp;

		//		linalg::BandedEigenSparseRowIterator jacCl(_globalJac, offsetC());
		//		linalg::BandedEigenSparseRowIterator jacCp(_globalJac, offsetCp(ParticleTypeIndex{ type }) + (_nParPoints[type] - 1) * strideParNode(type));

		//		for (unsigned int colNode = 0; colNode < _nPoints; colNode++)
		//		{
		//			for (unsigned int comp = 0; comp < _nComp; comp++, ++jacCp, ++jacCl) {
		//				// add Cl on Cl entries (added since these entries are also touched by bulk jacobian)
		//				// row: already at bulk phase. already at current node and component.
		//				// col: already at bulk phase. already at current node and component.
		//				if (!outliersOnly)
		//					jacCl[0] += static_cast<double>(filmDiff[comp]) * (1.0 - static_cast<double>(_colPorosity)) / static_cast<double>(_colPorosity)
		//					* _parGeomSurfToVol[type] / static_cast<double>(_parRadius[type])
		//					* _parTypeVolFrac[type + colNode * _nParType].getValue();
		//				// add Cl on Cp entries (added since these entries are also touched by bulk jacobian)
		//				// row: already at bulk phase. already at current node and component.
		//				// col: go to current particle phase entry.
		//				jacCl[jacCp.row() - jacCl.row()] = -static_cast<double>(filmDiff[comp]) * (1.0 - static_cast<double>(_colPorosity)) / static_cast<double>(_colPorosity)
		//					* _parGeomSurfToVol[type] / static_cast<double>(_parRadius[type])
		//					* _parTypeVolFrac[type + colNode * _nParType].getValue();

		//				// add Cp on Flux entries
		//				if (!_parExactInt[type]) {
		//					// row: already at particle. already at current node and liquid state.
		//					// col: already at particle. already at current node and liquid state.
		//					if (!outliersOnly) // Cp on Cb
		//						jacCp[0] += static_cast<double>(filmDiff[comp]) * 2.0 / static_cast<double>(_deltaR[_offsetMetric[type]]) * _parInvWeights[type][0] / static_cast<double>(_parPorosity[type]) / static_cast<double>(_poreAccessFactor[type * _nComp + comp]);
		//					// row: already at particle. already at current node and liquid state.
		//					// col: go to current bulk phase.
		//					jacCp[jacCl.row() - jacCp.row()] = -static_cast<double>(filmDiff[comp]) * 2.0 / static_cast<double>(_deltaR[_offsetMetric[type]]) * _parInvWeights[type][0] / static_cast<double>(_parPorosity[type]) / static_cast<double>(_poreAccessFactor[type * _nComp + comp]);
		//				}
		//				else {
		//					unsigned int entry = jacCp.row();
		//					for (int node = _parPolyDeg[type]; node >= 0; node--, jacCp -= strideParNode(type)) {
		//						// row: already at particle. Already at current node and liquid state.
		//						// col: original entry at outer node.
		//						if (!outliersOnly) // Cp on Cb
		//							jacCp[entry - jacCp.row()]
		//							+= static_cast<double>(filmDiff[comp]) * 2.0 / static_cast<double>(_deltaR[_offsetMetric[type]]) * _parInvMM[_offsetMetric[type] + _nParElem[type] - 1](node, _nParNode[type] - 1) * exIntLiftContribution / static_cast<double>(_parPorosity[type]) / static_cast<double>(_poreAccessFactor[type * _nComp + comp]);
		//						// row: already at particle. Already at current node and liquid state.
		//						// col: go to current bulk phase.
		//						jacCp[jacCl.row() - jacCp.row()]
		//							= -static_cast<double>(filmDiff[comp]) * 2.0 / static_cast<double>(_deltaR[_offsetMetric[type]]) * _parInvMM[_offsetMetric[type] + _nParElem[type] - 1](node, _nParNode[type] - 1) * exIntLiftContribution / static_cast<double>(_parPorosity[type]) / static_cast<double>(_poreAccessFactor[type * _nComp + comp]);
		//					}
		//					// set back iterator to first node as required by component loop
		//					jacCp += _nParNode[type] * strideParNode(type);
		//				}
		//			}
		//			if (colNode < _nPoints - 1) // execute iteration statement only when condition is true in next loop.
		//				jacCp += _strideBound[type] + (_nParPoints[type] - 1) * strideParNode(type);
		//		}
		//	}

		//	return 1;
		//}

	};

} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_PARTCICLEDIFFUSIONOPERATORDG_HPP_

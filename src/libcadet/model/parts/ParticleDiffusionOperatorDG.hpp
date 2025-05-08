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

		int residual(double t, unsigned int parType, unsigned int colNode, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, int const* const qsBinding, WithoutParamSensitivity);
		int residual(double t, unsigned int parType, unsigned int colNode, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithParamSensitivity);
		int residual(double t, unsigned int parType, unsigned int colNode, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithParamSensitivity);
		int residual(double t, unsigned int parType, unsigned int colNode, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithoutParamSensitivity);


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
		int residualImpl(double t, unsigned int parType, unsigned int colNode, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, int const* const qsBinding);

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
	};

} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_PARTCICLEDIFFUSIONOPERATORDG_HPP_

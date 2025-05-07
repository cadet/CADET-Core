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

#ifndef LIBCADET_PARTICLEDISPERSIONOPERATORDG_HPP_
#define LIBCADET_PARTICLEDISPERSIONOPERATORDG_HPP_

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
	class ParticleDispersionOperatorDG
	{
	public:

		ParticleDispersionOperatorDG();
		~ParticleDispersionOperatorDG() CADET_NOEXCEPT;

		bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int nParType, const int strideBulkComp);
		bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters);
		
		void setEquidistantRadialDisc(unsigned int parType);
		void setEquivolumeRadialDisc(unsigned int parType);
		void setUserdefinedRadialDisc(unsigned int parType);
		void updateRadialDisc();

		void clearParDepSurfDiffusion();

		bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active const* const filmDiff);

		int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yBulk, double const* yDot, double* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac);
		int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yBulk, double const* yDot, active* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac);
		int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yBulk, double const* yDot, double* res, WithoutParamSensitivity);
		int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yBulk, double const* yDot, active* res, WithParamSensitivity);
		int residual(const IModel& model, double t, unsigned int secIdx, active const* y, active const* yBulk, double const* yDot, active* res, WithParamSensitivity);
		int residual(const IModel& model, double t, unsigned int secIdx, active const* y, active const* yBulk, double const* yDot, active* res, WithoutParamSensitivity);


		/* Physical model parameters */

		std::vector<active> _parRadius; //!< Particle radius \f$ r_p \f$
		bool _singleParRadius;
		std::vector<active> _parCoreRadius; //!< Particle core radius \f$ r_c \f$
		bool _singleParCoreRadius;
		std::vector<active> _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$
		bool _singleParPorosity;
		std::vector<double> _parGeomSurfToVol; //!< Particle surface to volume ratio factor (i.e., 3.0 for spherical, 2.0 for cylindrical, 1.0 for hexahedral)
		std::vector<active> _parOuterSurfAreaPerVolume; //!< Particle shell outer sphere surface to volume ratio
		std::vector<active> _parInnerSurfAreaPerVolume; //!< Particle shell inner sphere surface to volume ratio

		std::vector<active> _filmDiffusion; //!< Particle diffusion coefficient \f$ D_p \f$
		MultiplexMode _filmDiffusionMode;
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
			 * Equidistant distribution of shell edges
			 */
			Equidistant,

			/**
			 * Volumes of shells are uniform
			 */
			Equivolume,

			/**
			 * Shell edges specified by user
			 */
			UserDefined
		};

		// todo when separating discretization and particle model
		//class Indexer
		//{
		//public:

		//	Indexer(const Discretization& disc) : _disc(disc) { }

		//	inline int strideColNode() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }

		//	protected:
		//		_disc;
		//};

		int _strideBulkComp;

		std::vector<ParticleDiscretizationMode> _parDiscMode; //!< Particle discretization mode

		std::vector<double> _parDiscVector; //!< Particle discretization shell edges

		unsigned int* _nParElem; //!< Array with number of radial elements in each particle type
		unsigned int* _nParPointsBeforeType; //!< Array with total number of radial points before a particle type (cumulative sum of nParPoints), additional last element contains total number of particle shells
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

		std::vector<active> _parShellSize; //!< Particle shell size
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

		template <typename StateType, typename ResidualType, typename ParamType, typename RowIteratorType, bool wantJac>
		int residualImpl(const IModel& model, double t, unsigned int secIdx, StateType const* y, StateType const* yBulk, double const* yDot, ResidualType* res, RowIteratorType jacBegin);

	};

} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_PARTICLEDISPERSIONOPERATORDG_HPP_

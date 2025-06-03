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
#include "model/ModelUtils.hpp"
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

	class IConfigHelper;

namespace model
{

class IDynamicReactionModel;
class IParameterStateDependence;

namespace parts
{
	namespace cell
	{
		struct CellParameters;
	}

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

		bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp);
		bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound, const int* reqBinding, const bool hasDynamicReactions);

		void setEquidistantRadialDisc();
		void setEquivolumeRadialDisc();
		void setUserdefinedRadialDisc();
		void updateRadialDisc();

		void clearParDepSurfDiffusion();

		bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active const* const filmDiff, active const* const poreAccessFactor);

		/**
		 * @brief Computes the residual of the transport equations
		 * @param [in] model Model that owns the operator
		 * @param [in] t Current time point
		 * @param [in] secIdx Index of the current section
		 * @param [in] yPar Pointer to particle phase entry in unit state vector
		 * @param [in] yBulk Pointer to corresponding bulk phase entry in unit state vector
		 * @param [in] yDotPar Pointer to particle phase derivative entry in unit state vector
		 * @param [out] resPar Pointer Pointer to particle phase entry in unit residual vector
		 * @param [out] colPos column position of the particle (particle coordinate zero)
		 * @param [in] jacIt Matrix iterator pointing to the particle phase entry in the unit Jacobian
		 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
		 */
		template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
		int residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, linalg::BandedEigenSparseRowIterator& jacBase);

		unsigned int _parTypeIdx; //!< Particle type index (wrt the unit operation that owns this particle model)

		/* Physical model parameters */

		active _parRadius; //!< Particle radius \f$ r_p \f$
		bool _singleParRadius;
		active _parCoreRadius; //!< Particle core radius \f$ r_c \f$
		bool _singleParCoreRadius;
		active _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$
		bool _singleParPorosity;
		std::vector<active> _invBetaP; //!< Ratio of solid to liquid particle volume
		double _parGeomSurfToVol; //!< Particle surface to volume ratio factor (i.e., 3.0 for spherical, 2.0 for cylindrical, 1.0 for hexahedral)
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
		IParameterStateDependence* _parDepSurfDiffusion; //!< Parameter dependencies for particle surface diffusion
		bool _singleParDepSurfDiffusion; //!< Determines whether a single parameter dependence for particle surface diffusion is used
		bool _hasParDepSurfDiffusion; //!< Determines whether particle surface diffusion parameter dependencies are present
		bool _hasSurfaceDiffusion; //!< Determines whether surface diffusion is present
		const int* _reqBinding; //!< Array of size @p _strideBound with flags whether binding is in rapid equilibrium for each bound state
		bool _hasDynamicReactions; //! Determines whether or not the binding has any dynamic reactions

		unsigned int _nComp; //!< Number of components

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

		double relativeCoordinate(const unsigned int nodeIdx) const CADET_NOEXCEPT
		{
			const unsigned int element = floor(nodeIdx / _nParNode);
			const unsigned int node = nodeIdx % _nParNode;
			// divide by particle radius to get relative position
			return static_cast<double>((_deltaR[element] * element + 0.5 * _deltaR[element] * (1 + _parNodes[node])) / (_parRadius - _parCoreRadius));
		}

		template<typename ParamType>
		ParamType surfaceToVolumeRatio() const CADET_NOEXCEPT
		{
			return _parGeomSurfToVol / static_cast<ParamType>(_parRadius);
		}

		int _strideBulkComp;
		inline int strideBulkComp() const CADET_NOEXCEPT { return _strideBulkComp; }
		inline int strideParComp() const CADET_NOEXCEPT { return 1; }
		inline int strideParLiquid() const CADET_NOEXCEPT { return static_cast<int>(_nComp); }
		inline int strideParBound() const CADET_NOEXCEPT { return static_cast<int>(_strideBound); }
		inline int strideParNode() const CADET_NOEXCEPT { return strideParLiquid() + strideParBound(); }
		inline int strideParElem() const CADET_NOEXCEPT { return strideParNode() * _nParNode; }
		inline int strideParBlock() const CADET_NOEXCEPT { return static_cast<int>(_nParPoints) * strideParNode(); }
		inline int offsetBoundComp(ComponentIndex comp) const CADET_NOEXCEPT { return _boundOffset[comp.value]; }

		ParticleDiscretizationMode _parDiscMode; //!< Particle discretization mode

		std::vector<double> _parDiscVector; //!< Particle discretization element boundary coodinates

		unsigned int _nParElem; //!< Array with number of radial elements
		unsigned int _parPolyDeg; //!< polynomial degree of particle elements
		unsigned int _nParNode; //!< Array with number of radial nodes per element
		unsigned int _nParPoints; //!< Array with number of radial nodes per element
		bool _parGSM; //!< specifies whether (single element) Galerkin spectral method should be used in particles
		unsigned int* _nBound; //!< Array with number of bound states for each component
		unsigned int* _boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
		unsigned int _strideBound; //!< Total number of bound states

		std::vector<active> _parElementSize; //!< Particle element size
		std::vector<active> _parCenterRadius; //!< Particle node-centered position for each particle node

		/* DG specific operators */

		active* _deltaR; //!< equidistant particle element spacing
		Eigen::VectorXd _parNodes; //!< Array with positions of nodes in radial reference element for each particle
		Eigen::MatrixXd _parPolyDerM; //!< Array with polynomial derivative Matrix for each particle
		Eigen::MatrixXd* _minus_InvMM_ST; //!< equals minus inverse mass matrix times transposed stiffness matrix.
		Eigen::VectorXd _parInvWeights; //!< Array with weights for LGL quadrature of size nNodes for each particle
		Eigen::MatrixXd* _parInvMM; //!< dense inverse mass matrix for exact integration of integrals with metrics, for each particle
		Eigen::MatrixXd _parInvMM_Leg; //!< dense inverse mass matrix (Legendre) for exact integration of integral without metric, for each particle
		Eigen::MatrixXd _secondOrderStiffnessM; //!< specific second order stiffness matrix
		Eigen::MatrixXd _minus_parInvMM_Ar; //!< inverse mass matrix times specific second order stiffness matrix
		Eigen::Vector<active, Dynamic>* _Ir; //!< metric part for each element

		Eigen::MatrixXd* _DGjacParDispBlocks; //!< particle dispersion blocks of DG jacobian

		Eigen::Vector<active, Dynamic> _g_p; //!< auxiliary variable g = dc_p / dr
		Eigen::Vector<active, Dynamic> _g_pSum; //!< auxiliary variable g = sum_{k \in p, s_i} dc_k / dr
		Eigen::Vector<active, Dynamic> _surfaceFluxParticle; //!< stores the surface flux values for each particle
		active* _localFlux; //!< stores the local (at respective particle) film diffusion flux

		void initializeDG();

		void initializeDGjac(const double parGeomSurfToVol);

		int calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool outliersOnly = false);

		int getParticleCoordinates(double* coords) const;

		typedef Eigen::Triplet<double> T;

		void calcParticleJacobianPattern(std::vector<T>& tripletList, unsigned int offsetPar, unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx);
		void parTimeDerJacPattern_GRM(std::vector<T>& tripletList, unsigned int offset, unsigned int colNode, unsigned int secIdx);
		void parBindingPattern_GRM(std::vector<T>& tripletList, const int offset, const unsigned int colNode);

		unsigned int calcParDiffNNZ();

		int calcStaticAnaParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac);

		bool setParameter(const ParameterId& pId, double value);
		bool setParameter(const ParameterId& pId, int value);
		bool setParameter(const ParameterId& pId, bool value);
		bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue);
		bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value);
		
	protected:

		int addSolidDGentries(const int secIdx, linalg::BandedEigenSparseRowIterator& jacBase, const int* const reqBinding);

		template<typename ResidualType, typename ParamType>
		void applyParInvMap(Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>>& state);

		template<typename StateType, typename ResidualType>
		void parGSMVolumeIntegral(Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer);

		template<typename StateType, typename ResidualType>
		void parVolumeIntegral(const bool aux, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer);

		template<typename StateType>
		void InterfaceFluxParticle(Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, const unsigned int strideCell, const unsigned int strideNode, const bool aux, const int comp, const bool addParDisc = false);

		template<typename StateType, typename ResidualType>
		void parSurfaceIntegral(Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state,
			Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer, unsigned const int strideCell, unsigned const int strideNode,
			const bool aux, const int comp = 0, const bool addParDisc = false);

		template<typename StateType>
		void solve_auxiliary_DG(Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<>>& conc, unsigned int strideCell, unsigned int strideNode, int comp);

		/*		DG particle Jacobian		*/

		Eigen::MatrixXd DGjacobianParDispBlock(unsigned int elemIdx,double parGeomSurfToVol);
		Eigen::MatrixXd GSMjacobianParDispBlock(double parGeomSurfToVol);
		Eigen::MatrixXd getParBMatrix(int element, double parGeomSurfToVol);
		Eigen::MatrixXd parAuxBlockGstar(unsigned int elemIdx, MatrixXd leftG, MatrixXd middleG, MatrixXd rightG);
		Eigen::MatrixXd getParGBlock(unsigned int elemIdx);

		void insertParJacBlock(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, const active* const parDiff, const active* const surfDiff, const active* const beta_p, const int* nonKinetic, unsigned int nBlocks, int offRowToCol);
		void addDiagonalSolidJacobianEntries(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, const int* const reqBinding, const active* const surfDiffPtr);

	};

} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_PARTCICLEDIFFUSIONOPERATORDG_HPP_

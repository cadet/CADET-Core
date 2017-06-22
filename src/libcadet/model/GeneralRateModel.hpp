// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2017: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines the general rate model (GRM).
 */

#ifndef LIBCADET_GENERALRATEMODEL_HPP_
#define LIBCADET_GENERALRATEMODEL_HPP_

#include "UnitOperation.hpp"
#include "model/BindingModel.hpp"
#include "cadet/SolutionExporter.hpp"
#include "AutoDiff.hpp"
#include "linalg/SparseMatrix.hpp"
#include "linalg/Gmres.hpp"
#include "MemoryPool.hpp"
#include "ParamIdUtil.hpp"
#include "Weno.hpp"
#include "model/ModelUtils.hpp"

#include <array>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <tuple>

#include "Benchmark.hpp"

namespace cadet
{

// Forward declarations
namespace linalg
{
	class BandMatrix;
	class FactorizableBandMatrix;
}

namespace model
{

/**
 * @brief General rate model of liquid column chromatography
 * @details See @cite Guiochon2006, @cite Gu1995, @cite Felinger2004
 * 
 * @f[\begin{align}
	\frac{\partial c_i}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax}} \frac{\partial^2 c_i}{\partial z^2} - \frac{1 - \varepsilon_c}{\varepsilon_c} \frac{3}{r_p} j_{f,i} \\
	\frac{\partial c_{p,i}}{\partial t} + \frac{1 - \varepsilon_p}{\varepsilon_p} \frac{\partial q_{i}}{\partial t} &= D_{p,i} \left( \frac{\partial^2 c_{p,i}}{\partial r^2} + \frac{2}{r} \frac{\partial c_{p,i}}{\partial r} \right) + D_{s,i} \frac{1 - \varepsilon_p}{\varepsilon_p} \left( \frac{\partial^2 q_{i}}{\partial r^2} + \frac{2}{r} \frac{\partial q_{i}}{\partial r} \right) \\
	a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c_p, q)
\end{align} @f]
@f[ \begin{align}
	j_{f,i} = k_{f,i} \left( c_i - c_{p,i} \left(\cdot, \cdot, r_p\right)\right)
\end{align} @f]
 * Danckwerts boundary conditions (see @cite Danckwerts1953)
@f[ \begin{align}
u c_{\text{in},i}(t) &= u c_i(t,0) - D_{\text{ax}} \frac{\partial c_i}{\partial z}(t,0) \\
\frac{\partial c_i}{\partial z}(t,L) &= 0 \\
\varepsilon_p D_{p,i} \frac{\partial c_{p,i}}{\partial r}(\cdot, \cdot, r_p) + (1-\varepsilon_p) D_{s,i} \frac{\partial q_{i}}{\partial r}(\cdot, \cdot, r_p) &= j_{f,i} \\
\frac{\partial c_{p,i}}{\partial r}(\cdot, \cdot, 0) &= 0
\end{align} @f]
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), @cite Puttmann2013 @cite Puttmann2016 (forward sensitivities, AD, band compression)
 */
class GeneralRateModel : public IUnitOperation
{
public:

	GeneralRateModel(UnitOpIdx unitOpIdx);
	virtual ~GeneralRateModel() CADET_NOEXCEPT;

	virtual unsigned int numDofs() const CADET_NOEXCEPT;
	virtual unsigned int numPureDofs() const CADET_NOEXCEPT;
	virtual bool usesAD() const CADET_NOEXCEPT;
	virtual unsigned int requiredADdirs() const CADET_NOEXCEPT;

	virtual UnitOpIdx unitOperationId() const CADET_NOEXCEPT { return _unitOpIdx; }
	virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
	virtual void setFlowRates(const active& in, const active& out) CADET_NOEXCEPT;
	virtual bool canAccumulate() const CADET_NOEXCEPT { return false; }

	static const char* identifier() { return "GENERAL_RATE_MODEL"; }
	virtual const char* unitOperationName() const CADET_NOEXCEPT { return "GENERAL_RATE_MODEL"; }

	virtual bool configure(IParameterProvider& paramProvider, IConfigHelper& helper);
	virtual bool reconfigure(IParameterProvider& paramProvider);
	virtual void notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset);

	virtual std::unordered_map<ParameterId, double> getAllParameterValues() const;
	virtual bool hasParameter(const ParameterId& pId) const;

	virtual bool setParameter(const ParameterId& pId, int value);
	virtual bool setParameter(const ParameterId& pId, double value);
	virtual bool setParameter(const ParameterId& pId, bool value);

	virtual bool setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue);
	virtual void setSensitiveParameterValue(const ParameterId& id, double value);

	virtual void clearSensParams();

	virtual void useAnalyticJacobian(const bool analyticJac);

	virtual void reportSolution(ISolutionRecorder& recorder, double const* const solution) const;
	virtual void reportSolutionStructure(ISolutionRecorder& recorder) const;

	virtual int residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res);

	virtual int residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int adDirOffset);

	virtual int residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
		double const* const y, double const* const yDot, active* const adRes);

	virtual int residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, active* const adRes, active* const adY, unsigned int adDirOffset);

	virtual int residualSensFwdCombine(const active& timeFactor, const std::vector<const double*>& yS, const std::vector<const double*>& ySdot,
		const std::vector<double*>& resS, active const* adRes, double* const tmp1, double* const tmp2, double* const tmp3);

	virtual int linearSolve(double t, double timeFactor, double alpha, double tol, double* const rhs, double const* const weight,
		double const* const y, double const* const yDot, double const* const res);

	virtual void prepareADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const;

	virtual void applyInitialCondition(double* const vecStateY, double* const vecStateYdot) { }
	virtual void applyInitialCondition(IParameterProvider& paramProvider, double* const vecStateY, double* const vecStateYdot);

	virtual void consistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol);
	virtual void consistentInitialTimeDerivative(double t, unsigned int secIdx, double timeFactor, double const* vecStateY, double* const vecStateYdot);

	virtual void consistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes);

	virtual void leanConsistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol);
	virtual void leanConsistentInitialTimeDerivative(double t, double timeFactor, double* const vecStateYdot, double* const res);

	virtual void leanConsistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes);

	virtual bool hasInlet() const CADET_NOEXCEPT { return true; }
	virtual bool hasOutlet() const CADET_NOEXCEPT { return true; }

	virtual unsigned int localOutletComponentIndex() const CADET_NOEXCEPT;
	virtual unsigned int localOutletComponentStride() const CADET_NOEXCEPT;
	virtual unsigned int localInletComponentIndex() const CADET_NOEXCEPT;
	virtual unsigned int localInletComponentStride() const CADET_NOEXCEPT;

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size);
	virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections) { }

	virtual void expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut);

	virtual void multiplyWithJacobian(double const* yS, double alpha, double beta, double* ret);
	virtual void multiplyWithDerivativeJacobian(double const* sDot, double* ret, double timeFactor);
	inline void multiplyWithJacobian(double const* yS, double* ret)
	{
		multiplyWithJacobian(yS, 1.0, 0.0, ret);
	}

#ifdef CADET_BENCHMARK_MODE
	virtual std::vector<double> benchmarkTimings() const
	{
		return std::vector<double>({
			_timerResidual.totalElapsedTime(),
			_timerResidualPar.totalElapsedTime(),
			_timerResidualSens.totalElapsedTime(),
			_timerResidualSensPar.totalElapsedTime(),
			_timerJacobianPar.totalElapsedTime(),
			_timerConsistentInit.totalElapsedTime(),
			_timerConsistentInitPar.totalElapsedTime(),
			_timerLinearSolve.totalElapsedTime(),
			_timerFactorize.totalElapsedTime(),
			_timerFactorizePar.totalElapsedTime(),
			_timerMatVec.totalElapsedTime(),
			_timerGmres.totalElapsedTime()
		});
	}

	virtual char const* const* benchmarkDescriptions() const
	{
		static const char* const desc[] = {
			"Residual",
			"ResidualPar",
			"ResidualSens",
			"ResidualSensPar",
			"JacobianPar",
			"ConsistentInit",
			"ConsistentInitPar",
			"LinearSolve",
			"Factorize",
			"FactorizePar",
			"MatVec",
			"Gmres"
		};
		return desc;
	}
#endif

protected:

	class Indexer;

	virtual int residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int adDirOffset, bool updateJacobian, bool paramSensitivity);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* const y, double const* const yDot, ResidualType* const res);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualBulk(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualBulkForwardsFlow(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualBulkBackwardsFlow(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualParticle(const ParamType& t, unsigned int colCell, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res);

	template <typename StateType, typename ResidualType, typename ParamType>
	int residualFlux(const ParamType& t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res);

	void assembleOffdiagJac(double t, unsigned int secIdx);
	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);
	void prepareBulkADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const;

	int schurComplementMatrixVector(double const* x, double* z) const;
	void assembleDiscretizedJacobianColumnBlock(unsigned int comp, double alpha, const Indexer& idxr, double timeFactor);
	void assembleDiscretizedJacobianParticleBlock(unsigned int pblk, double alpha, const Indexer& idxr, double timeFactor);
	
	void setEquidistantRadialDisc();
	void setEquivolumeRadialDisc();
	void setUserdefinedRadialDisc(const std::vector<double>& cellInterfaces);

	void addTimeDerivativeToJacobianColumnBlock(linalg::FactorizableBandMatrix& fbm, const Indexer& idxr, double alpha, double timeFactor);
	void addMobilePhaseTimeDerivativeToJacobianParticleBlock(linalg::FactorizableBandMatrix::RowIterator& jac, const Indexer& idxr, double alpha, double invBetaP, double timeFactor);
	void solveForFluxes(double* const vecState, const Indexer& idxr);

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	void checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	struct Discretization
	{
		unsigned int nComp; //!< Number of components
		unsigned int nCol; //!< Number of column cells
		unsigned int nPar; //!< Number of radial cells in each particle
		unsigned int* nBound; //!< Array with number of bound states for each component
		unsigned int* boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
		unsigned int strideBound; //!< Total number of bound states
	};

	UnitOpIdx _unitOpIdx; //!< Unit operation index
	Discretization _disc; //!< Discretization info
	IBindingModel* _binding; //!<  Binding model
//	IExternalFunction* _extFun; //!< External function (owned by library user)

	linalg::BandMatrix* _jacC; //!< Interstitial jacobian diagonal block
	linalg::BandMatrix* _jacP; //!< Particle jacobian diagonal blocks (all of them)

	linalg::DoubleSparseMatrix _jacCF; //!< Jacobian block connecting interstitial states and fluxes (interstitial transport equation)
	linalg::DoubleSparseMatrix _jacFC; //!< Jacobian block connecting fluxes and interstitial states (flux equation)
	linalg::DoubleSparseMatrix* _jacPF; //!< Jacobian blocks connecting particle states and fluxes (particle transport boundary condition)
	linalg::DoubleSparseMatrix* _jacFP; //!< Jacobian blocks connecting fluxes and particle states (flux equation)

	linalg::FactorizableBandMatrix* _jacCdisc; //!< Interstitial jacobian diagonal block with time derivatives from BDF method
	linalg::FactorizableBandMatrix* _jacPdisc; //!< Particle jacobian diagonal blocks (all of them) with time derivatives from BDF method

	linalg::DoubleSparseMatrix _jacInlet; //!< Jacobian inlet DOF block matrix connects inlet DOFs to first bulk cells

	active _colLength; //!< Column length \f$ L \f$
	active _colPorosity; //!< Column porosity (external porosity) \f$ \varepsilon_c \f$
	active _parRadius; //!< Particle radius \f$ r_p \f$
	active _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$
	active _crossSection; //!< Cross section area 

	// Section dependent parameters
	std::vector<active> _colDispersion; //!< Column dispersion (may be section dependent) \f$ D_{\text{ax}} \f$
	std::vector<active> _velocity; //!< Interstitial velocity (may be section dependent) \f$ u \f$
	active _curVelocity; //!< Current interstitial velocity \f$ u \f$ in this time section

	// Vectorial parameters
	std::vector<active> _filmDiffusion; //!< Film diffusion coefficient \f$ k_f \f$
	std::vector<active> _parDiffusion; //!< Particle diffusion coefficient \f$ D_p \f$
	std::vector<active> _parSurfDiffusion; //!< Particle surface diffusion coefficient \f$ D_s \f$

	std::unordered_map<ParameterId, active*> _parameters; //!< Provides access to all parameters
	bool _analyticJac; //!< Determines whether AD or analytic Jacobians are used

	ArrayPool _stencilMemory; //!< Provides memory for the stencil
	double* _wenoDerivatives; //!< Holds derivatives of the WENO scheme
	Weno _weno; //!< The WENO scheme implementation
	double _wenoEpsilon; //!< The @f$ \varepsilon @f$ of the WENO scheme (prevents division by zero)

	std::unordered_set<active*> _sensParams; //!< Holds all parameters with activated AD directions
	unsigned int _jacobianAdDirs; //!< Number of AD seed vectors required for Jacobian computation

	std::vector<double> _parCellSize; //!< Particle cell / shell size
	std::vector<double> _parCenterRadius; //!< Particle cell-centered position for each particle cell
	std::vector<double> _parOuterSurfAreaPerVolume;
	std::vector<double> _parInnerSurfAreaPerVolume;

	ArrayPool _discParFlux; //!< Storage for discretized @f$ k_f @f$ value

	bool _factorizeJacobian; //!< Determines whether the Jacobian needs to be factorized
	double* _tempState; //!< Temporary storage with the size of the state vector or nCol * nPar * _binding->consistentInitializationWorkspaceSize() whichever is larger
	linalg::Gmres _gmres; //!< GMRES algorithm for the Schur-complement in linearSolve()
	double _schurSafety; //!< Safety factor for Schur-complement solution

	BENCH_TIMER(_timerResidual)
	BENCH_TIMER(_timerResidualPar)
	BENCH_TIMER(_timerResidualSens)
	BENCH_TIMER(_timerResidualSensPar)
	BENCH_TIMER(_timerJacobianPar)
	BENCH_TIMER(_timerConsistentInit)
	BENCH_TIMER(_timerConsistentInitPar)
	BENCH_TIMER(_timerLinearSolve)
	BENCH_TIMER(_timerFactorize)
	BENCH_TIMER(_timerFactorizePar)
	BENCH_TIMER(_timerMatVec)
	BENCH_TIMER(_timerGmres)

	// Wrapper for calling the corresponding function in GeneralRateModel class
	friend int schurComplementMultiplierGRM(void* userData, double const* x, double* z);

	class Indexer
	{
	public:
		Indexer(const Discretization& disc) : _disc(disc) { }

		// Strides
		inline const int strideColCell() const CADET_NOEXCEPT { return 1; }
		inline const int strideColComp() const CADET_NOEXCEPT { return static_cast<int>(_disc.nCol); }

		inline const int strideParComp() const CADET_NOEXCEPT { return 1; }
		inline const int strideParLiquid() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline const int strideParBound() const CADET_NOEXCEPT { return static_cast<int>(_disc.strideBound); }
		inline const int strideParShell() const CADET_NOEXCEPT { return strideParLiquid() + strideParBound(); }
		inline const int strideParBlock() const CADET_NOEXCEPT { return static_cast<int>(_disc.nPar) * strideParShell(); }

		inline const int strideFluxCell() const CADET_NOEXCEPT { return 1; }
		inline const int strideFluxComp() const CADET_NOEXCEPT { return static_cast<int>(_disc.nCol); }

		// Offsets
		inline const int offsetC() const CADET_NOEXCEPT { return _disc.nComp; }
		inline const int offsetCp() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol + offsetC(); }
		inline const int offsetCp(unsigned int colCell) const CADET_NOEXCEPT { return offsetCp() + strideParBlock() * colCell; }
		inline const int offsetCp(unsigned int colCell, unsigned int parCell) const CADET_NOEXCEPT { return offsetCp() + strideParBlock() * colCell + strideParShell() * parCell; }
		inline const int offsetJf() const CADET_NOEXCEPT { return (_disc.nComp + strideParBlock()) * _disc.nCol + offsetC(); }
		inline const int offsetBoundComp(unsigned int comp) const CADET_NOEXCEPT { return _disc.boundOffset[comp]; }

		// Return pointer to first element of state variable in state vector
		template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
		template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

		template <typename real_t> inline real_t* cp(real_t* const data) const { return data + offsetCp(); }
		template <typename real_t> inline real_t const* cp(real_t const* const data) const { return data + offsetCp(); }

		template <typename real_t> inline real_t* q(real_t* const data) const { return data + offsetCp() + strideParLiquid(); }
		template <typename real_t> inline real_t const* q(real_t const* const data) const { return data + offsetCp() + strideParLiquid(); }

		template <typename real_t> inline real_t* jf(real_t* const data) const { return data + offsetJf(); }
		template <typename real_t> inline real_t const* jf(real_t const* const data) const { return data + offsetJf(); }

		// Return specific variable in state vector
		template <typename real_t> inline real_t& c(real_t* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp * strideColComp() + col]; }
		template <typename real_t> inline const real_t& c(real_t const* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp * strideColComp() + col]; }

		
		template <typename real_t> inline real_t& cp(real_t* const data, unsigned int col, unsigned int par, unsigned int comp) const
		{
			return data[offsetCp() + col * strideParBlock() + par * strideParShell() + comp];
		}
		template <typename real_t> inline const real_t& cp(real_t const* const data, unsigned int col, unsigned int par, unsigned int comp) const
		{
			return data[offsetCp() + col * strideParBlock() + par * strideParShell() + comp];
		}

		
		template <typename real_t> inline real_t& q(real_t* const data, unsigned int col, unsigned int par, unsigned int phase, unsigned int comp) const
		{
			return data[offsetCp() + col * strideParBlock() + par * strideParShell() + strideParLiquid() + offsetBoundComp(comp) + phase];
		}
		template <typename real_t> inline const real_t& q(real_t const* const data, unsigned int col, unsigned int par, unsigned int phase, unsigned int comp) const
		{
			return data[offsetCp() + col * strideParBlock() + par * strideParShell() + strideParLiquid() + offsetBoundComp(comp) + phase];
		}


		template <typename real_t> real_t& jf(real_t* const data, unsigned int col, unsigned int comp) const { return data[offsetJf() + col * strideFluxComp() + comp]; }
		template <typename real_t> const real_t& jf(real_t const* const data, unsigned int col, unsigned int comp) const { return data[offsetJf() + col * strideFluxComp() + comp]; }

		// Iterator-like access
		template <typename real_t> inline real_t& cNextComp(real_t* data) const { return *(data + strideColComp()); }
		template <typename real_t> inline const real_t& cNextComp(const real_t* data) const { return *(data + strideColComp()); }
		template <typename real_t> inline real_t& cNextCol(real_t* data) const { return *(data + strideColCell()); }
		template <typename real_t> inline const real_t& cNextCol(const real_t* data) const { return *(data + strideColCell()); }
		template <typename real_t> inline real_t& cPrevComp(real_t* data) const { return *(data - strideColComp()); }
		template <typename real_t> inline const real_t& cPrevComp(const real_t* data) const { return *(data - strideColComp()); }
		template <typename real_t> inline real_t& cPrevCol(real_t* data) const { return *(data - strideColCell()); }
		template <typename real_t> inline const real_t& cPrevCol(const real_t* data) const { return *(data - strideColCell()); }


		template <typename real_t> inline real_t& cpNextComp(real_t* data) const { return *(data + strideParComp()); }
		template <typename real_t> inline const real_t& cpNextComp(const real_t* data) const { return *(data + strideParComp()); }
		template <typename real_t> inline real_t& cpPrevComp(real_t* data) const { return *(data - strideParComp()); }
		template <typename real_t> inline const real_t& cpPrevComp(const real_t* data) const { return *(data - strideParComp()); }
		template <typename real_t> inline real_t& cpNextShell(real_t* data) const { return *(data + strideParShell()); }
		template <typename real_t> inline const real_t& cpNextShell(const real_t* data) const { return *(data + strideParShell()); }
		template <typename real_t> inline real_t& cpPrevShell(real_t* data) const { return *(data - strideParShell()); }
		template <typename real_t> inline const real_t& cpPrevShell(const real_t* data) const { return *(data - strideParShell()); }


		template <typename real_t> inline real_t& jfNextComp(real_t* data) const { return *(data + strideFluxComp()); }
		template <typename real_t> inline const real_t& jfNextComp(const real_t* data) const { return *(data + strideFluxComp()); }
		template <typename real_t> inline real_t& jfNextCol(real_t* data) const { return *(data + strideFluxCell()); }
		template <typename real_t> inline const real_t& jfNextCol(const real_t* data) const { return *(data + strideFluxCell()); }
		template <typename real_t> inline real_t& jfPrevComp(real_t* data) const { return *(data - strideFluxComp()); }
		template <typename real_t> inline const real_t& jfPrevComp(const real_t* data) const { return *(data - strideFluxComp()); }
		template <typename real_t> inline real_t& jfPrevCol(real_t* data) const { return *(data - strideFluxCell()); }
		template <typename real_t> inline const real_t& jfPrevCol(const real_t* data) const { return *(data - strideFluxCell()); }

	protected:
		const Discretization& _disc;
	};

	class Exporter : public ISolutionExporter
	{
	public:

		Exporter(const Discretization& disc, double const* data) : _disc(disc), _idx(disc), _data(data) { }
		Exporter(const Discretization&& disc, double const* data) = delete;

		virtual bool hasMultipleBoundStates() const CADET_NOEXCEPT { return cadet::model::hasMultipleBoundStates(_disc.nBound, _disc.nComp); }
		virtual bool hasNonBindingComponents() const CADET_NOEXCEPT { return cadet::model::hasNonBindingComponents(_disc.nBound, _disc.nComp); }
		virtual bool hasParticleFlux() const CADET_NOEXCEPT { return true; }
		virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT { return true; }

		virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
		virtual unsigned int numAxialCells() const CADET_NOEXCEPT { return _disc.nCol; }
		virtual unsigned int numRadialCells() const CADET_NOEXCEPT { return _disc.nPar; }
		virtual unsigned int numBoundStates() const CADET_NOEXCEPT { return _disc.strideBound; }
		virtual unsigned int const* numBoundStatesPerComponent() const CADET_NOEXCEPT { return _disc.nBound; }
		virtual unsigned int numBoundStates(unsigned int comp) const CADET_NOEXCEPT { return _disc.nBound[comp]; }
		virtual unsigned int numColumnDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol; }
		virtual unsigned int numParticleDofs() const CADET_NOEXCEPT { return (_disc.nComp + _disc.strideBound) * _disc.nPar * _disc.nCol; }
		virtual unsigned int numFluxDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol; }
		
		virtual double concentration(unsigned int component, unsigned int axialCell) const { return _idx.c(_data, axialCell, component); }
		virtual double flux(unsigned int component, unsigned int axialCell) const { return _idx.jf(_data, axialCell, component); }
		virtual double mobilePhase(unsigned int component, unsigned int axialCell, unsigned int radialCell) const { return _idx.cp(_data, axialCell, radialCell, component); }
		virtual double solidPhase(unsigned int component, unsigned int axialCell, unsigned int radialCell, unsigned int boundState) const
		{
			return _idx.q(_data, axialCell, radialCell, boundState, component);
		}
		
		virtual double const* concentration() const { return _idx.c(_data); }
		virtual double const* flux() const { return _idx.jf(_data); }
		virtual double const* mobilePhase() const { return _idx.cp(_data); }
		virtual double const* solidPhase() const { return _idx.q(_data); }
		virtual double const* inlet(unsigned int& stride) const
		{
			stride = _idx.strideColComp();
			return &_idx.c(_data, 0, 0);
		}
		virtual double const* outlet(unsigned int& stride) const
		{
			stride = _idx.strideColComp();
			return &_idx.c(_data, _disc.nCol - 1, 0);
		}

		virtual StateOrdering const* concentrationOrdering(unsigned int& len) const
		{
			len = _concentrationOrdering.size();
			return _concentrationOrdering.data();
		}

		virtual StateOrdering const* fluxOrdering(unsigned int& len) const
		{
			len = _fluxOrdering.size();
			return _fluxOrdering.data();
		}

		virtual StateOrdering const* mobilePhaseOrdering(unsigned int& len) const
		{
			len = _particleOrdering.size();
			return _particleOrdering.data();
		}

		virtual StateOrdering const* solidPhaseOrdering(unsigned int& len) const
		{
			len = _solidOrdering.size();
			return _solidOrdering.data();
		}

	protected:
		const Discretization& _disc;
		const Indexer _idx;
		double const* const _data;

		const std::array<StateOrdering, 2> _concentrationOrdering = { { StateOrdering::Component, StateOrdering::AxialCell } };
		const std::array<StateOrdering, 4> _particleOrdering = { { StateOrdering::AxialCell, StateOrdering::RadialCell, StateOrdering::Phase, StateOrdering::Component } };
		const std::array<StateOrdering, 2> _solidOrdering = { { StateOrdering::Component, StateOrdering::Phase } };
		const std::array<StateOrdering, 2> _fluxOrdering = { { StateOrdering::Component, StateOrdering::AxialCell } };
	};
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_GENERALRATEMODEL_HPP_

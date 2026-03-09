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

#include "model/parts/ConvectionDispersionOperatorDG.hpp"
#include "cadet/Exceptions.hpp"

#include "ParamReaderHelper.hpp"
#include "AdUtils.hpp"
#include "SimulationTypes.hpp"
//#include "model/parts/AxialConvectionDispersionKernelDG.hpp" // todo radial flow DG and outsource residual implementation to kernel
//#include "model/parts/RadialConvectionDispersionKernelDG.hpp"
#include "model/ParameterDependence.hpp"
#include "SensParamUtil.hpp"
#include "ConfigurationHelper.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>

namespace cadet
{

namespace model
{

namespace parts
{

/**
 * @brief Creates an AxialConvectionDispersionOperatorBaseDG
 */
AxialConvectionDispersionOperatorBaseDG::AxialConvectionDispersionOperatorBaseDG() :
	_dispersionDep(nullptr), _DGjacAxDispBlocks(nullptr), _auxState(nullptr), _subsState(nullptr)
{
}

AxialConvectionDispersionOperatorBaseDG::~AxialConvectionDispersionOperatorBaseDG() CADET_NOEXCEPT
{
	if (_dispersionDep)
		delete _dispersionDep;

	delete[] _DGjacAxDispBlocks;
	delete[] _auxState;
	delete[] _subsState;
}

/**
 * @brief Reads parameters and allocates memory
 * @details Has to be called once before the operator is used.
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [in] helper Configures parameter dependencies
 * @param [in] nComp Number of components
 * @param [in] exact_integration DG volume integral computation
 * @param [in] nElem Number of axial elements
 * @param [in] polyDeg Polynomial degree of DG approach
 * @param [in] strideNode node stride in state vector
 * @return @c true if configuration went fine, @c false otherwise
 */
bool AxialConvectionDispersionOperatorBaseDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, unsigned int nComp, int polynomial_integration_mode, unsigned int nElem, unsigned int polyDeg, unsigned int strideNode)
{
	const bool firstConfigCall = _auxState == nullptr; // used to not multiply allocate memory

	_nComp = nComp;
	_exactInt = static_cast<bool>(polynomial_integration_mode); // only integration mode 0 applies the inexact collocated diagonal LGL mass matrix
	_nElem = nElem;
	_polyDeg = polyDeg;
	_nNodes = _polyDeg + 1u;
	_nPoints = _nNodes * _nElem;
	_strideNode = strideNode;
	_strideElem = _nNodes * strideNode;

	/* Allocate space for DG discretization operations */
	_nodes.resize(_nNodes);
	_nodes.setZero();
	_invWeights.resize(_nNodes);
	_invWeights.setZero();
	_polyDerM.resize(_nNodes, _nNodes);
	_polyDerM.setZero();
	_invMM.resize(_nNodes, _nNodes);
	_invMM.setZero();

	if (firstConfigCall)
		_auxState = new active[_nPoints];
	if (firstConfigCall)
		_subsState = new active[_nPoints];
	for (int i = 0; i < _nPoints; i++) {
		_auxState[i] = 0.0;
		_subsState[i] = 0.0;
	}
	_boundary.setZero();
	_surfaceFlux.resize(_nElem + 1u);
	_surfaceFlux.setZero();

	_newStaticJac = true;

	dgtoolbox::lglNodesWeights(_polyDeg, _nodes, _invWeights, true);
	_invMM = dgtoolbox::invMMatrix(_polyDeg, _nodes);
	_polyDerM = dgtoolbox::derivativeMatrix(_polyDeg, _nodes);

	if(polynomial_integration_mode == 2) // use Gauss quadrature for exact integration
		_invMM = dgtoolbox::gaussQuadratureMMatrix(_nodes, _nNodes).inverse();

	if (paramProvider.exists("COL_DISPERSION_DEP"))
	{
		const std::string paramDepName = paramProvider.getString("COL_DISPERSION_DEP");
		_dispersionDep = helper.createParameterParameterDependence(paramDepName);
		if (!_dispersionDep)
			throw InvalidParameterException("Unknown parameter dependence " + paramDepName + " in COL_DISPERSION_DEP");

		_dispersionDep->configureModelDiscretization(paramProvider);
	}
	else
		_dispersionDep = helper.createParameterParameterDependence("CONSTANT_ONE");

	return true;
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] unitOpIdx Unit operation id of the owning unit operation model
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [out] parameters Map in which local parameters are inserted
 * @return @c true if configuration went fine, @c false otherwise
 */
bool AxialConvectionDispersionOperatorBaseDG::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters)
{
	const bool firstConfigCall = _DGjacAxDispBlocks == nullptr; // used to not multiply allocate memory

	// Read geometry parameters
	_colLength = paramProvider.getDouble("COL_LENGTH");
	_deltaZ = _colLength / _nElem;

	/* compute dispersion jacobian blocks(without parameters except element spacing, i.e. static entries) */
	// we only need unique dispersion blocks, which are given by elements 1, 2, nElem for inexact integration DG and by elements 1, 2, 3, nElem-1, nElem for eaxct integration DG
	if (firstConfigCall)
		_DGjacAxDispBlocks = new Eigen::MatrixXd[(_exactInt ? std::min(_nElem, 5u) : std::min(_nElem, 3u))];
	_DGjacAxDispBlocks[0] = DGjacobianDispBlock(1);
	if (_nElem > 1)
		_DGjacAxDispBlocks[1] = DGjacobianDispBlock(2);
	if (_nElem > 2 && _exactInt)
		_DGjacAxDispBlocks[2] = DGjacobianDispBlock(3);
	else if (_nElem > 2 && !_exactInt)
		_DGjacAxDispBlocks[2] = DGjacobianDispBlock(_nElem);
	if (_exactInt && _nElem > 3)
		_DGjacAxDispBlocks[3] = DGjacobianDispBlock(std::max(4u, _nElem - 1u));
	if (_exactInt && _nElem > 4)
		_DGjacAxDispBlocks[4] = DGjacobianDispBlock(_nElem);

	// note that convection jacobian block is computet in notifyDiscontinuousSectionTransition() since this block needs to be recomputed when flow direction changes

	// Read cross section area or set to -1
	_crossSection = -1.0;
	if (paramProvider.exists("CROSS_SECTION_AREA"))
	{
		_crossSection = paramProvider.getDouble("CROSS_SECTION_AREA");
	}

	// Read section dependent parameters (transport)

	// Read VELOCITY
	_velocity.clear();
	if (paramProvider.exists("VELOCITY"))
	{
		readScalarParameterOrArray(_velocity, paramProvider, "VELOCITY", 1);
	}
	_dir = 1;

	readScalarParameterOrArray(_colDispersion, paramProvider, "COL_DISPERSION", 1);
	if (paramProvider.exists("COL_DISPERSION_MULTIPLEX"))
	{
		const int mode = paramProvider.getInt("COL_DISPERSION_MULTIPLEX");
		if (mode == 0)
			// Comp-indep, sec-indep
			_dispersionCompIndep = true;
		else if (mode == 1)
			// Comp-dep, sec-indep
			_dispersionCompIndep = false;
		else if (mode == 2)
			// Comp-indep, sec-dep
			_dispersionCompIndep = true;
		else if (mode == 3)
			// Comp-dep, sec-dep
			_dispersionCompIndep = false;

		if (!_dispersionCompIndep && (_colDispersion.size() % _nComp != 0))
			throw InvalidParameterException("Number of elements in field COL_DISPERSION is not a positive multiple of NCOMP (" + std::to_string(_nComp) + ")");
		if ((mode == 0) && (_colDispersion.size() != 1))
			throw InvalidParameterException("Number of elements in field COL_DISPERSION inconsistent with COL_DISPERSION_MULTIPLEX (should be 1)");
		if ((mode == 1) && (_colDispersion.size() != _nComp))
			throw InvalidParameterException("Number of elements in field COL_DISPERSION inconsistent with COL_DISPERSION_MULTIPLEX (should be " + std::to_string(_nComp) + ")");
	}
	else
	{
		// Infer component dependence of COL_DISPERSION:
		//   size not divisible by NCOMP -> component independent
		_dispersionCompIndep = ((_colDispersion.size() % _nComp) != 0);
	}

	// Expand _colDispersion to make it component dependent
	if (_dispersionCompIndep)
	{
		std::vector<active> expanded(_colDispersion.size() * _nComp);
		for (std::size_t i = 0; i < _colDispersion.size(); ++i)
			std::fill(expanded.begin() + i * _nComp, expanded.begin() + (i + 1) * _nComp, _colDispersion[i]);

		_colDispersion = std::move(expanded);
	}

	if (_dispersionDep)
	{
		if (!_dispersionDep->configure(paramProvider, unitOpIdx, ParTypeIndep, BoundStateIndep, "COL_DISPERSION_DEP"))
			throw InvalidParameterException("Failed to configure dispersion parameter dependency (COL_DISPERSION_DEP)");
	}

	if (_velocity.empty() && (_crossSection <= 0.0))
	{
		throw InvalidParameterException("At least one of CROSS_SECTION_AREA and VELOCITY has to be set");
	}

	// Add parameters to map
	if (_dispersionCompIndep)
	{
		if (_colDispersion.size() > _nComp)
		{
			// Register only the first item in each section
			for (std::size_t i = 0; i < _colDispersion.size() / _nComp; ++i)
				parameters[makeParamId(hashString("COL_DISPERSION"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, i)] = &_colDispersion[i * _nComp];
		}
		else
		{
			// We have only one parameter
			parameters[makeParamId(hashString("COL_DISPERSION"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colDispersion[0];
		}
	}
	else
		registerParam2DArray(parameters, _colDispersion, [=](bool multi, unsigned int sec, unsigned int comp) { return makeParamId(hashString("COL_DISPERSION"), unitOpIdx, comp, ParTypeIndep, BoundStateIndep, ReactionIndep, multi ? sec : SectionIndep); }, _nComp);

	registerScalarSectionDependentParam(hashString("VELOCITY"), parameters, _velocity, unitOpIdx, ParTypeIndep);
	parameters[makeParamId(hashString("COL_LENGTH"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colLength;
	parameters[makeParamId(hashString("CROSS_SECTION_AREA"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_crossSection;

	return true;
}

/**
 * @brief Notifies the operator that a discontinuous section transition is in progress
 * @details In addition to changing flow direction internally, if necessary, the function returns whether
 *          the flow direction has changed.
 * @param [in] t Current time point
 * @param [in] secIdx Index of the new section that is about to be integrated
 * @return @c true if flow direction has changed, otherwise @c false
 */
bool AxialConvectionDispersionOperatorBaseDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, MatrixXd& jacInlet)
{

	_curSection = secIdx;
	_newStaticJac = true;

	// setFlowRates() was called before, so _curVelocity has direction dirOld
	const int dirOld = _dir;

	if (_crossSection <= 0.0)
	{
		// Use the provided _velocity (direction is also set), only update _dir
		_curVelocity = getSectionDependentScalar(_velocity, secIdx);
		_dir = (_curVelocity >= 0.0) ? 1 : -1;
	}
	else if (!_velocity.empty())
	{
		// Use network flow rate but take direction from _velocity
		_dir = (getSectionDependentScalar(_velocity, secIdx) >= 0.0) ? 1 : -1;

		// _curVelocity has correct magnitude but previous direction, so flip it if necessary
		if (dirOld * _dir < 0)
			_curVelocity *= -1.0;
	}

	// Remaining case: _velocity is empty and _crossSection <= 0.0
	// _curVelocity is goverend by network flow rate provided in setFlowRates().
	// Direction never changes (always forward, that is, _dir = 1)-
	// No action required.

	 // recompute convection jacobian block, which depends on flow direction
	_DGjacAxConvBlock = DGjacobianConvBlock();

	if (_curVelocity >= 0.0) { // forward flow upwind convection
		if (_exactInt)
			jacInlet = static_cast<double>(_curVelocity) * _DGjacAxConvBlock.col(0); // only first element depends on inlet concentration
		else
			jacInlet(0, 0) = static_cast<double>(_curVelocity) * _DGjacAxConvBlock(0, 0); // only first node depends on inlet concentration
	}
	else {  // backward flow upwind convection
		if (_exactInt)
			jacInlet = static_cast<double>(_curVelocity) * _DGjacAxConvBlock.col(_DGjacAxConvBlock.cols() - 1); // only last element depends on inlet concentration
		else
			jacInlet(0, 0) = static_cast<double>(_curVelocity) * _DGjacAxConvBlock(_DGjacAxConvBlock.rows() - 1, _DGjacAxConvBlock.cols() - 1); // only last node depends on inlet concentration
	}
	
	// Detect change in flow direction
	return (dirOld * _dir < 0);
}

/**
 * @brief Sets the flow rates for the current time section
 * @details The flow rates may change due to valve switches.
 * @param [in] in Total volumetric inlet flow rate
 * @param [in] out Total volumetric outlet flow rate
 * @param [in] colPorosity Porosity used for computing interstitial velocity from volumetric flow rate
 */
void AxialConvectionDispersionOperatorBaseDG::setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT
{
	// If we have cross section area, interstitial velocity is given by network flow rates
	if (_crossSection > 0.0)
		_curVelocity = _dir * in / (_crossSection * colPorosity);
}

/**
 * @brief Computes the residual of the transport equations
 * @param [in] model Model that owns the operator
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] y Pointer to unit operation's state vector
 * @param [in] yDot Pointer to unit operation's time derivative state vector
 * @param [out] res Pointer to unit operation's residual vector
 * @param [in] jac Matrix that holds the Jacobian
 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
 */
int AxialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac)
{
	//// Reset Jacobian but keep pattern
	//double* val = jac.valuePtr();
	//for (unsigned int entry = 0; entry < jac.nonZeros(); val++)
	//	val[0] = 0.0;

	linalg::BandedEigenSparseRowIterator jacIt(jac, 0);

	return residualImpl<double, double, double, linalg::BandedEigenSparseRowIterator, true>(model, t, secIdx, y, yDot, res, jacIt);
}

int AxialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, WithoutParamSensitivity)
{
	return residualImpl<double, double, double, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int AxialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithoutParamSensitivity)
{
	return residualImpl<active, active, double, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int AxialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac)
{
	// todo include jacobian in convDisp operator
	//// Reset Jacobian but keep pattern
	//double* val = jac.valuePtr();
	//for (unsigned int entry = 0; entry < jac.nonZeros(); val++)
	//	val[0] = 0.0;

	linalg::BandedEigenSparseRowIterator jacIt(jac, 0);

	return residualImpl<double, active, active, linalg::BandedEigenSparseRowIterator, true>(model, t, secIdx, y, yDot, res, jacIt);
}

int AxialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, WithParamSensitivity)
{
	return residualImpl<double, active, active, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int AxialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithParamSensitivity)
{
	return residualImpl<active, active, active, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

template <typename StateType, typename ResidualType, typename ParamType, typename RowIteratorType, bool wantJac>
int AxialConvectionDispersionOperatorBaseDG::residualImpl(const IModel& model, double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res, RowIteratorType jacBegin)
{
	for (unsigned int comp = 0; comp < _nComp; comp++) {

		// create Eigen objects
		Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> _C(y + offsetC() + comp, _nPoints, InnerStride<Dynamic>(_strideNode));
		Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> _resC(res + offsetC() + comp, _nPoints, InnerStride<Dynamic>(_strideNode));
		Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>> _h(reinterpret_cast<ResidualType*>(_subsState), _nPoints, InnerStride<>(1));
		Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> _g(reinterpret_cast<StateType*>(_auxState), _nPoints, InnerStride<>(1));

		// Add time derivative to bulk residual
		if (yDot)
		{
			Eigen::Map<const VectorXd, 0, InnerStride<Dynamic>> _Cdot(yDot + offsetC() + comp, _nPoints, InnerStride<Dynamic>(_strideNode));
			_resC = _Cdot.template cast<ResidualType>();
		}
		else
			_resC.setZero();

		const ParamType u = static_cast<ParamType>(_curVelocity);
		const ParamType d_ax = static_cast<ParamType>(getSectionDependentSlice(_colDispersion, _nComp, secIdx)[comp]);

		// ===================================//
		// reset cache                        //
		// ===================================//

		_h.setZero();
		_g.setZero();
		_boundary[0] = y[comp]; // copy inlet DOFs to ghost node

		// ======================================//
		// solve auxiliary system g = d c / d x  //
		// ======================================//

		 // DG volume integral in strong form
		volumeIntegral<StateType, StateType>(_C, _g);

		// calculate numerical flux values c*
		InterfaceFluxAuxiliary<StateType>(y + offsetC() + comp, _strideNode, _strideElem);

		// DG surface integral in strong form
		surfaceIntegral<StateType, StateType>(y + offsetC() + comp, &_g[0], _strideNode, _strideElem, 1u, _nNodes);

		// ======================================//
		// solve main equation RHS  d h / d x    //
		// ======================================//

		// calculate the substitute h(g(c), c) and apply inverse mapping jacobian (reference space)
		_h = 2.0 / static_cast<ParamType>(_deltaZ) * (-u * _C + d_ax * (-2.0 / static_cast<ParamType>(_deltaZ)) * _g).template cast<ResidualType>();

		// DG volume integral in strong form
		volumeIntegral<ResidualType, ResidualType>(_h, _resC);

		// update boundary values for auxiliary variable g (solid wall)
		calcBoundaryValues<StateType>();

		// calculate numerical flux values h*
		InterfaceFlux<StateType, ParamType>(y + offsetC() + comp, d_ax);

		// DG surface integral in strong form
		surfaceIntegral<ResidualType, ResidualType>(&_h[0], res + offsetC() + comp,
			1u, _nNodes, _strideNode, _strideElem);
	}

	return 0;
}
/**
* @brief analytically calculates the (static) state jacobian
* @return 1 if jacobain estimation fits the predefined pattern of the jacobian, 0 if not.
*/
int AxialConvectionDispersionOperatorBaseDG::calcTransportJacobian(Eigen::SparseMatrix<double, RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int bulkOffset) {

	// DG convection dispersion Jacobian
	if (_exactInt)
		calcConvDispExIntDGSEMJacobian(jacobian, jacInlet, bulkOffset);
	else
		calcConvDispCollocationDGSEMJacobian(jacobian, jacInlet, bulkOffset);

	if (!jacobian.isCompressed()) // if matrix lost its compressed storage, the pattern did not fit.
		return 0;

	return 1;
}
/**
 * @brief calculates the number of entris for the DG convection dispersion jacobian
 * @note only dispersion entries are relevant for jacobian NNZ as the convection entries are a subset of these
 */
unsigned int AxialConvectionDispersionOperatorBaseDG::nJacEntries(bool pureNNZ) {

	if (_exactInt) {
		if (pureNNZ) {
			return _nComp * ((3u * _nElem - 2u) * _nNodes * _nNodes + (2u * _nElem - 3u) * _nNodes); // dispersion entries
		}
		return _nComp * _nNodes * _nNodes + _nNodes // convection entries
			+ _nComp * ((3u * _nElem - 2u) * _nNodes * _nNodes + (2u * _nElem - 3u) * _nNodes); // dispersion entries
	}
	else {
		if (pureNNZ) {
			return _nComp * (_nElem * _nNodes * _nNodes + 8u * _nNodes); // dispersion entries
		}
		return _nComp * _nNodes * _nNodes + 1u // convection entries
			+ _nComp * (_nElem * _nNodes * _nNodes + 8u * _nNodes); // dispersion entries
	}
}
void model::parts::AxialConvectionDispersionOperatorBaseDG::convDispJacPattern(std::vector<T>& tripletList, const int bulkOffset)
{
	if (_exactInt)
		ConvDispExIntPattern(tripletList, bulkOffset);
	else
		ConvDispCollocationPattern(tripletList, bulkOffset);
}
/**
 * @brief Multiplies the time derivative Jacobian @f$ \frac{\partial F}{\partial \dot{y}}\left(t, y, \dot{y}\right) @f$ with a given vector
 * @details The operation @f$ z = \frac{\partial F}{\partial \dot{y}} x @f$ is performed.
 *          The matrix-vector multiplication is performed matrix-free (i.e., no matrix is explicitly formed).
 *          
 *          Note that this function only performs multiplication with the Jacobian of the (axial) transport equations.
 *          The input vectors are assumed to point to the beginning (including inlet DOFs) of the respective unit operation's arrays.
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in] sDot Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void AxialConvectionDispersionOperatorBaseDG::multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const
{
	double* localRet = ret + offsetC();
	double const* localSdot = sDot + offsetC();
	const int gapelement = strideColNode() - static_cast<int>(_nComp) * strideColComp();

	for (unsigned int i = 0; i < _nPoints; ++i, localRet += gapelement, localSdot += gapelement)
	{
		for (unsigned int j = 0; j < _nComp; ++j, ++localRet, ++localSdot)
		{
			*localRet = (*localSdot);
		}
	}
}

/**
 * @brief Adds the derivatives with respect to @f$ \dot{y} @f$ of @f$ F(t, y, \dot{y}) @f$ to the Jacobian
 * @details This functions computes 
 *          @f[ \begin{align*} \text{_jacCdisc} = \text{_jacCdisc} + \alpha \frac{\partial F}{\partial \dot{y}}. \end{align*} @f]
 *          The factor @f$ \alpha @f$ is useful when constructing the linear system in the time integration process.
 * @param [in] alpha Factor in front of @f$ \frac{\partial F}{\partial \dot{y}} @f$
 */
void AxialConvectionDispersionOperatorBaseDG::addTimeDerivativeToJacobian(double alpha, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacDisc, unsigned int blockOffset)
{
	const int gapelement = strideColNode() - static_cast<int>(_nComp) * strideColComp();
	linalg::BandedEigenSparseRowIterator jac(jacDisc, blockOffset);

	for (unsigned int point = 0; point < _nPoints; ++point, jac+=gapelement) {
		for (unsigned int comp = 0; comp < _nComp; ++comp, ++jac) {
			// dc_b / dt in transport equation
			jac[0] += alpha;
		}
	}
}

unsigned int AxialConvectionDispersionOperatorBaseDG::jacobianLowerBandwidth() const CADET_NOEXCEPT
{
	// @todo use more efficient seed vectors. currently, we treat the jacobian as banded, but the pattern is actually more sparse when multiple components are considered
	// (note that active type directions are limited)
	// We have different jacobian structure for exact integration and collocation DG scheme, i.e. we need different seed vectors
	// collocation DG: 2 * N_n * (N_c + N_q) + 1 = total bandwidth (main diagonal entries maximally depend on the next and last N_n liquid phase entries of same component)
	//    ex. int. DG: 4 * N_n * (N_c + N_q) + 1 = total bandwidth (main diagonal entries maximally depend on the next and last 2*N_n liquid phase entries of same component)

	return (_exactInt) ? 2 * _nNodes * strideColNode() : _nNodes * strideColNode();
}

unsigned int AxialConvectionDispersionOperatorBaseDG::jacobianUpperBandwidth() const CADET_NOEXCEPT
{
	return jacobianLowerBandwidth();
}

bool AxialConvectionDispersionOperatorBaseDG::setParameter(const ParameterId& pId, double value)
{
	// Check if parameter is in parameter dependence of column dispersion coefficient
	if (_dispersionDep)
	{
		if (_dispersionDep->hasParameter(pId))
		{
			_dispersionDep->setParameter(pId, value);
			return true;
		}
	}

	// We only need to do something if COL_DISPERSION is component independent
	if (!_dispersionCompIndep)
		return false;

	if ((pId.name != hashString("COL_DISPERSION")) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
		return false;

	if (_colDispersion.size() > _nComp)
	{
		// Section dependent
		if (pId.section == SectionIndep)
			return false;

		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[pId.section * _nComp + i].setValue(value);
	}
	else
	{
		// Section independent
		if (pId.section != SectionIndep)
			return false;

		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[i].setValue(value);
	}

	return true;
}

bool AxialConvectionDispersionOperatorBaseDG::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
{
	// Check if parameter is in parameter dependence of column dispersion coefficient
	if (_dispersionDep)
	{
		active* const param = _dispersionDep->getParameter(pId);
		if (param)
		{
			param->setValue(value);
			return true;
		}
	}

	// We only need to do something if COL_DISPERSION is component independent
	if (!_dispersionCompIndep)
		return false;

	if ((pId.name != hashString("COL_DISPERSION")) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
		return false;

	if (_colDispersion.size() > _nComp)
	{
		// Section dependent
		if (pId.section == SectionIndep)
			return false;

		if (!contains(sensParams, &_colDispersion[pId.section * _nComp]))
			return false;

		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[pId.section * _nComp + i].setValue(value);
	}
	else
	{
		// Section independent
		if (pId.section != SectionIndep)
			return false;

		if (!contains(sensParams, &_colDispersion[0]))
			return false;

		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[i].setValue(value);
	}

	return true;
}

bool AxialConvectionDispersionOperatorBaseDG::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
{
	// Check if parameter is in parameter dependence of column dispersion coefficient
	if (_dispersionDep)
	{
		active* const param = _dispersionDep->getParameter(pId);
		if (param)
		{
			param->setADValue(adDirection, adValue);
			return true;
		}
	}

	// We only need to do something if COL_DISPERSION is component independent
	if (!_dispersionCompIndep)
		return false;

	if ((pId.name != hashString("COL_DISPERSION")) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
		return false;

	if (_colDispersion.size() > _nComp)
	{
		// Section dependent
		if (pId.section == SectionIndep)
			return false;

		sensParams.insert(&_colDispersion[pId.section * _nComp]);
		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[pId.section * _nComp + i].setADValue(adDirection, adValue);
	}
	else
	{
		// Section independent
		if (pId.section != SectionIndep)
			return false;

		sensParams.insert(&_colDispersion[0]);
		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[i].setADValue(adDirection, adValue);
	}

	return true;
}

/**
 * @brief Creates an ConvectionDispersionOperatorDG
 */
template <typename Operator>
ConvectionDispersionOperatorDG<Operator>::ConvectionDispersionOperatorDG()
{
}

template <typename Operator>
ConvectionDispersionOperatorDG<Operator>::~ConvectionDispersionOperatorDG() CADET_NOEXCEPT
{
}

/**
 * @brief Returns the number of AD directions required for computing the Jacobian
 * @details Band compression is used to minimize the amount of AD directions.
 * @return Number of required AD directions
 */
template <typename Operator>
int ConvectionDispersionOperatorDG<Operator>::requiredADdirs() const CADET_NOEXCEPT
{
	return _baseOp.requiredADdirs();
}

/**
 * @brief Reads parameters and allocates memory
 * @details Has to be called once before the operator is used.
 * @param [in] unitOpIdx Unit operation id of the owning unit operation model
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [out] parameters Map in which local parameters are inserted
 * @param [in] nComp Number of components
 * @param [in] nElem Number of axial elements
 * @return @c true if configuration went fine, @c false otherwise
 */
template <typename Operator>
bool ConvectionDispersionOperatorDG<Operator>::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, unsigned int nComp, int polynomial_integration_mode, unsigned int nElem, unsigned int polyDeg, unsigned int strideNode)
{
	const bool retVal = _baseOp.configureModelDiscretization(paramProvider, helper, nComp, polynomial_integration_mode, nElem, polyDeg, strideNode);

	// todo: manage jacobians in convDispOp instead of unitOp ?
	//// Allocate memory
	//if (_disc.exactInt)
	//	_jacInlet.resize(_disc.nNodes, 1); // first element depends on inlet concentration (same for every component)
	//else
	//	_jacInlet.resize(1, 1); // first element depends on inlet concentration (same for every component)
	//_jac.resize((_disc.nComp + _disc.strideBound) * _disc.nPoints, (_disc.nComp + _disc.strideBound) * _disc.nPoints);
	//_jacDisc.resize((_disc.nComp + _disc.strideBound) * _disc.nPoints, (_disc.nComp + _disc.strideBound) * _disc.nPoints);

	return retVal;
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] unitOpIdx Unit operation id of the owning unit operation model
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [out] parameters Map in which local parameters are inserted
 * @return @c true if configuration went fine, @c false otherwise
 */
template <typename Operator>
bool ConvectionDispersionOperatorDG<Operator>::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters)
{
	return _baseOp.configure(unitOpIdx, paramProvider, parameters);
}

/**
 * @brief Notifies the operator that a discontinuous section transition is in progress
 * @details In addition to changing flow direction internally, if necessary, the function returns whether
 *          the flow direction has changed.
 * @param [in] t Current time point
 * @param [in] secIdx Index of the new section that is about to be integrated
 * @param [in,out] adJac Jacobian information for AD (AD vectors for residual and state, direction offset)
 * @return @c true if flow direction has changed, otherwise @c false
 */
template <typename Operator>
bool ConvectionDispersionOperatorDG<Operator>::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const AdJacobianParams& adJac)
{
	if (_baseOp.exactInt())
		_jacInlet.resize(_baseOp.nNodes(), 1); // first element depends on inlet concentration (same for every component)
	else
		_jacInlet.resize(1, 1); // first element depends on inlet concentration (same for every component)

	const bool hasChanged = _baseOp.notifyDiscontinuousSectionTransition(t, secIdx, _jacInlet);

	// Check whether flow direction has changed and we need to update AD vectors
	// Exit if we do not need to setup (secIdx == 0) or change (prevU and u differ in sign) matrices
	if ((secIdx != 0) && !hasChanged)
		return false;

	// todo: manage jacobians in convDispOp instead of unitOp ?
	//const unsigned int lb = _baseOp.jacobianLowerBandwidth();
	//const unsigned int ub = _baseOp.jacobianUpperBandwidth();
	//if (_baseOp.forwardFlow())
	//{
	//	// Forwards flow

	//	// Repartition column bulk Jacobians
	//	_jacC.repartition(lb, ub);
	//	_jacCdisc.repartition(lb, ub);
	//}
	//else
	//{
	//	// Backwards flow

	//	// Repartition column bulk Jacobians
	//	_jacC.repartition(ub, lb);
	//	_jacCdisc.repartition(ub, lb);
	//}

	// Update AD seed vectors since Jacobian structure has changed (bulk block bandwidths)
	prepareADvectors(adJac);

	return true;
}

/**
 * @brief Sets the AD seed vectors for the bulk transport variables
 * @param [in,out] adJac Jacobian information for AD (AD vectors for residual and state, direction offset)
 */
template <typename Operator>
void ConvectionDispersionOperatorDG<Operator>::prepareADvectors(const AdJacobianParams& adJac) const
{
	// Early out if AD is disabled
	if (!adJac.adY)
		return;

	// Get bandwidths of blocks
	const unsigned int lowerColBandwidth = _baseOp.exactInt() ? 2 * _baseOp.nNodes() * _baseOp.strideColNode() : _baseOp.nNodes() * _baseOp.strideColNode();
	const unsigned int upperColBandwidth = lowerColBandwidth;

	// Column block
	ad::prepareAdVectorSeedsForBandMatrix(adJac.adY + _baseOp.nComp(), adJac.adDirOffset, _baseOp.nComp() * _baseOp.nPoints(), lowerColBandwidth, upperColBandwidth, lowerColBandwidth);
}

/**
 * @brief Sets the flow rates for the current time section
 * @details The flow rates may change due to valve switches.
 * @param [in] in Total volumetric inlet flow rate
 * @param [in] out Total volumetric outlet flow rate
 * @param [in] colPorosity Porosity used for computing interstitial velocity from volumetric flow rate
 */
template <typename Operator>
void ConvectionDispersionOperatorDG<Operator>::setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT
{
	_baseOp.setFlowRates(in, out, colPorosity);
}

/**
 * @brief Computes the residual of the transport equations
 * @param [in] model Model that owns the operator
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] y Pointer to unit operation's state vector
 * @param [in] yDot Pointer to unit operation's time derivative state vector
 * @param [out] res Pointer to unit operation's residual vector
 * @param [in] wantJac Determines whether analytic Jacobian is computed
 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
 */
template <typename Operator>
int ConvectionDispersionOperatorDG<Operator>::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, bool wantJac, WithoutParamSensitivity)
{
	if (wantJac)
		return _baseOp.residual(model, t, secIdx, y, yDot, res, _jacC);
	else
		return _baseOp.residual(model, t, secIdx, y, yDot, res, WithoutParamSensitivity());
}

template <typename Operator>
int ConvectionDispersionOperatorDG<Operator>::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, bool wantJac, WithoutParamSensitivity)
{
	return _baseOp.residual(model, t, secIdx, y, yDot, res, WithoutParamSensitivity());
}

template <typename Operator>
int ConvectionDispersionOperatorDG<Operator>::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, bool wantJac, WithParamSensitivity)
{
	return _baseOp.residual(model, t, secIdx, y, yDot, res, WithParamSensitivity());
}

template <typename Operator>
int ConvectionDispersionOperatorDG<Operator>::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, bool wantJac, WithParamSensitivity)
{
	if (wantJac)
		return _baseOp.residual(model, t, secIdx, y, yDot, res, _jacC);
	else
		return _baseOp.residual(model, t, secIdx, y, yDot, res, WithParamSensitivity());
}

/**
 * @brief Multiplies the time derivative Jacobian @f$ \frac{\partial F}{\partial \dot{y}}\left(t, y, \dot{y}\right) @f$ with a given vector
 * @details The operation @f$ z = \frac{\partial F}{\partial \dot{y}} x @f$ is performed.
 *          The matrix-vector multiplication is performed matrix-free (i.e., no matrix is explicitly formed).
 *          
 *          Note that this function only performs multiplication with the Jacobian of the (axial) transport equations.
 *          The input vectors are assumed to point to the beginning (including inlet DOFs) of the respective unit operation's arrays.
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in] sDot Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [out] ret Vector @f$ z @f$ which stores the result of the operation
 */
template <typename Operator>
void ConvectionDispersionOperatorDG<Operator>::multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const
{
	_baseOp.multiplyWithDerivativeJacobian(simTime, sDot, ret);
}

/**
 * @brief Extracts the system Jacobian from band compressed AD seed vectors
 * @details The input vectors are assumed to point to the beginning (including inlet DOFs) of the respective unit operation's arrays.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
template <typename Operator>
void ConvectionDispersionOperatorDG<Operator>::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	const active* const adVec = adRes + offsetC();

	const int lowerBandwidth = _baseOp.jacobianLowerBandwidth();
	const int upperBandwidth = lowerBandwidth;
	const int stride = lowerBandwidth + 1 + upperBandwidth;

	int diagDir = lowerBandwidth;
	for (int eq = 0; eq < _jacC.rows(); ++eq)
	{
		// Start with lowest subdiagonal and stay in the range of the columns:
		// diagDir might be too big for the matrix and, hence, dir ranges between
		// diagDir - lowerBandwidth <= dir <= diagDir + upperBandwidth
		int dir = diagDir - lowerBandwidth + eq % stride;

		// Loop over diagonals
		for (int diag = 0; diag < stride; ++diag)
		{
			if (eq - lowerBandwidth + diag >= 0 && // left boundary
				eq - lowerBandwidth + diag < _jacC.cols() && // right boundary
				adVec[eq].getADValue(adDirOffset + dir) != 0.0 // keep pattern
				)
				_jacC.coeffRef(eq, eq - lowerBandwidth + diag) = adVec[eq].getADValue(adDirOffset + dir);

			// Wrap around at end of row and jump to lowest subdiagonal
			if (dir == diagDir + upperBandwidth)
				dir = diagDir - lowerBandwidth;
			else
				++dir;
		}
	}
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the corresponding band matrices.
 *          The input vectors are assumed to point to the beginning (including inlet DOFs) of the respective unit operation's arrays.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 * @return Maximum elementwise absolute difference between analytic and AD Jacobian
 */
template <typename Operator>
double ConvectionDispersionOperatorDG<Operator>::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	const int lowerBandwidth = _baseOp.jacobianLowerBandwidth();
	const int upperBandwidth = lowerBandwidth;
	const int diagDir = lowerBandwidth;

	// Column block
	const double maxDiffCol = ad::compareBandedEigenJacobianWithAd(adRes + offsetC(), adDirOffset, diagDir, lowerBandwidth, upperBandwidth, 0, _jacC.rows(), _jacC, 0);
	LOG(Debug) << "-> Col block diff: " << maxDiffCol;

	return maxDiffCol;
}

#endif

/**
 * @brief Assembles the axial transport Jacobian @f$ J_0 @f$ of the time-discretized equations
 * @details The system \f[ \left( \frac{\partial F}{\partial y} + \alpha \frac{\partial F}{\partial \dot{y}} \right) x = b \f]
 *          has to be solved. The system Jacobian of the original equations,
 *          \f[ \frac{\partial F}{\partial y}, \f]
 *          is already computed (by AD or manually in residualImpl() with @c wantJac = true). This function is responsible
 *          for adding
 *          \f[ \alpha \frac{\partial F}{\partial \dot{y}} \f]
 *          to the system Jacobian, which yields the Jacobian of the time-discretized equations
 *          \f[ F\left(t, y_0, \sum_{k=0}^N \alpha_k y_k \right) = 0 \f]
 *          when a BDF method is used. The time integrator needs to solve this equation for @f$ y_0 @f$, which requires
 *          the solution of the linear system mentioned above (@f$ \alpha_0 = \alpha @f$ given in @p alpha).
 *
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 */
template <typename Operator>
void ConvectionDispersionOperatorDG<Operator>::assembleDiscretizedJacobian(double alpha)
{
	// set to static jacobian entries
	_jacCdisc = _jacC;

	// Add time derivatives
	addTimeDerivativeToJacobian(alpha);
}


/**
 * @brief Adds the derivatives with respect to @f$ \dot{y} @f$ of @f$ F(t, y, \dot{y}) @f$ to the Jacobian
 * @details This functions computes 
 *          @f[ \begin{align*} \text{_jacCdisc} = \text{_jacCdisc} + \alpha \frac{\partial F}{\partial \dot{y}}. \end{align*} @f]
 *          The factor @f$ \alpha @f$ is useful when constructing the linear system in the time integration process.
 * @param [in] alpha Factor in front of @f$ \frac{\partial F}{\partial \dot{y}} @f$
 */
template <typename Operator>
void ConvectionDispersionOperatorDG<Operator>::addTimeDerivativeToJacobian(double alpha)
{
	_baseOp.addTimeDerivativeToJacobian(alpha, _jacCdisc);
}

/**
 * @brief Assembles and factorizes the time discretized Jacobian
 * @details See assembleDiscretizedJacobian() for assembly of the time discretized Jacobian.
 * @param [in] alpha Factor in front of @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @return @c true if factorization went fine, otherwise @c false
 */
template <typename Operator>
bool ConvectionDispersionOperatorDG<Operator>::assembleAndFactorizeDiscretizedJacobian(double alpha)
{
	// todo: this functionality is currently not needed for DG since we assemble and factorize a global jacobian, not blocks (as in FV)
	//assembleDiscretizedJacobian(alpha);
	//return _jacCdisc.factorize();
	return true;
}

/**
 * @brief Solves a (previously factorized) equation system
 * @details The (time discretized) Jacobian matrix has to be factorized before calling this function.
 *          Note that the given right hand side vector @p rhs is not shifted by the inlet DOFs. That
 *          is, it is assumed to point directly to the first axial DOF.
 * 
 * @param [in,out] rhs On entry, right hand side of the equation system. On exit, solution of the system.
 * @return @c true if the system was solved correctly, otherwise @c false
 */
template <typename Operator>
bool ConvectionDispersionOperatorDG<Operator>::solveDiscretizedJacobian(double* rhs) const
{
	// todo: this functionality is currently not needed for DG since we assemble and factorize a global jacobian, not blocks (as in FV)
	return true;// _jacCdisc.solve(rhs);
}

/**
 * @brief Solves a system with the time derivative Jacobian and given right hand side
 * @details Note that the given right hand side vector @p rhs is not shifted by the inlet DOFs. That
 *          is, it is assumed to point directly to the first axial DOF.
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in,out] rhs On entry, right hand side. On exit, solution of the system.
 * @return @c true if the system was solved correctly, @c false otherwise
 */
template <typename Operator>
bool ConvectionDispersionOperatorDG<Operator>::solveTimeDerivativeSystem(const SimulationTime& simTime, double* const rhs)
{
	// Assemble
	double* vals = _jacCdisc.valuePtr();
	for (int entry = 0; entry < _jacCdisc.nonZeros(); entry++)
		vals[entry] = 0.0;

	addTimeDerivativeToJacobian(1.0);

	Eigen::SparseLU<Eigen::SparseMatrix<double>> _linSolver;
	_linSolver.analyzePattern(_jacCdisc);
	_linSolver.factorize(_jacCdisc);

	if (_linSolver.info() != Success) {
		LOG(Error) << "factorization failed in sensitivity initialization";
	}

	Eigen::Map<Eigen::VectorXd> ret_vec(rhs, _jacCdisc.rows());
	ret_vec = _linSolver.solve(ret_vec);

	// Use the factors to solve the linear system 
	if (_linSolver.info() != Success) {
		LOG(Error) << "solve failed in sensitivity initialization";
	}

	return true;
}

// Template instantiations
template class ConvectionDispersionOperatorDG<AxialConvectionDispersionOperatorBaseDG>;
template class ConvectionDispersionOperatorDG<RadialConvectionDispersionOperatorBaseDG>;

/* ====================================================================================
 *  RADIAL CONVECTION DISPERSION OPERATOR DG IMPLEMENTATION
 * ==================================================================================== */

/**
 * @brief Creates a RadialConvectionDispersionOperatorBaseDG
 */
RadialConvectionDispersionOperatorBaseDG::RadialConvectionDispersionOperatorBaseDG() :
	_dispersionDep(nullptr), _auxState(nullptr), _subsState(nullptr),
	_variableDispersion(false), _overintegrate(false)
{
}

RadialConvectionDispersionOperatorBaseDG::~RadialConvectionDispersionOperatorBaseDG() CADET_NOEXCEPT
{
	if (_dispersionDep)
		delete _dispersionDep;

	delete[] _auxState;
	delete[] _subsState;
}

/**
 * @brief Reads parameters and allocates memory
 */
bool RadialConvectionDispersionOperatorBaseDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, unsigned int nComp, int polynomial_integration_mode, unsigned int nElem, unsigned int polyDeg, unsigned int strideNode)
{
	const bool firstConfigCall = _auxState == nullptr;

	_nComp = nComp;
	// Note: Radial DG uses exact integration only (polynomial_integration_mode ignored)
	_nElem = nElem;
	_polyDeg = polyDeg;
	_nNodes = _polyDeg + 1u;
	_nPoints = _nNodes * _nElem;
	_strideNode = strideNode;
	_strideElem = _nNodes * strideNode;

	// Allocate space for DG matrices
	_nodes.resize(_nNodes);
	_nodes.setZero();
	_invWeights.resize(_nNodes);
	_invWeights.setZero();
	_polyDerM.resize(_nNodes, _nNodes);
	_polyDerM.setZero();
	_invMM.resize(_nNodes, _nNodes);
	_invMM.setZero();
	_M00.resize(_nNodes, _nNodes);
	_M00.setZero();
	_M01.resize(_nNodes, _nNodes);
	_M01.setZero();

	// Auxiliary state arrays
	if (firstConfigCall)
		_auxState = new active[_nPoints];
	if (firstConfigCall)
		_subsState = new active[_nPoints];
	for (unsigned int i = 0; i < _nPoints; i++) {
		_auxState[i] = 0.0;
		_subsState[i] = 0.0;
	}
	_boundary.setZero();
	_surfaceFlux.resize(_nElem + 1u);
	_surfaceFlux.setZero();

	_newStaticJac = true;

	// Initialize standard DG matrices
	dgtoolbox::lglNodesWeights(_polyDeg, _nodes, _invWeights, true);
	_invMM = dgtoolbox::invMMatrix(_polyDeg, _nodes);
	_polyDerM = dgtoolbox::derivativeMatrix(_polyDeg, _nodes);

	// Standard and weighted mass matrices for radial integrals
	// M^{(0,0)} = integral of l_i * l_j (standard Lagrange mass matrix)
	// M^{(0,1)} = integral of l_i * l_j * (1+xi) (weighted Lagrange mass matrix)
	_M00 = dgtoolbox::mMatrix(_polyDeg, _nodes, 0.0, 0.0);
	_M01 = dgtoolbox::weightedMMatrix(_polyDeg, _nodes);

	if (paramProvider.exists("COL_DISPERSION_DEP"))
	{
		const std::string paramDepName = paramProvider.getString("COL_DISPERSION_DEP");
		_dispersionDep = helper.createParameterParameterDependence(paramDepName);
		if (!_dispersionDep)
			throw InvalidParameterException("Unknown parameter dependence " + paramDepName + " in COL_DISPERSION_DEP");

		_dispersionDep->configureModelDiscretization(paramProvider);

		// Check if dispersion dependence is non-trivial (not CONSTANT_ONE or IDENTITY)
		_variableDispersion = (paramDepName != "CONSTANT_ONE" && paramDepName != "IDENTITY" && paramDepName != "NONE");
	}
	else
	{
		_dispersionDep = helper.createParameterParameterDependence("CONSTANT_ONE");
		_variableDispersion = false;
	}

	// Read overintegration option (for nonlinear variable coefficients)
	_overintegrate = paramProvider.exists("OVERINTEGRATION") ? paramProvider.getBool("OVERINTEGRATION") : false;

	// Allocate storage for variable dispersion values
	_dispersionAtNodes.resize(_nElem);
	for (unsigned int elem = 0; elem < _nElem; ++elem)
		_dispersionAtNodes[elem].resize(_nNodes);
	_dispersionAtInterfaces.resize(_nElem + 1);

	return true;
}

/**
 * @brief Reads model parameters
 */
bool RadialConvectionDispersionOperatorBaseDG::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters)
{
	// Read geometry parameters
	_innerRadius = paramProvider.getDouble("COL_RADIUS_INNER");
	_outerRadius = paramProvider.getDouble("COL_RADIUS_OUTER");
	_deltaRho = (_outerRadius - _innerRadius) / static_cast<double>(_nElem);

	// COL_LENGTH is optional for radial models
	// For pure radial flow, use r_out - r_in as the characteristic length
	if (paramProvider.exists("COL_LENGTH"))
		_colLength = paramProvider.getDouble("COL_LENGTH");
	else
		_colLength = _outerRadius - _innerRadius;

	// Compute radial node coordinates
	computeRadialNodeCoordinates();

	// Compute cell-dependent matrices
	computeCellDependentMatrices();

	// Read cross section area (not directly used in radial but needed for interface)
	// In radial flow, cross section varies with radius

	// Allocate Jacobian blocks (cell-dependent for radial)
	_DGjacRadDispBlocks.resize(_nElem);
	_DGjacRadConvBlocks.resize(_nElem);

	// Read section-dependent velocity coefficient
	// For radial flow: v(r) = VELOCITY_COEFF / r
	_velocity.clear();
	if (paramProvider.exists("VELOCITY_COEFF"))
	{
		readScalarParameterOrArray(_velocity, paramProvider, "VELOCITY_COEFF", 1);
	}

	_dir = 1;

	// Register velocity coefficient
	registerScalarSectionDependentParam(hashString("VELOCITY_COEFF"), parameters, _velocity, unitOpIdx, ParTypeIndep);

	// Read dispersion coefficient
	readScalarParameterOrArray(_colDispersion, paramProvider, "COL_DISPERSION", 1);
	if (paramProvider.exists("COL_DISPERSION_MULTIPLEX"))
	{
		const int mode = paramProvider.getInt("COL_DISPERSION_MULTIPLEX");
		if (mode == 0)
			_dispersionCompIndep = true;
		else if (mode == 1)
			_dispersionCompIndep = false;
		else if (mode == 2)
			_dispersionCompIndep = true;
		else if (mode == 3)
			_dispersionCompIndep = false;
	}
	else
	{
		if (_colDispersion.size() == 1)
			_dispersionCompIndep = true;
		else
			_dispersionCompIndep = false;
	}

	if (_dispersionCompIndep)
	{
		if (_colDispersion.size() > 1)
		{
			// Section-dependent, component-independent
			for (unsigned int i = 0; i < _colDispersion.size(); ++i)
				registerParam1DArray(parameters, _colDispersion, [=](bool multi, unsigned int sec) { return makeParamId(hashString("COL_DISPERSION"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, (multi ? sec : SectionIndep)); });
		}
		else
		{
			// Section-independent, component-independent
			parameters[makeParamId(hashString("COL_DISPERSION"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colDispersion[0];
		}
	}
	else
		registerParam2DArray(parameters, _colDispersion, [=](bool multi, unsigned int comp, unsigned int sec) { return makeParamId(hashString("COL_DISPERSION"), unitOpIdx, comp, ParTypeIndep, BoundStateIndep, ReactionIndep, (multi ? sec : SectionIndep)); }, _nComp);

	// Configure dispersion dependence
	if (_dispersionDep)
	{
		if (!_dispersionDep->configure(paramProvider, unitOpIdx, ParTypeIndep, BoundStateIndep, "COL_DISPERSION_DEP"))
			throw InvalidParameterException("Failed to configure COL_DISPERSION_DEP");
	}

	return true;
}

/**
 * @brief Notifies the operator that a discontinuous section transition is in progress
 */
bool RadialConvectionDispersionOperatorBaseDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, Eigen::MatrixXd& jacInlet)
{
	_curSection = secIdx;

	// Set current velocity coefficient (u = Q / (2 * pi * L * epsilon))
	// For radial flow, this is later divided by radius to get local velocity
	if (_velocity.size() == 1)
		_curVelocity = _dir * _velocity[0];
	else
		_curVelocity = _dir * _velocity[secIdx];

	_newStaticJac = true;

	// u = velocity_coeff (NOT multiplied by rho, matches residualImpl)
	const double u = static_cast<double>(_curVelocity);

	for (unsigned int cell = 0; cell < _nElem; cell++)
	{
		_DGjacRadConvBlocks[cell] = u * DGjacobianConvBlockRadial(cell);
	}

	// Inlet Jacobian: jacInlet should have size (nNodes, 1) for exact integration
	// Similar to axial DG: jacInlet = velocity * convBlock.col(0) for forward flow
	jacInlet.resize(_nNodes, 1);

	if (_curVelocity >= 0.0)
	{
		// Forward flow (outward radial): inlet at inner radius, first element
		jacInlet = _DGjacRadConvBlocks[0].col(0);
	}
	else
	{
		// Backward flow (inward radial): inlet at outer radius, last element
		jacInlet = _DGjacRadConvBlocks[_nElem - 1].col(_DGjacRadConvBlocks[_nElem - 1].cols() - 1);
	}

	return true;
}

/**
 * @brief Computes radial node coordinates in physical space
 */
void RadialConvectionDispersionOperatorBaseDG::computeRadialNodeCoordinates()
{
	_rhoNodeCoords.resize(_nPoints);
	_rhoCellBounds.resize(_nElem + 1);

	const double deltaRho = static_cast<double>(_deltaRho);
	const double innerR = static_cast<double>(_innerRadius);

	// Compute cell boundaries
	for (unsigned int elem = 0; elem <= _nElem; ++elem)
		_rhoCellBounds[elem] = innerR + elem * deltaRho;

	// Compute node coordinates in physical space
	for (unsigned int elem = 0; elem < _nElem; ++elem)
	{
		const double rho_left = _rhoCellBounds[elem];
		for (unsigned int node = 0; node < _nNodes; ++node)
		{
			// Map reference node xi in [-1,1] to physical radial coordinate
			_rhoNodeCoords[elem * _nNodes + node] = rho_left + 0.5 * deltaRho * (1.0 + _nodes[node]);
		}
	}
}

/**
 * @brief Computes cell-dependent weighted mass matrices and stiffness matrices
 * @param [in] dispersion base dispersion coefficient (unused, for interface compatibility)
 */
void RadialConvectionDispersionOperatorBaseDG::computeCellDependentMatrices(double dispersion)
{
	_invMM_rho.resize(_nElem);
	_S_g.resize(_nElem);

	const double deltaRho = static_cast<double>(_deltaRho);

	for (unsigned int elem = 0; elem < _nElem; ++elem)
	{
		// Cell left boundary radial coordinate
		const double rho_i = _rhoCellBounds[elem];

		// Weighted mass matrix: M_rho = (deltaRho/2) * M^{(0,1)} + rho_i * M^{(0,0)}
		Eigen::MatrixXd M_rho = (deltaRho / 2.0) * _M01 + rho_i * _M00;
		_invMM_rho[elem] = M_rho.inverse();

		// Stiffness matrix for auxiliary equation: S_g = D^T * M_rho
		// Used in volume integral: dc/dt = ... - M_rho^{-1} * S_g * g
		// For constant dispersion, this is just D^T * M_rho (dispersion multiplied later)
		_S_g[elem] = _polyDerM.transpose() * M_rho;
	}
}

/**
 * @brief Updates dispersion values at nodes and interfaces for variable D(rho)
 * @detail Evaluates _dispersionDep at each DG node and cell interface
 */
void RadialConvectionDispersionOperatorBaseDG::updateDispersionValues(const IModel& model, unsigned int secIdx, unsigned int comp)
{
	if (!_variableDispersion)
		return;

	const active* const d_rad = currentDispersion(secIdx);
	const double baseDispersion = static_cast<double>(d_rad[comp]);

	// Evaluate dispersion at each node
	for (unsigned int elem = 0; elem < _nElem; ++elem)
	{
		for (unsigned int node = 0; node < _nNodes; ++node)
		{
			const double rho = _rhoNodeCoords[elem * _nNodes + node];
			// Normalize position to [0, 1] for parameter dependence evaluation
			const double relPos = (rho - static_cast<double>(_innerRadius)) /
			                      (static_cast<double>(_outerRadius) - static_cast<double>(_innerRadius));

			// Evaluate D(rho) = baseDispersion * dependence_factor(rho)
			const double depValue = _dispersionDep->getValue(model, ColumnPosition{relPos, 0.0, 0.0},
			                                                  comp, ParTypeIndep, BoundStateIndep, baseDispersion);
			_dispersionAtNodes[elem][node] = depValue;
		}
	}

	// Evaluate dispersion at interfaces
	for (unsigned int iface = 0; iface <= _nElem; ++iface)
	{
		const double rho = _rhoCellBounds[iface];
		const double relPos = (rho - static_cast<double>(_innerRadius)) /
		                      (static_cast<double>(_outerRadius) - static_cast<double>(_innerRadius));
		_dispersionAtInterfaces[iface] = _dispersionDep->getValue(model, ColumnPosition{relPos, 0.0, 0.0},
		                                                           comp, ParTypeIndep, BoundStateIndep, baseDispersion);
	}
}

/**
 * @brief Recomputes S_g matrices for variable dispersion using quadrature
 * @detail For variable D(rho): S_g[i,j] = ∫ dL_i/dξ * L_j * ρ(ξ) * D(ρ(ξ)) dξ
 */
void RadialConvectionDispersionOperatorBaseDG::recomputeDispersionMatrices()
{
	if (!_variableDispersion)
		return;

	const double deltaRho = static_cast<double>(_deltaRho);

	// Number of quadrature points: use 3p/2 for overintegration, p+2 otherwise
	const int nQuadPoints = _overintegrate ? static_cast<int>((3 * _polyDeg + 1) / 2 + 1) : static_cast<int>(_polyDeg + 2);

	for (unsigned int elem = 0; elem < _nElem; ++elem)
	{
		const double rho_left = _rhoCellBounds[elem];

		// Use the radialDispersionMatrix function from DGToolbox
		_S_g[elem] = dgtoolbox::radialDispersionMatrix(_polyDeg, _nodes, rho_left, deltaRho,
		                                                _dispersionAtNodes[elem], nQuadPoints);
	}
}

/**
 * @brief Returns the current velocity at a given relative position
 */
active RadialConvectionDispersionOperatorBaseDG::currentVelocity(double pos) const CADET_NOEXCEPT
{
	// pos is relative coordinate in [0, 1]
	const active radius = static_cast<double>(_innerRadius) + pos * (static_cast<double>(_outerRadius) - static_cast<double>(_innerRadius));
	return _curVelocity / radius;
}

/**
 * @brief Returns the Jacobian factor for the inlet boundary
 */
double RadialConvectionDispersionOperatorBaseDG::inletJacobianFactor() const CADET_NOEXCEPT
{
	const double rho_inlet = static_cast<double>(_innerRadius);
	const double denom = rho_inlet * static_cast<double>(_deltaRho);
	return static_cast<double>(_curVelocity) / denom;
}

unsigned int RadialConvectionDispersionOperatorBaseDG::jacobianLowerBandwidth() const CADET_NOEXCEPT
{
	// Radial DG always uses exact integration
	return 2 * _nNodes * strideColNode();
}

unsigned int RadialConvectionDispersionOperatorBaseDG::jacobianUpperBandwidth() const CADET_NOEXCEPT
{
	return jacobianLowerBandwidth();
}

bool RadialConvectionDispersionOperatorBaseDG::setParameter(const ParameterId& pId, double value)
{
	if (_dispersionDep)
	{
		if (_dispersionDep->hasParameter(pId))
		{
			_dispersionDep->setParameter(pId, value);
			return true;
		}
	}

	if (!_dispersionCompIndep)
		return false;

	if ((pId.name != hashString("COL_DISPERSION")) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
		return false;

	if (_colDispersion.size() > _nComp)
	{
		if (pId.section == SectionIndep)
			return false;

		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[pId.section * _nComp + i].setValue(value);
	}
	else
	{
		if (pId.section != SectionIndep)
			return false;

		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[i].setValue(value);
	}

	return true;
}

bool RadialConvectionDispersionOperatorBaseDG::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
{
	if (_dispersionDep)
	{
		active* const param = _dispersionDep->getParameter(pId);
		if (param)
		{
			param->setValue(value);
			return true;
		}
	}

	if (!_dispersionCompIndep)
		return false;

	if ((pId.name != hashString("COL_DISPERSION")) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
		return false;

	if (_colDispersion.size() > _nComp)
	{
		if (pId.section == SectionIndep)
			return false;

		if (!contains(sensParams, &_colDispersion[pId.section * _nComp]))
			return false;

		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[pId.section * _nComp + i].setValue(value);
	}
	else
	{
		if (pId.section != SectionIndep)
			return false;

		if (!contains(sensParams, &_colDispersion[0]))
			return false;

		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[i].setValue(value);
	}

	return true;
}

bool RadialConvectionDispersionOperatorBaseDG::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if (_dispersionDep)
	{
		active* const param = _dispersionDep->getParameter(pId);
		if (param)
		{
			param->setADValue(adDirection, adValue);
			return true;
		}
	}

	if (!_dispersionCompIndep)
		return false;

	if ((pId.name != hashString("COL_DISPERSION")) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
		return false;

	if (_colDispersion.size() > _nComp)
	{
		if (pId.section == SectionIndep)
			return false;

		sensParams.insert(&_colDispersion[pId.section * _nComp]);
		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[pId.section * _nComp + i].setADValue(adDirection, adValue);
	}
	else
	{
		if (pId.section != SectionIndep)
			return false;

		sensParams.insert(&_colDispersion[0]);
		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[i].setADValue(adDirection, adValue);
	}

	return true;
}

void RadialConvectionDispersionOperatorBaseDG::setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT
{
	const double pi = 3.1415926535897932384626434;

	if (_colLength > 0.0)
		_curVelocity = _dir * in / (2.0 * pi * _colLength * colPorosity);
}

void RadialConvectionDispersionOperatorBaseDG::multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const
{
	double const* localSdot = sDot + offsetC();
	double* localRet = ret + offsetC();

	for (unsigned int i = 0; i < _nPoints * _nComp; ++i)
		localRet[i] = localSdot[i];
}

void RadialConvectionDispersionOperatorBaseDG::addTimeDerivativeToJacobian(double alpha, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacDisc, unsigned int blockOffset)
{
	for (unsigned int i = 0; i < _nPoints * _nComp; ++i)
		jacDisc.coeffRef(blockOffset + i, blockOffset + i) += alpha;
}

unsigned int RadialConvectionDispersionOperatorBaseDG::nJacEntries(bool pureNNZ)
{
	// Radial DG always uses exact integration
	if (pureNNZ)
	{
		return _nComp * ((3u * _nElem - 2u) * _nNodes * _nNodes + (2u * _nElem - 3u) * _nNodes); // dispersion entries
	}
	return _nComp * _nNodes * _nNodes + _nNodes // convection entries
		+ _nComp * ((3u * _nElem - 2u) * _nNodes * _nNodes + (2u * _nElem - 3u) * _nNodes); // dispersion entries
}

unsigned int RadialConvectionDispersionOperatorBaseDG::nConvDispEntries(bool pureNNZ)
{
	// Similar pattern to axial but with radial geometry
	if (!pureNNZ)
		return _nComp * (2 * jacobianLowerBandwidth() + 1);

	// Pure non-zero count
	return _nComp * _nPoints * (2 * _nNodes + 1);
}

void RadialConvectionDispersionOperatorBaseDG::convDispJacPattern(std::vector<T>& tripletList, const int bulkOffset)
{
	// Similar to axial but adapted for radial geometry
	const int lowerBW = jacobianLowerBandwidth();
	const int upperBW = jacobianUpperBandwidth();

	for (unsigned int point = 0; point < _nPoints; ++point)
	{
		for (unsigned int comp = 0; comp < _nComp; ++comp)
		{
			const int row = bulkOffset + point * _strideNode + comp;
			for (int diag = -lowerBW; diag <= static_cast<int>(upperBW); ++diag)
			{
				const int col = row + diag;
				if (col >= bulkOffset && col < static_cast<int>(bulkOffset + _nPoints * _nComp))
					tripletList.push_back(T(row, col, 0.0));
			}
		}
	}
}

int RadialConvectionDispersionOperatorBaseDG::calcTransportJacobian(Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int bulkOffset)
{
	// Assemble Jacobian using exact integration (only scheme for radial DG)
	calcConvDispRadialDGSEMJacobian(jacobian, jacInlet, bulkOffset);

	_newStaticJac = false;

	if (!jacobian.isCompressed())
		return 0;

	return 1;
}

/* Note: The residual and Jacobian assembly methods would be implemented here.
 * Due to complexity, these are placeholders that will need full implementation.
 * The key structure follows the axial implementation but with:
 * 1. Cell-dependent mass matrices _invMM_rho[elem]
 * 2. Cell-dependent stiffness matrices _S_g[elem]
 * 3. Radial geometry factors (rho at boundaries)
 */

// Placeholder for residual methods - will implement the full DG scheme
int RadialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac)
{
	return residualImpl<double, double, double, linalg::BandedEigenSparseRowIterator, true>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator(jac, offsetC()));
}

int RadialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac)
{
	return residualImpl<double, active, active, linalg::BandedEigenSparseRowIterator, true>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator(jac, offsetC()));
}

int RadialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, WithoutParamSensitivity)
{
	return residualImpl<double, double, double, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int RadialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, WithParamSensitivity)
{
	return residualImpl<double, active, active, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int RadialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithParamSensitivity)
{
	return residualImpl<active, active, active, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int RadialConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithoutParamSensitivity)
{
	return residualImpl<active, active, active, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

// Jacobian block computation for radial DG (cell-dependent)
// Note: Radial DG uses exact integration only (no collocation scheme)
Eigen::MatrixXd RadialConvectionDispersionOperatorBaseDG::DGjacobianConvBlockRadial(unsigned int cellIdx)
{
	// Returns d(res_e)/dc / u for convection only (caller multiplies by u).
	// Block is N x (N+1):
	//   Forward flow: col 0 = inlet/prev-element last node, cols 1..N = current element
	//   Backward flow: cols 0..N-1 = current element, col N = outlet/next-element first node
	//
	// Derived from residualImpl: res -= u*invMM_rho*D^T*M00*c + surface flux corrections.
	// No (2/deltaRho) factor here -- coordinate transforms cancel in the volume integral.

	Eigen::MatrixXd convBlock = Eigen::MatrixXd::Zero(_nNodes, _nNodes + 1);

	const Eigen::MatrixXd volContrib = _invMM_rho[cellIdx] * (_polyDerM.transpose() * _M00);

	if (_curVelocity >= 0.0) // Forward (outward) flow
	{
		// Volume: cols 1..N map to current element nodes 0..N-1
		convBlock.block(0, 1, _nNodes, _nNodes) -= volContrib;

		if (cellIdx == 0)
		{
			// Danckwerts inlet BC: left_flux[0] = u*c[0,0] - 2*u*c_inlet
			// d(left_flux)/d(c_inlet)/u = -2  =>  col 0 = -2*invMM_rho[:,0]
			// d(left_flux)/d(c[0,0])/u = +1 added to col 1 (= current node 0)
			convBlock.block(0, 0, _nNodes, 1) = -2.0 * _invMM_rho[cellIdx].col(0);
			convBlock.block(0, 1, _nNodes, 1) += _invMM_rho[cellIdx].col(0);
		}
		else
		{
			// Upwind left flux from prev element's last node
			convBlock.block(0, 0, _nNodes, 1) = -_invMM_rho[cellIdx].col(0);
		}

		// Right surface: c_star[e+1] = c[e,N-1] (upwind), contributes +invMM_rho[:,N-1] at col N
		// Same for all forward-flow elements (including last: g_star[nElem]=0 so no dispersion)
		convBlock.block(0, _nNodes, _nNodes, 1) += _invMM_rho[cellIdx].col(_nNodes - 1);
	}
	else // Backward (inward) flow
	{
		// Volume: cols 0..N-1 map to current element nodes 0..N-1
		convBlock.block(0, 0, _nNodes, _nNodes) -= volContrib;

		// Left surface: c_star[e] = c[e,0] (upwind from right), contributes -invMM_rho[:,0] at col 0
		// Same for all backward-flow elements (including element 0: g_star[0]=0 so no dispersion)
		convBlock.block(0, 0, _nNodes, 1) -= _invMM_rho[cellIdx].col(0);

		if (cellIdx == _nElem - 1)
		{
			// Danckwerts inlet BC: right_flux[nElem] = 2*u*c_inlet - u*c[last,N-1]
			// d(right_flux)/d(c_inlet)/u = +2  =>  col N = +2*invMM_rho[:,N-1]
			// d(right_flux)/d(c[last,N-1])/u = -1 added to col N-1
			convBlock.block(0, _nNodes - 1, _nNodes, 1) -= _invMM_rho[cellIdx].col(_nNodes - 1);
			convBlock.block(0, _nNodes, _nNodes, 1) = 2.0 * _invMM_rho[cellIdx].col(_nNodes - 1);
		}
		else
		{
			// Upwind right flux to next element's first node
			convBlock.block(0, _nNodes, _nNodes, 1) = _invMM_rho[cellIdx].col(_nNodes - 1);
		}
	}

	return convBlock;
}

Eigen::MatrixXd RadialConvectionDispersionOperatorBaseDG::getGBlockRadial(unsigned int cellIdx)
{
	// Auxiliary G block for cell cellIdx (derivative of g w.r.t. c)
	// Radial DG uses exact integration only (no collocation)
	// NOTE: Uses standard inverse mass matrix _invMM (not _invMM_rho) because
	// the auxiliary equation g = dc/dρ doesn't have radial weighting

	Eigen::MatrixXd gBlock = Eigen::MatrixXd::Zero(_nNodes, _nNodes + 2);
	const double deltaRho = static_cast<double>(_deltaRho);

	// Interior derivative term: D * c
	gBlock.block(0, 1, _nNodes, _nNodes) = _polyDerM;

	// Interface flux corrections using standard inverse mass matrix
	if (cellIdx > 0 && cellIdx < _nElem - 1)
	{
		// Interior cell: corrections at both interfaces
		// Left interface correction
		gBlock.block(0, 0, _nNodes, 1) -= 0.5 * _invMM.col(0);
		gBlock.block(0, 1, _nNodes, 1) += 0.5 * _invMM.col(0);
		// Right interface correction
		gBlock.block(0, _nNodes, _nNodes, 1) -= 0.5 * _invMM.col(_nNodes - 1);
		gBlock.block(0, _nNodes + 1, _nNodes, 1) += 0.5 * _invMM.col(_nNodes - 1);
	}
	else if (cellIdx == 0)
	{
		// First cell: only right interface correction
		if (_nElem == 1)
			return gBlock * 2.0 / deltaRho;
		gBlock.block(0, _nNodes, _nNodes, 1) -= 0.5 * _invMM.col(_nNodes - 1);
		gBlock.block(0, _nNodes + 1, _nNodes, 1) += 0.5 * _invMM.col(_nNodes - 1);
	}
	else if (cellIdx == _nElem - 1)
	{
		// Last cell: only left interface correction
		gBlock.block(0, 0, _nNodes, 1) -= 0.5 * _invMM.col(0);
		gBlock.block(0, 1, _nNodes, 1) += 0.5 * _invMM.col(0);
	}

	gBlock *= 2.0 / deltaRho;

	return gBlock;
}

Eigen::MatrixXd RadialConvectionDispersionOperatorBaseDG::DGjacobianDispBlockRadial(unsigned int cellIdx)
{
	// Dispersion Jacobian block for cell cellIdx (d_comp factor applied by caller).
	// Block is N x (3N+2). Column j maps to physical offset (j - N - 1) from the
	// current element's first node (node 0):
	//   cols 0..N+1   : prev element nodes + ghost (prev-prev last / next first)
	//   cols N..2N+1  : current element nodes + ghosts (prev last / next first)
	//   cols 2N..3N+1 : next element nodes + ghost
	//
	// Derived from residualImpl:
	//   res -= d * invMM_rho * S_g * g_ref      (volume, no 2/deltaRho factor)
	//   res += invMM_rho[:,0]  * (-u*c* + rho_L*d*(2/dR)*g*_L)  (left surface)
	//   res += invMM_rho[:,N-1]* ( u*c* - rho_R*d*(2/dR)*g*_R)  (right surface)
	// where g_ref = auxiliary in reference coords, g* = 0.5*(g_left + g_right).
	//
	// getGBlockRadial(e) = -(2/deltaRho) * dg_ref_e/dc, so:
	//   dg_ref_e/dc = -(deltaRho/2) * G_e
	//
	// The (2/deltaRho) from surface fluxes and (deltaRho/2) from dg/dc cancel,
	// leaving clean factor 0.5 in surface terms.

	const double deltaRho = static_cast<double>(_deltaRho);
	const double halfDeltaRho = 0.5 * deltaRho;
	const double rho_left = _rhoCellBounds[cellIdx];
	const double rho_right = _rhoCellBounds[cellIdx + 1];

	Eigen::MatrixXd dispBlock = Eigen::MatrixXd::Zero(_nNodes, 3 * _nNodes + 2);

	Eigen::MatrixXd G_curr = getGBlockRadial(cellIdx);

	// Volume: -(deltaRho/2) * invMM_rho * S_g * G_curr
	// G_curr has N+2 cols; placed at dispBlock cols N..2N+1
	dispBlock.block(0, _nNodes, _nNodes, _nNodes + 2) =
	    -halfDeltaRho * (_invMM_rho[cellIdx] * _S_g[cellIdx] * G_curr);

	// Left surface: -0.5 * rho_left * invMM_rho[:,0] * (G_prev[N-1,:] + G_curr[0,:])
	// G_prev[N-1,:] placed at dispBlock cols 0..N+1
	// G_curr[0,:]   placed at dispBlock cols N..2N+1
	if (cellIdx > 0)
	{
		Eigen::MatrixXd G_prev = getGBlockRadial(cellIdx - 1);
		dispBlock.block(0, 0, _nNodes, _nNodes + 2) -=
		    0.5 * rho_left * _invMM_rho[cellIdx].col(0) * G_prev.row(_nNodes - 1);
		dispBlock.block(0, _nNodes, _nNodes, _nNodes + 2) -=
		    0.5 * rho_left * _invMM_rho[cellIdx].col(0) * G_curr.row(0);
	}

	// Right surface: +0.5 * rho_right * invMM_rho[:,N-1] * (G_curr[N-1,:] + G_next[0,:])
	// G_curr[N-1,:] placed at dispBlock cols N..2N+1
	// G_next[0,:]   placed at dispBlock cols 2N..3N+1
	if (cellIdx < _nElem - 1)
	{
		Eigen::MatrixXd G_next = getGBlockRadial(cellIdx + 1);
		dispBlock.block(0, _nNodes, _nNodes, _nNodes + 2) +=
		    0.5 * rho_right * _invMM_rho[cellIdx].col(_nNodes - 1) * G_curr.row(_nNodes - 1);
		dispBlock.block(0, 2 * _nNodes, _nNodes, _nNodes + 2) +=
		    0.5 * rho_right * _invMM_rho[cellIdx].col(_nNodes - 1) * G_next.row(0);
	}

	return dispBlock;
}

void RadialConvectionDispersionOperatorBaseDG::calcConvDispRadialDGSEMJacobian(Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int offC)
{
	// Assemble radial DG Jacobian (exact integration only)
	// Cell-dependent blocks - cannot batch like axial because each cell has unique matrices

	// u = velocity_coeff (NOT multiplied by rho, matches residualImpl)
	const double u = static_cast<double>(_curVelocity);

	const int strideNode = strideColNode();
	const int strideElem = strideColElement();
	const int strideColBound = strideNode - static_cast<int>(_nComp);

	// Pre-compute cell-dependent Jacobian blocks (already include u and d_rad factors from block computation)
	for (unsigned int cell = 0; cell < _nElem; cell++)
	{
		_DGjacRadConvBlocks[cell] = u * DGjacobianConvBlockRadial(cell);
		_DGjacRadDispBlocks[cell] = DGjacobianDispBlockRadial(cell);  // Don't multiply by d_rad yet, component-dependent
	}

	/*======================================================*/
	/*			Compute Dispersion Jacobian Block			*/
	/*======================================================*/

	// Each cell needs separate treatment because blocks are cell-dependent
	// Dispersion block structure: nNodes x (3*nNodes + 2)
	// The block column to physical node mapping is:
	//   - Col nNodes-1: prev element's last node (surface flux only)
	//   - Col nNodes: coupling via gBlock col 0 (prev element interface)
	//   - Cols nNodes+1..2*nNodes: current element's nNodes nodes
	//   - Col 2*nNodes+1: coupling via gBlock col nNodes+1 (next element interface)
	//
	// Physical node index for block col j, cell k:
	//   - j < nNodes: prev element, node (j) -- only col nNodes-1 is used
	//   - j = nNodes: prev element's last node (gBlock col 0 contribution)
	//   - j in [nNodes+1, 2*nNodes]: current element, node (j - nNodes - 1)
	//   - j = 2*nNodes+1: next element's first node (gBlock col nNodes+1)

	for (unsigned int cell = 0; cell < _nElem; cell++)
	{
		const Eigen::MatrixXd& dispBlock = _DGjacRadDispBlocks[cell];
		linalg::BandedEigenSparseRowIterator cellJac(jacobian, offC + cell * strideElem);

		for (unsigned int i = 0; i < _nNodes; i++, cellJac += strideColBound)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++, ++cellJac)
			{
				const double d_comp = static_cast<double>(currentDispersion(_curSection)[comp]);

				// Iterate all 3N+2 columns of the dispersion block.
				// Column j maps to physical offset (j - N - 1) from current element's first node.
				// relOffset from row node i = (j - N - 1) - i.
				// Zero entries (at domain boundaries) are skipped by the abs check.
				for (unsigned int j = 0; j < 3 * _nNodes + 2; j++)
				{
					const double val = dispBlock(i, j) * d_comp;
					if (std::abs(val) > 1e-15)
					{
						const int relOffset = static_cast<int>(j) - static_cast<int>(_nNodes) - 1 - static_cast<int>(i);
						cellJac[relOffset * strideNode] += val;
					}
				}
			}
		}
	}

	/*======================================================*/
	/*			Compute Convection Jacobian Block			*/
	/*======================================================*/

	// Convection block structure: nNodes x (nNodes + 1)
	// Forward flow: col 0 = upwind (inlet/prev element), cols 1..nNodes = current element
	// Backward flow: cols 0..nNodes-1 = current element, col nNodes = downwind (inlet/next element)

	linalg::BandedEigenSparseRowIterator jacConv(jacobian, offC);

	if (_curVelocity >= 0.0)  // Forward flow (outward radial flow)
	{
		// First element: inlet contribution (col 0 of block)
		jacInlet = _DGjacRadConvBlocks[0].col(0);

		// Insert first element block (cols 1..nNodes, current element only)
		const Eigen::MatrixXd& firstBlock = _DGjacRadConvBlocks[0];
		for (unsigned int i = 0; i < _nNodes; i++, jacConv += strideColBound)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++, ++jacConv)
			{
				for (unsigned int j = 1; j <= _nNodes; j++)
				{
					// j-1 is the node index in current element, i is row's node index
					int nodeOffset = static_cast<int>(j - 1) - static_cast<int>(i);
					jacConv[nodeOffset * strideNode] += firstBlock(i, j);
				}
			}
		}

		// Remaining elements: full block (col 0 = prev element's last node)
		for (unsigned int cell = 1; cell < _nElem; cell++)
		{
			const Eigen::MatrixXd& convBlock = _DGjacRadConvBlocks[cell];
			for (unsigned int i = 0; i < _nNodes; i++, jacConv += strideColBound)
			{
				for (unsigned int comp = 0; comp < _nComp; comp++, ++jacConv)
				{
					for (unsigned int j = 0; j < convBlock.cols(); j++)
					{
						// j=0 is prev element's last node (offset -1 from current first)
						// j=1..nNodes is current element (offset 0..nNodes-1)
						int nodeOffset = static_cast<int>(j) - 1 - static_cast<int>(i);
						jacConv[nodeOffset * strideNode] += convBlock(i, j);
					}
				}
			}
		}
	}
	else  // Backward flow (inward radial flow)
	{
		// Non-inlet elements first (all except last)
		for (unsigned int cell = 0; cell < _nElem - 1; cell++)
		{
			const Eigen::MatrixXd& convBlock = _DGjacRadConvBlocks[cell];
			for (unsigned int i = 0; i < _nNodes; i++, jacConv += strideColBound)
			{
				for (unsigned int comp = 0; comp < _nComp; comp++, ++jacConv)
				{
					for (unsigned int j = 0; j < convBlock.cols(); j++)
					{
						// cols 0..nNodes-1 = current element, col nNodes = next element's first
						int nodeOffset = static_cast<int>(j) - static_cast<int>(i);
						jacConv[nodeOffset * strideNode] += convBlock(i, j);
					}
				}
			}
		}

		// Last element: inlet contribution (last column of block)
		const Eigen::MatrixXd& lastBlock = _DGjacRadConvBlocks[_nElem - 1];
		jacInlet = lastBlock.col(lastBlock.cols() - 1);

		// Insert last element block (cols 0..nNodes-1, current element only)
		for (unsigned int i = 0; i < _nNodes; i++, jacConv += strideColBound)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++, ++jacConv)
			{
				for (unsigned int j = 0; j < _nNodes; j++)
				{
					int nodeOffset = static_cast<int>(j) - static_cast<int>(i);
					jacConv[nodeOffset * strideNode] += lastBlock(i, j);
				}
			}
		}
	}
}

}  // namespace parts

}  // namespace model

}  // namespace cadet

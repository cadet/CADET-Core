// SPDX-License-Identifier: AGPL-3.0-or-later
// =================================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Affero General Public
//  License v3.0 (or, at your option, any later version).
// =================================================================================

#include "model/parts/ConvectionDispersionOperatorDG.hpp"
#include "cadet/Exceptions.hpp"

#include "ParamReaderHelper.hpp"
#include "AdUtils.hpp"
#include "SimulationTypes.hpp"
#include "model/ParameterDependence.hpp"
#include "SensParamUtil.hpp"
#include "ConfigurationHelper.hpp"
#include "ParamReaderHelper.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <cstring>

namespace cadet
{

namespace model
{

namespace parts
{

/**
 * @brief Creates an AxialConvectionDispersionOperatorBaseCollocationDG
 */
AxialConvectionDispersionOperatorBaseCollocationDG::AxialConvectionDispersionOperatorBaseCollocationDG() :
	_DGjacAxDispBlocks(nullptr), _auxState(nullptr), _subsState(nullptr), _dispersionDep(nullptr)
{
}

AxialConvectionDispersionOperatorBaseCollocationDG::~AxialConvectionDispersionOperatorBaseCollocationDG() CADET_NOEXCEPT
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
 * @param [in] nElem Number of axial elements
 * @param [in] polyDeg Polynomial degree of DG approach
 * @param [in] strideNode node stride in state vector
 * @return @c true if configuration went fine, @c false otherwise
 */
bool AxialConvectionDispersionOperatorBaseCollocationDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, unsigned int nComp, unsigned int nElem, unsigned int polyDeg, unsigned int strideNode)
{
	const bool firstConfigCall = _auxState == nullptr; // used to not multiply allocate memory

	_nComp = nComp;
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
	_polyDerM = dgtoolbox::derivativeMatrix(_polyDeg, _nodes);

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
bool AxialConvectionDispersionOperatorBaseCollocationDG::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters)
{
	const bool firstConfigCall = _DGjacAxDispBlocks == nullptr; // used to not multiply allocate memory

	// Read geometry parameters
	_bedLength = paramProvider.getDouble("BED_LENGTH");
	_deltaZ = _bedLength / _nElem;

	/* compute dispersion jacobian blocks(without parameters except element spacing, i.e. static entries) */
	// we only need unique dispersion blocks, which are given by elements 1, 2, nElem for inexact integration DG and by elements 1, 2, 3, nElem-1, nElem for eaxct integration DG
	// note that convection jacobian block is computed in notifyDiscontinuousSectionTransition() since this block needs to be recomputed when flow direction changes
	if (firstConfigCall)
		_DGjacAxDispBlocks = new Eigen::MatrixXd[std::min(_nElem, 3u)];
	_DGjacAxDispBlocks[0] = DGjacobianDispBlock(1);
	if (_nElem > 1)
		_DGjacAxDispBlocks[1] = DGjacobianDispBlock(2);
	if (_nElem > 2)
		_DGjacAxDispBlocks[2] = DGjacobianDispBlock(_nElem);

	_crossSection = paramProvider.getDouble("CROSS_SECTION_AREA");

	// Read section dependent parameters (transport)
	_dir.clear();
	const std::vector<double> fwdFlow = paramProvider.getDoubleArray("FORWARD_FLOW");
	for (std::size_t i = 0; i < fwdFlow.size(); ++i)
		_dir.push_back(fwdFlow[i] ? 1.0 : -1.0);

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

	parameters[makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_curVelocity;
	parameters[makeParamId(hashString("BED_LENGTH"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_bedLength;
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
bool AxialConvectionDispersionOperatorBaseCollocationDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, MatrixXd& jacInlet)
{
	_curSection = secIdx;
	_newStaticJac = true;

	 // recompute convection jacobian block, which depends on flow direction
	_DGjacAxConvBlock = DGjacobianConvBlock();

	if (jacInlet.size() == 0) // initialize inlet Jacobian if not already done
		jacInlet = Eigen::MatrixXd::Zero(1, 1);

	if (_curVelocity >= 0.0) // forward flow upwind convection
		jacInlet(0, 0) = static_cast<double>(_curVelocity) * _DGjacAxConvBlock(0, 0); // only first node depends on inlet concentration
	else // backward flow upwind convection
		jacInlet(0, 0) = static_cast<double>(_curVelocity) * _DGjacAxConvBlock(_DGjacAxConvBlock.rows() - 1, _DGjacAxConvBlock.cols() - 1); // only last node depends on inlet concentration

	// Detect change in flow direction
	const double curDir = getSectionDependentScalar(_dir, secIdx);
	const bool changedDirection = secIdx > 0 ? (getSectionDependentScalar(_dir, secIdx - 1) * curDir < 0) : false;
	if (changedDirection)
		_curVelocity *= -1.0;

	return changedDirection;
}

/**
 * @brief Sets the flow rates for the current time section
 * @details The flow rates may change due to valve switches.
 * @param [in] in Total volumetric inlet flow rate
 * @param [in] out Total volumetric outlet flow rate
 * @param [in] colPorosity Porosity used for computing interstitial velocity from volumetric flow rate
 */
void AxialConvectionDispersionOperatorBaseCollocationDG::setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT
{
	_curVelocity = in / (_crossSection * colPorosity);
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
int AxialConvectionDispersionOperatorBaseCollocationDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac)
{
	//// Reset Jacobian but keep pattern
	//double* val = jac.valuePtr();
	//for (unsigned int entry = 0; entry < jac.nonZeros(); val++)
	//	val[0] = 0.0;

	linalg::BandedEigenSparseRowIterator jacIt(jac, 0);

	return residualImpl<double, double, double, linalg::BandedEigenSparseRowIterator, true>(model, t, secIdx, y, yDot, res, jacIt);
}

int AxialConvectionDispersionOperatorBaseCollocationDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, WithoutParamSensitivity)
{
	return residualImpl<double, double, double, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int AxialConvectionDispersionOperatorBaseCollocationDG::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithoutParamSensitivity)
{
	return residualImpl<active, active, double, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int AxialConvectionDispersionOperatorBaseCollocationDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac)
{
	// todo include jacobian in convDisp operator
	//// Reset Jacobian but keep pattern
	//double* val = jac.valuePtr();
	//for (unsigned int entry = 0; entry < jac.nonZeros(); val++)
	//	val[0] = 0.0;

	linalg::BandedEigenSparseRowIterator jacIt(jac, 0);

	return residualImpl<double, active, active, linalg::BandedEigenSparseRowIterator, true>(model, t, secIdx, y, yDot, res, jacIt);
}

int AxialConvectionDispersionOperatorBaseCollocationDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, WithParamSensitivity)
{
	return residualImpl<double, active, active, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int AxialConvectionDispersionOperatorBaseCollocationDG::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithParamSensitivity)
{
	return residualImpl<active, active, active, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

template <typename StateType, typename ResidualType, typename ParamType, typename RowIteratorType, bool wantJac>
int AxialConvectionDispersionOperatorBaseCollocationDG::residualImpl(const IModel& model, double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res, RowIteratorType jacBegin)
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
template <typename StateType>
int AxialConvectionDispersionOperatorBaseCollocationDG::calcTransportJacobian(const IModel&, double t, unsigned int secIdx, Eigen::SparseMatrix<double, RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int bulkOffset, const StateType* const y)
{
	calcConvDispCollocationDGSEMJacobian(jacobian, jacInlet, bulkOffset);

	return jacobian.isCompressed();
}

int AxialConvectionDispersionOperatorBaseCollocationDG::calcTransportJacobian(const IModel& model, double t, unsigned int secIdx, Eigen::SparseMatrix<double, RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, int bulkOffset, double const* y)
{
	return calcTransportJacobian<double>(model, t, secIdx, jacobian, jacInlet, bulkOffset, y);
}

int AxialConvectionDispersionOperatorBaseCollocationDG::calcTransportJacobian(const IModel& model, double t, unsigned int secIdx, Eigen::SparseMatrix<double, RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, int bulkOffset, active const* y)
{
	return calcTransportJacobian<active>(model, t, secIdx, jacobian, jacInlet, bulkOffset, y);
}

/**
 * @brief calculates the number of entris for the DG convection dispersion jacobian
 * @note only dispersion entries are relevant for jacobian NNZ as the convection entries are a subset of these
 */
unsigned int AxialConvectionDispersionOperatorBaseCollocationDG::nJacEntries(bool pureNNZ) const CADET_NOEXCEPT
{
	if (pureNNZ)
		return _nComp * (_nElem * _nNodes * _nNodes + 8u * _nNodes); // dispersion entries
	
	return _nComp * _nNodes * _nNodes + 1u // convection entries
		+ _nComp * (_nElem * _nNodes * _nNodes + 8u * _nNodes); // dispersion entries
}
void model::parts::AxialConvectionDispersionOperatorBaseCollocationDG::convDispJacPattern(std::vector<T>& tripletList, const int bulkOffset) const
{
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
void AxialConvectionDispersionOperatorBaseCollocationDG::multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const
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
void AxialConvectionDispersionOperatorBaseCollocationDG::addTimeDerivativeToJacobian(double alpha, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacDisc, unsigned int blockOffset) const
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

unsigned int AxialConvectionDispersionOperatorBaseCollocationDG::jacobianLowerBandwidth() const CADET_NOEXCEPT
{
	// @todo use more efficient seed vectors. currently, we treat the jacobian as banded, but the pattern is actually more sparse when multiple components are considered
	// (note that active type directions are limited)
	// We have a specific jacobian structure for the collocation DG scheme:
	// 2 * N_n * (N_c + N_q) + 1 = total bandwidth (main diagonal entries maximally depend on the next and last N_n liquid phase entries of same component)
	return _nNodes * strideColNode();
}

unsigned int AxialConvectionDispersionOperatorBaseCollocationDG::jacobianUpperBandwidth() const CADET_NOEXCEPT
{
	return jacobianLowerBandwidth();
}

bool AxialConvectionDispersionOperatorBaseCollocationDG::setParameter(const ParameterId& pId, double value)
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

bool AxialConvectionDispersionOperatorBaseCollocationDG::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
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

bool AxialConvectionDispersionOperatorBaseCollocationDG::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
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

VariableCrossSectionConvectionDispersionOperatorBaseDG::VariableCrossSectionConvectionDispersionOperatorBaseDG()
	: _auxState(nullptr), _subsState(nullptr), _dispersionDep(nullptr)
{
}

VariableCrossSectionConvectionDispersionOperatorBaseDG::~VariableCrossSectionConvectionDispersionOperatorBaseDG() CADET_NOEXCEPT
{
	if (_dispersionDep)
		delete _dispersionDep;
	delete[] _auxState;
	delete[] _subsState;
}

bool VariableCrossSectionConvectionDispersionOperatorBaseDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, unsigned int nComp, unsigned int nElem, unsigned int polyDeg, unsigned int strideNode)
{
	const bool firstConfigCall = (_auxState == nullptr);

	_nComp = nComp;
	_nElem = nElem;
	_polyDeg = polyDeg;
	_nNodes = _polyDeg + 1u;
	_nPoints = _nNodes * _nElem;
	_strideNode = strideNode;
	_strideElem = _nNodes * strideNode;

	// Read frustum geometry parameters
	_geometryType = geometryTypeFromString(paramProvider.getString("GEOMETRY"));
	_flowFraction = 1.0;

	switch (_geometryType)
	{
	case GeometryType::AxialFlowCylinder:
	{
		const double crossSection = paramProvider.getDouble("CROSS_SECTION_AREA");
		const double pi = 3.14159265358979323846;
		_radiusXStart = std::sqrt(crossSection / pi);
		_radiusXEnd = _radiusXStart;
		_bedLength = paramProvider.getDouble("BED_LENGTH");
		_colHeight = _bedLength;

		if (!(crossSection > 0.0 && static_cast<double>(_bedLength) > 0.0))
			throw InvalidParameterException("Geometry parameters for AXIAL_FLOW_CYLINDER must satisfy CROSS_SECTION_AREA > 0.0, BED_LENGTH > 0.0");

		break;
	}
	case GeometryType::RadialFlowCylinderShell:
	{
		if (paramProvider.getString("GEOMETRY") == "RADIAL_FLOW_CYLINDER_SHELL_WEDGE")
			_flowFraction = paramProvider.getDouble("CIRCLE_FRACTION");

		const double areaOuter = paramProvider.getDouble("CROSS_SECTION_AREA_OUTER");
		_colHeight = paramProvider.getDouble("CYLINDER_HEIGHT");
		_bedLength = paramProvider.getDouble("BED_LENGTH");
		// A = 2 pi r h => r = A / (2 pi h)
		const double pi = 3.14159265358979323846;
		_radiusXEnd = areaOuter / (2.0 * pi * _colHeight);
		_radiusXStart = _radiusXEnd - _bedLength;

		if (paramProvider.exists("CROSS_SECTION_AREA_INNER")) // check this field for user convenience
		{
			const double areaInner = paramProvider.getDouble("CROSS_SECTION_AREA_INNER");
			const double areaInnerDerived = 2.0 * pi * static_cast<double>(_radiusXStart) * static_cast<double>(_colHeight);
			if (std::abs(areaInner - areaInnerDerived) > 1e-12 * std::max(1.0, std::abs(areaInnerDerived)))
				throw InvalidParameterException("For radial flow cylinder shell, CROSS_SECTION_AREA_INNER must equal the area derived from CROSS_SECTION_AREA_OUTER, CYLINDER_HEIGHT, and BED_LENGTH");
		}

		if (!(static_cast<double>(_radiusXStart) > 0.0 && static_cast<double>(_radiusXEnd) >= static_cast<double>(_radiusXStart) && static_cast<double>(_colHeight) > 0.0))
			throw InvalidParameterException("Geometry parameters for RADIAL_FLOW_CYLINDER_SHELL must satisfy CROSS_SECTION_AREA_OUTER > 0.0, CYLINDER_HEIGHT > 0.0, 0 < BED_LENGTH <= outer radius");

		break;
	}
	case GeometryType::AxialFlowFrustum:
	{
		const double areaLarge = paramProvider.getDouble("CROSS_SECTION_AREA_LARGE_END");
		const double areaSmall = paramProvider.getDouble("CROSS_SECTION_AREA_SMALL_END");
		const double pi = 3.14159265358979323846;
		// A = pi r^2 => r = sqrt(A / pi)
		_radiusXStart = sqrt(areaLarge / pi);
		_radiusXEnd = sqrt(areaSmall / pi);
		_bedLength = paramProvider.getDouble("BED_LENGTH");
		_colHeight = _bedLength;

		if (!(static_cast<double>(_radiusXStart) > 0.0 &&
			static_cast<double>(_radiusXStart) >= static_cast<double>(_radiusXEnd) &&
			static_cast<double>(_bedLength) > 0.0))
			throw InvalidParameterException("Geometry parameters for AXIAL_FLOW_FRUSTUM must satisfy 0 < CROSS_SECTION_AREA_SMALL_END <= CROSS_SECTION_AREA_LARGE_END, BED_LENGTH > 0.0");

		break;
	}
	default:
		throw InvalidParameterException("Unsupported geometry type " + paramProvider.getString("COL_GEOMETRY"));
		break;
	}

	_deltaX = static_cast<double>(_bedLength) / static_cast<double>(_nElem);

	// Allocate DG matrices
	_nodes.resize(_nNodes);
	_nodes.setZero();
	_invWeights.resize(_nNodes);
	_invWeights.setZero();
	_polyDerM.resize(_nNodes, _nNodes);
	_polyDerM.setZero();
	_M00.resize(_nNodes, _nNodes);
	_M00.setZero();

	if (firstConfigCall)
	{
		_auxState = new active[_nPoints];
		_subsState = new active[_nPoints];
	}
	for (unsigned int i = 0; i < _nPoints; i++) {
		_auxState[i] = 0.0;
		_subsState[i] = 0.0;
	}

	_surfaceFluxC.resize(_nElem + 1u);
	_surfaceFluxC.setZero();
	_surfaceFluxG.resize(_nElem + 1u);
	_surfaceFluxG.setZero();

	_newStaticJac = true;

	// LGL nodes and DG operators (except geometry dependent ones)
	dgtoolbox::lglNodesWeights(_polyDeg, _nodes, _invWeights, /*invertWeights=*/true);
	_polyDerM = dgtoolbox::derivativeMatrix(_polyDeg, _nodes);
	_M00 = dgtoolbox::mMatrix(_polyDeg, _nodes, 0.0, 0.0);
	_invM00 = _M00.inverse();

	// Dispersion parameter dependence (optional)
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

	_variableDispersion = std::strcmp(_dispersionDep->name(), "CONSTANT_ONE") != 0;
	if (_variableDispersion)
	{
		if (!paramProvider.exists("DISPERSION_SPATIAL_DEPENDENCE_POLYDEG"))
			throw InvalidParameterException("Missing DISPERSION_SPATIAL_DEPENDENCE_POLYDEG parameter for variable dispersion, see COL_DISPERSION_DEP");

		_axDispQuadDeg = paramProvider.getInt("DISPERSION_SPATIAL_DEPENDENCE_POLYDEG");
	}
	else _axDispQuadDeg = 0;

	// Allocate per-cell Jacobian blocks
	_DGjacDispBlocks.resize(_nComp);
	for (auto& comps : _DGjacDispBlocks)
		comps.resize(_nElem);
	_DGjacConvBlocks.resize(_nElem);

	return true;
}

bool VariableCrossSectionConvectionDispersionOperatorBaseDG::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters)
{
	readScalarParameterOrArray(_colDispersion, paramProvider, "COL_DISPERSION", 1);

	readScalarParameterOrArray(_forwardFlow, paramProvider, "FORWARD_FLOW", 1);
	_curFwdFlow = static_cast<bool>(_forwardFlow[0]);

	if (paramProvider.exists("COL_DISPERSION_MULTIPLEX"))
	{
		const int mode = paramProvider.getInt("COL_DISPERSION_MULTIPLEX");
		if (mode == 0 || mode == 2)
			_dispersionCompIndep = true;
		else
			_dispersionCompIndep = false;

		if (!_dispersionCompIndep && (_colDispersion.size() % _nComp != 0))
			throw InvalidParameterException("Number of elements in field COL_DISPERSION is not a positive multiple of NCOMP (" + std::to_string(_nComp) + ")");
	}
	else
	{
		_dispersionCompIndep = ((_colDispersion.size() % _nComp) != 0);
	}

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

	// Register parameters
	if (_dispersionCompIndep)
	{
		if (_colDispersion.size() > _nComp)
		{
			for (std::size_t i = 0; i < _colDispersion.size() / _nComp; ++i)
				parameters[makeParamId(hashString("COL_DISPERSION"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, static_cast<int>(i))] = &_colDispersion[i * _nComp];
		}
		else
		{
			parameters[makeParamId(hashString("COL_DISPERSION"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colDispersion[0];
		}
	}
	else
		registerParam2DArray(parameters, _colDispersion, [=](bool multi, unsigned int sec, unsigned int comp) {
		return makeParamId(hashString("COL_DISPERSION"), unitOpIdx, comp, ParTypeIndep, BoundStateIndep, ReactionIndep, multi ? static_cast<int>(sec) : SectionIndep);
			}, _nComp);

	parameters[makeParamId(hashString("BED_LENGTH"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_bedLength;
	parameters[makeParamId(hashString("COL_RADIUS_SMALL_END"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_radiusXStart;
	parameters[makeParamId(hashString("COL_RADIUS_LARGE_END"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_radiusXEnd;

	// Geometry dependent operators
	computeGeometry();
	computeOperators(0);

	return true;
}

bool VariableCrossSectionConvectionDispersionOperatorBaseDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, Eigen::MatrixXd& jacInlet)
{
	_curSection = secIdx;
	_newStaticJac = true;

	// note: flow direction is set by setFlowRates() before notifyDiscontinuousSectionTransition() is called
	_curFwdFlow = static_cast<bool>(getSectionDependentScalar(_forwardFlow, secIdx));
	const bool changedDirection = secIdx > 0 ? (static_cast<bool>(getSectionDependentScalar(_forwardFlow, secIdx - 1)) && _curFwdFlow) : false;

	// Recompute Jacobian blocks
	for (unsigned int elem = 0; elem < _nElem; elem++)
	{
		_DGjacConvBlocks[elem] = DGjacobianConvBlock(elem);
		for (unsigned int comp = 0; comp < _nComp; comp++)
			_DGjacDispBlocks[comp][elem] = DGjacobianDispBlock(elem, comp);
	}

	if (_curFwdFlow)
		jacInlet = _DGjacConvBlocks[0].col(0);
	else
		jacInlet = _DGjacConvBlocks[_nElem - 1].col(_DGjacConvBlocks[_nElem - 1].cols() - 1);

	// some model parameters are baked into the DG operators but may be updated at section transitions
	computeOperators(secIdx);

	return changedDirection;
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT
{
	_flowRate = in;
	_QOverEps = _curFwdFlow ? _flowRate / colPorosity / _flowFraction : -_flowRate / colPorosity / _flowFraction;
}

unsigned int VariableCrossSectionConvectionDispersionOperatorBaseDG::jacobianLowerBandwidth() const CADET_NOEXCEPT
{
	return 2 * _nNodes * strideColNode();
}

unsigned int VariableCrossSectionConvectionDispersionOperatorBaseDG::jacobianUpperBandwidth() const CADET_NOEXCEPT
{
	return jacobianLowerBandwidth();
}

unsigned int VariableCrossSectionConvectionDispersionOperatorBaseDG::nJacEntries(bool pureNNZ) const CADET_NOEXCEPT
{
	if (pureNNZ)
		return _nComp * ((3u * _nElem - 2u) * _nNodes * _nNodes + (2u * _nElem - 3u) * _nNodes);
	return _nComp * _nNodes * _nNodes + _nNodes
		+ _nComp * ((3u * _nElem - 2u) * _nNodes * _nNodes + (2u * _nElem - 3u) * _nNodes);
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::convDispJacPattern(
	std::vector<T>& tripletList, const int bulkOffset) const
{
	const int lowerBW = static_cast<int>(jacobianLowerBandwidth());
	const int upperBW = static_cast<int>(jacobianUpperBandwidth());

	for (unsigned int point = 0; point < _nPoints; ++point)
	{
		for (unsigned int comp = 0; comp < _nComp; ++comp)
		{
			const int row = bulkOffset + static_cast<int>(point) * strideColNode() + static_cast<int>(comp);
			for (int diag = -lowerBW; diag <= upperBW; ++diag)
			{
				const int col = row + diag;
				if (col >= bulkOffset && col < static_cast<int>(bulkOffset + _nPoints * strideColNode()))
					tripletList.push_back(T(row, col, 0.0));
			}
		}
	}
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const
{
	double const* localSdot = sDot + offsetC();
	double* localRet = ret + offsetC();
	const int gapElem = strideColNode() - static_cast<int>(_nComp) * strideColComp();

	for (unsigned int i = 0; i < _nPoints; ++i, localRet += gapElem, localSdot += gapElem)
		for (unsigned int j = 0; j < _nComp; ++j, ++localRet, ++localSdot)
			*localRet = *localSdot;
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::addTimeDerivativeToJacobian(double alpha, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacDisc, unsigned int blockOffset) const
{
	const int gapElem = strideColNode() - static_cast<int>(_nComp) * strideColComp();
	linalg::BandedEigenSparseRowIterator jac(jacDisc, blockOffset);

	for (unsigned int point = 0; point < _nPoints; ++point, jac += gapElem)
		for (unsigned int comp = 0; comp < _nComp; ++comp, ++jac)
			jac[0] += alpha;
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::computeGeometry()
{
	switch (_geometryType)
	{
	case GeometryType::AxialFlowCylinder:
		computeGeometryAxial();
		break;
	case GeometryType::RadialFlowCylinderShell:
		computeGeometryRadial();
		break;
	case GeometryType::AxialFlowFrustum:
		computeGeometryFrustum();
		break;
	default:
		break;
	}
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::computeOperators(const unsigned int secIdx)
{
	switch (_geometryType)
	{
	case GeometryType::AxialFlowCylinder:
		computeOperatorsAxial();
		break;
	case GeometryType::RadialFlowCylinderShell:
		computeOperatorsRadial(secIdx);
		break;
	case GeometryType::AxialFlowFrustum:
		computeOperatorsFrustum(secIdx);
		break;
	default:
		break;
	}
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::computeGeometryAxial()
{
	_crossSectionArea.resize(_nPoints);
	const double pi = 3.14159265358979323846;
	const double radius = static_cast<double>(_radiusXStart);

	for (int elem = 0; elem < _nElem; ++elem)
	{
		for (int node = 0; node < _nNodes; ++node)
		{
			_crossSectionArea[elem * _nNodes + node] = pi * radius * radius;
		}
	}
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::computeOperatorsAxial()
{
	const double pi = 3.14159265358979323846;
	const double radius = static_cast<double>(_radiusXStart);
	const double A = pi * radius * radius;   // constant cross-sectional area
	const Eigen::MatrixXd M_A = A * _M00;    // area-weighted mass matrix

	_invMM_A.resize(_nElem);
	_invMM_A_times_ST_AD.resize(_nComp);
	for (auto& comps : _invMM_A_times_ST_AD)
		comps.resize(_nElem);
	_invMM_A_times_DT_timesM00.resize(_nElem);

	for (unsigned int elem = 0; elem < _nElem; ++elem)
	{
		_invMM_A[elem] = M_A.inverse();  // = (1/A) * M00^{-1}

		for (int comp = 0; comp < _nComp; comp++)
		{
			// (M^A)^{-1} * D^T * d * M^A = M00^{-1} * D^T * d * M00
			_invMM_A_times_ST_AD[comp][elem] = _invMM_A[elem] * _polyDerM.transpose()
				* static_cast<double>(_colDispersion[comp]) * M_A;
		}
		// (M^A)^{-1} * D^T * M00 = (1/A) * M00^{-1} * D^T * M00
		_invMM_A_times_DT_timesM00[elem] = _invMM_A[elem] * _polyDerM.transpose() * _M00;
	}
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::computeGeometryRadial()
{
	const double cylinderHeight = static_cast<double>(_colHeight);
	const double drho = static_cast<double>(_deltaX);
	const double rhoInner = static_cast<double>(_radiusXStart);
	const double pi = 3.14159265358979323846;

	_crossSectionArea.resize(_nPoints);

	for (int elem = 0; elem < _nElem; ++elem)
	{
		for (int node = 0; node < _nNodes; ++node)
		{
			const double rho = rhoInner + elem * drho + 0.5 * drho * (1.0 + _nodes[node]);

			_crossSectionArea[elem * _nNodes + node] = 2.0 * pi * rho * cylinderHeight;
		}
	}
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::computeOperatorsRadial(const unsigned int secIdx)
{
	const double pi = 3.14159265358979323846;
	Eigen::MatrixXd M01 = dgtoolbox::mMatrix(_polyDeg, _nodes, 0.0, 1.0);

	_invMM_A.resize(_nElem);
	_invMM_A_times_ST_AD.resize(_nComp);
	for (auto& comps : _invMM_A_times_ST_AD)
		comps.resize(_nElem);
	_invMM_A_times_DT_timesM00.resize(_nElem);

	// number of nodes for exact gauss quadrature of the integral with dispersion weight and geometric factor -> _axDispQuadDeg + 1
	const int nQuadNodes = std::ceil((_axDispQuadDeg + 1 + 2 * _polyDeg + 1) / 2);
	// geometric weight factors as ( \sum_n (1 - \xi)^n \alpha_n + \sum_m (1 + \xi)^m \beta_m + gamma) = A(x^e(\xi)) = A((\xi + 1) \DeltaX_i / 2 + x_i)
	// with A = 2 * \pi *  H * [ (\xi + 1) \DeltaX_i / 2 + x_i ]
	const std::vector<double> alpha = { };
	std::vector<double> beta = { 0.0 };

	for (unsigned int elem = 0; elem < _nElem; ++elem)
	{
		const double x_L = static_cast<double>(_radiusXStart) + elem * static_cast<double>(_deltaX);
		beta[0] = 2.0 * pi * static_cast<double>(_colHeight) * static_cast<double>(_deltaX) / 2.0;
		const double gamma = 2.0 * pi * static_cast<double>(_colHeight) * x_L;

		// M^A = 2 pi H (x_i M^(0,0) + \frac{\Delta x_i}{2} M^(0,1))
		const Eigen::MatrixXd M_A = 2.0 * pi * static_cast<double>(_colHeight) * (x_L * _M00 + 0.5 * static_cast<double>(_deltaX) * M01);
		_invMM_A[elem] = M_A.inverse();

		for (int comp = 0; comp < _nComp; comp++)
		{
			const double baseDispersion = static_cast<double>(currentDispersion(secIdx)[comp]);

			if (!_variableDispersion) // we can compute the integral exactly
			{
				// Matrix product (M^A)^-1 * D^T * \tilde{S}^(AD) = (M^A)^-1 * D^T * \tilde{M}^(AD)
				_invMM_A_times_ST_AD[comp][elem] = _invMM_A[elem] * _polyDerM.transpose() * baseDispersion * M_A;
			}
			else
			{
				// Evaluate dispersion at each quadrature node
				Eigen::VectorXd dispAtQNodes(nQuadNodes);
				Eigen::VectorXd quadNodes = Eigen::VectorXd::Zero(nQuadNodes);
				Eigen::VectorXd quadWeights = Eigen::VectorXd::Zero(nQuadNodes);
				dgtoolbox::lgNodesWeights(nQuadNodes - 1, quadNodes, quadWeights, false);

				for (unsigned int node = 0; node < nQuadNodes; ++node)
				{
					const double physicalPos = x_L + 0.5 * static_cast<double>(_deltaX) * (1.0 + quadNodes[node]);
					// Normalize position to [0, 1] for parameter dependence evaluation
					const double relPos = physicalPos / static_cast<double>(_bedLength);

					// Evaluate D(rho) = baseDispersion * dependence_factor(rho)
					dispAtQNodes[node] = _dispersionDep->getValue(ColumnPosition{ relPos, 0.0, 0.0 },
						comp, ParTypeIndep, BoundStateIndep, baseDispersion);
				}

				// \tilde{M}^(AD)_{i,j} = \int_{-1}^1 \ell_i(\xi) \frac{\partial \ell_j(\xi)}{\partial \xi} w(\xi) d\xi
				Eigen::MatrixXd gaussMM_AD = dgtoolbox::weightedQuadMassMatrix(
					_nodes,
					dispAtQNodes,
					alpha,
					beta,
					gamma,
					quadNodes,
					quadWeights);

				// Matrix product (M^A)^-1 * D^T * \tilde{S}^(AD) = (M^A)^-1 * D^T * \tilde{M}^(AD)
				_invMM_A_times_ST_AD[comp][elem] = _invMM_A[elem] * _polyDerM.transpose() * gaussMM_AD;
			}
		}
		// Matrix product (M^A)^-1 D^T * M00
		_invMM_A_times_DT_timesM00[elem] = _invMM_A[elem] * _polyDerM.transpose() * _M00;
	}
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::computeGeometryFrustum()
{
	const double H = static_cast<double>(_bedLength);
	const double r0 = static_cast<double>(_radiusXStart);
	const double rH = static_cast<double>(_radiusXEnd);
	const double dx = static_cast<double>(_deltaX);
	const double pi = 3.14159265358979323846;

	_crossSectionArea.resize(_nPoints);

	for (int elem = 0; elem < _nElem; ++elem)
	{
		for (int node = 0; node < _nNodes; ++node)
		{
			const double x = elem * dx + 0.5 * dx * (1.0 + _nodes[node]);
			const double r = r0 + (x / H) * (rH - r0);

			_crossSectionArea[elem * _nNodes + node] = pi * r * r;
		}
	}
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::computeOperatorsFrustum(const unsigned int secIdx)
{
	const double pi = 3.14159265358979323846;
	const double H = static_cast<double>(_bedLength);
	const double dx = static_cast<double>(_deltaX);
	const double r0 = static_cast<double>(_radiusXStart);
	const double rH = static_cast<double>(_radiusXEnd);
	const double rDiff = rH - r0;

	Eigen::MatrixXd M01 = dgtoolbox::mMatrix(_polyDeg, _nodes, 0.0, 1.0);
	Eigen::MatrixXd M02 = dgtoolbox::mMatrix(_polyDeg, _nodes, 0.0, 2.0);

	_invMM_A.resize(_nElem);
	_invMM_A_times_ST_AD.resize(_nComp);
	for (auto& comps : _invMM_A_times_ST_AD)
		comps.resize(_nElem);
	_invMM_A_times_DT_timesM00.resize(_nElem);

	// number of nodes for exact gauss quadrature of the integral with dispersion weight
	const int frustumGeomFactorDegree = 2; // we have r(x)^2 = (r0 + x/H * (rH - r0))^2, which is a polynomial of degree 2 in x
	const int nQuadNodes = std::ceil((_axDispQuadDeg + frustumGeomFactorDegree + 2 * _polyDeg + 1) / 2);
	// geometric weight factors as ( \sum_n (1 - \xi)^n \alpha_n + \sum_m (1 + \xi)^m \beta_m + gamma) = A(x^e(\xi)) = A((\xi + 1) \DeltaX_i / 2 + x_i)
	// with A = \pi r(x)^2 and r(x) = r_0 + \frac{x}{H} \left( r_{L^\mathrm{b}} - r_0 \right)
	const std::vector<double> alpha = { };
	std::vector<double> beta = { 0.0, 0.0 };

	for (unsigned int elem = 0; elem < _nElem; ++elem)
	{
		const double x_L = elem * dx;
		const double beta1 = (r0 * rDiff * dx / H + rDiff * rDiff / H / H * x_L * dx);
		const double beta2 = (rDiff * rDiff / H / H * dx * dx / 2.0 / 2.0);
		beta[0] = beta1; beta[1] = beta2;
		const double gamma = (r0 * r0 + 2.0 * r0 * x_L / H * rDiff + x_L * x_L / H / H * rDiff * rDiff);

		Eigen::MatrixXd M_A = gamma * _M00;
		M_A += beta1 * M01;
		M_A += beta2 * M02;
		M_A *= pi;
		_invMM_A[elem] = M_A.inverse();

		for (int comp = 0; comp < _nComp; comp++)
		{
			const double baseDispersion = static_cast<double>(currentDispersion(secIdx)[comp]);

			if (!_variableDispersion)
			{
				// we can compute the integral exactly
				Eigen::MatrixXd gaussMM_AD = baseDispersion * M_A;
				// Matrix product (M^A)^-1 * D^T * \tilde{S}^(AD) = (M^A)^-1 * D^T * \tilde{M}^(AD)
				_invMM_A_times_ST_AD[comp][elem] = _invMM_A[elem] * _polyDerM.transpose() * gaussMM_AD;
			}
			else
			{
				// Evaluate dispersion at each quadrature node
				Eigen::VectorXd dispAtQNodes(nQuadNodes);
				Eigen::VectorXd quadNodes = Eigen::VectorXd::Zero(nQuadNodes);
				Eigen::VectorXd quadWeights = Eigen::VectorXd::Zero(nQuadNodes);
				dgtoolbox::lgNodesWeights(nQuadNodes - 1, quadNodes, quadWeights, false);

				for (unsigned int node = 0; node < nQuadNodes; ++node)
				{
					const double physicalPos = x_L + 0.5 * dx * (1.0 + quadNodes[node]);
					// Normalize position to [0, 1] for parameter dependence evaluation
					const double relPos = physicalPos / H;

					// Evaluate D(rho) = baseDispersion * dependence_factor(rho)
					dispAtQNodes[node] = _dispersionDep->getValue(ColumnPosition{ relPos, 0.0, 0.0 },
						comp, ParTypeIndep, BoundStateIndep, baseDispersion);
				}

				// \tilde{M}^(AD)_{i,j} = \int_{-1}^1 \ell_i(\xi) \frac{\partial \ell_j(\xi)}{\partial \xi} w(\xi) d\xi
				Eigen::MatrixXd gaussMM_AD = dgtoolbox::weightedQuadMassMatrix(
					_nodes,
					dispAtQNodes,
					alpha,
					beta,
					gamma,
					quadNodes,
					quadWeights);

				// Matrix product (M^A)^-1 * D^T * \tilde{S}^(AD) = (M^A)^-1 * D^T * \tilde{M}^(AD)
				_invMM_A_times_ST_AD[comp][elem] = _invMM_A[elem] * _polyDerM.transpose() * gaussMM_AD;
			}
		}
		// Matrix product (M^A)^-1 D^T * M00
		_invMM_A_times_DT_timesM00[elem] = _invMM_A[elem] * _polyDerM.transpose() * _M00;
	}
}

// ============================================================================
//  Jacobian block helpers
// ============================================================================

Eigen::MatrixXd VariableCrossSectionConvectionDispersionOperatorBaseDG::getGBlock(unsigned int elemIdx)
{
	// Identical to radial getGBlockRadial — auxiliary equation uses standard M^{(0,0)} inverse.
	// g = (2/deltaX) * [D*c + M^{-1} * B * (c* - c)]
	Eigen::MatrixXd gBlock = Eigen::MatrixXd::Zero(_nNodes, _nNodes + 2);
	const double deltaX = static_cast<double>(_deltaX);

	gBlock.block(0, 1, _nNodes, _nNodes) = _polyDerM;

	if (elemIdx > 0 && elemIdx < _nElem - 1)
	{
		gBlock.block(0, 0, _nNodes, 1) -= 0.5 * _invM00.col(0);
		gBlock.block(0, 1, _nNodes, 1) += 0.5 * _invM00.col(0);
		gBlock.block(0, _nNodes, _nNodes, 1) -= 0.5 * _invM00.col(_nNodes - 1);
		gBlock.block(0, _nNodes + 1, _nNodes, 1) += 0.5 * _invM00.col(_nNodes - 1);
	}
	else if (elemIdx == 0)
	{
		if (_nElem == 1)
			return gBlock * 2.0 / deltaX;
		gBlock.block(0, _nNodes, _nNodes, 1) -= 0.5 * _invM00.col(_nNodes - 1);
		gBlock.block(0, _nNodes + 1, _nNodes, 1) += 0.5 * _invM00.col(_nNodes - 1);
	}
	else // last cell
	{
		gBlock.block(0, 0, _nNodes, 1) -= 0.5 * _invM00.col(0);
		gBlock.block(0, 1, _nNodes, 1) += 0.5 * _invM00.col(0);
	}

	return gBlock * (2.0 / deltaX);
}

Eigen::MatrixXd VariableCrossSectionConvectionDispersionOperatorBaseDG::DGjacobianConvBlock(unsigned int elemIdx)
{
	const double invHalfDeltaX = 2.0 / static_cast<double>(_deltaX);
	const double QOverEps = static_cast<double>(_QOverEps);

	Eigen::MatrixXd convBlock = Eigen::MatrixXd::Zero(_nNodes, _nNodes + 1);

	if (_curFwdFlow)
	{
		convBlock.block(0, 1, _nNodes, _nNodes) -= QOverEps * _invMM_A_times_DT_timesM00[elemIdx];

		convBlock.block(0, 0, _nNodes, 1) = -QOverEps * _invMM_A[elemIdx].col(0);

		convBlock.block(0, _nNodes, _nNodes, 1) += QOverEps * _invMM_A[elemIdx].col(_nNodes - 1);
	}
	else // backward flow
	{
		convBlock.block(0, 0, _nNodes, _nNodes) -= QOverEps * _invMM_A_times_DT_timesM00[elemIdx];

		convBlock.block(0, 0, _nNodes, 1) -= QOverEps * _invMM_A[elemIdx].col(0);

		convBlock.block(0, _nNodes, _nNodes, 1) = QOverEps * _invMM_A[elemIdx].col(_nNodes - 1);
	}

	return invHalfDeltaX * convBlock;
}

Eigen::MatrixXd VariableCrossSectionConvectionDispersionOperatorBaseDG::DGjacobianDispBlock(unsigned int elemIdx, unsigned int comp)
{
	const double invHalfDeltaX = 2.0 / static_cast<double>(_deltaX);
	const double d_comp = static_cast<double>(getSectionDependentSlice(_colDispersion, _nComp, _curSection)[comp]);

	Eigen::MatrixXd dispBlock = Eigen::MatrixXd::Zero(_nNodes, 3 * _nNodes + 2);
	Eigen::MatrixXd G_curr = getGBlock(elemIdx);

	// Volume: invMA * S_g * G_curr
	dispBlock.block(0, _nNodes, _nNodes, _nNodes + 2) = invHalfDeltaX * (_invMM_A_times_ST_AD[comp][elemIdx] * G_curr);

	const double Aleft = _crossSectionArea[elemIdx * _nNodes];
	const double Aright = _crossSectionArea[(elemIdx + 1) * _nNodes - 1];

	// Left surface
	if (elemIdx > 0)
	{
		Eigen::MatrixXd G_prev = getGBlock(elemIdx - 1);
		dispBlock.block(0, 0, _nNodes, _nNodes + 2) +=
			invHalfDeltaX * 0.5 * Aleft * d_comp * _invMM_A[elemIdx].col(0) * G_prev.row(_nNodes - 1);
		dispBlock.block(0, _nNodes, _nNodes, _nNodes + 2) +=
			invHalfDeltaX * 0.5 * Aleft * d_comp * _invMM_A[elemIdx].col(0) * G_curr.row(0);
	}

	// Right surface
	if (elemIdx < _nElem - 1)
	{
		Eigen::MatrixXd G_next = getGBlock(elemIdx + 1);
		dispBlock.block(0, _nNodes, _nNodes, _nNodes + 2) -=
			invHalfDeltaX * 0.5 * Aright * d_comp * _invMM_A[elemIdx].col(_nNodes - 1) * G_curr.row(_nNodes - 1);
		dispBlock.block(0, 2 * _nNodes, _nNodes, _nNodes + 2) -=
			invHalfDeltaX * 0.5 * Aright * d_comp * _invMM_A[elemIdx].col(_nNodes - 1) * G_next.row(0);
	}

	return dispBlock;
}

void VariableCrossSectionConvectionDispersionOperatorBaseDG::calcConvDispDGSEMJacobian(Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int offC)
{
	const int strideNode = strideColNode();
	const int strideElem = strideColElement();
	const int strideColBound = strideNode - static_cast<int>(_nComp);

	/*======================================================*/
	/*         Dispersion Jacobian                          */
	/*======================================================*/

	for (unsigned int elem = 0; elem < _nElem; elem++)
	{
		linalg::BandedEigenSparseRowIterator elemJac(jacobian, offC + elem * strideElem);

		for (unsigned int i = 0; i < _nNodes; i++, elemJac += strideColBound)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++, ++elemJac)
			{
				for (unsigned int j = 0; j < 3 * _nNodes + 2; j++)
				{
					const double val = _DGjacDispBlocks[comp][elem](i, j);
					if (std::abs(val) > 1e-15)
					{
						const int relOffset = static_cast<int>(j) - static_cast<int>(_nNodes) - 1 - static_cast<int>(i);
						elemJac[relOffset * strideNode] += val;
					}
				}
			}
		}
	}

	/*======================================================*/
	/*         Convection Jacobian                          */
	/*======================================================*/

	linalg::BandedEigenSparseRowIterator jacConv(jacobian, offC);

	if (_curFwdFlow)
	{
		jacInlet = _DGjacConvBlocks[0].col(0);

		for (unsigned int i = 0; i < _nNodes; i++, jacConv += strideColBound)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++, ++jacConv)
			{
				for (unsigned int j = 1; j <= _nNodes; j++)
				{
					const int nodeOffset = static_cast<int>(j - 1) - static_cast<int>(i);
					jacConv[nodeOffset * strideNode] += _DGjacConvBlocks[0](i, j);
				}
			}
		}

		for (unsigned int elem = 1; elem < _nElem; elem++)
		{
			for (unsigned int i = 0; i < _nNodes; i++, jacConv += strideColBound)
			{
				for (unsigned int comp = 0; comp < _nComp; comp++, ++jacConv)
				{
					for (unsigned int j = 0; j < static_cast<unsigned int>(_DGjacConvBlocks[elem].cols()); j++)
					{
						const int nodeOffset = static_cast<int>(j) - 1 - static_cast<int>(i);
						jacConv[nodeOffset * strideNode] += _DGjacConvBlocks[elem](i, j);
					}
				}
			}
		}
	}
	else // backward flow
	{
		for (unsigned int elem = 0; elem < _nElem - 1; elem++)
		{
			for (unsigned int i = 0; i < _nNodes; i++, jacConv += strideColBound)
			{
				for (unsigned int comp = 0; comp < _nComp; comp++, ++jacConv)
				{
					for (unsigned int j = 0; j < static_cast<unsigned int>(_DGjacConvBlocks[elem].cols()); j++)
					{
						const int nodeOffset = static_cast<int>(j) - static_cast<int>(i);
						jacConv[nodeOffset * strideNode] += _DGjacConvBlocks[elem](i, j);
					}
				}
			}
		}

		jacInlet = _DGjacConvBlocks[_nElem - 1].col(_DGjacConvBlocks[_nElem - 1].cols() - 1);

		for (unsigned int i = 0; i < _nNodes; i++, jacConv += strideColBound)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++, ++jacConv)
			{
				for (unsigned int j = 0; j < _nNodes; j++)
				{
					const int nodeOffset = static_cast<int>(j) - static_cast<int>(i);
					jacConv[nodeOffset * strideNode] += _DGjacConvBlocks[_nElem - 1](i, j);
				}
			}
		}
	}

	_newStaticJac = false;
}

// ============================================================================
//  calcTransportJacobian
// ============================================================================

template <typename StateType>
int VariableCrossSectionConvectionDispersionOperatorBaseDG::calcTransportJacobian(const IModel&, double t, unsigned int secIdx, Eigen::SparseMatrix<double, RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int bulkOffset, const StateType* const /*y*/)
{
	calcConvDispDGSEMJacobian(jacobian, jacInlet, bulkOffset);
	return jacobian.isCompressed();
}

int VariableCrossSectionConvectionDispersionOperatorBaseDG::calcTransportJacobian(
	const IModel& model, double t, unsigned int secIdx,
	Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian,
	Eigen::MatrixXd& jacInlet, int bulkOffset, double const* y)
{
	return calcTransportJacobian<double>(model, t, secIdx, jacobian, jacInlet, bulkOffset, y);
}

int VariableCrossSectionConvectionDispersionOperatorBaseDG::calcTransportJacobian(
	const IModel& model, double t, unsigned int secIdx,
	Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian,
	Eigen::MatrixXd& jacInlet, int bulkOffset, active const* y)
{
	return calcTransportJacobian<active>(model, t, secIdx, jacobian, jacInlet, bulkOffset, y);
}

// ============================================================================
//  residual overloads
// ============================================================================

int VariableCrossSectionConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac)
{
	return residualImpl<double, double, double, linalg::BandedEigenSparseRowIterator, true>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator(jac, offsetC()));
}

int VariableCrossSectionConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, Eigen::SparseMatrix<double, Eigen::RowMajor>& jac)
{
	return residualImpl<double, active, active, linalg::BandedEigenSparseRowIterator, true>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator(jac, offsetC()));
}

int VariableCrossSectionConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, WithoutParamSensitivity)
{
	return residualImpl<double, double, double, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int VariableCrossSectionConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, WithParamSensitivity)
{
	return residualImpl<double, active, active, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int VariableCrossSectionConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithParamSensitivity)
{
	return residualImpl<active, active, active, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

int VariableCrossSectionConvectionDispersionOperatorBaseDG::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithoutParamSensitivity)
{
	return residualImpl<active, active, active, linalg::BandedEigenSparseRowIterator, false>(model, t, secIdx, y, yDot, res, linalg::BandedEigenSparseRowIterator());
}

template <typename StateType, typename ResidualType, typename ParamType, typename RowIteratorType, bool wantJac>
int VariableCrossSectionConvectionDispersionOperatorBaseDG::residualImpl(const IModel& model, double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res, RowIteratorType /*jacBegin*/)
{

	for (unsigned int comp = 0; comp < _nComp; comp++)
	{
		StateType const* yBulk = y + offsetC() + comp;
		ResidualType* resBulk = res + offsetC() + comp;

		// Initialise residual with dc/dt
		if (yDot)
		{
			for (unsigned int i = 0; i < _nPoints; i++)
				resBulk[i * _strideNode] = static_cast<ResidualType>(yDot[offsetC() + comp + i * _strideNode]);
		}
		else
		{
			for (unsigned int i = 0; i < _nPoints; i++)
				resBulk[i * _strideNode] = ResidualType(0.0);
		}

		// Copy over inlet ghost node
		_inletC = y[comp];

		// ----------------------------------------------------------------
		// Step 1: auxiliary variable g = (2/deltaX) * dc/dx
		// ----------------------------------------------------------------

		StateType* g = reinterpret_cast<StateType*>(_auxState);
		int strideNode_g = 1u;
		int strideElem_g = _nNodes;

		for (int i = 0; i < _nPoints; i++)
			g[i] = StateType(0.0);

		Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> yMap(yBulk, _nPoints, InnerStride<Dynamic>(_strideNode));
		Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> gMap(g, _nPoints, InnerStride<>(strideNode_g));

		// g = - D c
		volumeIntegralAuxImpl<StateType, StateType>(yMap, gMap);

		// c^{*,d} = 0.5 (c^- + c^+)
		InterfaceFluxAuxiliaryImpl<StateType>(yBulk, _strideNode, _strideElem);

		// g = - D c - M^{(0,0),-1} * B * [c - c*]
		surfaceIntegralAuxImpl<StateType, StateType>(yBulk, g, _strideNode, _strideElem, strideNode_g, strideElem_g);

		// g = 2/dx_i * (- D c - M^{(0,0),-1} * B * [c - c*])
		gMap *= static_cast<StateType>(2.0 / _deltaX);

		// ----------------------------------------------------------------
		// Step 2: numerical fluxes c^{*,a} and g^{*,d}
		// ----------------------------------------------------------------

		// c^{*,a} = c^-, g^{*,d} = 0.5 (g^- + g^+)
		computeNumericalFluxes<StateType>(yBulk);

		// ----------------------------------------------------------------
		// Step 3: main equation volume and surface integrals
		// ----------------------------------------------------------------

		surfaceIntegralMainImpl<StateType, ResidualType>(resBulk, comp);

		volumeIntegralMainImpl<StateType, ResidualType>(yBulk, resBulk, g, comp);
	}

	return 0;
}

bool VariableCrossSectionConvectionDispersionOperatorBaseDG::setParameter(const ParameterId& pId, double value)
{
	if (_dispersionDep && _dispersionDep->hasParameter(pId))
	{
		_dispersionDep->setParameter(pId, value);
		return true;
	}

	if (!_dispersionCompIndep)
		return false;

	if ((pId.name != hashString("COL_DISPERSION")) || (pId.component != CompIndep) ||
		(pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) ||
		(pId.particleType != ParTypeIndep))
		return false;

	if (_colDispersion.size() > _nComp)
	{
		if (pId.section == SectionIndep) return false;
		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[pId.section * _nComp + i].setValue(value);
	}
	else
	{
		if (pId.section != SectionIndep) return false;
		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[i].setValue(value);
	}
	return true;
}

bool VariableCrossSectionConvectionDispersionOperatorBaseDG::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
{
	if (_dispersionDep)
	{
		active* const param = _dispersionDep->getParameter(pId);
		if (param) { param->setValue(value); return true; }
				}

	if (!_dispersionCompIndep)
		return false;

	if ((pId.name != hashString("COL_DISPERSION")) || (pId.component != CompIndep) ||
		(pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) ||
		(pId.particleType != ParTypeIndep))
		return false;

	if (_colDispersion.size() > _nComp)
	{
		if (pId.section == SectionIndep) return false;
		if (!contains(sensParams, &_colDispersion[pId.section * _nComp])) return false;
		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[pId.section * _nComp + i].setValue(value);
	}
	else
	{
		if (pId.section != SectionIndep) return false;
		if (!contains(sensParams, &_colDispersion[0])) return false;
		for (unsigned int i = 0; i < _nComp; ++i)
			_colDispersion[i].setValue(value);
		}
	return true;
	}

bool VariableCrossSectionConvectionDispersionOperatorBaseDG::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
{
	// No geometric and axial diffusion sensitivities available since that would dispatch active type DG operators which is currently not implemented.
	// Only report and reject pId if it actually names one of these parameters; any other parameter must fall through so it can be handled elsewhere.
	auto checkAvailabilityOfParameter = [&](const std::string& pName) -> bool
		{
			if (pId.name != hashStringRuntime(pName))
				return false;

			if (_geometryType == GeometryType::AxialFlowCylinder)
				LOG(Error) << "Sensitivities for " << pName <<
					" not supported for exact integration DG, please change to collocation DG or FV discretization.";
			else
				LOG(Error) << "Sensitivities for " << pName <<
					" in variable cross section column geometries are not supported.";

			return true;
		};

	if (checkAvailabilityOfParameter("COL_DISPERSION")) return false;
	if (checkAvailabilityOfParameter("CROSS_SECTION_AREA_LARGE_END")) return false;
	if (checkAvailabilityOfParameter("CROSS_SECTION_AREA_SMALL_END")) return false;
	if (checkAvailabilityOfParameter("BED_LENGTH")) return false;
	if (checkAvailabilityOfParameter("CYLINDER_HEIGHT")) return false;
	if (checkAvailabilityOfParameter("CROSS_SECTION_AREA_INNER")) return false;
	if (checkAvailabilityOfParameter("CROSS_SECTION_AREA_OUTER")) return false;
	if (checkAvailabilityOfParameter("CROSS_SECTION_AREA")) return false;

	return false;
}

template int AxialConvectionDispersionOperatorBaseCollocationDG::calcTransportJacobian<double>(const IModel&, double t, unsigned int secIdx, Eigen::SparseMatrix<double, RowMajor>&, Eigen::MatrixXd&, const int, const double* const);
template int AxialConvectionDispersionOperatorBaseCollocationDG::calcTransportJacobian<active>(const IModel&, double t, unsigned int secIdx, Eigen::SparseMatrix<double, RowMajor>&, Eigen::MatrixXd&, const int, const active* const);
template int VariableCrossSectionConvectionDispersionOperatorBaseDG::calcTransportJacobian<double>(const IModel&, double, unsigned int, Eigen::SparseMatrix<double, RowMajor>&, Eigen::MatrixXd&, const int, const double* const);
template int VariableCrossSectionConvectionDispersionOperatorBaseDG::calcTransportJacobian<active>(const IModel&, double, unsigned int, Eigen::SparseMatrix<double, RowMajor>&, Eigen::MatrixXd&, const int, const active* const);

}  // namespace parts

}  // namespace model

}  // namespace cadet

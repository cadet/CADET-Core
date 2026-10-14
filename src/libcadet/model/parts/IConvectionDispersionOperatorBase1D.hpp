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

	/**
	* @file
	* Defines the interface for one-dimensional convection dispersion transport operators.
	*/

	#ifndef LIBCADET_ICONVECTIONDISPERSIONOPERATOR1D_HPP_
	#define LIBCADET_ICONVECTIONDISPERSIONOPERATOR1D_HPP_

	#include "ParamIdUtil.hpp"
	#include "AutoDiff.hpp"
	#include "SimulationTypes.hpp"

	#include <Eigen/Dense>
	#include <Eigen/Sparse>

	#include <unordered_map>
	#include <unordered_set>
	#include <vector>

	namespace cadet
	{

	class IParameterProvider;
	class IConfigHelper;
	struct SimulationTime;
	class IModel;

	namespace model
	{
	namespace parts
	{

		/**
		 * @brief Interface for one-dimensional convection-dispersion transport base operators
		 * @details Abstracts the common API shared by all 1D convection-dispersion operator
		 *          implementations (DG and FV, axial and radial). Covers configuration,
		 *          residual evaluation, Jacobian assembly, parameter handling, and solution export.
		 */
		class IConvectionDispersionOperatorBase1D
		{
		public:

			virtual ~IConvectionDispersionOperatorBase1D() CADET_NOEXCEPT = default;

			/**
			* @brief Configures the spatial discretization
			* @param [in] paramProvider Parameter provider
			* @param [in] helper Configuration helper
			* @param [in] nComp Number of components
			* @param [in] nElem Number of elements (DG) or cells (FV)
			* @param [in] polyDeg Polynomial degree; pass 0 for FV
			* @param [in] strideNode Node stride in the global state vector
			* @return @c true on success
			*/
			virtual bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, unsigned int nComp, unsigned int nElem, unsigned int polyDeg, unsigned int strideNode) = 0;

			/**
			* @brief Configures section-dependent model parameters
			* @param [in] unitOpIdx Index of the owning unit operation
			* @param [in] paramProvider Parameter provider
			* @param [in,out] parameters Map in which sensitive parameters are registered
			* @return @c true on success
			*/
			virtual bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters) = 0;

			/**
			* @brief Notifies the operator of a discontinuous section transition
			* @param [in] t Current time
			* @param [in] secIdx New section index
			* @param [in,out] jacInlet Dense Jacobian block connecting inlet DOFs to the first bulk cells
			* @return @c true if the Jacobian sparsity pattern has changed
			*/
			virtual bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, Eigen::MatrixXd& jacInlet) = 0;

			/**
			* @brief Sets the volumetric flow rates for the current time section
			* @param [in] in  Inlet volumetric flow rate
			* @param [in] out Outlet volumetric flow rate
			* @param [in] colPorosity Column porosity \f$ \varepsilon_c \f$
			*/
			virtual void setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT = 0;

			virtual int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, WithoutParamSensitivity) = 0;

			virtual int residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, WithParamSensitivity) = 0;

			virtual int residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithParamSensitivity) = 0;

			virtual int residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, WithoutParamSensitivity) = 0;

			/**
			* @brief Populates the sparsity pattern of the convection-dispersion Jacobian block
			* @param [in,out] tripletList Triplet list to append entries to
			* @param [in] bulkOffset Row/column offset to the first bulk DOF in the global state vector
			*/
			virtual void convDispJacPattern(std::vector<Eigen::Triplet<double>>& tripletList,
				int bulkOffset) const = 0;

			/**
			* @brief Returns the number of non-zero Jacobian entries contributed by this operator
			* @param [in] pureNNZ If @c true, returns only structural NNZ (dispersion stencil);
			*                     if @c false, includes convection entries (a subset of the dispersion stencil)
			*/
			virtual unsigned int nJacEntries(bool pureNNZ) const CADET_NOEXCEPT = 0;

			/**
			* @brief Analytically assembles the static transport Jacobian (double state)
			* @details Non-template virtual override; concrete implementations delegate to the
			*          corresponding template instantiation calcTransportJacobian<double>.
			* @param [in]     model     Owning model (used by some FV implementations)
			* @param [in]     t         Current time
			* @param [in]     secIdx    Current section index
			* @param [in,out] jacobian  Sparse row-major system Jacobian to fill
			* @param [in,out] jacInlet  Dense inlet Jacobian block
			* @param [in]     bulkOffset Offset to the first bulk DOF
			* @param [in]     y         Current state vector (double)
			* @return @c true if the filled pattern matches the stored sparsity pattern
			*/
			virtual int calcTransportJacobian(const IModel& model, double t, unsigned int secIdx, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, int bulkOffset, double const* y) = 0;

			/**
			* @brief Analytically assembles the static transport Jacobian (active/AD state)
			* @details Non-template virtual override; concrete implementations delegate to the
			*          corresponding template instantiation calcTransportJacobian<active>.
			*/
			virtual int calcTransportJacobian(const IModel& model, double t, unsigned int secIdx, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, int bulkOffset, active const* y) = 0;

			/**
			* @brief Multiplies the time-derivative Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$
			*        with the given vector
			* @param [in]  simTime Simulation time point and section index
			* @param [in]  sDot    Input vector
			* @param [out] ret     Result vector
			*/
			virtual void multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const = 0;

			/**
			* @brief Adds @f$ \alpha \frac{\partial F}{\partial \dot{y}} @f$ to the discretized Jacobian
			* @details Used during BDF time integration to assemble the full linear system Jacobian.
			* @param [in]     alpha       BDF time-integration factor
			* @param [in,out] jacDisc     Sparse Jacobian to which the time-derivative contribution is added
			* @param [in]     blockOffset Row/column offset to the bulk block in the global Jacobian
			*/
			virtual void addTimeDerivativeToJacobian(double alpha, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacDisc, unsigned int blockOffset) const = 0;

			/**
			* @brief Returns the relative coordinate in [0, 1] of a global discrete point
			* @param [in] idx Zero-based global discrete point (node) index
			*/
			virtual double relativeCoordinate(unsigned int idx) const CADET_NOEXCEPT = 0;

			/**
			* @brief Returns @c true if the current interstitial velocity is non-negative (forward flow)
			*/
			virtual bool forwardFlow() const CADET_NOEXCEPT = 0;

			/**
			* @brief Returns the number of AD seed directions required for full Jacobian computation
			*/
			virtual unsigned int requiredADdirs() const CADET_NOEXCEPT = 0;

			virtual bool hasSmoothnessIndicator() const CADET_NOEXCEPT = 0;
			virtual int writeSmoothnessIndicator(double* buffer) const CADET_NOEXCEPT = 0;

			/**
			* @brief Writes the physical primary coordinates (axial or radial) of all discrete points
			* @param [out] coords Output buffer of length nPoints
			* @return Number of coordinates written
			*/
			virtual int writeCoordinates(double* coords) const CADET_NOEXCEPT = 0;

			virtual bool setParameter(const ParameterId& pId, double value) = 0;

			virtual bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue) = 0;

			virtual bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& id, double value) = 0;
		};

	} // namespace parts
	} // namespace model
	} // namespace cadet

	#endif  // LIBCADET_ICONVECTIONDISPERSIONOPERATOR1D_HPP_
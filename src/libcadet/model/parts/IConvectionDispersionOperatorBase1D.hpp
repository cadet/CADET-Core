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

			/**
			* @brief Hands an axially varying bulk porosity profile to the operator
			* @details The transport operator sees the porosity only through the flow through area
			*          @f$ w(x) = \varepsilon_c(x) A(x) @f$. A constant porosity cancels from every
			*          ratio the operators are built from and is handed over through setFlowRates()
			*          instead; a profile has to be baked into the operators, which requires an
			*          operator whose weight is not a closed-form analytic expression. Operators
			*          that cannot represent a varying weight return @c false, which the owning
			*          unit operation turns into a configuration error.
			*          Has to be called after configureModelDiscretization() and before configure(),
			*          since configure() assembles the geometry dependent operators.
			* @param [in] porosity Porosity profile with one entry per discrete bulk point, owned by
			*             the caller, which has to outlive the operator. Pass @c nullptr for a
			*             constant porosity.
			* @param [in] nPoints Number of entries in @p porosity, i.e. @c 1 for a constant porosity
			*             and the number of discrete bulk points for a profile
			* @return @c true if the operator can represent the given profile
			*/
			virtual bool setPorosityProfile(active const* porosity, unsigned int nPoints)
			{
				// By default, an operator only supports a constant porosity
				return (porosity == nullptr) || (nPoints <= 1u);
			}

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

			/**
			* @brief Number of discrete points that the bulk side of the film transfer term is spread over
			* @details The weak form multiplies the bulk balance by the flow through area
			*          \f$ A^\varepsilon = \varepsilon^{\mathrm{b}} A \f$, so that the time derivative carries the
			*          weight \f$ A^\varepsilon \f$ while the film transfer source carries the weight
			*          \f$ A^\varepsilon \left(1 - \varepsilon^{\mathrm{b}}\right) k^{\mathrm{f}} / \varepsilon^{\mathrm{b}} \f$.
			*          The two are proportional, and the source is therefore integrated node locally, unless the
			*          porosity or the film transfer coefficient varies within an element.
			*          Returns one as long as this operator cannot produce anything but a node local source,
			*          whether because of its geometry or because its operators are not set up yet.
			* @return Upper bound on the number of discrete points of the stencil
			*/
			virtual int filmTransferStencilSize() const CADET_NOEXCEPT { return 1; }

			/**
			* @brief Assembles the bulk side operator of the film transfer term
			* @details Inverting the mass matrix of the time derivative leaves the element-wise operator
			*
			*              T^(e) = (M^{A^eps})^{-1} Mtilde^{A^eps (1 - eps) kf / eps},
			*
			*          which maps the nodal values of \f$ k^{\mathrm{f}} \left(c^{\mathrm{b}} - c^{\mathrm{p}}\right) \f$
			*          onto the bulk residual. The weights returned here are divided by the nodal film transfer
			*          coefficient, so that they apply to the film transfer flux itself and reduce to
			*          \f$ \left(1 - \varepsilon^{\mathrm{b}}\right) / \varepsilon^{\mathrm{b}} \f$ for a porosity and a
			*          film transfer coefficient that are constant on an element.
			* @param [in] colPorosity Bulk porosity at each of the nPoints discrete points
			* @param [in] filmDiffusion Film transfer coefficient at each discrete point and component, point major
			* @param [in] nComp Number of components
			* @param [out] weights nPoints * nComp * filmTransferStencilSize() weights, point major, then component, then stencil
			* @param [out] offsets nPoints * filmTransferStencilSize() signed offsets in discrete points
			* @return @c false if the source is node local for the given porosity and film transfer coefficient,
			*         in which case nothing is written and the caller keeps using the plain prefactor
			*/
			virtual bool assembleFilmTransferOperator(active const* const colPorosity, active const* const filmDiffusion, const unsigned int nComp, active* weights, int* offsets) const { return false; }

			virtual bool setParameter(const ParameterId& pId, double value) = 0;

			virtual bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue) = 0;

			virtual bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& id, double value) = 0;
		};

	} // namespace parts
	} // namespace model
	} // namespace cadet

	#endif  // LIBCADET_ICONVECTIONDISPERSIONOPERATOR1D_HPP_
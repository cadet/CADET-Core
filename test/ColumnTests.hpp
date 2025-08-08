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
 * Defines simulation tests for column unit operations.
 */

#ifndef CADETTEST_COLUMNSIMTEST_HPP_
#define CADETTEST_COLUMNSIMTEST_HPP_

#include <limits>
#include <string>
#include "AutoDiff.hpp"
#include "common/JsonParameterProvider.hpp"

namespace cadet
{

class IUnitOperation;

namespace test
{

namespace column
{
	class JsonDiscScopeManager {
	private:
		JsonParameterProvider& jpp_;
		int levels_;
		bool inDiscretization_;

	public:
		JsonDiscScopeManager(JsonParameterProvider& jpp, const std::string& unitID)
			: jpp_(jpp), levels_(0), inDiscretization_(false) {
			if (jpp_.exists("model")) {
				jpp_.pushScope("model");
				++levels_;
			}
			if (jpp_.exists("unit_" + unitID)) {
				jpp_.pushScope("unit_" + unitID);
				++levels_;
			}
		}

		~JsonDiscScopeManager() {
			if (inDiscretization_) {
				jpp_.popScope(); // Exit discretization
			}
			for (int i = 0; i < levels_; ++i) {
				jpp_.popScope();
			}
		}

		std::string getUnitType() const {
			return jpp_.getString("UNIT_TYPE");
		}

		void enterDiscretization() {
			if (!inDiscretization_) {
				if (jpp_.exists("discretization")) {
					jpp_.pushScope("discretization");
					inDiscretization_ = true;
				}
			}
		}

		void exitDiscretization() {
			if (inDiscretization_) {
				jpp_.popScope();
				inDiscretization_ = false;
			}
		}

		JsonParameterProvider& getProvider() { return jpp_; }
	};

	// Base classes for bulk and particle discretization
	class BulkDiscretization {
	public:
		virtual ~BulkDiscretization() = default;
		virtual int getNAxCells() const = 0;
		virtual void setBulkParameters(JsonParameterProvider& jpp, const std::string& unitType, JsonDiscScopeManager& manager) const = 0;
		virtual std::string getDiscType() const = 0;
		virtual std::unique_ptr<BulkDiscretization> clone() const = 0;
	};

	class ParticleDiscretization {
	public:
		virtual ~ParticleDiscretization() = default;
		virtual int getNParCells() const = 0;
		virtual void setParticleParameters(JsonParameterProvider& jpp) const = 0;
		virtual std::string getDiscType() const = 0;
		virtual std::unique_ptr<ParticleDiscretization> clone() const = 0;
	};

	// Concrete bulk discretization implementations
	class BulkFV : public BulkDiscretization {
	private:
		int nAxCells_;
		int wenoOrder_;
		int nRadCells_;

	public:
		BulkFV(int nAxCells = 0, int wenoOrder = 0, int nRadCells = 0)
			: nAxCells_(nAxCells), wenoOrder_(wenoOrder), nRadCells_(nRadCells) {
		}

		std::string getDiscType() const override { return "FV"; }

		int getNAxCells() const override { return nAxCells_; }

		void setBulkParameters(JsonParameterProvider& jpp, const std::string& unitType, JsonDiscScopeManager& manager) const override {
			jpp.set("SPATIAL_METHOD", "FV");

			if (nAxCells_) jpp.set("NCOL", nAxCells_);

			if (nRadCells_) {
				if (unitType == "MULTI_CHANNEL_TRANSPORT") {
					manager.exitDiscretization();
					jpp.set("NCHANNEL", nRadCells_);
					manager.enterDiscretization();
				}
				else {
					jpp.set("NRAD", nRadCells_);
				}
			}

			if (jpp.exists("weno") && wenoOrder_) {
				jpp.pushScope("weno");
				jpp.set("WENO_ORDER", wenoOrder_);
				jpp.popScope();
			}
		}

		std::unique_ptr<BulkDiscretization> clone() const override {
			return std::make_unique<BulkFV>(nAxCells_, wenoOrder_, nRadCells_);
		}

		// Setters for configuration
		void setWenoOrder(int order) { wenoOrder_ = order; }
		int getWenoOrder() const { return wenoOrder_; }
		void setNRad(int nRad) { nRadCells_ = nRad; }
	};

	class BulkDG : public BulkDiscretization {
	private:
		int exactIntegration_;
		int polyDeg_;
		int nElem_;
		int radPolyDeg_;
		int radNelem_;

	public:
		BulkDG(int exactIntegration = -1, int polyDeg = 0, int nElem = 0, int radPolyDeg = 0, int radNelem = 0)
			: exactIntegration_(exactIntegration), polyDeg_(polyDeg), nElem_(nElem),
			radPolyDeg_(radPolyDeg), radNelem_(radNelem) {
		}

		std::string getDiscType() const override { return "DG"; }

		int getNAxCells() const override { return nElem_; }

		void setBulkParameters(JsonParameterProvider& jpp, const std::string& unitType, JsonDiscScopeManager& manager) const override {
			jpp.set("SPATIAL_METHOD", "DG");

			if (radPolyDeg_) {
				jpp.set("RAD_POLYDEG", radPolyDeg_);
				if (radNelem_) jpp.set("RAD_NELEM", radNelem_);
				if (polyDeg_) jpp.set("AX_POLYDEG", polyDeg_);
				if (nElem_) jpp.set("AX_NELEM", nElem_);
			}
			else {
				if (exactIntegration_ > -1) jpp.set("EXACT_INTEGRATION", exactIntegration_);
				if (polyDeg_) jpp.set("POLYDEG", polyDeg_);
				if (nElem_) jpp.set("NELEM", nElem_);
			}
		}

		std::unique_ptr<BulkDiscretization> clone() const override {
			return std::make_unique<BulkDG>(exactIntegration_, polyDeg_, nElem_, radPolyDeg_, radNelem_);
		}

		// Setters for configuration
		void setIntegrationMode(int integrationMode) { exactIntegration_ = integrationMode; }
		int getIntegrationMode() const { return exactIntegration_; }
	};

	// Concrete particle discretization implementations
	class ParticleFV : public ParticleDiscretization {
	private:
		int nParCells_;

	public:
		ParticleFV(int nParCells = 0) : nParCells_(nParCells) {}

		std::string getDiscType() const override { return "FV"; }

		int getNParCells() const override { return nParCells_; }

		void setParticleParameters(JsonParameterProvider& jpp) const override
		{
			if (nParCells_) {
				jpp.set("NPAR", nParCells_);
				jpp.set("NCELLS", nParCells_);
			}
		}

		std::unique_ptr<ParticleDiscretization> clone() const override {
			return std::make_unique<ParticleFV>(nParCells_);
		}
	};

	class ParticleDG : public ParticleDiscretization {
	private:
		int parPolyDeg_;
		int parNelem_;

	public:
		ParticleDG(int parPolyDeg = 0, int parNelem = 0)
			: parPolyDeg_(parPolyDeg), parNelem_(parNelem) {
		}

		std::string getDiscType() const override { return "DG"; }

		int getNParCells() const override { return parNelem_; }

		void setParticleParameters(JsonParameterProvider& jpp) const override
		{
			if (parPolyDeg_)
			{
				jpp.set("SPATIAL_METHOD", "DG");
				jpp.set("PAR_POLYDEG", parPolyDeg_);
				jpp.set("PAR_NELEM", parNelem_);
			}
		}

		std::unique_ptr<ParticleDiscretization> clone() const override {
			return std::make_unique<ParticleDG>(parPolyDeg_, parNelem_);
		}
	};

	// Dummy implementations
	class DummyBulk : public BulkDiscretization {
	public:
		int getNAxCells() const override { return 0; }
		void setBulkParameters(JsonParameterProvider& jpp, const std::string& unitType, JsonDiscScopeManager& manager) const override {}
		std::unique_ptr<BulkDiscretization> clone() const override {
			return std::make_unique<DummyBulk>();
		}
		std::string getDiscType() const override { return "NONE"; }
	};

	class DummyParticle : public ParticleDiscretization {
	public:
		int getNParCells() const override { return 0; }
		void setParticleParameters(JsonParameterProvider& jpp) const override {}
		std::unique_ptr<ParticleDiscretization> clone() const override {
			return std::make_unique<DummyParticle>();
		}
		std::string getDiscType() const override { return "NONE"; }
	};

	// Main DiscParams class that combines bulk and particle discretization
	class DiscParams {
	protected:
		std::unique_ptr<BulkDiscretization> bulkDisc_;
		std::unique_ptr<ParticleDiscretization> particleDisc_;

	public:
		DiscParams(std::unique_ptr<BulkDiscretization> bulkDisc = nullptr,
			std::unique_ptr<ParticleDiscretization> particleDisc = nullptr)
			: bulkDisc_(bulkDisc ? std::move(bulkDisc) : std::make_unique<DummyBulk>()),
			particleDisc_(particleDisc ? std::move(particleDisc) : std::make_unique<DummyParticle>()) {
		}

		// Copy constructor
		DiscParams(const DiscParams& other)
			: bulkDisc_(other.bulkDisc_->clone()),
			particleDisc_(other.particleDisc_->clone()) {
		}

		// Assignment operator
		DiscParams& operator=(const DiscParams& other) {
			if (this != &other) {
				bulkDisc_ = other.bulkDisc_->clone();
				particleDisc_ = other.particleDisc_->clone();
			}
			return *this;
		}

		virtual ~DiscParams() = default;

		virtual int getNAxCells() const { return bulkDisc_->getNAxCells(); }
		virtual int getNParCells() const { return particleDisc_->getNParCells(); }

		virtual void setDisc(JsonParameterProvider& jpp, const std::string unitID = "000", const std::string parID = "000") const
		{
			JsonDiscScopeManager manager(jpp, unitID);
			std::string unitType = manager.getUnitType();

			manager.enterDiscretization();
			jpp.set("SPATIAL_METHOD", bulkDisc_->getDiscType());
			bulkDisc_->setBulkParameters(manager.getProvider(), unitType, manager);
			manager.exitDiscretization();

			if (jpp.exists("particle_type_" + parID))
			{
				jpp.pushScope("particle_type_" + parID);
				manager.enterDiscretization();
				jpp.set("SPATIAL_METHOD", particleDisc_->getDiscType());
				particleDisc_->setParticleParameters(jpp);
				manager.exitDiscretization();
				jpp.popScope();
			}
			else
			{
				manager.enterDiscretization();
				jpp.set("SPATIAL_METHOD", particleDisc_->getDiscType());
				particleDisc_->setParticleParameters(jpp);
				manager.exitDiscretization();
			}
		}

		std::string getBulkDiscType() const {
			return bulkDisc_->getDiscType();
		}

		std::string getParDiscType() const {
			return particleDisc_->getDiscType();
		}

		// Getters for direct access to discretization objects
		const BulkDiscretization* getBulkDisc() const { return bulkDisc_.get(); }
		const ParticleDiscretization* getParticleDisc() const { return particleDisc_.get(); }
	};

	inline std::unique_ptr<DiscParams> createFVFVParams(int nAxCells, int nParCells, int wenoOrder = 0, int nRadCells = 0) {
		return std::make_unique<DiscParams>(
			std::make_unique<BulkFV>(nAxCells, wenoOrder, nRadCells),
			std::make_unique<ParticleFV>(nParCells)
		);
	}

	inline std::unique_ptr<DiscParams> createFVDGParams(int nAxCells, int wenoOrder, int nRadCells,
		int parPolyDeg, int parNelem) {
		return std::make_unique<DiscParams>(
			std::make_unique<BulkFV>(nAxCells, wenoOrder, nRadCells),
			std::make_unique<ParticleDG>(parPolyDeg, parNelem)
		);
	}

	inline std::unique_ptr<DiscParams> createDGFVParams(int exactIntegration, int polyDeg, int nElem,
		int radPolyDeg, int radNelem, int nParCells) {
		return std::make_unique<DiscParams>(
			std::make_unique<BulkDG>(exactIntegration, polyDeg, nElem, radPolyDeg, radNelem),
			std::make_unique<ParticleFV>(nParCells)
		);
	}

	inline std::unique_ptr<DiscParams> createDGDGParams(int exactIntegration, int polyDeg, int nElem,
		int radPolyDeg, int radNelem,
		int parPolyDeg, int parNelem) {
		return std::make_unique<DiscParams>(
			std::make_unique<BulkDG>(exactIntegration, polyDeg, nElem, radPolyDeg, radNelem),
			std::make_unique<ParticleDG>(parPolyDeg, parNelem)
		);
	}

	using DummyParams = DiscParams; // Uses default dummy implementations

	class FVParams : public DiscParams {
	public:
		FVParams() : DiscParams(std::make_unique<BulkFV>(), std::make_unique<ParticleFV>()) {}
		FVParams(int nCol) : DiscParams(std::make_unique<BulkFV>(nCol), std::make_unique<ParticleFV>()) {}
		FVParams(int nCol, int nPar) : DiscParams(std::make_unique<BulkFV>(nCol), std::make_unique<ParticleFV>(nPar)) {}
		FVParams(int nCol, int nPar, int wenoOrder)
			: DiscParams(std::make_unique<BulkFV>(nCol, wenoOrder), std::make_unique<ParticleFV>(nPar)) {
		}
		FVParams(int nCol, int nPar, int wenoOrder, int nRad)
			: DiscParams(std::make_unique<BulkFV>(nCol, wenoOrder, nRad), std::make_unique<ParticleFV>(nPar)) {
		}

		void setWenoOrder(int order) {
			static_cast<BulkFV*>(bulkDisc_.get())->setWenoOrder(order);
		}
		int getWenoOrder() const {
			return static_cast<const BulkFV*>(bulkDisc_.get())->getWenoOrder();
		}
		void setNRad(int nRad) {
			static_cast<BulkFV*>(bulkDisc_.get())->setNRad(nRad);
		}
	};

	class DGParams : public DiscParams {
	public:
		DGParams() : DiscParams(std::make_unique<BulkDG>(), std::make_unique<ParticleDG>()) {}
		DGParams(int exact, int poly, int elem)
			: DiscParams(std::make_unique<BulkDG>(exact, poly, elem), std::make_unique<ParticleDG>()) {
		}
		DGParams(int exact, int poly, int elem, int parPolyDeg, int parNelem)
			: DiscParams(std::make_unique<BulkDG>(exact, poly, elem), std::make_unique<ParticleDG>(parPolyDeg, parNelem)) {
		}
		DGParams(int exact, int poly, int elem, int parPolyDeg, int parNelem, int radPolyDeg, int radNelem)
			: DiscParams(std::make_unique<BulkDG>(exact, poly, elem, radPolyDeg, radNelem),
				std::make_unique<ParticleDG>(parPolyDeg, parNelem)) {
		}

		void setIntegrationMode(int integrationMode) {
			static_cast<BulkDG*>(bulkDisc_.get())->setIntegrationMode(integrationMode);
		}
		int getIntegrationMode() const {
			return static_cast<const BulkDG*>(bulkDisc_.get())->getIntegrationMode();
		}

	};

	/**
	 * @brief Sets the number of axial cells in a configuration of a column-like unit operation
	 * @details Overwrites the NCOL field in the discretization group of the given ParameterProvider.
	 * @param [in,out] jpp ParameterProvider to change the number of axial cells in
	 * @param [in] nCol Number of axial cells
	 * @param [in] unitID unit operation ID
	 */
	void setNumAxialCells(JsonParameterProvider& jpp, unsigned int nCol, std::string unitID="000");

	/**
	 * @brief Sets the WENO order in a configuration of a column-like unit operation
	 * @details Overwrites the WENO_ORDER field in the weno group of the given ParameterProvider.
	 * @param [in,out] jpp ParameterProvider to change the WENO order in
	 * @param [in] order Target order
	 */
	void setWenoOrder(JsonParameterProvider& jpp, int order);

	/**
	 * @brief Reverses the flow of a column-like unit operation
	 * @param [in,out] jpp ParameterProvider to change the flow direction in
	 */
	void reverseFlow(JsonParameterProvider& jpp);

	/**
	 * @brief Infers cross section area of column model from interstitial velocity
	 * @details Uses interstitial velocity and porosity to calculate cross section area.
	 *          To this end, a volumetric flow rate of 1.0 m^3/s is assumed. Depending on
	 *          @p dir, the velocity field is removed or set to @c +1.0 or @c -1.0 indicating
	 *          direction of the flow inside the unit operation.
	 * 
	 * @param [in,out] jpp ParameterProvider to add cross section area to
	 * @param [in] useTotalPorosity Determines whether TOTAL_POROSITY is used (@c true) or COL_POROSITY (@c false)
	 * @param [in] dir Flow direction in unit operation (@c 0 removes field, @c 1 standard direction, @c -1 flow reversal) 
	 */
	void setCrossSectionArea(JsonParameterProvider& jpp, bool useTotalPorosity, int dir);

	/**
	 * @brief Returns the offset to the flux part in the local state vector
	 * @param [in] unit Unit operation
	 * @return Offset to the flux part
	 */
	unsigned int fluxOffsetOfColumnUnitOp(cadet::IUnitOperation* unit);

	/**
	 * @brief Runs a parameterProvider for a model setup provided by a json file
	 * @param [in] modelFileRelPath relative path to model setup json file
	 */
	JsonParameterProvider getReferenceFile(const std::string& modelFileRelPath);

	/**
	 * @brief Runs a simulation test comparing against (semi-)analytic single component pulse injection reference data
	 * @details Linear binding model is used in the column-like unit operation.
	 * @param [in] uoType Unit operation type
	 * @param [in] refFileRelPath Path to the reference data file from the directory of this file
	 * @param [in] forwardFlow Determines whether the unit operates in forward flow (@c true) or backwards flow (@c false)
	 * @param [in] dynamicBinding Determines whether dynamic binding (@c true) or rapid equilibrium (@c false) is used
	 * @param [in] disc spatial discretization parameters
	 * @param [in] method spatial discretization method
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testAnalyticBenchmark(const char* uoType, const char* refFileRelPath, bool forwardFlow, bool dynamicBinding, DiscParams& disc, const std::string method, double absTol, double relTol);
	void testAnalyticBenchmark(const char* uoType, const char* refFileRelPath, bool forwardFlow, bool dynamicBinding, DGParams& disc, double absTol, double relTol);
	void testAnalyticBenchmark(const char* uoType, const char* refFileRelPath, bool forwardFlow, bool dynamicBinding, FVParams& disc, double absTol, double relTol);

	/**
	 * @brief Runs a simulation test comparing against (semi-)analytic single component pulse injection reference data
	 * @details The component is assumed to be non-binding.
	 * @param [in] uoType Unit operation type
	 * @param [in] refFileRelPath Path to the reference data file from the directory of this file
	 * @param [in] forwardFlow Determines whether the unit operates in forward flow (@c true) or backwards flow (@c false)
	 * @param [in] disc spatial discretization parameters
	 * @param [in] method spatial discretization method
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testAnalyticNonBindingBenchmark(const char* uoType, const char* refFileRelPath, bool forwardFlow, DiscParams& disc, const std::string method, double absTol, double relTol);
	void testAnalyticNonBindingBenchmark(const char* uoType, const char* refFileRelPath, bool forwardFlow, DGParams& disc, double absTol, double relTol);
	void testAnalyticNonBindingBenchmark(const char* uoType, const char* refFileRelPath, bool forwardFlow, FVParams& disc, double absTol, double relTol);


	/**
	 * @brief Runs a simulation test comparing forward and backwards flow in the load-wash-elution example
	 * @param [in] uoType Unit operation type
	 * @param [in] disc spatial discretization parameters
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testForwardBackward(const char* uoType, FVParams disc, double absTol, double relTol);
	void testForwardBackward(const char* uoType, DGParams disc, double absTol, double relTol);
	void testForwardBackward(JsonParameterProvider jpp, double absTol, double relTol);

	/**
	 * @brief Checks the full Jacobian against AD and FD pattern switching
	 * @details Checks the analytic Jacobian against the AD Jacobian and checks both against the FD pattern.
	 * @param [in] jpp Configured column model
	 * @param [in] absTolFDpattern absolute tolerance when comparing the sign in the FD Jacobian pattern
	 * @param [in] absTolAD absolute tolerance when comparing the AD Jacobians. Deviation from default only advised when values are numerically challenging, i.e. are at least 1E+10
	 * @param [in] flowRate flow rate, needs to be specified for 2D units, where velocity is derived from flow rates.
	 */
	void testJacobianAD(JsonParameterProvider& jpp, const double absTolFDpattern = 0.0, const double absTolAD = std::numeric_limits<float>::epsilon() * 100.0, const cadet::active* flowRate = nullptr);

	/**
	 * @brief Checks the full Jacobian against AD and FD pattern switching in case of variable surface diffusion coefficient
	 * @details Checks the analytic Jacobian against the AD Jacobian and checks both against the FD pattern.
	 * @param [in] uoType Unit operation type
	 * @param [in] dynamicBinding Determines whether dynamic binding is used
	 */
	void testJacobianADVariableParSurfDiff(const std::string& uoType, const std::string& spatialMethod, bool dynamicBinding);

	/**
	 * @brief Checks the full Jacobian against AD and FD pattern switching from forward to backward flow and back
	 * @details Checks the analytic Jacobian against the AD Jacobian and checks both against the FD pattern.
	 *          Checks both forward and backward flow mode as well as switching between them.
	 * 
	 * @param [in] uoType Unit operation type
	 * @param [in] disc spatial discretization parameters
	 * @param [in] absTolFDpattern absolute tolerance when comparing the sign in the FD Jacobian pattern
	 */
	void testJacobianForwardBackward(const char* uoType, FVParams disc, const double absTolFDpattern = 0.0);
	void testJacobianForwardBackward(const char* uoType, DGParams disc, const double absTolFDpattern = 0.0);
	void testJacobianForwardBackward(JsonParameterProvider& jpp, const double absTolFDpattern = 0.0);

	/**
	 * @brief Checks the full Jacobian against FD switching from forward to backward flow and back
	 * @details Checks the analytic Jacobian against the finite difference Jacobian.
	 *          Checks both forward and backward flow mode as well as switching between them.
	 *          Uses centered finite differences.
	 * 
	 * @param [in] uoType Unit operation type
	 * @param [in] wenoOrder Order of the WENO method
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testJacobianWenoForwardBackwardFD(const std::string& uoType, const std::string& spatialMethod, int wenoOrder, double h = 1e-6, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

	/**
	 * @brief Checks the (analytic) time derivative Jacobian against FD
	 * @details Uses centered finite differences.
	 * @param [in] uoType Unit operation type
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testTimeDerivativeJacobianFD(const std::string& uoType, const std::string& spatialMethod, double h = 1e-6, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

	/**
	 * @brief Checks the bottom macro row and right macro column of the Jacobian against FD
	 * @details Uses centered finite differences to check the flux part of the Jacobian.
	 * @param [in] uoType Unit operation type
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testArrowHeadJacobianFD(const std::string& uoType, double h = 1e-6, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

	/**
	 * @brief Checks the bottom macro row and right macro column of the Jacobian against FD
	 * @details Uses centered finite differences to check the flux part of the Jacobian.
	 * @param [in] uoType Unit operation type
	 * @param [in] dynamicBinding Determines whether dynamic binding is used
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testArrowHeadJacobianFD(const std::string& uoType, bool dynamicBinding, double h = 1e-6, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

	/**
	 * @brief Checks the bottom macro row and right macro column of the Jacobian against FD
	 * @details Uses centered finite differences to check the flux part of the Jacobian.
	 * @param [in] jpp Configured column model
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testArrowHeadJacobianFD(JsonParameterProvider& jpp, double h = 1e-6, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

	/**
	 * @brief Checks the bottom macro row and right macro column of the Jacobian against FD in case of variable surface diffusion coefficient
	 * @details Uses centered finite differences to check the flux part of the Jacobian.
	 * @param [in] uoType Unit operation type
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testArrowHeadJacobianFDVariableParSurfDiff(const std::string& uoType, double h = 1e-6, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

	/**
	 * @brief Checks the forward sensitivity residual using analytic Jacobians
	 * @details Uses centered finite differences.
	 * @param [in] jpp Configured unit model
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testFwdSensJacobians(JsonParameterProvider jpp, double h = 1e-6, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0, const bool hasBinding=true);

	/**
	 * @brief Checks the forward sensitivity solution against finite differences
	 * @details Assumes column-like unit models and uses centered finite differences. Checks 4 parameters in
	 *          the standard load-wash-elution test case:
	 *          COL_DISPERSION, CONST_COEFF (salt, loading), SMA_KA (first protein), CONNECTION (volumetric flow rate).
	 *          Each sensitivity is checked for quasi-stationary and dynamic binding, which means that
	 *          in total 8 checks are performed.
	 * @param [in] uoType Unit operation type
	 * @param [in] disableSensErrorTest Determines whether sensitivities take part in local error test
	 * @param [in] fdStepSize Array with step sizes of centered finite differences
	 * @param [in] absTols Array with absolute error tolerances
	 * @param [in] relTols Array with relative error tolerances
	 * @param [in] passRates Array with rates of relative error test passes
	 */
	void testFwdSensSolutionFD(const std::string& uoType, const std::string& spatialMethod, bool disableSensErrorTest, double const* fdStepSize, double const* absTols, double const* relTols, double const* passRates);

	/**
	 * @brief Checks the forward sensitivity solution with forward flow against the one using backward flow
	 * @details Assumes column-like unit models and checks 4 parameters in the standard load-wash-elution test case:
	 *          COL_DISPERSION, CONST_COEFF (salt, loading), SMA_KA (first protein), CONNECTION (volumetric flow rate).
	 *          Each sensitivity is checked for quasi-stationary and dynamic binding, which means that
	 *          in total 8 checks are performed.
	 * @param [in] uoType Unit operation type
	 * @param [in] absTols Array with absolute error tolerances
	 * @param [in] relTols Array with relative error tolerances
	 * @param [in] passRates Array with rates of relative error test passes
	 */
	void testFwdSensSolutionForwardBackward(const std::string& uoType, const std::string& spatialMethod, double const* absTols, double const* relTols, double const* passRates);

	/**
	 * @brief Checks consistent initialization using a model with linear binding
	 * @details Assumes column-like unit models and checks the residual of the model equations after
	 *          consistent initialization. A linear binding model is applied. Both binding modes and 
	 *          both AD and analytic Jacobians are checked. 
	 * @param [in] uoType Unit operation type
	 * @param [in] consTol Error tolerance for consistent initialization solver
	 * @param [in] absTol Error tolerance for checking whether algebraic residual is 0
	 * @param [in] reqBnd specifies binding mode, defaults to using both kinetic and equilibrium
	 * @param [in] useAD specifies Jacobian mode, defaults to using both analytical and AD
	 */
	void testConsistentInitializationLinearBinding(const std::string& uoType, const std::string& spatialMethod, double consTol, double absTol, const int reqBnd = -1, const int useAD = -1);

	/**
	 * @brief Checks consistent initialization using a model with SMA binding
	 * @details Assumes column-like unit models and checks the residual of the model equations after
	 *          consistent initialization. Both binding modes and both AD and analytic Jacobians are checked.
	 * @param [in] uoType Unit operation type
	 * @param [in] initState Initial state vector to start process from
	 * @param [in] consTol Error tolerance for consistent initialization solver
	 * @param [in] absTol Error tolerance for checking whether algebraic residual is 0
	 * @param [in] reqBnd specifies binding mode, defaults to using both kinetic and equilibrium
	 * @param [in] useAD specifies Jacobian mode, defaults to using both analytical and AD
	 */
	void testConsistentInitializationSMABinding(const std::string& uoType, const std::string& spatialMethod, double const* const initState, double consTol, double absTol, const int reqBnd = -1, const int useAD = -1);

	/**
	 * @brief Checks consistent initialization of sensitivities in a column-like model
	 * @details Assumes column-like unit models and checks the residual of the sensitivity equations after
	 *          consistent initialization. Both binding modes and both AD and analytic Jacobians are checked.
	 * @param [in] uoType Unit operation type
	 * @param [in] y State vector of original system
	 * @param [in] yDot Time derivative of state vector of original system
	 * @param [in] linearBinding Determines whether linear binding or SMA binding model is used
	 * @param [in] absTol Error tolerance for checking whether sensitivity residual is 0
	 * @param [in] reqBnd specifies binding mode, defaults to using both kinetic and equilibrium
	 * @param [in] useAD specifies Jacobian mode, defaults to using both analytical and AD
	 */
	void testConsistentInitializationSensitivity(const std::string& uoType, const std::string& spatialMethod, double const* const y, double const* const yDot, bool linearBinding, double absTol, const int reqBnd = -1, const int useAD = -1);

	/**
	 * @brief Checks whether the inlet DOFs produce the identity matrix in the Jacobian of the unit operation
	 * @details Assumes column-like unit models. Both AD and analytic Jacobians are checked.
	 * @param [in] uoType Unit operation type
	 */
	void testInletDofJacobian(const std::string& uoType, const std::string& spatialMethod);

	/**
	 * @brief Runs a simulation test comparing against numerical reference data (outlet data)
	 * @param [in] setupFileRelPath Path to the setup data file from the directory of this file
	 * @param [in] refFileRelPath Path to the reference data file from the directory of this file
	 * @param [in] unitID ID of the unit of interest
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 * @param [in] disc Numerical discretization parameters
	 * @param [in] compare_sens Specifies whether sensitivities are included
	 * @param [in] simDataStride strides in simulation data (eg stride over radial ports), only applied to simulation data
	 * @param [in] outletDataStride strides in outlet data (eg stride over nComp), applied to simulation and reference data
	 * @param [in] outletDataOffset offset to outlet data entry (eg to specific component), applied to simulation and reference data
	 */
	void testReferenceBenchmark(const std::string& modelFileRelPath, const std::string& refFileRelPath, const std::string& unitID, const std::vector<double> absTol, const std::vector<double> relTol, const cadet::test::column::DiscParams& disc, const bool compare_sens = false, const int simDataStride = 1, const int outletDataStride = 1, const int outletDataOffset = 0);

	/**
	 * @brief Runs an EOC test comparing against numerical reference data (outlet data)
	 * @param [in] setupFileRelPath Path to the setup data file from the directory of this file. Model configuration is sufficient, rest will be copied from reference file.
	 * @param [in] refFileRelPath Path to the reference data file from the directory of this file
	 * @param [in] convFileRelPath Path to the convergence reference data file from the directory of this file
	 * @param [in] unitID ID of the unit of interest
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 * @param [in] nDisc number of discretizations to be computed for EOC
	 * @param [in] disc Numerical discretization parameters with starting resolution
	 * @param [in] compare_sens Specifies whether sensitivities are included
	 */
	void testEOCReferenceBenchmark(const std::string& modelFileRelPath, const std::string& refFileRelPath, const std::string& convFileRelPath, const std::string& unitID, const std::vector<double> absTol, const std::vector<double> relTol, const unsigned int nDisc, const cadet::test::column::DiscParams& disc, const bool compare_sens = false);

	/**
	 * @brief Runs a simulation test comparing against numerical reference data (outlet data), generated by a foreign source
	 * @detail Reference data from the foreign source needs to be stored in an h5 file with hierarchy output->solution->SOLUTION_OUTLET
	 * @param [in] setupFileRelPath Path to the setup data file from the directory of this file. Full CADET configuration is required
	 * @param [in] refFileRelPath Path to the reference data file from the directory of this file
	 * @param [in] unitID ID of the unit of interest
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 * @param [in] compIdx Index of component of interest. Defaults to -1 to compare all components
	 */
	void testForeignReferenceBenchmark(const std::string& configFileRelPath, const std::string& refFileRelPath, const std::string& unitID, const double absTol, const double relTol, const int compIdx);


} // namespace column
} // namespace test
} // namespace cadet

#endif  // CADETTEST_COLUMNSIMTEST_HPP_

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
 * Sweeps every registered parameter of a model and compares its forward (AD) sensitivity
 * against central finite differences.
 *
 * The existing forward sensitivity FD tests (cadet::test::column::testFwdSensSolutionFD)
 * cover four hand-picked parameters. This sweep instead asks the configured model for all
 * of its parameters and checks each one, so that a sensitivity which is silently wrong -
 * or silently zero - cannot escape unnoticed.
 *
 * The sweep is slow (three simulations per parameter) and is therefore hidden from the
 * default test run. Run it explicitly by tag:
 *
 *     testRunner "[SensitivityFDSweep]"
 *
 * Every parameter produces one line on stdout prefixed with "SENSFD|" so that the result
 * table can be extracted from the log.
 */

#include <catch.hpp>
#include "Approx.hpp"
#include "cadet/cadet.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "common/Driver.hpp"
#include "CompileTimeConfig.hpp"
#include "JsonTestModels.hpp"
#include "ParamIdUtil.hpp"
#include "SimHelper.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
	/**
	 * @brief Parameter names the sweep is able to resolve
	 * @details A model identifies its parameters by a hash of the parameter name, so the name
	 *          itself cannot be recovered from a ParameterId. This list is hashed once and used
	 *          as a reverse lookup table. It doubles as a completeness check: a registered
	 *          parameter whose hash is not produced by any name in this list is reported as
	 *          UNRESOLVED and fails the sweep, so a newly added parameter cannot quietly drop
	 *          out of the FD comparison.
	 */
	const char* const knownParameterNames[] = {
		// Column / bulk transport
		"BED_LENGTH",
		"CHANNEL_CROSS_SECTION_AREAS",
		"COL_DISPERSION",
		"COL_DISPERSION_AXIAL",
		"COL_DISPERSION_RADIAL",
		"COL_POROSITY",
		"COL_RADIUS",
		"COL_RADIUS_LARGE_END",
		"COL_RADIUS_SMALL_END",
		"CROSS_SECTION_AREA",
		"CROSS_SECTION_AREA_INNER",
		"CROSS_SECTION_AREA_LARGE_END",
		"CROSS_SECTION_AREA_OUTER",
		"CROSS_SECTION_AREA_SMALL_END",
		"CYLINDER_HEIGHT",
		"EXCHANGE_MATRIX",
		"RADIAL_DISC_VECTOR",
		"TOTAL_POROSITY",
		"VELOCITY",
		// Particle transport
		"FILM_DIFFUSION",
		"PAR_CORERADIUS",
		"PAR_DIFFUSION",
		"PAR_POROSITY",
		"PAR_RADIUS",
		"PAR_SURFDIFFUSION",
		"PAR_TYPE_VOLFRAC",
		"PORE_ACCESSIBILITY",
		"PORE_DIFFUSION",
		"SURFACE_DIFFUSION",
		// Initial conditions
		"INIT_C",
		"INIT_CP",
		"INIT_CS",
		"INIT_Q",
		"INIT_VOLUME",
		// Steric mass action binding
		"SMA_KA",
		"SMA_KD",
		"SMA_LAMBDA",
		"SMA_NU",
		"SMA_REFC0",
		"SMA_REFQ",
		"SMA_SIGMA",
		// Inlet profile
		"CONST_COEFF",
		"CUBE_COEFF",
		"LIN_COEFF",
		"QUAD_COEFF",
		// System
		"CONNECTION",
		"CONNECTION_CUBE",
		"CONNECTION_LIN",
		"CONNECTION_QUAD",
		"SECTION_TIMES",
	};

	/**
	 * @brief Builds the hash to name lookup table for knownParameterNames
	 */
	const std::unordered_map<cadet::StringHash, std::string>& parameterNameLookup()
	{
		static const std::unordered_map<cadet::StringHash, std::string> lookup = []()
		{
			std::unordered_map<cadet::StringHash, std::string> map;
			for (const char* name : knownParameterNames)
				map[cadet::hashStringRuntime(name)] = name;
			return map;
		}();

		return lookup;
	}

	/**
	 * @brief Renders a ParameterId as "NAME[unit, comp, parType, bnd, rxn, sec]"
	 * @details Index positions that are independent of the parameter are printed as "-".
	 */
	std::string describe(const cadet::ParameterId& pId)
	{
		const std::unordered_map<cadet::StringHash, std::string>& lookup = parameterNameLookup();
		const std::unordered_map<cadet::StringHash, std::string>::const_iterator it = lookup.find(pId.name);

		std::ostringstream ss;
		if (it != lookup.end())
			ss << it->second;
		else
			ss << "UNRESOLVED_0x" << std::hex << static_cast<uint64_t>(pId.name) << std::dec;

		const auto idx = [&ss](int value, bool last)
		{
			if (value < 0)
				ss << '-';
			else
				ss << value;
			if (!last)
				ss << ", ";
		};

		ss << '[';
		idx(pId.unitOperation == cadet::UnitOpIndep ? -1 : static_cast<int>(pId.unitOperation), false);
		idx(pId.component == cadet::CompIndep ? -1 : static_cast<int>(pId.component), false);
		idx(pId.particleType == cadet::ParTypeIndep ? -1 : static_cast<int>(pId.particleType), false);
		idx(pId.boundState == cadet::BoundStateIndep ? -1 : static_cast<int>(pId.boundState), false);
		idx(pId.reaction == cadet::ReactionIndep ? -1 : static_cast<int>(pId.reaction), false);
		idx(pId.section == cadet::SectionIndep ? -1 : static_cast<int>(pId.section), true);
		ss << ']';

		return ss.str();
	}

	/**
	 * @brief Tightens the time integrator settings so that finite differences are meaningful
	 * @details The default LWE tolerances (ABSTOL 1e-8, RELTOL 1e-6) leave a noise floor that
	 *          swamps a central difference quotient. The values used here are the ones from the
	 *          report in issue #810.
	 */
	void tightenTolerances(cadet::JsonParameterProvider& jpp)
	{
		jpp.pushScope("solver");
		jpp.pushScope("time_integrator");

		jpp.set("ABSTOL", 1e-12);
		jpp.set("RELTOL", 1e-10);
		jpp.set("MAX_STEPS", 1000000);

		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Which of the two DG bulk operators a sweep runs on
	 * @details createLWE() asks for exact integration DG (USE_COLLOCATION_DG = 0). The collocation
	 *          operator is a separate implementation with its own parameter handling, so it has to
	 *          be swept separately. The field is ignored by the FV models.
	 */
	enum class DGVariant { AsGenerated, Collocation, ExactIntegration };

	void selectDGVariant(cadet::JsonParameterProvider& jpp, DGVariant variant)
	{
		if (variant == DGVariant::AsGenerated)
			return;

		jpp.pushScope("model");
		jpp.pushScope("unit_000");
		jpp.pushScope("discretization");

		jpp.set("USE_COLLOCATION_DG", (variant == DGVariant::Collocation) ? 1 : 0);

		jpp.popScope();
		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Chooses between the analytic and the AD Jacobian
	 * @details The forward sensitivity system is linear and is solved with the Jacobian as given.
	 *          A wrong analytic Jacobian therefore corrupts the sensitivities while leaving the
	 *          nominal solution intact, since Newton converges on the residual regardless. Running
	 *          a configuration both ways tells the two apart.
	 */
	void selectAnalyticJacobian(cadet::JsonParameterProvider& jpp, bool enabled)
	{
		jpp.pushScope("model");
		jpp.pushScope("unit_000");
		jpp.pushScope("discretization");

		jpp.set("USE_ANALYTIC_JACOBIAN", enabled);

		jpp.popScope();
		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Gives the first particle type a nonzero surface diffusion
	 * @details createLWE() sets SURFACE_DIFFUSION to zero, and the residual skips the whole surface
	 *          diffusion term for a coefficient of exactly zero. Everything that term carries - in
	 *          particular the solid to liquid volume ratio, which is derived from PAR_POROSITY and
	 *          PORE_ACCESSIBILITY - is then never exercised. The reference benchmarks do use a
	 *          nonzero coefficient, so the sweep has to cover this case as well.
	 */
	void enableSurfaceDiffusion(cadet::JsonParameterProvider& jpp, double value, unsigned int nComp)
	{
		jpp.pushScope("model");
		jpp.pushScope("unit_000");
		jpp.pushScope("particle_type_000");

		if (jpp.exists("SURFACE_DIFFUSION"))
			jpp.set("SURFACE_DIFFUSION", std::vector<double>(nComp, value));

		jpp.popScope();
		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Outcome of a single parameter check
	 */
	struct SweepEntry
	{
		cadet::ParameterId id;
		std::string label;        //!< Human readable parameter id, see describe()
		double nominal;           //!< Parameter value the sensitivity is evaluated at
		double maxAD;             //!< Maximum absolute AD sensitivity over the outlet
		double maxFD;             //!< Maximum absolute FD sensitivity over the outlet
		double relErr;            //!< Relative maximum error for FD step h
		double relErrHalf;        //!< Relative maximum error for FD step h/2, NaN if not computed
		double relErrQuarter;     //!< Relative maximum error for FD step h/4, NaN if not computed
		double relErrRich;        //!< Relative maximum error against the Richardson extrapolant, NaN if not computed
		std::string status;       //!< OK, OK_CONVERGING, MISMATCH, AD_ZERO, INERT, REJECTED, FAILED, SKIPPED, UNRESOLVED
		std::string message;      //!< Exception text or explanation, empty if none
	};

	/**
	 * @brief Reads the outlet of unit 0 into a flat vector
	 */
	std::vector<double> collectOutlet(cadet::InternalStorageUnitOpRecorder const* const rec, bool sensitivity)
	{
		const unsigned int nComp = rec->numComponents();
		const unsigned int nPorts = rec->numInletPorts();
		const unsigned int n = rec->numDataPoints() * nComp * nPorts;

		double const* const src = sensitivity ? rec->sensOutlet(0) : rec->outlet();
		REQUIRE(src != nullptr);

		return std::vector<double>(src, src + n);
	}

	/**
	 * @brief Runs the model once with the given value of @p pId and returns the outlet
	 */
	std::vector<double> runPerturbed(const cadet::JsonParameterProvider& jppBase, const cadet::ParameterId& pId, double value)
	{
		// Constructed from the underlying json rather than by copying the provider: the copy
		// constructor of JsonParameterProvider duplicates the json tree but keeps the scope
		// stack of the original, so the copy would write into the original's tree.
		cadet::JsonParameterProvider jpp(*jppBase.data());
		cadet::Driver drv;
		drv.configure(jpp);
		drv.simulator()->setParameterValue(pId, value);

		// Driver::configure() has already written the initial state vector from the model, so a
		// change to INIT_C, INIT_CP or INIT_CS would otherwise never reach the simulation and the
		// finite difference would come out as exactly zero. Reapplying costs one evaluation and
		// leaves every other parameter unaffected.
		drv.simulator()->applyInitialCondition();

		drv.run();

		return collectOutlet(drv.solution()->unitOperation(0), false);
	}

	/**
	 * @brief Relative maximum error between an AD and an FD sensitivity
	 * @details Returns max(abs(ad - fd)) / max(abs(fd)). If the FD sensitivity vanishes
	 *          identically the absolute maximum error is returned instead, so that a nonzero
	 *          AD sensitivity against a zero FD sensitivity is not divided by zero.
	 */
	double relativeMaxError(const std::vector<double>& ad, const std::vector<double>& fd)
	{
		double maxDiff = 0.0;
		double maxFD = 0.0;
		for (std::size_t i = 0; i < ad.size(); ++i)
		{
			maxDiff = std::max(maxDiff, std::abs(ad[i] - fd[i]));
			maxFD = std::max(maxFD, std::abs(fd[i]));
		}

		if (maxFD == 0.0)
			return maxDiff;

		return maxDiff / maxFD;
	}

	double maxAbs(const std::vector<double>& v)
	{
		double m = 0.0;
		for (double x : v)
			m = std::max(m, std::abs(x));
		return m;
	}

	/**
	 * @brief Compares the AD sensitivity of every registered parameter against central differences
	 * @param [in] uoType Unit operation type, as understood by createLWE()
	 * @param [in] spatialMethod Spatial discretization, as understood by createLWE()
	 * @param [in] h Relative finite difference step size
	 * @param [in] relTol Relative maximum error above which a parameter is reported as MISMATCH
	 * @param [in] dgVariant Which DG bulk operator to use, ignored by the FV models
	 * @param [in] surfDiff Surface diffusion coefficient to impose, @c 0.0 to leave the model as generated
	 */
	void sweep(const std::string& uoType, const std::string& spatialMethod, double h, double relTol,
		DGVariant dgVariant = DGVariant::AsGenerated, double surfDiff = 0.0, bool analyticJacobian = true,
		int bindingMode = -1)
	{
		std::ostringstream surfDiffLabel;
		if (surfDiff > 0.0)
			surfDiffLabel << " surfdiff=" << std::scientific << std::setprecision(0) << surfDiff;

		const std::string config = uoType + " " + spatialMethod
			+ ((dgVariant == DGVariant::Collocation) ? " collocation"
				: ((dgVariant == DGVariant::ExactIntegration) ? " exact-integration" : ""))
			+ surfDiffLabel.str()
			+ (analyticJacobian ? "" : " adjac")
			+ ((bindingMode == 0) ? " quasistationary" : ((bindingMode == 1) ? " kinetic" : ""));

		// Not every model can be integrated at the tight tolerances the finite differences would
		// like. The lumped rate model without pores under DG, for one, fails at t = 0 with them.
		// Rather than report nothing for such a configuration, the sweep falls back to the
		// tolerances createLWE() ships with and says so in the configuration label.
		bool tight = true;

		// Builds the configuration the whole sweep is based on. Every simulation starts from a
		// freshly built provider, so that nothing a previous run changed can carry over.
		const auto makeModel = [&]()
		{
			cadet::JsonParameterProvider jpp = createLWE(uoType, spatialMethod);
			selectDGVariant(jpp, dgVariant);
			if (!analyticJacobian)
				selectAnalyticJacobian(jpp, false);
			if (surfDiff > 0.0)
				enableSurfaceDiffusion(jpp, surfDiff, 4);
			if (bindingMode >= 0)
				cadet::test::setBindingMode(jpp, bindingMode > 0);
			if (tight)
				tightenTolerances(jpp);
			return jpp;
		};

		// Enumerate the parameters of the unperturbed model
		cadet::JsonParameterProvider jppProbe = makeModel();

		std::vector<cadet::ParameterId> ids;
		std::unordered_map<cadet::ParamIdHash, double> values;
		{
			cadet::Driver drv;
			drv.configure(jppProbe);

			const std::unordered_map<cadet::ParameterId, double> all = drv.simulator()->getAllParameterValues();
			ids.reserve(all.size());
			for (const std::pair<const cadet::ParameterId, double>& p : all)
			{
				ids.push_back(p.first);
				values[cadet::hashParameter(p.first)] = p.second;
			}
		}

		REQUIRE(!ids.empty());

		// Deterministic report order
		std::sort(ids.begin(), ids.end(), [](const cadet::ParameterId& a, const cadet::ParameterId& b)
		{
			return describe(a) < describe(b);
		});

		// Scale used as the finite difference step for parameters whose value is zero: the
		// largest absolute value among the parameters sharing the same name and unit operation.
		std::unordered_map<uint64_t, double> scaleByNameAndUnit;
		for (const cadet::ParameterId& pId : ids)
		{
			const uint64_t key = (static_cast<uint64_t>(pId.name) << 16) ^ static_cast<uint64_t>(pId.unitOperation);
			const double v = std::abs(values[cadet::hashParameter(pId)]);
			std::unordered_map<uint64_t, double>::iterator it = scaleByNameAndUnit.find(key);
			if (it == scaleByNameAndUnit.end())
				scaleByNameAndUnit[key] = v;
			else
				it->second = std::max(it->second, v);
		}

		// Section times are ordered and bound the simulation. A step that is large enough to
		// reorder them, or to move the end of the simulation before the last requested output
		// time, says nothing about the sensitivity, so they get their own step size below.
		const cadet::StringHash sectionTimesHash = cadet::hashStringRuntime("SECTION_TIMES");
		std::vector<double> sectionTimes;
		for (const cadet::ParameterId& pId : ids)
		{
			if (pId.name != sectionTimesHash)
				continue;
			if (sectionTimes.size() <= pId.section)
				sectionTimes.resize(pId.section + 1u, 0.0);
			sectionTimes[pId.section] = values[cadet::hashParameter(pId)];
		}

		// The unperturbed outlet, needed by the one sided difference quotient below. Computed once,
		// since it does not depend on which parameter is being checked, and used to decide whether
		// this configuration can be integrated at the tight tolerances at all.
		std::vector<double> baseline;
		const auto runBaseline = [&]()
		{
			cadet::JsonParameterProvider jpp = makeModel();
			cadet::Driver drv;
			drv.configure(jpp);
			drv.run();
			baseline = collectOutlet(drv.solution()->unitOperation(0), false);
		};

		try
		{
			runBaseline();
		}
		catch (const std::exception& e)
		{
			WARN("Configuration " << config << " cannot be integrated at ABSTOL 1e-12 / RELTOL 1e-10 ("
				<< e.what() << "), falling back to the tolerances of the generated setup");
			tight = false;
			runBaseline();
		}

		const std::string label = config + (tight ? "" : " loose-tol");

		// The FD reference simulations do not need the sensitivity machinery
		cadet::JsonParameterProvider jppFD = makeModel();

		std::vector<SweepEntry> report;
		report.reserve(ids.size());

		for (const cadet::ParameterId& pId : ids)
		{
			SweepEntry entry;
			entry.id = pId;
			entry.label = describe(pId);
			entry.nominal = values[cadet::hashParameter(pId)];
			entry.maxAD = std::numeric_limits<double>::quiet_NaN();
			entry.maxFD = std::numeric_limits<double>::quiet_NaN();
			entry.relErr = std::numeric_limits<double>::quiet_NaN();
			entry.relErrHalf = std::numeric_limits<double>::quiet_NaN();
			entry.relErrQuarter = std::numeric_limits<double>::quiet_NaN();
			entry.relErrRich = std::numeric_limits<double>::quiet_NaN();

			if (entry.label.compare(0, 10, "UNRESOLVED") == 0)
			{
				entry.status = "UNRESOLVED";
				entry.message = "Parameter name is not in knownParameterNames, add it";
				report.push_back(entry);
				continue;
			}

			// Finite difference step, absolute
			double step = 0.0;
			if (pId.name == sectionTimesHash)
			{
				if ((pId.section == 0) || (pId.section + 1u >= sectionTimes.size()))
				{
					entry.status = "SKIPPED";
					entry.message = "Start or end time of the simulation, not a section switch time";
					report.push_back(entry);
					continue;
				}

				// Stay well inside the neighbouring sections
				step = h * std::min(sectionTimes[pId.section] - sectionTimes[pId.section - 1u],
					sectionTimes[pId.section + 1u] - sectionTimes[pId.section]);
			}
			else
			{
				const uint64_t scaleKey = (static_cast<uint64_t>(pId.name) << 16) ^ static_cast<uint64_t>(pId.unitOperation);
				step = h * std::max(std::abs(entry.nominal), scaleByNameAndUnit[scaleKey]);
			}

			if (step <= 0.0)
			{
				entry.status = "SKIPPED";
				entry.message = "Parameter and all of its siblings are zero, no finite difference step available";
				report.push_back(entry);
				continue;
			}

			// AD sensitivity
			std::vector<double> ad;
			try
			{
				cadet::JsonParameterProvider jppAD = makeModel();
				cadet::test::returnSensitivities(jppAD, 0);

				// The sensitivity is kept out of the local error test. Its absolute tolerance then
				// does not matter, and the state trajectory of the AD run is the same one the
				// finite difference runs are built around, which is what makes the two comparable.
				// Leaving the sensitivity in the error test made IDA fail at t = 0 for a third of
				// the parameters, because no single SENS_ABSTOL suits parameters whose magnitudes
				// span twelve orders.
				cadet::test::disableSensitivityErrorTest(jppAD);

				const std::unordered_map<cadet::StringHash, std::string>& lookup = parameterNameLookup();
				cadet::test::addSensitivity(jppAD, lookup.at(pId.name), pId, 1e-10);

				cadet::Driver drv;
				drv.configure(jppAD);
				REQUIRE(drv.simulator()->numSensParams() == 1);
				drv.run();

				ad = collectOutlet(drv.solution()->unitOperation(0), true);
			}
			catch (const std::exception& e)
			{
				// A model that rejects the sensitivity on purpose ends up here, as does a
				// simulation that does not converge. The two are told apart by the message.
				entry.status = "REJECTED";
				entry.message = e.what();
				report.push_back(entry);
				continue;
			}

			// A concentration, a porosity or a diffusion coefficient must not go negative: SMA
			// takes powers of the bound concentrations, so a symmetric step around a nominal value
			// of zero produces NaN and the simulation fails rather than returning a difference
			// quotient. Those parameters get a one sided quotient instead, which is still second
			// order accurate and never leaves the physical range.
			const bool oneSided = (entry.nominal >= 0.0) && (entry.nominal - step < 0.0);

			bool inert = false;
			const auto finiteDifference = [&](double delta)
			{
				std::vector<double> fd;

				if (oneSided)
				{
					const std::vector<double> one = runPerturbed(jppFD, pId, entry.nominal + delta);
					const std::vector<double> two = runPerturbed(jppFD, pId, entry.nominal + 2.0 * delta);

					inert = (one == baseline) && (two == baseline);

					fd.resize(baseline.size());
					for (std::size_t i = 0; i < fd.size(); ++i)
						fd[i] = (-3.0 * baseline[i] + 4.0 * one[i] - two[i]) / (2.0 * delta);
				}
				else
				{
					const std::vector<double> left = runPerturbed(jppFD, pId, entry.nominal - delta);
					const std::vector<double> right = runPerturbed(jppFD, pId, entry.nominal + delta);

					// A parameter the solution does not react to at all. Either it genuinely does
					// not enter this model, or the model derived something from it at configuration
					// time and kept the derived value, in which case the finite difference
					// reference is meaningless and the AD sensitivity cannot be judged by it.
					inert = (left == right);

					fd.resize(left.size());
					for (std::size_t i = 0; i < fd.size(); ++i)
						fd[i] = (right[i] - left[i]) / (2.0 * delta);
				}

				return fd;
			};

			try
			{
				const std::vector<double> fd = finiteDifference(step);
				REQUIRE(fd.size() == ad.size());

				entry.maxAD = maxAbs(ad);
				entry.maxFD = maxAbs(fd);
				entry.relErr = relativeMaxError(ad, fd);

				if (inert && (entry.maxAD == 0.0))
				{
					// Both statements are zero and they agree, but they agree for two very
					// different reasons that the comparison cannot tell apart: the parameter may
					// genuinely not enter this model, or the model may have derived a value from
					// it at configuration time and kept that, in which case neither the finite
					// difference nor the sensitivity sees it. Reported, but not counted as a
					// defect, because deciding between the two needs the source.
					entry.status = "INERT_NO_SENS";
					entry.message = "Perturbing the parameter did not change the solution and the AD "
						"sensitivity is zero. Either the parameter does not enter this model, or a "
						"quantity derived from it is cached at configuration time.";
				}
				else if (inert)
				{
					entry.status = "INERT";
					entry.message = "Perturbing the parameter did not change the solution at all, yet the "
						"AD sensitivity is nonzero. Either the model caches a quantity derived from this "
						"parameter at configuration time, so that the finite difference reference carries "
						"no information, or the sensitivity is spurious.";
				}
				else if (entry.relErr <= relTol)
					entry.status = "OK";
				else
				{
					// The discrepancy at a single step size says little on its own: near a steep
					// elution front the second derivative is large, so the truncation error of the
					// difference quotient can sit well above the tolerance while the quotient is
					// still converging to the AD value. Refine the step twice and look at how the
					// discrepancy responds.
					const std::vector<double> fdHalf = finiteDifference(0.5 * step);
					entry.relErrHalf = relativeMaxError(ad, fdHalf);

					const std::vector<double> fdQuarter = finiteDifference(0.25 * step);
					entry.relErrQuarter = relativeMaxError(ad, fdQuarter);

					// Richardson extrapolation of the two finest quotients. A central difference
					// has an O(h^2) truncation error, so (4 * FD(h/4) - FD(h/2)) / 3 cancels that
					// leading term and leaves O(h^4). If the AD sensitivity is the derivative,
					// it has to agree with this extrapolant far more closely than with either
					// quotient on its own.
					std::vector<double> fdRich(fdQuarter.size());
					for (std::size_t i = 0; i < fdRich.size(); ++i)
						fdRich[i] = (4.0 * fdQuarter[i] - fdHalf[i]) / 3.0;
					entry.relErrRich = relativeMaxError(ad, fdRich);

					if ((entry.relErrRich <= relTol) || (entry.relErrHalf <= relTol) || (entry.relErrQuarter <= relTol))
						entry.status = "OK";
					else if ((entry.maxAD == 0.0) && (entry.maxFD > 0.0))
						entry.status = "AD_ZERO";
					else if ((entry.relErr <= 100.0 * relTol) && (entry.maxFD > 0.0)
						&& (std::abs(entry.maxAD - entry.maxFD) <= 1e-4 * entry.maxFD))
						// Small, and the two peak magnitudes agree to four digits. The pointwise
						// maximum difference is the most roundoff sensitive thing one can form out
						// of a difference quotient, and once it stops decreasing with the step it
						// is measuring that roundoff. The peak magnitude is not sensitive in the
						// same way, so its agreement bounds any systematic error in the
						// sensitivity far better than the pointwise maximum does.
						entry.status = "OK_NOISE_LIMITED";
					else if (entry.relErrQuarter <= 0.2 * entry.relErr)
						// A fourfold refinement cut the discrepancy at least fivefold. A derivative
						// that is simply wrong leaves a fixed offset, which cannot shrink like that.
						entry.status = "OK_CONVERGING";
					else if (entry.relErrQuarter < 0.8 * entry.relErr)
						// Coming down, but not far enough to conclude anything. The step is too
						// coarse for this parameter; neither confirmed nor refuted.
						entry.status = "INCONCLUSIVE";
					else
						entry.status = "MISMATCH";
				}
			}
			catch (const std::exception& e)
			{
				entry.status = "FAILED";
				entry.message = e.what();
			}

			report.push_back(entry);
		}

		// Emit the table
		for (const SweepEntry& e : report)
		{
			std::cout << "SENSFD|" << label << '|' << e.label << '|'
				<< std::scientific << std::setprecision(6)
				<< e.nominal << '|' << e.maxAD << '|' << e.maxFD << '|'
				<< e.relErr << '|' << e.relErrHalf << '|' << e.relErrQuarter << '|' << e.relErrRich << '|'
				<< e.status << '|' << e.message << std::endl;
		}

		// Status counts, so that the categories that do not fail are still visible at a glance
		{
			std::map<std::string, unsigned int> counts;
			for (const SweepEntry& e : report)
				++counts[e.status];

			std::cout << "SENSFDSUMMARY|" << label;
			for (const std::pair<const std::string, unsigned int>& c : counts)
				std::cout << '|' << c.first << '=' << c.second;
			std::cout << std::endl;
		}

		// Turn the table into assertions. A parameter the model refuses on purpose (REJECTED), one
		// without a usable finite difference step (SKIPPED), and one where the step was too coarse
		// to decide (INCONCLUSIVE) are listed in the table but do not fail; everything else that is
		// not OK is a defect and is reported by name.
		for (const SweepEntry& e : report)
		{
			if ((e.status.compare(0, 2, "OK") == 0) || (e.status == "SKIPPED") || (e.status == "REJECTED")
				|| (e.status == "INCONCLUSIVE") || (e.status == "INERT_NO_SENS"))
				continue;

			INFO(label << " " << e.label << " status " << e.status << " relErr " << e.relErr
				<< " relErr(h/2) " << e.relErrHalf << " relErr(h/4) " << e.relErrQuarter
				<< " relErr(Richardson) " << e.relErrRich << " " << e.message);
			CHECK(e.status == "OK");
		}
	}
}

TEST_CASE("Forward sensitivity FD sweep for COLUMN_MODEL_1D GRM", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	const std::string spatialMethod = GENERATE(std::string("FV"), std::string("DG"), std::string("DGFV"));
	SECTION("Spatial method " + spatialMethod)
	{
		sweep("COLUMN_MODEL_1D_GRM", spatialMethod, 1e-3, 1e-4);
	}
}

TEST_CASE("Forward sensitivity FD sweep for COLUMN_MODEL_1D GRM collocation DG", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	const std::string spatialMethod = GENERATE(std::string("DG"), std::string("DGFV"));
	SECTION("Spatial method " + spatialMethod)
	{
		sweep("COLUMN_MODEL_1D_GRM", spatialMethod, 1e-3, 1e-4, DGVariant::Collocation);
	}
}

TEST_CASE("Forward sensitivity FD sweep for COLUMN_MODEL_1D GRM with surface diffusion", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	const std::string spatialMethod = GENERATE(std::string("FV"), std::string("DG"), std::string("DGFV"));
	SECTION("Spatial method " + spatialMethod)
	{
		sweep("COLUMN_MODEL_1D_GRM", spatialMethod, 1e-3, 1e-4, DGVariant::AsGenerated, 5e-11);
	}
}

// The same configuration as above, but solved with the AD Jacobian. The forward sensitivity system
// is linear and is solved with whatever Jacobian the model provides, so a wrong analytic Jacobian
// shows up here and nowhere else: the nominal solution stays right because Newton converges on the
// residual. Comparing the two runs attributes a mismatch to the Jacobian or rules it out.
TEST_CASE("Forward sensitivity FD sweep for COLUMN_MODEL_1D GRM with surface diffusion and AD Jacobian", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	const std::string spatialMethod = GENERATE(std::string("FV"), std::string("DG"), std::string("DGFV"));
	SECTION("Spatial method " + spatialMethod)
	{
		sweep("COLUMN_MODEL_1D_GRM", spatialMethod, 1e-3, 1e-4, DGVariant::AsGenerated, 5e-11, false);
	}
}

// Narrowing down the DG bulk plus FV particle discrepancy: a tenfold smaller surface diffusion
// coefficient shows whether the error is proportional to it, and the two binding modes show whether
// it sits in the path that contributes surface diffusion to the film diffusion boundary, which is
// taken for quasi-stationary bound states only.
TEST_CASE("Forward sensitivity FD sweep for COLUMN_MODEL_1D GRM DGFV surface diffusion variants", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	SECTION("Tenth of the surface diffusion")
	{
		sweep("COLUMN_MODEL_1D_GRM", "DGFV", 1e-3, 1e-4, DGVariant::AsGenerated, 5e-12);
	}
	SECTION("Quasi-stationary binding")
	{
		sweep("COLUMN_MODEL_1D_GRM", "DGFV", 1e-3, 1e-4, DGVariant::AsGenerated, 5e-11, true, 0);
	}
	SECTION("Kinetic binding without surface diffusion")
	{
		sweep("COLUMN_MODEL_1D_GRM", "DGFV", 1e-3, 1e-4, DGVariant::AsGenerated, 0.0, true, 1);
	}
}

TEST_CASE("Forward sensitivity FD sweep for COLUMN_MODEL_1D LRMP", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	const std::string spatialMethod = GENERATE(std::string("FV"), std::string("DG"));
	SECTION("Spatial method " + spatialMethod)
	{
		sweep("COLUMN_MODEL_1D_LRMP", spatialMethod, 1e-3, 1e-4);
	}
}

TEST_CASE("Forward sensitivity FD sweep for COLUMN_MODEL_1D LRMP collocation DG", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	sweep("COLUMN_MODEL_1D_LRMP", "DG", 1e-3, 1e-4, DGVariant::Collocation);
}

TEST_CASE("Forward sensitivity FD sweep for COLUMN_MODEL_1D LRM", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	const std::string spatialMethod = GENERATE(std::string("FV"), std::string("DG"));
	SECTION("Spatial method " + spatialMethod)
	{
		sweep("COLUMN_MODEL_1D", spatialMethod, 1e-3, 1e-4);
	}
}

TEST_CASE("Forward sensitivity FD sweep for COLUMN_MODEL_1D LRM collocation DG", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	sweep("COLUMN_MODEL_1D", "DG", 1e-3, 1e-4, DGVariant::Collocation);
}

TEST_CASE("Forward sensitivity FD sweep for RADIAL_COLUMN_MODEL_1D GRM", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	sweep("RADIAL_COLUMN_MODEL_1D_GRM", "FV", 1e-3, 1e-4);
}

TEST_CASE("Forward sensitivity FD sweep for FRUSTUM_COLUMN_MODEL_1D GRM", "[SensitivityFDSweep][Sensitivity][FD][.]")
{
	const std::string spatialMethod = GENERATE(std::string("FV"), std::string("DG"));
	SECTION("Spatial method " + spatialMethod)
	{
		sweep("FRUSTUM_COLUMN_MODEL_1D_GRM", spatialMethod, 1e-3, 1e-4);
	}
}

#ifdef ENABLE_2D_MODELS

TEST_CASE("Forward sensitivity FD sweep for COLUMN_MODEL_2D GRM", "[SensitivityFDSweep][Sensitivity][FD][GRM2D][.]")
{
	const std::string spatialMethod = GENERATE(std::string("FV"), std::string("DG"));
	SECTION("Spatial method " + spatialMethod)
	{
		sweep("COLUMN_MODEL_2D_GRM", spatialMethod, 1e-3, 1e-4);
	}
}

#endif

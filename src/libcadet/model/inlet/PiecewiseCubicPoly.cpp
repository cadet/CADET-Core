// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2016: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Provides a piecewise cubic polynomial inlet profiles.
 */

#include "cadet/InletProfile.hpp"
#include "cadet/ParameterProvider.hpp"
#include "common/CompilerSpecific.hpp"

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <limits>
#include <string>
#include <sstream>
#include <iomanip>

namespace cadet
{

namespace model
{

/**
 * @brief A piecewise cubic polynomial inlet profile
 * @details On each section, a cubic polynomial is given as profile of each component.
 *          Since continuity or smoothness is not enforced, this ansatz is more general
 *          than a cubic spline.
 *          
 *          Note that the section times of the simulator have to match the breaks of the
 *          piecewise cubic polynomial.
 */
class PiecewiseCubicPolyInlet : public cadet::IInletProfile
{
public:
#if CADET_COMPILETIME_HASH
	PiecewiseCubicPolyInlet() { }
#else
	PiecewiseCubicPolyInlet() : _hashCons(hashString("CONST_COEFF")), _hashLin(hashString("LIN_COEFF")),
		_hashQuad(hashString("QUAD_COEFF")), _hashCub(hashString("CUBE_COEFF")), _hashSectionTimes(hashString("SECTION_TIMES"))
	{
	}
#endif

	virtual ~PiecewiseCubicPolyInlet() CADET_NOEXCEPT { }

	static const char* identifier() { return "PIECEWISE_CUBIC_POLY"; }
	virtual const char* name() const CADET_NOEXCEPT { return PiecewiseCubicPolyInlet::identifier(); }

	virtual std::vector<cadet::ParameterId> availableParameters(unsigned int unitOpIdx) CADET_NOEXCEPT
	{
		std::vector<cadet::ParameterId> params;
		params.reserve((_sectionTimes.size() - 1) * _nComp * 4);
		for (unsigned int sec = 0; sec < _sectionTimes.size() - 1; ++sec)
		{
			for (unsigned int comp = 0; comp < _nComp; ++comp)
			{
				params.push_back(cadet::makeParamId(_hashCons, unitOpIdx, comp, cadet::BoundPhaseIndep, cadet::ReactionIndep, sec));
				params.push_back(cadet::makeParamId(_hashLin, unitOpIdx, comp, cadet::BoundPhaseIndep, cadet::ReactionIndep, sec));
				params.push_back(cadet::makeParamId(_hashQuad, unitOpIdx, comp, cadet::BoundPhaseIndep, cadet::ReactionIndep, sec));
				params.push_back(cadet::makeParamId(_hashCub, unitOpIdx, comp, cadet::BoundPhaseIndep, cadet::ReactionIndep, sec));
			}
			params.push_back(cadet::makeParamId(_hashSectionTimes, cadet::UnitOpIndep, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, sec));
		}

		return params;
	}

	virtual void inletConcentration(double t, unsigned int sec, double* inletConc)
	{
		// This function evaluates a piecewise cubic polynomial given on some intervals
		// called sections. On each section a polynomial of degree 3 is evaluated:
		// 
		//   p_i(t) = CUBE_COEFF[i] * (t - t_i)^3 + QUAD_COEFF[i] * (t - t_i)^2 + LIN_COEFF[i] * (t - t_i) + CONST_COEFF[i],
		//   
		// where p_i is the polynomial on section i given by the interval [t_i, t_{i+1}].

		cadet_assert(sec < _sectionTimes.size());

		const double tShift = t - _sectionTimes[sec];
		double const* const con = _const.data() + sec * _nComp;
		double const* const lin = _lin.data() + sec * _nComp;
		double const* const quad = _quad.data() + sec * _nComp;
		double const* const cub = _cub.data() + sec * _nComp;

		// Evaluate polynomial using Horner's scheme
		for (unsigned int comp = 0; comp < _nComp; ++comp)
			inletConc[comp] = con[comp] + tShift * (lin[comp] + tShift * (quad[comp] + tShift * cub[comp]));
	}

	virtual void parameterDerivative(double t, unsigned int sec, const cadet::ParameterId& pId, double* paramDeriv)
	{
		cadet_assert(sec < _sectionTimes.size());

		std::fill(paramDeriv, paramDeriv + _nComp, 0.0);

		if (pId.section != sec)
			// Wrong section leads to zero derivative
			return;

		const double tShift = t - _sectionTimes[sec];

		// SECTION_TIMES is global and, thus, has no associated unitOp
		if ((pId.name == _hashSectionTimes) && (pId.unitOperation == cadet::UnitOpIndep) && (pId.reaction == cadet::ReactionIndep) &&
			(pId.component == cadet::CompIndep) && (pId.boundPhase == cadet::BoundPhaseIndep))
		{
			double const* const lin = _lin.data() + sec * _nComp;
			double const* const quad = _quad.data() + sec * _nComp;
			double const* const cub = _cub.data() + sec * _nComp;

			// p_i(t) = CUBE_COEFF[i] * (t - t_i)^3 + QUAD_COEFF[i] * (t - t_i)^2 + LIN_COEFF[i] * (t - t_i) + CONST_COEFF[i]
			// Now suppose t_i = t_i(q) and t = t(q), then
			// dp_i / dt_i = -3 * CUBE_COEFF[i] * (t - t_i)^2 - 2 * QUAD_COEFF[i] * (t - t_i) - LIN_COEFF[i]
			for (unsigned int i = 0; i < _nComp; ++i)
				paramDeriv[i] = -lin[i] - tShift * (2.0 * quad[i] + tShift * 3.0 * cub[i]);

			return;
		}

#if CADET_COMPILETIME_HASH
		// Assume that t' = dt / dq = 0 for all handled parameters q
		switch (pId.name)
		{
			case _hashCons:
				paramDeriv[pId.component] = 1.0;
				break;
			case _hashLin:
				paramDeriv[pId.component] = tShift;
				break;
			case _hashQuad:
				paramDeriv[pId.component] = tShift * tShift;
				break;
			case _hashCub:
				paramDeriv[pId.component] = tShift * tShift * tShift;
				break;
		}
#else
		// Assume that t' = dt / dq = 0 for all handled parameters q
		if (pId.name == _hashCons)
			paramDeriv[pId.component] = 1.0;
		else if (pId.name == _hashLin)
			paramDeriv[pId.component] = tShift;
		else if (pId.name == _hashQuad)
			paramDeriv[pId.component] = tShift * tShift;
		else if (pId.name == _hashCub)
			paramDeriv[pId.component] = tShift * tShift * tShift;
#endif
	}

	virtual void timeDerivative(double t, unsigned int sec, double* timeDerivative)
	{
		cadet_assert(sec < _sectionTimes.size());
		const double tShift = t - _sectionTimes[sec];
		double const* const lin = _lin.data() + sec * _nComp;
		double const* const quad = _quad.data() + sec * _nComp;
		double const* const cub = _cub.data() + sec * _nComp;

		// Evaluate polynomial using Horner's scheme
		for (unsigned int comp = 0; comp < _nComp; ++comp)
			timeDerivative[comp] = lin[comp] + tShift * (2.0 * quad[comp] + tShift * 3.0 * cub[comp]);
	}

	virtual void setParameterValue(const cadet::ParameterId& pId, double value)
	{
#if CADET_COMPILETIME_HASH
		switch (pId.name)
		{
			case _hashCons:
				_const[pId.component + pId.section * _nComp] = value;
				break;
			case _hashLin:
				_lin[pId.component + pId.section * _nComp] = value;
				break;
			case _hashQuad:
				_quad[pId.component + pId.section * _nComp] = value;
				break;
			case _hashCub:
				_cub[pId.component + pId.section * _nComp] = value;
				break;
			case _hashSectionTimes:
				_sectionTimes[pId.section] = value;
				break;
		}
#else
		if (pId.name == _hashCons)
			_const[pId.component + pId.section * _nComp] = value;
		else if (pId.name == _hashLin)
			_lin[pId.component + pId.section * _nComp] = value;
		else if (pId.name == _hashQuad)
			_quad[pId.component + pId.section * _nComp] = value;
		else if (pId.name == _hashCub)
			_cub[pId.component + pId.section * _nComp] = value;
		else if (pId.name == _hashSectionTimes)
			_sectionTimes[pId.section] = value;
#endif
	}

	virtual double getParameterValue(const cadet::ParameterId& pId)
	{
#if CADET_COMPILETIME_HASH
		switch (pId.name)
		{
			case _hashCons:
				return _const[pId.component + pId.section * _nComp];
			case _hashLin:
				return _lin[pId.component + pId.section * _nComp];
			case _hashQuad:
				return _quad[pId.component + pId.section * _nComp];
			case _hashCub:
				return _cub[pId.component + pId.section * _nComp];
			case _hashSectionTimes:
				return _sectionTimes[pId.section];
		}
#else
		if (pId.name == _hashCons)
			return _const[pId.component + pId.section * _nComp];
		else if (pId.name == _hashLin)
			return _lin[pId.component + pId.section * _nComp];
		else if (pId.name == _hashQuad)
			return _quad[pId.component + pId.section * _nComp];
		else if (pId.name == _hashCub)
			return _cub[pId.component + pId.section * _nComp];
		else if (pId.name == _hashSectionTimes)
			return _sectionTimes[pId.section];
#endif
		return std::numeric_limits<double>::quiet_NaN();
	}

	virtual void numComponents(unsigned int nComp) CADET_NOEXCEPT { _nComp = nComp; }

	inline const std::vector<double>& sectionTimes() const CADET_NOEXCEPT { return _sectionTimes; }
	inline void sectionTimes(const std::vector<double>& sectionTimes) CADET_NOEXCEPT { _sectionTimes = sectionTimes; }

	inline const std::vector<double>& constantCoeff() const CADET_NOEXCEPT { return _const; }
	inline void constantCoeff(const std::vector<double>& cons) CADET_NOEXCEPT { _const = cons; }

	inline const std::vector<double>& linearCoeff() const CADET_NOEXCEPT { return _lin; }
	inline void linearCoeff(const std::vector<double>& lin) CADET_NOEXCEPT { _lin = lin; }

	inline const std::vector<double>& quadraticCoeff() const CADET_NOEXCEPT { return _quad; }
	inline void quadraticCoeff(const std::vector<double>& quad) CADET_NOEXCEPT { _quad = quad; }

	inline const std::vector<double>& cubicCoeff() const CADET_NOEXCEPT { return _cub; }
	inline void cubicCoeff(const std::vector<double>& cub) CADET_NOEXCEPT { _cub = cub; }

	virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections) CADET_NOEXCEPT
	{
		_sectionTimes = std::vector<double>(secTimes, secTimes + nSections + 1);
	}

	virtual bool configure(IParameterProvider* paramProvider, unsigned int nComp)
	{
		_nComp = nComp;

		_const.clear();
		_lin.clear();
		_quad.clear();
		_cub.clear();		

		if (!paramProvider)
			return false;

		unsigned int i = 0;
		std::ostringstream oss;
		oss << "sec_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
		while (paramProvider->exists(oss.str()))
		{
			paramProvider->pushScope(oss.str());

			const std::vector<double> cons = paramProvider->getDoubleArray("CONST_COEFF");
			const std::vector<double> lin = paramProvider->getDoubleArray("LIN_COEFF");
			const std::vector<double> quad = paramProvider->getDoubleArray("QUAD_COEFF");
			const std::vector<double> cub = paramProvider->getDoubleArray("CUBE_COEFF");

			// Preallocate memory if possible
			if ((i == 0) && (_sectionTimes.size() > 0))
			{
				const unsigned int nComp = cons.size();

				_const.reserve(nComp * _sectionTimes.size());
				_lin.reserve(nComp * _sectionTimes.size());
				_quad.reserve(nComp * _sectionTimes.size());
				_cub.reserve(nComp * _sectionTimes.size());
			}

			_const.insert(_const.end(), cons.begin(), cons.end());
			_lin.insert(_lin.end(), lin.begin(), lin.end());
			_quad.insert(_quad.end(), quad.begin(), quad.end());
			_cub.insert(_cub.end(), cub.begin(), cub.end());

			paramProvider->popScope();

			++i;
			oss.str("");
			oss << "sec_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
		}
		return true;
	}

private:

#if CADET_COMPILETIME_HASH
	static const cadet::StringHash _hashCons = "CONST_COEFF"_hash;
	static const cadet::StringHash _hashLin = "LIN_COEFF"_hash;
	static const cadet::StringHash _hashQuad = "QUAD_COEFF"_hash;
	static const cadet::StringHash _hashCub = "CUBE_COEFF"_hash;
	static const cadet::StringHash _hashSectionTimes = "SECTION_TIMES"_hash;
#else
	const cadet::StringHash _hashCons;
	const cadet::StringHash _hashLin;
	const cadet::StringHash _hashQuad;
	const cadet::StringHash _hashCub;
	const cadet::StringHash _hashSectionTimes;
#endif

	std::vector<double> _sectionTimes;
	unsigned int _nComp;

	std::vector<double> _const;
	std::vector<double> _lin;
	std::vector<double> _quad;
	std::vector<double> _cub;
};

namespace inlet
{
	void registerPiecewiseCubicPoly(std::unordered_map<std::string, std::function<IInletProfile*()>>& inlets)
	{
		inlets[PiecewiseCubicPolyInlet::identifier()] = []() { return new PiecewiseCubicPolyInlet(); };
	}
} // namespace inlet

} // namespace model
} // namespace cadet

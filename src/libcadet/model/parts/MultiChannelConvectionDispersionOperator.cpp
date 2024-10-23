// =============================================================================
//  CADET
//
//  Copyright © The CADET Authors
//            Please see the CONTRIBUTORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/parts/MultiChannelConvectionDispersionOperator.hpp"
#include "cadet/Exceptions.hpp"

#include "Stencil.hpp"
#include "ParamReaderHelper.hpp"
#include "AdUtils.hpp"
#include "SensParamUtil.hpp"
#include "SimulationTypes.hpp"
#include "linalg/CompressedSparseMatrix.hpp"
#include "model/parts/AxialConvectionDispersionKernel.hpp"
#include "model/ParameterDependence.hpp"

#ifdef SUPERLU_FOUND
	#include "linalg/SuperLUSparseMatrix.hpp"
#endif
#ifdef UMFPACK_FOUND
	#include "linalg/UMFPackSparseMatrix.hpp"
#endif

#include "linalg/BandMatrix.hpp"
#include "linalg/Gmres.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>

namespace
{

cadet::model::MultiplexMode readAndRegisterMultiplexParam(cadet::IParameterProvider& paramProvider, std::unordered_map<cadet::ParameterId, cadet::active*>& parameters, std::vector<cadet::active>& values, const std::string& name, unsigned int nComp, unsigned int nChannel, cadet::UnitOpIdx uoi)
{
	cadet::model::MultiplexMode mode = cadet::model::MultiplexMode::Independent;
	readParameterMatrix(values, paramProvider, name, nComp * nChannel, 1);
	unsigned int nSec = 1;
	if (paramProvider.exists(name + "_MULTIPLEX"))
	{
		const int modeConfig = paramProvider.getInt(name + "_MULTIPLEX");
		if (modeConfig == 0)
		{
			mode = cadet::model::MultiplexMode::Independent;
			if (values.size() > 1)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be 1)");
		}
		else if (modeConfig == 1)
		{
			mode = cadet::model::MultiplexMode::Radial;
			if (values.size() != nChannel)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nChannel) + ")");
		}
		else if (modeConfig == 2)
		{
			mode = cadet::model::MultiplexMode::Component;
			if (values.size() != nComp)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nComp) + ")");
		}
		else if (modeConfig == 3)
		{
			mode = cadet::model::MultiplexMode::ComponentRadial;
			if (values.size() != nComp * nChannel)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nComp * nChannel) + ")");
		}
		else if (modeConfig == 4)
		{
			mode = cadet::model::MultiplexMode::Section;
			nSec = values.size();
		}
		else if (modeConfig == 5)
		{
			mode = cadet::model::MultiplexMode::RadialSection;
			if (values.size() % nChannel != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " is not a positive multiple of NCHANNEL (" + std::to_string(nChannel) + ")");

			nSec = values.size() / nChannel;
		}
		else if (modeConfig == 6)
		{
			mode = cadet::model::MultiplexMode::ComponentSection;
			if (values.size() % nComp != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " is not a positive multiple of NCOMP (" + std::to_string(nComp) + ")");

			nSec = values.size() / nComp;
		}
		else if (modeConfig == 7)
		{
			mode = cadet::model::MultiplexMode::ComponentRadialSection;
			if (values.size() % (nComp * nChannel) != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " is not a positive multiple of NCOMP * NCHANNEL (" + std::to_string(nComp * nChannel) + ")");

			nSec = values.size() / (nComp * nChannel);
		}
	}
	else
	{
		if (values.size() == 1)
			mode = cadet::model::MultiplexMode::Independent;
		else if (values.size() == nComp)
			mode = cadet::model::MultiplexMode::Component;
		else if (values.size() == nChannel)
			mode = cadet::model::MultiplexMode::Radial;
		else if (values.size() == nChannel * nComp)
			mode = cadet::model::MultiplexMode::ComponentRadial;
		else if (values.size() % nComp == 0)
		{
			mode = cadet::model::MultiplexMode::ComponentSection;
			nSec = values.size() / nComp;
		}
		else if (values.size() % nChannel == 0)
		{
			mode = cadet::model::MultiplexMode::RadialSection;
			nSec = values.size() / nChannel;
		}
		else if (values.size() % (nChannel * nComp) == 0)
		{
			mode = cadet::model::MultiplexMode::ComponentRadialSection;
			nSec = values.size() / (nComp * nChannel);
		}
		else
			throw cadet::InvalidParameterException("Could not infer multiplex mode of field " + name + ", set " + name + "_MULTIPLEX or change number of elements");

		// Do not infer cadet::model::MultiplexMode::Section in case of no matches (might hide specification errors)
	}

	const cadet::StringHash nameHash = cadet::hashStringRuntime(name);
	switch (mode)
	{
		case cadet::model::MultiplexMode::Independent:
		case cadet::model::MultiplexMode::Section:
			{
				std::vector<cadet::active> p(nComp * nChannel * nSec);
				for (unsigned int s = 0; s < nSec; ++s)
					std::fill(p.begin() + s * nChannel * nComp, p.begin() + (s+1) * nChannel * nComp, values[s]);

				values = std::move(p);

				for (unsigned int s = 0; s < nSec; ++s)
					parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, (mode == cadet::model::MultiplexMode::Independent) ? cadet::SectionIndep : s)] = &values[s * nChannel * nComp];
			}
			break;
		case cadet::model::MultiplexMode::Component:
		case cadet::model::MultiplexMode::ComponentSection:
			{
				std::vector<cadet::active> p(nComp * nChannel * nSec);
				for (unsigned int s = 0; s < nSec; ++s)
				{
					for (unsigned int i = 0; i < nComp; ++i)
						std::copy(values.begin() + s * nComp, values.begin() + (s+1) * nComp, p.begin() + i * nComp + s * nComp * nChannel);
				}

				values = std::move(p);

				for (unsigned int s = 0; s < nSec; ++s)
				{
					for (unsigned int i = 0; i < nComp; ++i)
						parameters[cadet::makeParamId(nameHash, uoi, i, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, (mode == cadet::model::MultiplexMode::Component) ? cadet::SectionIndep : s)] = &values[s * nChannel * nComp + i];
				}
			}
			break;
		case cadet::model::MultiplexMode::Radial:
		case cadet::model::MultiplexMode::RadialSection:
			{
				std::vector<cadet::active> p(nComp * nChannel * nSec);
				for (unsigned int i = 0; i < nChannel * nSec; ++i)
					std::fill(p.begin() + i * nComp, p.begin() + (i+1) * nComp, values[i]);

				values = std::move(p);

				for (unsigned int s = 0; s < nSec; ++s)
				{
					for (unsigned int i = 0; i < nChannel; ++i)
						parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, i, cadet::BoundStateIndep, cadet::ReactionIndep, (mode == cadet::model::MultiplexMode::Radial) ? cadet::SectionIndep : s)] = &values[s * nChannel * nComp + i * nComp];
				}
			}
			break;
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
			cadet::registerParam3DArray(parameters, values, [=](bool multi, unsigned int sec, unsigned int compartment, unsigned int comp) { return cadet::makeParamId(nameHash, uoi, comp, compartment, cadet::BoundStateIndep, cadet::ReactionIndep, multi ? sec : cadet::SectionIndep); }, nComp, nChannel);
			break;
		case cadet::model::MultiplexMode::Axial:
		case cadet::model::MultiplexMode::AxialRadial:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			cadet_assert(false);
			break;
	}

	return mode;
}

bool multiplexParameterValue(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nComp, unsigned int nChannel, double value, std::unordered_set<cadet::active*> const* sensParams)
{
	if (pId.name != nameHash)
		return false;

	switch (mode)
	{
		case cadet::model::MultiplexMode::Independent:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[0]))
					return false;

				for (std::size_t i = 0; i < data.size(); ++i)
					data[i].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::Section:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.section * nComp * nChannel]))
					return false;

				for (unsigned int i = 0; i < nComp * nChannel; ++i)
					data[i + pId.section * nComp * nChannel].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::Component:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component]))
					return false;

				for (unsigned int i = 0; i < nChannel; ++i)
					data[i * nComp + pId.component].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component + pId.section * nComp * nChannel]))
					return false;

				for (unsigned int i = 0; i < nChannel; ++i)
					data[i * nComp + pId.component + pId.section * nComp * nChannel].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::Radial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.particleType * nComp]))
					return false;

				for (unsigned int i = 0; i < nComp; ++i)
					data[i + pId.particleType * nComp].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::RadialSection:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.particleType * nComp + pId.section * nComp * nChannel]))
					return false;

				for (unsigned int i = 0; i < nComp; ++i)
					data[i + pId.particleType * nComp + pId.section * nComp * nChannel].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentRadial:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component + pId.particleType * nComp]))
					return false;

				data[pId.component + pId.particleType * nComp].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentRadialSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component + pId.particleType * nComp + pId.section * nComp * nChannel]))
					return false;

				data[pId.component + pId.particleType * nComp + pId.section * nComp * nChannel].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::Axial:
		case cadet::model::MultiplexMode::AxialRadial:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			cadet_assert(false);
			break;
	}

	return false;
}

bool multiplexParameterAD(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nComp, unsigned int nChannel, unsigned int adDirection, double adValue, std::unordered_set<cadet::active*>& sensParams)
{
	if (pId.name != nameHash)
		return false;

	switch (mode)
	{
		case cadet::model::MultiplexMode::Independent:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[0]);

				for (std::size_t i = 0; i < data.size(); ++i)
					data[i].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::Section:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.section * nComp * nChannel]);

				for (unsigned int i = 0; i < nComp * nChannel; ++i)
					data[i + pId.section * nComp * nChannel].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::Component:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.component]);

				for (unsigned int i = 0; i < nChannel; ++i)
					data[i * nComp + pId.component].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.component + pId.section * nComp * nChannel]);

				for (unsigned int i = 0; i < nChannel; ++i)
					data[i * nComp + pId.component + pId.section * nComp * nChannel].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::Radial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.particleType * nComp]);

				for (unsigned int i = 0; i < nComp; ++i)
					data[i + pId.particleType * nComp].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::RadialSection:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.particleType * nComp + pId.section * nComp * nChannel]);

				for (unsigned int i = 0; i < nComp; ++i)
					data[i + pId.particleType * nComp + pId.section * nComp * nChannel].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentRadial:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.component + pId.particleType * nComp]);

				data[pId.component + pId.particleType * nComp].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentRadialSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.component + pId.particleType * nComp + pId.section * nComp * nChannel]);

				data[pId.component + pId.particleType * nComp + pId.section * nComp * nChannel].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::Axial:
		case cadet::model::MultiplexMode::AxialRadial:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			cadet_assert(false);
			break;
	}

	return false;
}

}  // namespace

namespace cadet
{

namespace model
{

namespace parts
{

class MultiChannelConvectionDispersionOperator::LinearSolver
{
public:

	virtual ~LinearSolver() CADET_NOEXCEPT { }

	virtual bool initialize(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int nChannel, const Weno& weno) = 0;
	virtual void setSparsityPattern(const cadet::linalg::SparsityPattern& pattern) = 0;
	virtual void assembleDiscretizedJacobian(double alpha) = 0;
	virtual bool factorize() = 0;
	virtual bool solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const = 0;
};

int matrixMultiplierMultiChannelCDO(void* userData, double const* x, double* z);

class MultiChannelConvectionDispersionOperator::GmresSolver : public MultiChannelConvectionDispersionOperator::LinearSolver
{
public:

	GmresSolver(linalg::CompressedSparseMatrix const* jacC) : _jacC(jacC) { }
	virtual ~GmresSolver() CADET_NOEXCEPT { }

	virtual bool initialize(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int nChannel, const Weno& weno)
	{
		_gmres.initialize(nCol * nComp * nChannel, 0, linalg::toOrthogonalization(1), 0);
		_gmres.matrixVectorMultiplier(&matrixMultiplierMultiChannelCDO, this);
		_cache.resize(nCol * nComp * nChannel, 0.0);

		return true;
	}

	virtual void setSparsityPattern(const linalg::SparsityPattern& pattern) { }

	virtual void assembleDiscretizedJacobian(double alpha)
	{
		_alpha = alpha;
	}

	virtual bool factorize() { return true; }

	virtual bool solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const
	{
		if (init)
			std::copy(init, init + _cache.size(), _cache.begin());

		const int gmresResult = _gmres.solve(0.05 * outerTol * std::sqrt(_jacC->rows()), weight, rhs, _cache.data());
		std::copy(_cache.begin(), _cache.end(), rhs);

		return gmresResult == 0;
	}

protected:
	linalg::CompressedSparseMatrix const* const _jacC;
	double _alpha;
	mutable linalg::Gmres _gmres; //!< GMRES algorithm for the Schur-complement in linearSolve()
	mutable std::vector<double> _cache; //!< GMRES cache for result

	int matrixVectorMultiply(double const* x, double* z) const
	{
		std::fill(z, z + _jacC->rows(), _alpha);
		_jacC->multiplyVector(x, 1.0, 1.0, z);
		return 0;
	}

	// Wrapper for calling the corresponding function in this class
	friend int matrixMultiplierMultiChannelCDO(void* userData, double const* x, double* z);
};

int matrixMultiplierMultiChannelCDO(void* userData, double const* x, double* z)
{
	MultiChannelConvectionDispersionOperator::GmresSolver* const cdo = static_cast<MultiChannelConvectionDispersionOperator::GmresSolver*>(userData);
	return cdo->matrixVectorMultiply(x, z);
}

#if defined(UMFPACK_FOUND) || defined(SUPERLU_FOUND)

	template <typename sparse_t>
	class MultiChannelConvectionDispersionOperator::SparseDirectSolver : public MultiChannelConvectionDispersionOperator::LinearSolver
	{
	public:

		SparseDirectSolver(linalg::CompressedSparseMatrix const* jacC) : _jacC(jacC) { }
		virtual ~SparseDirectSolver() CADET_NOEXCEPT { }

		virtual bool initialize(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int nChannel, const Weno& weno)
		{
			return true;
		}

		virtual void setSparsityPattern(const linalg::SparsityPattern& pattern)
		{
			_jacCdisc.assignPattern(pattern);
			_jacCdisc.prepare();
		}

		virtual void assembleDiscretizedJacobian(double alpha)
		{
			// Copy normal matrix over to factorizable matrix
			_jacCdisc.copyFromSamePattern(*_jacC);

			for (int i = 0; i < _jacC->rows(); ++i)
				_jacCdisc.centered(i, 0) += alpha;
		}

		virtual bool factorize()
		{
			return _jacCdisc.factorize();
		}

		virtual bool solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const
		{
			return _jacCdisc.solve(rhs);
		}

	protected:
		linalg::CompressedSparseMatrix const* const _jacC;
		sparse_t _jacCdisc;
	};

#endif

class MultiChannelConvectionDispersionOperator::DenseDirectSolver : public MultiChannelConvectionDispersionOperator::LinearSolver
{
public:

	DenseDirectSolver(linalg::CompressedSparseMatrix const* jacC) : _jacC(jacC) { }
	virtual ~DenseDirectSolver() CADET_NOEXCEPT { }

	virtual bool initialize(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int nChannel, const Weno& weno)
	{
		// Note that we have to increase the lower bandwidth by 1 because the WENO stencil is applied to the
		// right cell face (lower + 1 + upper) and to the left cell face (shift the stencil by -1 because influx of cell i
		// is outflux of cell i-1)
		// We also have to make sure that there's at least one sub and super diagonal for the dispersion term
		const unsigned int lb = std::max(weno.lowerBandwidth() + 1u, 1u) * nComp * nChannel;

		// We have to make sure that there's at least one sub and super diagonal for the dispersion term
		const unsigned int ub = std::max(weno.upperBandwidth(), 1u) * nComp * nChannel;

		// When flow direction is changed, the bandwidths of the Jacobian swap.
		// Hence, we have to reserve memory such that the swapped Jacobian can fit into the matrix.
		const unsigned int mb = std::max(lb, ub);

		// Allocate matrices such that bandwidths can be switched (backwards flow support)
		_jacCdisc.resize(nCol * nComp * nChannel, mb, mb);
		return true;
	}

	virtual void setSparsityPattern(const linalg::SparsityPattern& pattern) { }

	virtual void assembleDiscretizedJacobian(double alpha)
	{
		// Copy normal matrix over to factorizable matrix
		_jacCdisc.setAll(0.0);

		linalg::FactorizableBandMatrix::RowIterator jac = _jacCdisc.row(0);
		for (int i = 0; i < _jacC->rows(); ++i, ++jac)
		{
			linalg::sparse_int_t const* const colIdx = _jacC->columnIndicesOfRow(i);
			double const* const vals = _jacC->valuesOfRow(i);
			const int nnz = _jacC->numNonZerosInRow(i);

			// Copy row from sparse matrix to banded matrix
			for (int c = 0; c < nnz; ++c)
			{
				const linalg::sparse_int_t diag = colIdx[c] - i;
				jac[diag] = vals[c];
			}

			// Add time derivative
			jac[0] += alpha;
		}
	}

	virtual bool factorize()
	{
		return _jacCdisc.factorize();
	}

	virtual bool solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const
	{
		return _jacCdisc.solve(rhs);
	}

protected:
	linalg::CompressedSparseMatrix const* const _jacC;
	linalg::FactorizableBandMatrix _jacCdisc;
};


/**
 * @brief Creates a MultiChannelConvectionDispersionOperator
 */
MultiChannelConvectionDispersionOperator::MultiChannelConvectionDispersionOperator() : _dir(0), _stencilMemory(sizeof(active) * Weno::maxStencilSize()),
	_wenoDerivatives(new double[Weno::maxStencilSize()]), _weno(), _linearSolver(nullptr)
{
}

MultiChannelConvectionDispersionOperator::~MultiChannelConvectionDispersionOperator() CADET_NOEXCEPT
{
	delete[] _wenoDerivatives;
	delete _linearSolver;
}

/**
 * @brief Reads parameters and allocates memory
 * @details Has to be called once before the operator is used.
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [in] nComp Number of components
 * @param [in] nCol Number of axial cells
 * @param [in] dynamicReactions Determines whether the sparsity pattern accounts for dynamic reactions
 * @return @c true if configuration went fine, @c false otherwise
 */
bool MultiChannelConvectionDispersionOperator::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, unsigned int nComp, unsigned int nCol, unsigned int nChannel, bool dynamicReactions)
{
	_nComp = nComp;
	_nCol = nCol;
	_nChannel = nChannel;
	_hasDynamicReactions = dynamicReactions;
	paramProvider.pushScope("discretization");

	_dispersionDep = helper.createParameterParameterDependence("CONSTANT_ONE");

	// Read WENO settings and apply them
	paramProvider.pushScope("weno");
	_weno.order(paramProvider.getInt("WENO_ORDER"));
	_weno.boundaryTreatment(paramProvider.getInt("BOUNDARY_MODEL"));
	_weno.epsilon(paramProvider.getDouble("WENO_EPS"));
	paramProvider.popScope();

	// Read solver settings
	if (paramProvider.exists("LINEAR_SOLVER_BULK"))
	{
		const std::string sol = paramProvider.getString("LINEAR_SOLVER_BULK");
		if (sol == "DENSE")
			_linearSolver = new DenseDirectSolver(&_jacC);
		else if (sol == "GMRES")
			_linearSolver = new GmresSolver(&_jacC);
#ifdef UMFPACK_FOUND
		else if (sol == "UMFPACK")
			_linearSolver = new SparseDirectSolver<linalg::UMFPackSparseMatrix>(&_jacC);
#endif
#ifdef SUPERLU_FOUND
		else if (sol == "SUPERLU")
			_linearSolver = new SparseDirectSolver<linalg::SuperLUSparseMatrix>(&_jacC);
#endif
		else
			throw InvalidParameterException("Unknown linear solver " + sol + " in field LINEAR_SOLVER_BULK");
	}

	// Default to sparse solver if available (preferably UMFPACK), fall back to dense
	if (!_linearSolver)
	{
#if defined(UMFPACK_FOUND)
		_linearSolver = new SparseDirectSolver<linalg::UMFPackSparseMatrix>(&_jacC);
#elif defined(SUPERLU_FOUND)
		_linearSolver = new SparseDirectSolver<linalg::SuperLUSparseMatrix>(&_jacC);
#else
		_linearSolver = new DenseDirectSolver(&_jacC);
#endif
//		LOG(Info) << "Default to dense banded linear solver due to invalid or missing LINEAR_SOLVER_BULK setting";
	}

	paramProvider.popScope();

	_crossSections.resize(nChannel);
	_curVelocity.resize(nChannel);

	return _linearSolver->initialize(paramProvider, nComp, nCol, nChannel, _weno);
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] unitOpIdx Unit operation id of the owning unit operation model
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [out] parameters Map in which local parameters are inserted
 * @return @c true if configuration went fine, @c false otherwise
 */
bool MultiChannelConvectionDispersionOperator::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters)
{
	// Read geometry parameters
	_colLength = paramProvider.getDouble("COL_LENGTH");
	const std::vector<double> initCrossSections = paramProvider.getDoubleArray("CHANNEL_CROSS_SECTION_AREAS");
	ad::copyToAd(initCrossSections.data(), _crossSections.data(), _nChannel);
	readParameterMatrix(_exchangeMatrix, paramProvider, "EXCHANGE_MATRIX", _nChannel * _nChannel * _nComp, 1);

	// Read section dependent parameters (transport)

	// Read VELOCITY
	_velocity.clear();
	if (paramProvider.exists("VELOCITY"))
	{
		readScalarParameterOrArray(_velocity, paramProvider, "VELOCITY", 1);

		if (paramProvider.exists("VELOCITY_MULTIPLEX"))
		{
			const int mode = paramProvider.getInt("VELOCITY_MULTIPLEX");
			if (mode == 0)
				// Rad-indep, sec-indep
				_singleVelocity = true;
			else if (mode == 1)
				// Rad-dep, sec-indep
				_singleVelocity = false;
			else if (mode == 2)
				// Rad-indep, sec-dep
				_singleVelocity = true;
			else if (mode == 3)
				// Rad-dep, sec-dep
				_singleVelocity = false;

			if (!_singleVelocity && (_velocity.size() % _nChannel != 0))
				throw InvalidParameterException("Number of elements in field VELOCITY is not a positive multiple of NCHANNEL (" + std::to_string(_nChannel) + ")");
			if ((mode == 0) && (_velocity.size() != 1))
				throw InvalidParameterException("Number of elements in field VELOCITY inconsistent with VELOCITY_MULTIPLEX (should be 1)");
			if ((mode == 1) && (_velocity.size() != _nChannel))
				throw InvalidParameterException("Number of elements in field VELOCITY inconsistent with VELOCITY_MULTIPLEX (should be " + std::to_string(_nChannel) + ")");
		}
		else
		{
			// Infer radial dependence of VELOCITY:
			//   size not divisible by NCHANNEL -> radial independent
			_singleVelocity = ((_velocity.size() % _nChannel) != 0);
		}

		// Expand _velocity to make it component dependent
		if (_singleVelocity)
		{
			std::vector<active> expanded(_velocity.size() * _nChannel);
			for (std::size_t i = 0; i < _velocity.size(); ++i)
				std::fill(expanded.begin() + i * _nChannel, expanded.begin() + (i + 1) * _nChannel, _velocity[i]);

			_velocity = std::move(expanded);
		}
	}
	else
	{
		_singleVelocity = false;
		_velocity.resize(_nChannel, 1.0);
	}

	// Register VELOCITY
	if (_singleVelocity)
	{
		if (_velocity.size() > _nChannel)
		{
			// Register only the first item in each section
			for (std::size_t i = 0; i < _velocity.size() / _nChannel; ++i)
				parameters[makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, i)] = &_velocity[i * _nChannel];
		}
		else
		{
			// We have only one parameter
			parameters[makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_velocity[0];
		}
	}
	else
		registerParam2DArray(parameters, _velocity, [=](bool multi, unsigned int sec, unsigned int compartment) { return makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, compartment, BoundStateIndep, ReactionIndep, multi ? sec : SectionIndep); }, _nChannel);

	_dir = std::vector<int>(_nChannel, 1);

	_axialDispersionMode = readAndRegisterMultiplexParam(paramProvider, parameters, _axialDispersion, "COL_DISPERSION", _nComp, _nChannel, unitOpIdx);

	// Add parameters to map
	parameters[makeParamId(hashString("COL_LENGTH"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colLength;
	registerParam3DArray(parameters, _exchangeMatrix, [=](bool multi, unsigned int channelSrc, unsigned int channelDest, unsigned comp) { return makeParamId(hashString("EXCHANGE_MATRIX"), unitOpIdx, multi ? comp : CompIndep, channelDest, channelSrc, ReactionIndep, SectionIndep); }, _nComp, _nChannel);

	setSparsityPattern();

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
bool MultiChannelConvectionDispersionOperator::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx)
{
	bool hasChanged = false;

	if (!_velocity.empty())
	{
		// _curVelocity has already been set to the network flow rate in setFlowRates()
		// the direction of the flow (i.e., sign of _curVelocity) is given by _velocity
		active const* const dirNew = getSectionDependentSlice(_velocity, _nChannel, secIdx);

		for (unsigned int i = 0; i < _nChannel; ++i)
		{
			const int newDir = (dirNew[i] >= 0) ? 1 : -1;
			if (_dir[i] * newDir < 0)
			{
				hasChanged = true;
				_curVelocity[i] *= -1;
			}
			_dir[i] = newDir;
		}
	}


	// Change the sparsity pattern if necessary
	if ((secIdx == 0) || hasChanged)
		setSparsityPattern();

	return hasChanged || (secIdx == 0);
}

/**
 * @brief Sets the AD seed vectors for the bulk transport variables
 * @param [in,out] adJac Jacobian information for AD (AD vectors for residual and state, direction offset)
 */
void MultiChannelConvectionDispersionOperator::prepareADvectors(const AdJacobianParams& adJac) const
{
	// Early out if AD is disabled
	if (!adJac.adY)
		return;

	// todo: implement compressed AD seeding
	active* const adVec = adJac.adY + _nChannel * _nComp;
	const int adOff = adJac.adDirOffset;

	const int nDof = _nComp * _nCol * _nChannel;

	for (int eq = 0; eq < nDof; ++eq)
	{
		// Clear previously set directions
		adVec[eq].fillADValue(adOff, 0.0);
		// Set directions
		adVec[eq].setADValue(adOff + eq, 1.0);
	}
}

/**
 * @brief Sets the flow rates for the current time section
 * @details The flow rates may change due to valve switches.
 * @param [in] compartment Index of the compartment
 * @param [in] in Total volumetric inlet flow rate
 * @param [in] out Total volumetric outlet flow rate
 */
void MultiChannelConvectionDispersionOperator::setFlowRates(int compartment, const active& in, const active& out) CADET_NOEXCEPT
{
	_curVelocity[compartment] = _dir[compartment] * in / (_crossSections[compartment]);
}

void MultiChannelConvectionDispersionOperator::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	for (unsigned int compartment = 0; compartment < _nChannel; ++compartment)
		_curVelocity[compartment] = _dir[compartment] * in[compartment] / (_crossSections[compartment]);
}

double MultiChannelConvectionDispersionOperator::inletFactor(unsigned int idxSec, int idxRad) const CADET_NOEXCEPT
{
	const double h = static_cast<double>(_colLength) / static_cast<double>(_nCol);
	return -std::abs(static_cast<double>(_curVelocity[idxRad])) / h;
}

const active& MultiChannelConvectionDispersionOperator::axialDispersion(unsigned int idxSec, int idxRad, int idxComp) const CADET_NOEXCEPT
{
	return *(getSectionDependentSlice(_axialDispersion, _nChannel * _nComp, idxSec) + idxRad * _nComp + idxComp);
}


/**
 * @brief Computes the residual of the transport equations
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] y Pointer to unit operation's state vector
 * @param [in] yDot Pointer to unit operation's time derivative state vector
 * @param [out] res Pointer to unit operation's residual vector
 * @param [in] wantJac Determines whether the Jacobian is computed or not
 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
 */
int MultiChannelConvectionDispersionOperator::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, bool wantJac, WithoutParamSensitivity)
{
	if (wantJac)
		return residualImpl<double, double, double, true>(model, t, secIdx, y, yDot, res);
	else
		return residualImpl<double, double, double, false>(model, t, secIdx, y, yDot, res);
}

int MultiChannelConvectionDispersionOperator::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, bool wantJac, WithoutParamSensitivity)
{
	if (wantJac)
		return residualImpl<active, active, double, true>(model, t, secIdx, y, yDot, res);
	else
		return residualImpl<active, active, double, false>(model, t, secIdx, y, yDot, res);
}

int MultiChannelConvectionDispersionOperator::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, bool wantJac, WithParamSensitivity)
{
	if (wantJac)
		return residualImpl<double, active, active, true>(model, t, secIdx, y, yDot, res);
	else
		return residualImpl<double, active, active, false>(model, t, secIdx, y, yDot, res);
}

int MultiChannelConvectionDispersionOperator::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, bool wantJac, WithParamSensitivity)
{
	if (wantJac)
		return residualImpl<active, active, active, true>(model, t, secIdx, y, yDot, res);
	else
		return residualImpl<active, active, active, false>(model, t, secIdx, y, yDot, res);
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int MultiChannelConvectionDispersionOperator::residualImpl(const IModel& model, double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res)
{
	if (wantJac)
	{
		// Reset Jacobian
		_jacC.setAll(0.0);
	}

	// Handle convection, axial dispersion (WENO)
	const ParamType h = static_cast<ParamType>(_colLength) / static_cast<double>(_nCol);
	for (unsigned int i = 0; i < _nChannel; ++i)
	{
		active const* const d_c = getSectionDependentSlice(_axialDispersion, _nChannel * _nComp, secIdx) + i * _nComp;

		convdisp::AxialFlowParameters<ParamType, cadet::Weno> fp{
			static_cast<ParamType>(_curVelocity[i]),
			d_c,
			h,
			_wenoDerivatives,
			&_weno,
			&_stencilMemory,
			static_cast<int>(_nComp * _nChannel),  // Stride between two cells
			_nComp,
			_nCol,
			_nComp * i,                        // Offset to the first component of the inlet DOFs in the local state vector
			_nComp * (_nChannel + i),          // Offset to the first component of the first bulk cell in the local state vector
			_dispersionDep,
			model
		};

		if (wantJac)
			convdisp::residualKernelAxial<StateType, ResidualType, ParamType, cadet::Weno, linalg::BandedSparseRowIterator, true>(SimulationTime{t, secIdx}, y, yDot, res, _jacC.row(i * _nComp), fp);
		else
			convdisp::residualKernelAxial<StateType, ResidualType, ParamType, cadet::Weno, linalg::BandedSparseRowIterator, false>(SimulationTime{t, secIdx}, y, yDot, res, _jacC.row(i * _nComp), fp);
	}

	// Handle inter-channel transport
	if (cadet_unlikely(_nChannel <= 1))
		return 0;

	const unsigned int offsetC = _nChannel * _nComp;
	for (unsigned int col = 0; col < _nCol; ++col)
	{
		const unsigned int offsetColBlock = col * _nChannel * _nComp;
		ResidualType* const resColBlock = res + offsetC + offsetColBlock;
		StateType const* const yColBlock = y + offsetC + offsetColBlock;

		for (unsigned int rad_orig = 0; rad_orig < _nChannel; ++rad_orig)
		{
			const unsigned int offsetToRadOrigBlock = rad_orig * _nComp;
			const unsigned int offsetColRadOrigBlock = offsetColBlock + offsetToRadOrigBlock;
			ResidualType* const resColRadOrigBlock = resColBlock + offsetToRadOrigBlock;
			StateType const* const yColRadOrigBlock = yColBlock + offsetToRadOrigBlock;

            for (unsigned int rad_dest = 0; rad_dest < _nChannel; ++rad_dest)
            {
				if (rad_orig == rad_dest)
					continue;

                const unsigned int offsetToRadDestBlock = rad_dest * _nComp;
                const unsigned int offsetColRadDestBlock = offsetColBlock + offsetToRadDestBlock;
                ResidualType* const resColRadDestBlock = resColBlock + offsetToRadDestBlock;
                // StateType const* const yColRadDestBlock = yColBlock + offsetToRadDestBlock;

                for (unsigned int comp = 0; comp < _nComp; ++comp)
                {
                    const unsigned int offsetCur_orig = offsetColRadOrigBlock + comp;
                    const unsigned int offsetCur_dest = offsetColRadDestBlock + comp;
                    StateType const* const yCur_orig = yColRadOrigBlock + comp;
                    // StateType const* const yCur_dest = yColRadDestBlock + comp;
                    ResidualType* const resCur_orig = resColRadOrigBlock + comp;
                    ResidualType* const resCur_dest = resColRadDestBlock + comp;

					const ParamType exchange_orig_dest_comp = static_cast<ParamType>(_exchangeMatrix[rad_orig * _nChannel * _nComp + rad_dest * _nComp + comp]);
					if (cadet_likely(exchange_orig_dest_comp > 0.0))
					{
						*resCur_orig += exchange_orig_dest_comp * yCur_orig[0];
						*resCur_dest -= exchange_orig_dest_comp * yCur_orig[0] * static_cast<ParamType>(_crossSections[rad_orig]) / static_cast<ParamType>(_crossSections[rad_dest]);

						if (wantJac)
						{
							_jacC.centered(offsetCur_orig, 0) += static_cast<double>(exchange_orig_dest_comp);
							_jacC.centered(offsetCur_dest, static_cast<int>(offsetCur_orig) - static_cast<int>(offsetCur_dest)) -= static_cast<double>(exchange_orig_dest_comp);
						}
					}
                }

            }

		}
	}

	return 0;
}

void MultiChannelConvectionDispersionOperator::setSparsityPattern()
{
	// Note that we have to increase the lower non-zeros by 1 because the WENO stencil is applied to the
	// right cell face (lower + 1 + upper) and to the left cell face (shift the stencil by -1 because influx of cell i
	// is outflux of cell i-1)
	// We also have to make sure that there's at least one sub and super diagonal for the dispersion term
	const unsigned int lowerNonZeros = std::max(_weno.lowerBandwidth() + 1u, 1u);
	const unsigned int upperNonZeros = std::max(_weno.upperBandwidth(), 1u);
	// Total number of non-zeros per row is WENO stencil (lowerNonZeros + 1u + upperNonZeros) + radial dispersion (2)
	cadet::linalg::SparsityPattern pattern(_nComp * _nCol * _nChannel, lowerNonZeros + 1u + upperNonZeros + 2u);

	// Handle convection, axial dispersion (WENO)
	for (unsigned int i = 0; i < _nChannel; ++i)
		cadet::model::parts::convdisp::sparsityPatternAxial(pattern.row(i * _nComp), _nComp, _nCol, _nComp * _nChannel, static_cast<double>(_curVelocity[i]), _weno);

	if (_nChannel > 1)
	{
		for (unsigned int col = 0; col < _nCol; ++col)
		{
			// Axial column element
			const unsigned int idxColBlock = col * _nChannel * _nComp;

			// Connecting from all compartments (orig) to all compartments (dest)
			for (unsigned int rad_orig = 0; rad_orig < _nChannel ; ++rad_orig)
			{
				const unsigned int idxColRadBlock_orig = idxColBlock + rad_orig * _nComp;

				for (unsigned int rad_dest = 0; rad_dest < _nChannel ; ++rad_dest)
				{
					// Don't connect compartments with themselves
					if (rad_orig == rad_dest)
						continue;

					const unsigned int idxColRadBlock_dest = idxColBlock + rad_dest * _nComp;

					for (unsigned int comp = 0; comp < _nComp; ++comp)
					{
						const active& exchange_orig_dest_comp = _exchangeMatrix[rad_orig * _nChannel * _nComp + rad_dest * _nComp + comp];
						if (exchange_orig_dest_comp == 0.0)
							continue;

						const unsigned int idxCur_orig = idxColRadBlock_orig + comp;
						const unsigned int idxCur_dest = idxColRadBlock_dest + comp;
						pattern.add(idxCur_orig, idxCur_dest);
					}
				}
			}
		}
	}

	// Add space for dynamic reactions
	if (_hasDynamicReactions)
	{
		// Add nComp x nComp diagonal blocks (everything can react with everything)
		for (unsigned int col = 0; col < _nCol; ++col)
		{
			const unsigned int idxColBlock = col * _nChannel * _nComp;

			for (unsigned int rad = 0; rad < _nChannel; ++rad)
			{
				const unsigned int idxColRadBlock = idxColBlock + rad * _nComp;

				for (unsigned int comp = 0; comp < _nComp; ++comp)
				{
					const unsigned int idxCur = idxColRadBlock + comp;
					for (unsigned int comp2 = 0; comp2 < _nComp; ++comp2)
						pattern.add(idxCur, idxColRadBlock + comp2);
				}
			}
		}
	}

	_jacC.assignPattern(pattern);
	_linearSolver->setSparsityPattern(pattern);
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
void MultiChannelConvectionDispersionOperator::multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const
{
	double* localRet = ret + _nComp * _nChannel;
	double const* localSdot = sDot + _nComp * _nChannel;
	for (unsigned int i = 0; i < _nCol * _nComp * _nChannel; ++i)
		localRet[i] = localSdot[i];
}

/**
 * @brief Extracts the system Jacobian from band compressed AD seed vectors
 * @details The input vectors are assumed to point to the beginning (including inlet DOFs) of the respective unit operation's arrays.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void MultiChannelConvectionDispersionOperator::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	// todo: implement compressed AD

	const int Dofs = _nComp * _nChannel * _nCol;
	active const* const adVec = adRes + adDirOffset;

	for (int row = 0; row < Dofs; row++)
		for (int col = 0; col < Dofs; col++)
			_jacC.native(row, col) = adVec[row].getADValue(col);
}

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
void MultiChannelConvectionDispersionOperator::assembleDiscretizedJacobian(double alpha)
{
	_linearSolver->assembleDiscretizedJacobian(alpha);
}

/**
 * @brief Assembles and factorizes the time discretized Jacobian
 * @details See assembleDiscretizedJacobian() for assembly of the time discretized Jacobian.
 * @param [in] alpha Factor in front of @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @return @c true if factorization went fine, otherwise @c false
 */
bool MultiChannelConvectionDispersionOperator::assembleAndFactorizeDiscretizedJacobian(double alpha)
{
	assembleDiscretizedJacobian(alpha);
	return _linearSolver->factorize();
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
bool MultiChannelConvectionDispersionOperator::solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const
{
	return _linearSolver->solveDiscretizedJacobian(rhs, weight, init, outerTol);
}

/**
 * @brief Solves a system with the time derivative Jacobian and given right hand side
 * @details Note that the given right hand side vector @p rhs is not shifted by the inlet DOFs. That
 *          is, it is assumed to point directly to the first axial DOF.
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in,out] rhs On entry, right hand side. On exit, solution of the system.
 * @return @c true if the system was solved correctly, @c false otherwise
 */
bool MultiChannelConvectionDispersionOperator::solveTimeDerivativeSystem(const SimulationTime& simTime, double* const rhs)
{
	return true;
}

bool MultiChannelConvectionDispersionOperator::setParameter(const ParameterId& pId, double value)
{

	if (_singleVelocity && (pId.name == hashString("VELOCITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep) && (pId.reaction == ReactionIndep))
	{
		if (_velocity.size() > _nChannel)
		{
			// Section dependent
			if (pId.section == SectionIndep)
				return false;

			for (unsigned int i = 0; i < _nChannel; ++i)
				_velocity[pId.section * _nChannel + i].setValue(value);
		}
		else
		{
			// Section independent
			if (pId.section != SectionIndep)
				return false;

			for (unsigned int i = 0; i < _nChannel; ++i)
				_velocity[i].setValue(value);
		}
	}

	const bool ad = multiplexParameterValue(pId, hashString("COL_DISPERSION"), _axialDispersionMode, _axialDispersion, _nComp, _nChannel, value, nullptr);
	if (ad)
		return true;

	return false;
}

bool MultiChannelConvectionDispersionOperator::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
{

	if (_singleVelocity && (pId.name == hashString("VELOCITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep) && (pId.reaction == ReactionIndep))
	{
		if (_velocity.size() > _nChannel)
		{
			// Section dependent
			if ((pId.section == SectionIndep) || !contains(sensParams, &_velocity[pId.section * _nChannel]))
				return false;

			for (unsigned int i = 0; i < _nChannel; ++i)
				_velocity[pId.section * _nChannel + i].setValue(value);
		}
		else
		{
			// Section independent
			if ((pId.section != SectionIndep) || !contains(sensParams, &_velocity[0]))
				return false;

			for (unsigned int i = 0; i < _nChannel; ++i)
				_velocity[i].setValue(value);
		}
	}

	const bool ad = multiplexParameterValue(pId, hashString("COL_DISPERSION"), _axialDispersionMode, _axialDispersion, _nComp, _nChannel, value, &sensParams);
	if (ad)
		return true;

	return false;
}

bool MultiChannelConvectionDispersionOperator::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
{

	if (_singleVelocity && (pId.name == hashString("VELOCITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep) && (pId.reaction == ReactionIndep))
	{
		if (_velocity.size() > _nChannel)
		{
			// Section dependent
			if (pId.section == SectionIndep)
				return false;

			sensParams.insert(&_velocity[pId.section * _nChannel]);
			for (unsigned int i = 0; i < _nChannel; ++i)
				_velocity[pId.section * _nChannel + i].setADValue(adDirection, adValue);
		}
		else
		{
			// Section independent
			if (pId.section != SectionIndep)
				return false;

			sensParams.insert(&_velocity[0]);
			for (unsigned int i = 0; i < _nChannel; ++i)
				_velocity[i].setADValue(adDirection, adValue);
		}
	}

	const bool ad = multiplexParameterAD(pId, hashString("COL_DISPERSION"), _axialDispersionMode, _axialDispersion, _nComp, _nChannel, adDirection, adValue, sensParams);
	if (ad)
		return true;

	return false;
}

}  // namespace parts

}  // namespace model

}  // namespace cadet

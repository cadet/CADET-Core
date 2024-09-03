// =============================================================================
//  CADET
//
//  Copyright © 2008-2022: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/parts/TwoDimensionalConvectionDispersionOperatorDG.hpp"
#include "model/parts/DGToolbox.hpp"
#include "cadet/Exceptions.hpp"
#include "ParamReaderHelper.hpp"
#include "AdUtils.hpp"
#include "SensParamUtil.hpp"
#include "SimulationTypes.hpp"
#include "model/parts/AxialConvectionDispersionKernel.hpp"
#include "model/ParameterDependence.hpp"
#include "ConfigurationHelper.hpp"
#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <unsupported/Eigen/KroneckerProduct>

using namespace Eigen;

namespace
{

cadet::model::MultiplexMode readAndRegisterMultiplexParam(cadet::IParameterProvider& paramProvider, std::unordered_map<cadet::ParameterId, cadet::active*>& parameters, std::vector<cadet::active>& values, const std::string& name, unsigned int nComp, unsigned int radNPoints, cadet::UnitOpIdx uoi)
{
	cadet::model::MultiplexMode mode = cadet::model::MultiplexMode::Independent;
	readParameterMatrix(values, paramProvider, name, nComp * radNPoints, 1);
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
			if (values.size() != radNPoints)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(radNPoints) + ")");
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
			if (values.size() != nComp * radNPoints)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nComp * radNPoints) + ")");
		}
		else if (modeConfig == 4)
		{
			mode = cadet::model::MultiplexMode::Section;
			nSec = values.size();
		}
		else if (modeConfig == 5)
		{
			mode = cadet::model::MultiplexMode::RadialSection;
			if (values.size() % radNPoints != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " is not a positive multiple of radNPoints (" + std::to_string(radNPoints) + ")");

			nSec = values.size() / radNPoints;
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
			if (values.size() % (nComp * radNPoints) != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " is not a positive multiple of NCOMP * radNPoints (" + std::to_string(nComp * radNPoints) + ")");

			nSec = values.size() / (nComp * radNPoints);
		}
	}
	else
	{
		if (values.size() == 1)
			mode = cadet::model::MultiplexMode::Independent;
		else if (values.size() == nComp)
			mode = cadet::model::MultiplexMode::Component;
		else if (values.size() == radNPoints)
			mode = cadet::model::MultiplexMode::Radial;
		else if (values.size() == radNPoints * nComp)
			mode = cadet::model::MultiplexMode::ComponentRadial;
		else if (values.size() % nComp == 0)
		{
			mode = cadet::model::MultiplexMode::ComponentSection;
			nSec = values.size() / nComp;
		}
		else if (values.size() % radNPoints == 0)
		{
			mode = cadet::model::MultiplexMode::RadialSection;
			nSec = values.size() / radNPoints;
		}
		else if (values.size() % (radNPoints * nComp) == 0)
		{
			mode = cadet::model::MultiplexMode::ComponentRadialSection;
			nSec = values.size() / (nComp * radNPoints);
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
				std::vector<cadet::active> p(nComp * radNPoints * nSec);
				for (unsigned int s = 0; s < nSec; ++s)
					std::fill(p.begin() + s * radNPoints * nComp, p.begin() + (s+1) * radNPoints * nComp, values[s]);

				values = std::move(p);

				for (unsigned int s = 0; s < nSec; ++s)
					parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, (mode == cadet::model::MultiplexMode::Independent) ? cadet::SectionIndep : s)] = &values[s * radNPoints * nComp];
			}
			break;
		case cadet::model::MultiplexMode::Component:
		case cadet::model::MultiplexMode::ComponentSection:
			{
				std::vector<cadet::active> p(nComp * radNPoints * nSec);
				for (unsigned int s = 0; s < nSec; ++s)
				{
					for (unsigned int i = 0; i < nComp; ++i)
						std::copy(values.begin() + s * nComp, values.begin() + (s+1) * nComp, p.begin() + i * nComp + s * nComp * radNPoints);
				}

				values = std::move(p);

				for (unsigned int s = 0; s < nSec; ++s)
				{
					for (unsigned int i = 0; i < nComp; ++i)
						parameters[cadet::makeParamId(nameHash, uoi, i, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, (mode == cadet::model::MultiplexMode::Component) ? cadet::SectionIndep : s)] = &values[s * radNPoints * nComp + i];
				}
			}
			break;
		case cadet::model::MultiplexMode::Radial:
		case cadet::model::MultiplexMode::RadialSection:
			{
				std::vector<cadet::active> p(nComp * radNPoints * nSec);
				for (unsigned int i = 0; i < radNPoints * nSec; ++i)
					std::fill(p.begin() + i * nComp, p.begin() + (i+1) * nComp, values[i]);

				values = std::move(p);

				for (unsigned int s = 0; s < nSec; ++s)
				{
					for (unsigned int i = 0; i < radNPoints; ++i)
						parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, i, cadet::BoundStateIndep, cadet::ReactionIndep, (mode == cadet::model::MultiplexMode::Radial) ? cadet::SectionIndep : s)] = &values[s * radNPoints * nComp + i * nComp];
				}
			}
			break;
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
			cadet::registerParam3DArray(parameters, values, [=](bool multi, unsigned int sec, unsigned int compartment, unsigned int comp) { return cadet::makeParamId(nameHash, uoi, comp, compartment, cadet::BoundStateIndep, cadet::ReactionIndep, multi ? sec : cadet::SectionIndep); }, nComp, radNPoints);
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

bool multiplexParameterValue(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nComp, unsigned int radNPoints, double value, std::unordered_set<cadet::active*> const* sensParams)
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

				if (sensParams && !cadet::contains(*sensParams, &data[pId.section * nComp * radNPoints]))
					return false;

				for (unsigned int i = 0; i < nComp * radNPoints; ++i)
					data[i + pId.section * nComp * radNPoints].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::Component:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component]))
					return false;

				for (unsigned int i = 0; i < radNPoints; ++i)
					data[i * nComp + pId.component].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component + pId.section * nComp * radNPoints]))
					return false;

				for (unsigned int i = 0; i < radNPoints; ++i)
					data[i * nComp + pId.component + pId.section * nComp * radNPoints].setValue(value);

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

				if (sensParams && !cadet::contains(*sensParams, &data[pId.particleType * nComp + pId.section * nComp * radNPoints]))
					return false;

				for (unsigned int i = 0; i < nComp; ++i)
					data[i + pId.particleType * nComp + pId.section * nComp * radNPoints].setValue(value);

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

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component + pId.particleType * nComp + pId.section * nComp * radNPoints]))
					return false;

				data[pId.component + pId.particleType * nComp + pId.section * nComp * radNPoints].setValue(value);

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

bool multiplexParameterAD(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nComp, unsigned int radNPoints, unsigned int adDirection, double adValue, std::unordered_set<cadet::active*>& sensParams)
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

				sensParams.insert(&data[pId.section * nComp * radNPoints]);

				for (unsigned int i = 0; i < nComp * radNPoints; ++i)
					data[i + pId.section * nComp * radNPoints].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::Component:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.component]);

				for (unsigned int i = 0; i < radNPoints; ++i)
					data[i * nComp + pId.component].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.component + pId.section * nComp * radNPoints]);

				for (unsigned int i = 0; i < radNPoints; ++i)
					data[i * nComp + pId.component + pId.section * nComp * radNPoints].setADValue(adDirection, adValue);

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

				sensParams.insert(&data[pId.particleType * nComp + pId.section * nComp * radNPoints]);

				for (unsigned int i = 0; i < nComp; ++i)
					data[i + pId.particleType * nComp + pId.section * nComp * radNPoints].setADValue(adDirection, adValue);

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

				sensParams.insert(&data[pId.component + pId.particleType * nComp + pId.section * nComp * radNPoints]);

				data[pId.component + pId.particleType * nComp + pId.section * nComp * radNPoints].setADValue(adDirection, adValue);

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

/**
 * @brief Creates a TwoDimensionalConvectionDispersionOperatorDG
 */
TwoDimensionalConvectionDispersionOperatorDG::TwoDimensionalConvectionDispersionOperatorDG() : _colPorosities(0), _dir(0), _dispersionDep(nullptr)
{
}

TwoDimensionalConvectionDispersionOperatorDG::~TwoDimensionalConvectionDispersionOperatorDG() CADET_NOEXCEPT
{
	delete _dispersionDep;
}

MatrixXd TwoDimensionalConvectionDispersionOperatorDG::calcTildeMr(const unsigned int elemIdx, const active* const dispersion)
{
	MatrixXd tildeMr = MatrixXd::Zero(_radNNodes, _qNNodes);
	MatrixXd ell = MatrixXd::Zero(_radNNodes, _qNNodes);

	for (unsigned int j = 0; j < _radNNodes; j++)
		ell.row(j) = dgtoolbox::evalLagrangeBasis(j, _radNodes, _qNodes);

	for (unsigned int j = 0; j < _radNNodes; j++)
	{
		for (unsigned int k = 0; k < _qNNodes; k++)
		{
			tildeMr(j, k) = static_cast<double>(dgtoolbox::mapRefToPhys<active>(_radDelta, elemIdx, _qNodes[k])) * static_cast<double>(dispersion[elemIdx * _qNNodes + k]) * ell(j, k) * _qWeights[k];
		}
	}

	return tildeMr;
}

MatrixXd TwoDimensionalConvectionDispersionOperatorDG::calcTildeMrDash(const unsigned int elemIdx, const active * const dispersion)
{
	MatrixXd ellEll = MatrixXd::Zero(_radNNodes, _qNNodes);
	MatrixXd tildeMrDash = MatrixXd::Zero(_radNNodes, _radNNodes * _qNNodes);

	for (unsigned int j = 0; j < _radNNodes; j++)
		ellEll.row(j) = dgtoolbox::evalLagrangeBasis(j, _radNodes, _qNodes);

	for (unsigned int j = 0; j < _radNNodes; j++)
	{
		for (unsigned int m = 0; m < _radNNodes; m++)
		{
			for (unsigned int k = 0; k < _qNNodes; k++)
			{
				tildeMrDash(j, m * _qNNodes + k) = static_cast<double>(dgtoolbox::mapRefToPhys<active>(_radDelta, elemIdx, _qNodes[k])) * static_cast<double>(dispersion[elemIdx * _qNNodes + k]) * ellEll(m, k) * ellEll(j, k) * _qWeights[k];
			}
		}
	}

	return tildeMrDash;
}

MatrixXd TwoDimensionalConvectionDispersionOperatorDG::calcTildeSrDash(const unsigned int elemIdx, const active* const dispersion)
{
	MatrixXd ellEll = MatrixXd::Zero(_radNNodes, _qNNodes);
	MatrixXd tildeSrDash = MatrixXd::Zero(_radNNodes, _radNNodes * _qNNodes);

	for (unsigned int j = 0; j < _radNNodes; j++)
		ellEll.row(j) = dgtoolbox::evalLagrangeBasis(j, _radNodes, _qNodes);

	MatrixXd tildeDr = dgtoolbox::derivativeMatrix(_quadratureOrder, _qNodes);
	MatrixXd tildeSr = _qWeights.asDiagonal() * tildeDr;

	for (unsigned int j = 0; j < _radNNodes; j++)
	{
		for (unsigned int m = 0; m < _radNNodes; m++)
		{
			for (unsigned int k = 0; k < _qNNodes; k++)
			{
				tildeSrDash(j, m * _qNNodes + k) = static_cast<double>(dgtoolbox::mapRefToPhys<active>(_radDelta, elemIdx, _qNodes[k])) * static_cast<double>(dispersion[elemIdx * _qNNodes + k]) * tildeSr(j, k) * ellEll(m, k);
			}
		}
	}

	return tildeSrDash;
}

void TwoDimensionalConvectionDispersionOperatorDG::writeAxialCoordinates(double* coords) const
{
	double* leftElemBndries = new double[_axNElem];
	for (int elem = 0; elem < _axNElem; elem++)
		leftElemBndries[elem] = static_cast<double>(_axDelta) * elem;
	dgtoolbox::writeDGCoordinates(coords, _axNElem, _axNNodes, static_cast<const double*>(&_axNodes[0]), static_cast<double>(_colLength), leftElemBndries);
	delete[] leftElemBndries;
}
void TwoDimensionalConvectionDispersionOperatorDG::writeRadialCoordinates(double* coords) const
{
	double* leftElemBndries = new double[_radNElem];
	for (int elem = 0; elem < _radNElem; elem++)
		leftElemBndries[elem] = static_cast<double>(_radDelta[elem]) * elem;
	dgtoolbox::writeDGCoordinates(coords, _radNElem, _radNNodes, static_cast<const double*>(&_radNodes[0]), static_cast<double>(_colLength), leftElemBndries);
	delete[] leftElemBndries;
}

void TwoDimensionalConvectionDispersionOperatorDG::initializeDG()
{
	_axInvWeights = VectorXd::Zero(_axNNodes);
	_axNodes = VectorXd::Zero(_axNNodes);
	_radInvWeights = VectorXd::Zero(_radNNodes);
	_radNodes = VectorXd::Zero(_radNNodes);
	_qNodes = VectorXd::Zero(_qNNodes);
	_qWeights = VectorXd::Zero(_qNNodes);
	// LGL nodes and weights for axial and radial reference element and radial quadrature
	dgtoolbox::lglNodesWeights(_axPolyDeg, _axNodes, _axInvWeights, true);
	dgtoolbox::lglNodesWeights(_radPolyDeg, _radNodes, _radInvWeights, true);
	dgtoolbox::lglNodesWeights(_quadratureOrder, _qNodes, _qWeights, false);

	// auxiliary variables
	_axAuxStateG.resize(_axNPoints * _radNPoints);
	_radAuxStateG.resize(_axNPoints * _radNPoints);
	_axAuxStateGTilde.resize(_axNPoints * _radNElem * _qNNodes);
	_radAuxStateGTilde.resize(_axNPoints * _radNElem * _qNNodes);

	// numerical fluxes
	_fStarAux1.resize(2 * _radNNodes);
	_fStarAux2.resize(_axNNodes * 2);
	_fStarConv.resize(2 * _radNNodes);
	_gZStarDispTilde.resize(2 * _qNNodes);
	_gRStarDisp.resize(_axNNodes * 2);

	// operators
	_axInvMM = dgtoolbox::invMMatrix(_axPolyDeg, _axNodes, 0.0, 0.0); // todo collocation option?
	_radInvTransMM = dgtoolbox::mMatrix(_radPolyDeg, _radNodes, 0.0, 0.0).transpose().inverse();
	_axStiffM = dgtoolbox::stiffnessMatrix(_axPolyDeg, _axNodes, 0.0, 0.0); // todo collocation option?
	_axTransStiffM = dgtoolbox::stiffnessMatrix(_axPolyDeg, _axNodes, 0.0, 0.0).transpose();
	_radStiffM = dgtoolbox::stiffnessMatrix(_radPolyDeg, _radNodes, 0.0, 0.0);
	_axLiftM = dgtoolbox::liftingMatrix(_axNNodes);
	_radLiftM = dgtoolbox::liftingMatrix(_radNNodes).transpose();
	_matrixProductCache.resize(_axNNodes * _radNNodes); // operator matrix product cache

	VectorXd radBaryWeights = dgtoolbox::barycentricWeights(_quadratureOrder, _qNodes);
	_radInterpolationM = dgtoolbox::polynomialInterpolationMatrix(_qNodes, _radNodes, radBaryWeights);

	_transMrCyl = new MatrixXd[_radNElem];
	_invTransMrCyl = new MatrixXd[_radNElem];
	_transTildeMr = new MatrixXd[_radNElem];
	_transTildeMrDash = new MatrixXd[_radNElem];
	_transTildeSrDash = new MatrixXd[_radNElem];
	_SrCyl = new MatrixXd[_radNElem];
	_radLiftMCyl = new MatrixXd[_radNElem];

	// Jacobian blocks
	const int uAxElem = std::min(5, static_cast<int>(_axNElem)); // number of unique axial Jacobian blocks (per radial element)
	_jacConvection = new MatrixXd[_radNElem];
	_jacAxDispersion = new MatrixXd[_radNElem * uAxElem];
	_jacRadDispersion = new MatrixXd[_radNElem];

	for (int rElem = 0; rElem < _radNElem; rElem++)
	{
		_jacConvection[rElem] = MatrixXd::Zero(_elemNPoints, 5 * _elemNPoints);
		for (int i = 0; i < uAxElem; i++)
			_jacAxDispersion[rElem * uAxElem + i] = MatrixXd::Zero(_elemNPoints, 5 * _elemNPoints);
		_jacRadDispersion[rElem] = MatrixXd::Zero(_elemNPoints, 5 * _elemNPoints);
	}
}

/**
 * @brief Reads parameters and allocates memory
 * @details Has to be called once before the operator is used.
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [in] nComp Number of components
 * @param [in] axNPoints Number of axial cells
 * @param [in] dynamicReactions Determines whether the sparsity pattern accounts for dynamic reactions
 * @return @c true if configuration went fine, @c false otherwise
 */
bool TwoDimensionalConvectionDispersionOperatorDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const unsigned int nComp, const unsigned int radNodeStride)
{
	_nComp = nComp;
	_strideBound = radNodeStride - nComp;
	//_hasDynamicReactions = dynamicReactions; // todo needed in pattern? see FV operator

	// TODO: Add support for parameter dependent dispersion
	_dispersionDep = helper.createParameterParameterDependence("CONSTANT_ONE");

	paramProvider.pushScope("discretization");

	if (paramProvider.exists("AX_POLYDEG"))
		_axPolyDeg = paramProvider.getInt("AX_POLYDEG");
	else
		_axPolyDeg = 4u; // default value
	if (paramProvider.getInt("AX_POLYDEG") < 1)
		throw InvalidParameterException("Axial polynomial degree must be at least 1!");
	else if (_axPolyDeg < 3u)
		LOG(Warning) << "Polynomial degree > 2 in axial bulk discretization (cf. AX_POLYDEG) is always recommended for performance reasons.";

	_axNNodes = _axPolyDeg + 1u;

	if (paramProvider.exists("AX_NELEM"))
		_axNElem = paramProvider.getInt("AX_NELEM");
	else if (paramProvider.exists("NCOL"))
		_axNElem = std::max(1u, paramProvider.getInt("NCOL") / _axNNodes); // number of elements is rounded down
	else
		throw InvalidParameterException("Specify field AX_NELEM (or NCOL)");

	if (_axNElem < 1)
		throw InvalidParameterException("Number of column elements must be at least 1!");

	_axNPoints = _axNNodes * _axNElem;

	if (paramProvider.exists("RAD_POLYDEG"))
		_radPolyDeg = paramProvider.getInt("RAD_POLYDEG");
	else
		_radPolyDeg = 4u; // default value
	if (paramProvider.getInt("RAD_POLYDEG") < 1)
		throw InvalidParameterException("Radial polynomial degree must be at least 1!");

	_radNNodes = _radPolyDeg + 1;

	if (paramProvider.exists("RAD_NELEM"))
		_radNElem = paramProvider.getInt("RAD_NELEM");
	else if (paramProvider.exists("NRAD"))
		_radNElem = std::max(1u, paramProvider.getInt("NRAD") / _radNNodes); // number of elements is rounded down
	else
		throw InvalidParameterException("Specify field RAD_NELEM (or NRAD)");

	if (_radNElem < 1)
		throw InvalidParameterException("Number of column elements must be at least 1!");

	_radNPoints = _radNNodes * _radNElem;
	_elemNPoints = _axNNodes * _radNNodes;

	if (paramProvider.exists("QUADRATURE_RULE"))
	{
		const std::string quadratureRule = paramProvider.getString("QUADRATURE_RULE");
		if (quadratureRule == "LOBATTO")
			_quadratureRule = 0;
		else if (quadratureRule == "GAUSS")
			_quadratureRule = 1;
		else
			throw InvalidParameterException("Unknown quadrature rule " + quadratureRule);

		_quadratureOrder = paramProvider.exists("QUADRATURE_ORDER") ? paramProvider.getInt("QUADRATURE_ORDER") : _radPolyDeg + 1; // todo or other default?
	}
	else
	{
		_quadratureRule = 0;
		_quadratureOrder = paramProvider.exists("QUADRATURE_ORDER") ? paramProvider.getInt("QUADRATURE_ORDER") : _radPolyDeg; // todo or nNodes?
	}
	_qNNodes = _quadratureOrder + 1;
	paramProvider.popScope();

	_axNodeStride = radNodeStride * _radNPoints;
	_axElemStride = _axNNodes * _axNodeStride;
	_radNodeStride = radNodeStride;
	_radElemStride = _radNNodes * _radNodeStride;

	//_radialCoordinates.resize(_radNPoints + 1); // todo not needed, delete?
	_radialElemInterfaces.resize(_radNElem + 1);
	_radDelta.resize(_radNElem);
	_nodalCrossSections.resize(_radNPoints);
	_curVelocity.resize(_radNPoints);

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
bool TwoDimensionalConvectionDispersionOperatorDG::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters)
{
	// Read geometry parameters
	_colLength = paramProvider.getDouble("COL_LENGTH");
	_colRadius = paramProvider.getDouble("COL_RADIUS");
	readScalarParameterOrArray(_colPorosities, paramProvider, "COL_POROSITY", 1);
	_axDelta = _colLength / _axNElem; // deltaR is treated in updateRadialDisc

	if ((_colPorosities.size() != 1) && (_colPorosities.size() != _radNPoints))
		throw InvalidParameterException("Number of elements in field COL_POROSITY is neither 1 nor radNPoints (" + std::to_string(_radNPoints) + ")");

	_singlePorosity = (_colPorosities.size() == 1);
	if (_singlePorosity)
		_colPorosities = std::vector<active>(_radNPoints, _colPorosities[0]);

	// Read radial discretization mode and default to "EQUIDISTANT"
	paramProvider.pushScope("discretization");
	const std::string rdt = paramProvider.getString("RADIAL_DISC_TYPE");
	if (rdt == "EQUIVOLUME")
		_radialDiscretizationMode = RadialDiscretizationMode::Equivolume;
	else if (rdt == "USER_DEFINED")
	{
		_radialDiscretizationMode = RadialDiscretizationMode::UserDefined;
		readScalarParameterOrArray(_radialElemInterfaces, paramProvider, "RADIAL_COMPARTMENTS", 1);

		if (_radialElemInterfaces.size() < _radNElem + 1)
			throw InvalidParameterException("Number of elements in field RADIAL_COMPARTMENTS is less than radNElem + 1 (" + std::to_string(_radNElem + 1) + ")");

		registerParam1DArray(parameters, _radialElemInterfaces, [=](bool multi, unsigned int i) { return makeParamId(hashString("RADIAL_COMPARTMENTS"), unitOpIdx, CompIndep, i, BoundStateIndep, ReactionIndep, SectionIndep); });
	}
	else
		_radialDiscretizationMode = RadialDiscretizationMode::Equidistant;
	paramProvider.popScope();

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

			if (!_singleVelocity && (_velocity.size() % _radNPoints != 0))
				throw InvalidParameterException("Number of elements in field VELOCITY is not a positive multiple of radNPoints (" + std::to_string(_radNPoints) + ")");
			if ((mode == 0) && (_velocity.size() != 1))
				throw InvalidParameterException("Number of elements in field VELOCITY inconsistent with VELOCITY_MULTIPLEX (should be 1)");
			if ((mode == 1) && (_velocity.size() != _radNPoints))
				throw InvalidParameterException("Number of elements in field VELOCITY inconsistent with VELOCITY_MULTIPLEX (should be " + std::to_string(_radNPoints) + ")");
		}
		else
		{
			// Infer radial dependence of VELOCITY:
			//   size not divisible by radNPoints -> radial independent
			_singleVelocity = ((_velocity.size() % _radNPoints) != 0);
		}

		// Expand _velocity to make it component dependent
		if (_singleVelocity)
		{
			std::vector<active> expanded(_velocity.size() * _radNPoints);
			for (std::size_t i = 0; i < _velocity.size(); ++i)
				std::fill(expanded.begin() + i * _radNPoints, expanded.begin() + (i + 1) * _radNPoints, _velocity[i]);

			_velocity = std::move(expanded);
		}
	}
	else
	{
		_singleVelocity = false;
		_velocity.resize(_radNPoints, 1.0);
	}

	// Register VELOCITY
	if (_singleVelocity)
	{
		if (_velocity.size() > _radNPoints)
		{
			// Register only the first item in each section
			for (std::size_t i = 0; i < _velocity.size() / _radNPoints; ++i)
				parameters[makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, i)] = &_velocity[i * _radNPoints];
		}
		else
		{
			// We have only one parameter
			parameters[makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_velocity[0];
		}
	}
	else
		registerParam2DArray(parameters, _velocity, [=](bool multi, unsigned int sec, unsigned int compartment) { return makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, compartment, BoundStateIndep, ReactionIndep, multi ? sec : SectionIndep); }, _radNPoints);

	_dir = std::vector<int>(_radNPoints, 1);

	_axialDispersionMode = readAndRegisterMultiplexParam(paramProvider, parameters, _axialDispersion, "COL_DISPERSION_AXIAL", _nComp, _radNPoints, unitOpIdx);
	_radialDispersionMode = readAndRegisterMultiplexParam(paramProvider, parameters, _radialDispersion, "COL_DISPERSION_RADIAL", _nComp, _radNPoints, unitOpIdx);

	// Add parameters to map
	parameters[makeParamId(hashString("COL_LENGTH"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colLength;
	parameters[makeParamId(hashString("COL_RADIUS"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colRadius;
	registerParam1DArray(parameters, _colPorosities, [=](bool multi, unsigned int i) { return makeParamId(hashString("COL_POROSITY"), unitOpIdx, CompIndep, multi ? i : ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep); });

	// compute standard DG operators
	initializeDG();

	// interpolate radial dispersion to quadrature nodes
	_curRadialDispersionTilde = std::vector<active>(_radNElem * _qNNodes, 0.0);
	_curAxialDispersionTilde = std::vector<active>(_radNElem * _qNNodes, 0.0);
	// todo component dependence! getSectionDependentSlice gives dispersion parameter in radial position major
	const active* const curRadialDispersion = getSectionDependentSlice(_radialDispersion, _radNPoints * _nComp, 0);
	const active* const curAxialDispersion = getSectionDependentSlice(_axialDispersion, _radNPoints * _nComp, 0);
	for (unsigned int i = 0; i < _radInterpolationM.rows(); i++) {
		for (unsigned int j = 0; j < _radInterpolationM.cols(); j++) {
			_curRadialDispersionTilde[i] += _radInterpolationM(i, j) * curRadialDispersion[j];
			_curAxialDispersionTilde[i] += _radInterpolationM(i, j) * curAxialDispersion[j];
		}
	}

	// compute nodal cross section areas and radial element interfaces
	updateRadialDisc();

	// compute DG operators that depend on radial geometry and dispersion after updateRadialDisc
	const active* const d_rad = getSectionDependentSlice(_radialDispersion, _radNPoints * _nComp, 0);
	// todo component dependence of radial dispersion
	const int comp = 0;
	for (unsigned int rElem = 0; rElem < _radNElem; rElem++)
	{
		// todo ? use active types for column radius sensitivity
		_radLiftMCyl[rElem] = MatrixXd::Zero(2, _radNNodes);
		_radLiftMCyl[rElem](0, 0) = -static_cast<double>(_radialElemInterfaces[rElem]);
		_radLiftMCyl[rElem](1, _radNNodes - 1) = static_cast<double>(_radialElemInterfaces[rElem + 1]);

		_transMrCyl[rElem] = (static_cast<double>(_radialElemInterfaces[rElem]) * dgtoolbox::mMatrix(_radPolyDeg, _radNodes, 0.0, 0.0) + static_cast<double>(_radDelta[rElem]) / 2.0 * dgtoolbox::mMatrix(_radPolyDeg, _radNodes, 0.0, 1.0)).transpose();
		_invTransMrCyl[rElem] = _transMrCyl[rElem].inverse();

		_transTildeMr[rElem] = calcTildeMr(rElem, &_curAxialDispersionTilde[0]).transpose(); // todo still needed in new derivation?
		_transTildeMrDash[rElem] = calcTildeMrDash(rElem, &_curAxialDispersionTilde[0]).transpose(); // todo still needed in new derivation?
		_transTildeSrDash[rElem] = calcTildeSrDash(rElem, &_curRadialDispersionTilde[0]).transpose(); // todo still needed in new derivation?

		_SrCyl[rElem] = _transMrCyl[rElem].transpose() * dgtoolbox::derivativeMatrix(_radPolyDeg, _radNodes);
	}

	computeConvDispJacobianBlocks();

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
bool TwoDimensionalConvectionDispersionOperatorDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx)
{
	bool hasChanged = false;

	// todo update operators for section dependent parameters
	const active* const curRadialDispersion = getSectionDependentSlice(_radialDispersion, _radNPoints * _nComp, 0);
	for (unsigned int i = 0; i < _radInterpolationM.rows(); i++) {
		for (unsigned int j = 0; j < _radInterpolationM.cols(); j++) {
			_curRadialDispersionTilde[i] += _radInterpolationM(i, j) * curRadialDispersion[j];
		}
	}

	if (!_velocity.empty())
	{
		// _curVelocity has already been set to the network flow rate in setFlowRates()
		// the direction of the flow (i.e., sign of _curVelocity) is given by _velocity
		active const* const dirNew = getSectionDependentSlice(_velocity, _radNPoints, secIdx);

		for (unsigned int i = 0; i < _radNPoints; ++i)
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

	// todo: recompute operators that involve section dependent parameters, if (secIdx > 0)
	
	computeConvDispJacobianBlocks();

	//// todo backward flow

	return hasChanged || (secIdx == 0);
}

/**
 * @brief Sets the flow rates for the current time section
 * @details The flow rates may change due to valve switches.
 * @param [in] compartment Index of the compartment
 * @param [in] in Total volumetric inlet flow rate
 * @param [in] out Total volumetric outlet flow rate
 */
void TwoDimensionalConvectionDispersionOperatorDG::setFlowRates(int compartment, const active& in, const active& out) CADET_NOEXCEPT
{
	_curVelocity[compartment] = _dir[compartment] * in / (_nodalCrossSections[compartment] * _colPorosities[compartment]);
}

void TwoDimensionalConvectionDispersionOperatorDG::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	for (unsigned int compartment = 0; compartment < _radNPoints; ++compartment)
		_curVelocity[compartment] = in[compartment] / (_nodalCrossSections[compartment] * _colPorosities[compartment]);
}

double TwoDimensionalConvectionDispersionOperatorDG::inletFactor(unsigned int idxSec, int idxRad) const CADET_NOEXCEPT // todo what for?
{
	const double h = static_cast<double>(_colLength) / static_cast<double>(_axNPoints);
	return -std::abs(static_cast<double>(_curVelocity[idxRad])) / h;
}

const active& TwoDimensionalConvectionDispersionOperatorDG::axialDispersion(unsigned int idxSec, int idxRad, int idxComp) const CADET_NOEXCEPT
{
	return *(getSectionDependentSlice(_axialDispersion, _radNPoints * _nComp, idxSec) + idxRad * _nComp + idxComp);
}

const active& TwoDimensionalConvectionDispersionOperatorDG::radialDispersion(unsigned int idxSec, int idxRad, int idxComp) const CADET_NOEXCEPT
{
	return *(getSectionDependentSlice(_radialDispersion, _radNPoints * _nComp, idxSec) + idxRad * _nComp + idxComp);
}

/**
 * @brief Computes the residual of the transport equations
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] y Pointer to unit operation's state vector
 * @param [in] yDot Pointer to unit operation's time derivative state vector
 * @param [out] res Pointer to unit operation's residual vector
 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
 */
int TwoDimensionalConvectionDispersionOperatorDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, double* res, bool wantJac, WithoutParamSensitivity)
{
	return residualImpl<double, double, double>(model, t, secIdx, y, yDot, res);
}

int TwoDimensionalConvectionDispersionOperatorDG::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, bool wantJac, WithoutParamSensitivity)
{
	return residualImpl<active, active, double>(model, t, secIdx, y, yDot, res);
}

int TwoDimensionalConvectionDispersionOperatorDG::residual(const IModel& model, double t, unsigned int secIdx, double const* y, double const* yDot, active* res, bool wantJac, WithParamSensitivity)
{
	return residualImpl<double, active, active>(model, t, secIdx, y, yDot, res);
}

int TwoDimensionalConvectionDispersionOperatorDG::residual(const IModel& model, double t, unsigned int secIdx, active const* y, double const* yDot, active* res, bool wantJac, WithParamSensitivity)
{
	return residualImpl<active, active, active>(model, t, secIdx, y, yDot, res);
}

template <typename StateType, typename ResidualType, typename ParamType>
int TwoDimensionalConvectionDispersionOperatorDG::residualImpl(const IModel& model, double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res)
{
	const unsigned int offsetC = _radNPoints * _nComp;

	const int auxRadElemStride = _radNNodes;
	const int auxAxNodeStride = _radNElem * auxRadElemStride;
	const int auxAxElemStride = _axNNodes * auxAxNodeStride;
	const int auxTildeAxNodeStride = _radNElem * _qNNodes;
	const int auxTildeAxElemStride = _axNNodes * auxTildeAxNodeStride;

	const active* const curAxialDispersion = getSectionDependentSlice(_axialDispersion, _radNPoints * _nComp, 0);
	const active* const curRadialDispersion = getSectionDependentSlice(_radialDispersion, _radNPoints * _nComp, 0);

	for (unsigned int comp = 0; comp < _nComp; comp++)
	{
		// todo radial position dependent dispersion
		//Eigen::Map<const Vector<ParamType, Dynamic>, 0, InnerStride<Dynamic>> D_ax(curAxialDispersion, _radNPoints, InnerStride<Dynamic>(_nComp));
		const double D_ax = static_cast<double>(curAxialDispersion[0]);
		const double D_rad = static_cast<double>(curRadialDispersion[0]);

		/*	auxiliary equations	*/

		for (unsigned int zEidx = 0; zEidx < _axNElem; zEidx++)
		{
			for (unsigned int rEidx = 0; rEidx < _radNElem; rEidx++)
			{
				// todo allow ParamType : performance should be improved by adding this only in the main equation, not the auxiliary already
				const double DeltaZ = static_cast<double>(_axDelta);
				const double DeltaR = static_cast<double>(_radDelta[rEidx]);

				const int elemOffset = zEidx * _axElemStride + rEidx * _radElemStride;
				const int auxElemOffset = zEidx * auxAxElemStride + rEidx * auxRadElemStride;

				ConstMatrixMap<StateType> _C(y + offsetC + comp + elemOffset, _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(_radNodeStride, _axNodeStride));

				MatrixMap<StateType> _Gz(reinterpret_cast<StateType*>(&_axAuxStateG[0]) + auxElemOffset, _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(1, _radNPoints));
				MatrixMap<StateType> _Gr(reinterpret_cast<StateType*>(&_radAuxStateG[0]) + auxElemOffset, _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(1, _radNPoints));

				MatrixMap<StateType> _fAux1(reinterpret_cast<StateType*>(&_fStarAux1[0]), 2, _radNNodes, Stride<Dynamic, Dynamic>(1, _radNNodes));
				MatrixMap<StateType> _fAux2(reinterpret_cast<StateType*>(&_fStarAux2[0]), _axNNodes, 2, Stride<Dynamic, Dynamic>(1, 2));

				// axial auxiliary flux: central flux
				_fAux1.row(0) = _C.row(0);
				if (zEidx != 0) // else inlet auxiliary boundary condition fulfilled
				{
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> prevNodeZ(y + offsetC + comp + elemOffset - _axNodeStride, _radNNodes, InnerStride<Dynamic>(_radNodeStride));
					_fAux1.row(0) += prevNodeZ;
					_fAux1.row(0) *= 0.5;
				}

				_fAux1.row(1) = _C.row(_axPolyDeg);
				if (zEidx != _axNElem - 1) // else already correct as per BC
				{
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> nextNodeZ(y + offsetC + comp + elemOffset + _axElemStride, _radNNodes, InnerStride<Dynamic>(_radNodeStride));
					_fAux1.row(1) += nextNodeZ;
					_fAux1.row(1) *= 0.5;
				}

				// radial auxiliary flux: central flux
				if (rEidx == 0) // else already set from last iteration, see below
					_fAux2.col(0) = _C.col(0); // auxiliary bc

				_fAux2.col(1) = _C.col(_radPolyDeg);
				if (rEidx != _radNElem - 1) // else auxiliary bc already fulfilled as per last line
				{
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> nextNodeR(y + offsetC + comp + elemOffset + _radElemStride, _axNNodes, InnerStride<Dynamic>(_axNodeStride));
					_fAux2.col(1) += nextNodeR;
					_fAux2.col(1) *= 0.5;
				}

				// auxiliary equation g^z // todo block of G as its global. Or make it all lokal?
				// todo radial position dependent D_ax
				_Gz = 2.0 / DeltaZ * D_ax * _axInvMM.template cast<StateType>() * (_axLiftM.template cast<StateType>() * _fAux1 - _axTransStiffM.template cast<StateType>() * _C);

				// auxiliary equation g^r // todo block of G as its global. Or make it all lokal?
				 // todo radial position dependent D_rad
				_Gr = 2.0 / DeltaR * D_rad * (_fAux2 * _radLiftM.template cast<StateType>() - _C * _radStiffM.template cast<StateType>()) * _radInvTransMM.template cast<StateType>();

				// reuse "right" radial flux as "left" radial flux in next iteration
				_fAux2.col(0) = _fAux2.col(1);
			}
		}

		/*	main equation	*/

		// Add time derivative to bulk residual
		if (yDot)
		{
			Eigen::Map<const VectorXd, 0, InnerStride<Dynamic>> _cDot(yDot + offsetC + comp, _axNPoints * _radNPoints, InnerStride<Dynamic>(_radNodeStride));
			Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> _cRes(res + offsetC + comp, _axNPoints * _radNPoints, InnerStride<Dynamic>(_radNodeStride));
			_cRes = _cDot.template cast<ResidualType>();
		}
		else
		{
			Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> _cRes(res + offsetC + comp, _axNPoints * _radNPoints, InnerStride<Dynamic>(_radNodeStride));
			_cRes.setZero();
		}

		// note: auxiliary states have been computed without prefactor _delta[elem]

		// todo: reuse memory of fStarAux for gStarDisp

		// todo interpolate auxiliary variable and num. flux to quadrature nodes. For now, we assume these points to be collocated, ie LGL quadrature
		{
			ConstMatrixMap<StateType> _GzGlobal(reinterpret_cast<StateType*>(&_axAuxStateG[0]), _axNPoints, _radNPoints, Stride<Dynamic, Dynamic>(1, _radNPoints));
			ConstMatrixMap<StateType> _GrGlobal(reinterpret_cast<StateType*>(&_radAuxStateG[0]), _axNPoints, _radNPoints, Stride<Dynamic, Dynamic>(1, _radNPoints));

			MatrixMap<StateType> _GzInterpolatedGlobal(reinterpret_cast<StateType*>(&_axAuxStateGTilde[0]), _axNPoints, _qNNodes * _radNElem, Stride<Dynamic, Dynamic>(1, _qNNodes * _radNElem));
			MatrixMap<StateType> _GrInterpolatedGlobal(reinterpret_cast<StateType*>(&_radAuxStateGTilde[0]), _axNPoints, _qNNodes * _radNElem, Stride<Dynamic, Dynamic>(1, _qNNodes * _radNElem));

			_GzInterpolatedGlobal = _GzGlobal;
			_GrInterpolatedGlobal = _GrGlobal;

		}

		for (unsigned int zEidx = 0; zEidx < _axNElem; zEidx++)
		{
			for (unsigned int rEidx = 0; rEidx < _radNElem; rEidx++)
			{
				const int elemOffset = zEidx * _axElemStride + rEidx * _radElemStride;
				const int auxElemOffset = zEidx * auxAxElemStride + rEidx * _radNNodes;
				const int auxTildeElemOffset = zEidx * auxTildeAxElemStride + rEidx * _qNNodes;

				ConstMatrixMap<StateType> _C(y + offsetC + comp + elemOffset, _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(_radNodeStride, _axNodeStride));
				MatrixMap<ResidualType> _Res(res + offsetC + comp + elemOffset, _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(_radNodeStride, _axNodeStride));

				// Note: _GzTilde is actually just the first, non repeating block of _GzTilde from the manuscript. Thats why we perform the loop for the axial dispersion integral ( see axial dispersion step 1).
				// Otherwise, we'd have to construct a matrix with repeating values
				ConstMatrixMap<StateType> _GzTilde(reinterpret_cast<StateType*>(&_axAuxStateGTilde[0]) + auxTildeElemOffset, _axNNodes, _qNNodes, Stride<Dynamic, Dynamic>(1, _qNNodes * _radNElem));
				
				// Note: Same as for _GzTilde
				ConstMatrixMap<StateType> _GrTilde(reinterpret_cast<StateType*>(&_radAuxStateGTilde[0]) + auxTildeElemOffset, _axNNodes, _qNNodes, Stride<Dynamic, Dynamic>(1, _qNNodes * _radNElem));

				MatrixMap<ResidualType> _fStarConvZ(reinterpret_cast<ResidualType*>(&_fStarConv[0]), 2, _radNNodes, Stride<Dynamic, Dynamic>(1, _radNNodes));
				MatrixMap<StateType> _gStarDispTildeZ(reinterpret_cast<StateType*>(&_gZStarDispTilde[0]), 2, _qNNodes, Stride<Dynamic, Dynamic>(1, _qNNodes));
				MatrixMap<StateType> _gStarDispR(reinterpret_cast<StateType*>(&_gRStarDisp[0]), _axNNodes, 2, Stride<Dynamic, Dynamic>(1, 2));
				
				/*	numerical fluxes	*/

				// radial dispersion (central) flux
				{
					ConstMatrixMap<StateType> _Gr(reinterpret_cast<StateType*>(&_radAuxStateG[0]) + auxElemOffset, _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(1, _radNPoints));

					if (rEidx == 0) // else already set from last iteration, see below
						_gStarDispR.col(0).setZero(); // solid wall BC

					if (rEidx == _radNElem - 1)
						_gStarDispR.col(1).setZero();
					else
					{
						Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> GRnextNode(reinterpret_cast<StateType*>(&_radAuxStateG[0]) + auxElemOffset + auxRadElemStride, _axNNodes, InnerStride<Dynamic>(auxAxNodeStride));
						_gStarDispR.col(1) = 0.5 * (_Gr.col(_radPolyDeg) + GRnextNode);
					}

				}

				// axial dispersion (central) flux and axial convection (upwind) flux
				if (zEidx == 0) // Danckwerts inlet BC
				{
					_gStarDispTildeZ.row(0).setZero();
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> inlet(y + rEidx * _radElemStride + comp, _radNNodes, InnerStride<Dynamic>(_radNodeStride));

					//VectorXd inlet = VectorXd::Ones(_radNNodes); // todo delete

					_fStarConvZ.row(0) = -static_cast<ParamType>(_curVelocity[0]) * inlet.template cast<StateType>(); // todo radially dependent velocity
				}
				else
				{
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> GZTildePrevNode(reinterpret_cast<StateType*>(&_axAuxStateGTilde[0]) + auxTildeElemOffset - auxTildeAxNodeStride, _qNNodes, InnerStride<Dynamic>(1));
					_gStarDispTildeZ.row(0) = 0.5 * (_GzTilde.row(0) + GZTildePrevNode.transpose());

					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> prevNodeZ(y + offsetC + comp + elemOffset - _axNodeStride, _radNNodes, InnerStride<Dynamic>(_radNodeStride));
					_fStarConvZ.row(0) = -static_cast<ParamType>(_curVelocity[0]) * prevNodeZ; // todo radially dependent velocity
				}

				_fStarConvZ.row(1) = -static_cast<ParamType>(_curVelocity[0]) * _C.block(_axPolyDeg, 0, 1, _radNNodes); // todo radially dependent velocity
				if (zEidx != _axNElem - 1)
				{
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> GZTildeNextNode(reinterpret_cast<StateType*>(&_axAuxStateGTilde[0]) + auxTildeElemOffset + auxTildeAxElemStride, _qNNodes, InnerStride<Dynamic>(1));
					_gStarDispTildeZ.row(1) = 0.5 * (_GzTilde.row(_axPolyDeg) + GZTildeNextNode.transpose());
				}
				else // Danckwert outflow boundary condition
					_gStarDispTildeZ.row(1).setZero();

				// todo velocity radial position dependence!

				// transport residual (i.e. LHS - RHS)

				// Axial convection
				_Res -= 2.0 / static_cast<ParamType>(_axDelta) * _axInvMM.template cast<ResidualType>() * (
					// axial surface integral
					_axLiftM.template cast<ResidualType>() * _fStarConvZ
					// axial volume integral
					- _axTransStiffM.template cast<ResidualType>() * (-static_cast<ParamType>(_curVelocity[0]) * _C.template cast<ResidualType>())
					);

				// Axial dispersion
				
				// TODO old code for radial position dependent parameters : delete or fix
				//// Step 1: To get _GzTilde * _transTildeMrDash * _invTransMrCyl, we need to multiply a matrix A of dimensionality N x M with M x M blocks of a matrix B of dimensionality M*K x M.
				////		   Additionally Matrix A might have actives while B has doubles and the resulting matrix is immediately multiplied by another double matrix C with dimensionality N x M.
				////		   We thus first multiply B and C, then cast the result to A's scalar type and block-wise multiply by the result of B * C which again has dimensionality M*K x M.
				//MatrixMap<ResidualType> matrixCache(reinterpret_cast<ResidualType*>(&_matrixProductCache[0]), _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(1, _radNNodes));
				//matrixCache.setZero();
				//
				//for (int qNode = 0; qNode < _qNNodes; qNode++)
				//{
				//	matrixCache += _GzTilde.template cast<ResidualType>() * (_transTildeMrDash[rEidx].block(qNode * _radNNodes, 0, _radNNodes, _radNNodes) * _invTransMrCyl[rEidx]).template cast<ResidualType>();
				//}

				_Res -= 2.0 / static_cast<ParamType>(_axDelta) * _axInvMM.template cast<ResidualType>() * (
						// surface integral
					    _axLiftM.template cast<ResidualType>() * (
						_gStarDispTildeZ.template cast<ResidualType>()
						)
					    // volume integral
						-_axTransStiffM.template cast<ResidualType>() * (
							_GzTilde.template cast<ResidualType>()
							)
						);

				// Radial dispersion
				_Res -= 2.0 / static_cast<ParamType>(_radDelta[rEidx]) * (
					_gStarDispR.template cast<ResidualType>() * _radLiftMCyl[rEidx].template cast<ResidualType>()
					- _GrTilde.template cast<ResidualType>() * _SrCyl[rEidx].template cast<ResidualType>()
					) * _invTransMrCyl[rEidx].template cast<ResidualType>();

				// reuse "right" radial flux as "left" radial flux in next iteration
				_gStarDispR.col(0) = _gStarDispR.col(1);

			}
		}
	}

	return 0;
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
void TwoDimensionalConvectionDispersionOperatorDG::multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const
{
	double* localRet = ret + _nComp * _radNPoints;
	double const* localSdot = sDot + _nComp * _radNPoints;
	for (unsigned int i = 0; i < _axNPoints * _nComp * _radNPoints; ++i)
		localRet[i] = localSdot[i];
}

/**
 * @brief Adds the derivatives with respect to @f$ \dot{y} @f$ of @f$ F(t, y, \dot{y}) @f$ to the Jacobian
 * @details This functions computes
 *          @f[ \begin{align*} \text{_jacCdisc} = \text{_jacCdisc} + \alpha \frac{\partial F}{\partial \dot{y}}. \end{align*} @f]
 *          The factor @f$ \alpha @f$ is useful when constructing the linear system in the time integration process.
 * @param [in] alpha Factor in front of @f$ \frac{\partial F}{\partial \dot{y}} @f$
 */
void TwoDimensionalConvectionDispersionOperatorDG::addTimeDerivativeToJacobian(double alpha, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacDisc, unsigned int blockOffset)
{
	const int gapCell = _radNodeStride - static_cast<int>(_nComp); // only != 0 for LRM2D
	linalg::BandedEigenSparseRowIterator jac(jacDisc, blockOffset);

	for (unsigned int point = 0; point < _axNPoints * _radNPoints; ++point, jac += gapCell) {
		for (unsigned int comp = 0; comp < _nComp; ++comp, ++jac) {
			// dc_b / dt in transport equation
			jac[0] += alpha;
		}
	}
}

int TwoDimensionalConvectionDispersionOperatorDG::nJacEntries()
{
	int nAxDepBlocks = 0;
	int nRadDepBlocks = 0;

	for (int zElem = 0; zElem < _axNElem; zElem++)
	{
		const int nRightAxElemDep = std::min(2, static_cast<int>(_axNElem) - 1 - zElem); //<! number of right axial elements, this element depends on
		const int nLeftAxElemDep = std::min(2, zElem); //<! number of left axial elements, this element depends on

		nAxDepBlocks += nLeftAxElemDep + 1 + nRightAxElemDep;
	}

	for (int rElem = 0; rElem < _radNElem; rElem++)
	{
		const int nRightRadElemDep = std::min(2, static_cast<int>(_radNElem) - 1 - rElem); //<! number of right radial elements, this element depends on
		const int nLeftRadElemDep = std::min(2, rElem); //<! number of left radial elements, this element depends on

		nRadDepBlocks += nLeftRadElemDep + 1 + nRightRadElemDep;
	}

	// reaction block entries, which we inclde in the convDisp pattern so that its always of the same shape
	const int reacEntries = _radNodeStride;

	return _radNPoints * _axNPoints * _nComp * (nAxDepBlocks * _axNNodes + nRadDepBlocks * _radNNodes) + _radNPoints * _axNPoints * reacEntries * reacEntries;
}

/**
 * @brief Computes the individual transport Jacobian blocks
 */
bool TwoDimensionalConvectionDispersionOperatorDG::computeConvDispJacobianBlocks()
{
	// todo : should we just allocate the auxiliary/intermediate matrices here since this function is only called every time section?
	const int Np = _elemNPoints; //<! number of points per 2D element

	MatrixXd cStarDer = MatrixXd::Zero(Np, 5 * Np);
	cStarDer.block(0, Np + Np - _radNNodes, _radNNodes, _radNNodes) = MatrixXd::Identity(_radNNodes, _radNNodes);
	cStarDer.block(Np - _radNNodes, 2 * Np + Np - _radNNodes, _radNNodes, _radNNodes) = MatrixXd::Identity(_radNNodes, _radNNodes);

	MatrixXd cDer = MatrixXd::Zero(Np, 5 * Np);
	cDer.block(0, 2 * Np, Np, Np) = MatrixXd::Identity(Np, Np);

	// define lambda function for a concise notation, since the KroneckerProduct implementation does not accept expression such as A.inverse()
	auto kroneckerProduct = [](const MatrixXd& A, const MatrixXd& B, MatrixXd& C) {
		KroneckerProduct<MatrixXd, MatrixXd> kroneckerProductObj(A, B);
		kroneckerProductObj.evalTo(C);
		};

	const double deltaZ = static_cast<double>(_axDelta);

	/* auxiliary block axial dispersion */
	MatrixXd* GzDer = new MatrixXd[3]; // we have three unique auxiliary blocks
	MatrixXd fStarAux1 = MatrixXd::Zero(Np, 3 * Np);

	MatrixXd _radMM = dgtoolbox::mMatrix(_radPolyDeg, _radNodes, 0.0, 0.0);
	MatrixXd MzKronMrInv = MatrixXd::Zero(Np, Np);
	kroneckerProduct(_axInvMM.inverse(), _radMM, MzKronMrInv);
	MzKronMrInv = MzKronMrInv.inverse();

	MatrixXd BzKronMr = MatrixXd::Zero(Np, Np);
	kroneckerProduct(dgtoolbox::liftingMatrixQuadratic(_axNNodes), _radMM, BzKronMr);

	MatrixXd SzTKronMr = MatrixXd::Zero(Np, Np);
	kroneckerProduct(_axTransStiffM, _radMM, SzTKronMr);

	for (int i = 0; i < 3; i++) // we have three unique auxiliary blocks
	{
		if (i == 0)
			fStarAux1.block(0, Np, _radNNodes, _radNNodes) = MatrixXd::Identity(_radNNodes, _radNNodes);
		else
		{
			fStarAux1.block(0, Np, _radNNodes, _radNNodes) = 0.5 * MatrixXd::Identity(_radNNodes, _radNNodes);
			fStarAux1.block(0, Np - _radNNodes, _radNNodes, _radNNodes) = 0.5 * MatrixXd::Identity(_radNNodes, _radNNodes);
		}
		if (i == 2 || _axNElem == 1) // that is, in the special case of having one axial element, we store the auxiliary block at GzDer[0]
			fStarAux1.block(Np - _radNNodes, 2 * Np - _radNNodes, _radNNodes, _radNNodes) = MatrixXd::Identity(_radNNodes, _radNNodes);
		else
		{
			fStarAux1.block(Np - _radNNodes, 2 * Np - _radNNodes, _radNNodes, _radNNodes) = 0.5 * MatrixXd::Identity(_radNNodes, _radNNodes);
			fStarAux1.block(Np - _radNNodes, 2 * Np, _radNNodes, _radNNodes) = 0.5 * MatrixXd::Identity(_radNNodes, _radNNodes);
		}

		GzDer[i] = 2.0 / deltaZ * MzKronMrInv * (BzKronMr * fStarAux1 - SzTKronMr * cDer.block(0, Np, Np, 3 * Np));

		fStarAux1.setZero();
	}

	/* auxiliary block radial dispersion */

	MatrixXd* GrDer = new MatrixXd[3]; // we have three unique auxiliary blocks
	MatrixXd fStarAux2 = MatrixXd::Zero(Np, 3 * Np);

	MatrixXd MzKronBr = MatrixXd::Zero(Np, Np);
	MatrixXd _axMM = dgtoolbox::mMatrix(_axPolyDeg, _axNodes, 0.0, 0.0);
	kroneckerProduct(_axMM, dgtoolbox::liftingMatrixQuadratic(_radNNodes), MzKronBr);

	MatrixXd MzKronSrT = MatrixXd::Zero(Np, Np);
	kroneckerProduct(_axMM, _radStiffM.transpose(), MzKronSrT);

	{
		MatrixXd auxCder = MatrixXd::Zero(Np, 3 * Np);
		MatrixXd auxCderSubBlock = MatrixXd::Zero(_radNNodes, 3 * _radNNodes);
		auxCderSubBlock.block(0, _radNNodes, _radNNodes, _radNNodes) = MatrixXd::Identity(_radNNodes, _radNNodes);
		for (int j = 0; j < _axNNodes; j++)
			auxCder.block(j * _radNNodes, j * 3 * _radNNodes, _radNNodes, 3 * _radNNodes) = auxCderSubBlock;

		auto Yblock = [](const double a, const double b, const double c, const double d, const int size) {

			MatrixXd A = MatrixXd::Zero(size, 3 * size);
			A(0, size - 1) = a;
			A(0, size) = b;
			A(size - 1, 2 * size - 1) = c;
			A(size - 1, 2 * size) = d;

			return A;
			};

		for (int i = 0; i < 3; i++) // we have three unique auxiliary blocks
		{
			MatrixXd Y1 = Yblock(0.5, 0.5, 0.5, 0.5, _radNNodes);

			if (_radNElem == 1)
				Y1 = Yblock(0.0, 1.0, 1.0, 0.0, _radNNodes);
			else if (i == 0)
				Y1 = Yblock(0.0, 1.0, 0.5, 0.5, _radNNodes);
			else if (i == 2)
				Y1 = Yblock(0.5, 0.5, 1.0, 0.0, _radNNodes);

			for (int j = 0; j < _axNNodes; j++)
				fStarAux2.block(j * _radNNodes, j * 3 * _radNNodes, _radNNodes, 3 * _radNNodes) = Y1;

			// note: without deltaR, which will be added in final Jacobian block, so that we only have three unique auxiliary blocks here
			GrDer[i] = MzKronMrInv * (MzKronBr * fStarAux2 - MzKronSrT * auxCder);

			fStarAux2.setZero();
		}
	}

	const int uAxElem = std::min(5, static_cast<int>(_axNElem)); // number of unique axial Jacobian blocks (per radial element)

	/* Actual Jacobian blocks */
	for (int rElem = 0; rElem < _radNElem; rElem++)
	{
		MatrixXd MzKronMrCylInv = MatrixXd::Zero(Np, Np);
		kroneckerProduct(_axInvMM.inverse(), _transMrCyl[rElem].transpose(), MzKronMrCylInv);
		MzKronMrCylInv = MzKronMrCylInv.inverse();

		MatrixXd BzKronMrCyl = MatrixXd::Zero(Np, Np);
		kroneckerProduct(dgtoolbox::liftingMatrixQuadratic(_axNNodes), _transMrCyl[rElem].transpose(), BzKronMrCyl);

		MatrixXd SzTKronMrCyl = MatrixXd::Zero(Np, Np);
		kroneckerProduct(_axTransStiffM, _transMrCyl[rElem].transpose(), SzTKronMrCyl);

		/* convection block */
		const double u = static_cast<double>(_curVelocity[0]); // @todo

		_jacConvection[rElem] = MzKronMrCylInv * 2.0 / deltaZ * u * (SzTKronMrCyl * cDer - BzKronMrCyl * cStarDer);

		// filter numerical noise
		_jacConvection[rElem] = _jacConvection[rElem].unaryExpr([](double val) {
			return (std::abs(val) < 1e-14) ? 0.0 : val;
			});

		/* axial dispersion block */
		const active* const curAxialDispersion = getSectionDependentSlice(_axialDispersion, _radNPoints * _nComp, 0);
		const double axDisp = static_cast<double>(curAxialDispersion[0]);

		MatrixXd gStarZ = MatrixXd::Zero(Np, 5 * Np);

		for (int i = 0; i < uAxElem; i++)
		{
			// first, we need to create numerical flux block gStarZ
			// note: we have three unique auxiliary block indices and need to find the indices wrt the currently considered element
			int auxIdx;
			if (i == 0)
				auxIdx = 0; // note that if _axNElem == 0, the corresponding auxiliary block is also stored at index 0
			else
				auxIdx = i == uAxElem - 1 ? 2 : 1;

			if (i > 0) // left boundary condition -> zero
			{
				const int leftAuxIdx = (i == 1) ? 0 : 1; // if the left neighbour is the left boundary element, set index to 0, else to 1
				gStarZ.block(0, 0, _radNNodes, 3 * Np) += GzDer[leftAuxIdx].block(Np - _radNNodes, 0, _radNNodes, 3 * Np);
				gStarZ.block(0, Np, _radNNodes, 3 * Np) += GzDer[auxIdx].block(0, 0, _radNNodes, 3 * Np);
			}

			if (i != uAxElem - 1) // right boundary condition -> zero
			{
				// if the right neighbour is the right boundary element, set index to 2, else to 1
				const int rightAuxIdx = (i + 1 == uAxElem - 1) ? 2 : 1;
				gStarZ.block(Np - _radNNodes, Np, _radNNodes, 3 * Np) += GzDer[auxIdx].block(Np - _radNNodes, 0, _radNNodes, 3 * Np);
				gStarZ.block(Np - _radNNodes, 2 * Np, _radNNodes, 3 * Np) += GzDer[rightAuxIdx].block(0, 0, _radNNodes, 3 * Np);
			}

			gStarZ *= 0.5;

			_jacAxDispersion[rElem * uAxElem + i] = MzKronMrCylInv * 2.0 / deltaZ * axDisp * (BzKronMrCyl * gStarZ);
			_jacAxDispersion[rElem * uAxElem + i].block(0, Np, Np, 3 * Np) -= MzKronMrCylInv * 2.0 / deltaZ * axDisp * (SzTKronMrCyl * GzDer[auxIdx]);

			// filter numerical noise
			_jacAxDispersion[rElem * uAxElem + i] = _jacAxDispersion[rElem * uAxElem + i].unaryExpr([](double val) {
				return (std::abs(val) < 1e-14) ? 0.0 : val;
				});

			gStarZ.setZero();
		}

		/* radial dispersion block */
		const active* const curRadialDispersion = getSectionDependentSlice(_radialDispersion, _radNPoints * _nComp, 0);
		const double radDisp = static_cast<double>(curRadialDispersion[0]);

		MatrixXd gStarR = MatrixXd::Zero(Np, 5 * Np);

		for (int i = 0; i < _axNNodes; i++)
		{
			if (rElem > 0) // else zero as per boundary condition
			{
				const int ll = rElem == 1 ? 0 : 1;
				const int lr = rElem == _radNElem - 1 ? 2 : 1;

				gStarR.block(i * _radNNodes, i * 5 * _radNNodes, 1, 3 * _radNNodes) = 0.5 * GrDer[ll].block((i + 1) * _radNNodes - 1, i * 3 * _radNNodes, 1, 3 * _radNNodes);
				gStarR.block(i * _radNNodes, i * 5 * _radNNodes + _radNNodes, 1, 3 * _radNNodes) += 0.5 * GrDer[lr].block(i * _radNNodes, i * 3 * _radNNodes, 1, 3 * _radNNodes);
			}
			if (rElem < _radNElem - 1) // else zero as per boundary condition
			{
				const int rl = rElem == 0 ? 0 : 1;
				const int rr = rElem + 1 == _radNElem - 1 ? 2 : 1;
				gStarR.block(i * _radNNodes + _radPolyDeg, i * 5 * _radNNodes + _radNNodes, 1, 3 * _radNNodes) = 0.5 * GrDer[rl].block((i + 1) * _radNNodes - 1, i * 3 * _radNNodes, 1, 3 * _radNNodes);
				gStarR.block(i * _radNNodes + _radPolyDeg, i * 5 * _radNNodes + 2 * _radNNodes, 1, 3 * _radNNodes) += 0.5 * GrDer[rr].block(i * _radNNodes, i * 3 * _radNNodes, 1, 3 * _radNNodes);
			}
		}

		MatrixXd _BrCyl = dgtoolbox::liftingMatrixQuadratic(_radNNodes);
		_BrCyl(0, 0) = -static_cast<double>(_radialElemInterfaces[rElem]);
		_BrCyl(_radPolyDeg, _radPolyDeg) = static_cast<double>(_radialElemInterfaces[rElem + 1]);

		MatrixXd MzKronBrCyl = MatrixXd::Zero(Np, Np);
		kroneckerProduct(_axMM, _BrCyl, MzKronBrCyl);

		MatrixXd MzKronSrTCyl = MatrixXd::Zero(Np, Np);
		kroneckerProduct(_axMM, _SrCyl[rElem].transpose(), MzKronSrTCyl);

		int auxIdx;
		if (rElem == 0)
			auxIdx = 0; // note that if _radNElem == 0, the corresponding auxiliary block is also stored at index 0
		else
			auxIdx = rElem == _radNElem - 1 ? 2 : 1;

		// note: accounted for deltaR twice since it was not included in auxiliary block
		_jacRadDispersion[rElem] = 2.0 / static_cast<double>(_radDelta[rElem]) * 2.0 / static_cast<double>(_radDelta[rElem]) * radDisp * MzKronMrCylInv * (MzKronBrCyl * gStarR);
		for (int zNode = 0; zNode < _axNNodes; zNode++) // todo can this be done more elegantly?
			_jacRadDispersion[rElem].block(0, (5 * zNode + 1) * _radNNodes, Np, 3 * _radNNodes) -= 2.0 / static_cast<double>(_radDelta[rElem]) * 2.0 / static_cast<double>(_radDelta[rElem]) * radDisp * MzKronMrCylInv * (MzKronSrTCyl * GrDer[auxIdx].block(0, zNode * 3 * _radNNodes, Np, 3 * _radNNodes));

		// filter numerical noise
		_jacRadDispersion[rElem] = _jacRadDispersion[rElem].unaryExpr([](double val) {
			return (std::abs(val) < 1e-14) ? 0.0 : val;
			});

		gStarR.setZero();
	}

	return 1;
}

typedef Eigen::Triplet<double> T;

/**
 * @brief adds an element block for all components to the system Jacobian
 * @detail the element block has axNNodes * radNNodes rows and arbitrary number of columns
 * @param [in] block to be added
 * @param [in] offRow offset to first considered row
 * @param [in] offColumn column to row offset
 * @param [in] depElem number of elements the element block element depends on
 * @param [in] addEntry function to add entry to Jacobian or sparsity list
 */
template <typename Action>
void TwoDimensionalConvectionDispersionOperatorDG::addAxElemBlockToJac(const Eigen::MatrixXd& block, const int offRow, const int offColumn, const int depElem, Action addEntry)
{
	const int depBlockStride = _axElemStride;
	int jac = offRow;

	for (unsigned int zNode = 0; zNode < _axNNodes; zNode++, jac += _axNodeStride - _radElemStride) // move the iterator to the next axial node, but account for row changes within the loop
	{
		for (unsigned int rNode = 0; rNode < _radNNodes; rNode++, jac += _strideBound) // move the iterator to the next radial node, but account for row changes within the loop
		{
			// inside the loop, the iterator is moved over all components
			for (unsigned int comp = 0; comp < _nComp; comp++, ++jac) // insert the entries for all components
			{
				for (unsigned int depBlock = 0; depBlock < depElem; depBlock++) // iterate over all axial block dependencies
				{
					for (unsigned int i = 0; i < _axNNodes; i++) // iterate over all axial node dependencies
					{
						for (unsigned int j = 0; j < _radNNodes; j++) // iterate over all radial node dependencies
						{
							const double entry = block(zNode * _radNNodes + rNode, depBlock * _elemNPoints + i * _radNNodes + j);
							if(std::abs(entry) > 1e-14)
								// row: at current node and component
								// col: add offset, go to current element, jump to current node
								addEntry(jac, offRow + offColumn + depBlock * depBlockStride + j * _radNodeStride + i * _axNodeStride + comp, entry);
						}
					}
				}
			}
		}
	}
}
/**
 * @brief adds an element block for all components to the system Jacobian
 * @detail the element block has axNNodes * radNNodes rows and arbitrary number of columns
 * @param [in] block to be added
 * @param [in] jacobian jacobian matrix
 * @param [in] offRow offset to first considered row
 * @param [in] offColumn column to row offset
 * @param [in] depElem number of elements the element block element depends on
 */
void TwoDimensionalConvectionDispersionOperatorDG::addAxElemBlockToJac(const Eigen::MatrixXd& block, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, const int offRow, const int offColumn, const int depElem)
{
	auto addEntry = [&jacobian](int row, int col, double value) {
		jacobian.coeffRef(row, col) += value;
		};

	addAxElemBlockToJac(block, offRow, offColumn, depElem, addEntry);
}
/**
 * @brief adds an element block for all components to the system Jacobian
 * @detail the element block has axNNodes * radNNodes rows and arbitrary number of columns
 * @param [in] block to be added
 * @param [in] jacobian sparsity triplet list
 * @param [in] offRow offset to first considered row
 * @param [in] offColumn column to row offset
 * @param [in] depElem number of elements the element block element depends on
 */
void TwoDimensionalConvectionDispersionOperatorDG::addAxElemBlockToJac(const Eigen::MatrixXd& block, std::vector<T>& tripletList, const int offRow, const int offColumn, const int depElem)
{
	auto addEntry = [&tripletList](int row, int col, double value) {
		//if (std::abs(value) > 1e-14)
			tripletList.push_back(T(row, col, 0.0));
		};

	addAxElemBlockToJac(block, offRow, offColumn, depElem, addEntry);
}
/**
 * @brief adds an element block for all components to the system Jacobian
 * @detail the element block has axNNodes * radNNodes rows and arbitrary number of columns
 * @param [in] block to be added
 * @param [in] offRow offset to first considered row
 * @param [in] nLeftRadElemDep number of left radial elements this block depends on
 * @param [in] depElem number of elements the element block element depends on
 * @param [in] addEntry function to add entry to Jacobian or sparsity list
 */
template <typename Action>
void TwoDimensionalConvectionDispersionOperatorDG::addRadElemBlockToJac(const Eigen::MatrixXd& block, const int offRow, const int nLeftRadElemDep, const int depElem, Action addEntry)
{
	const int depBlockStride = _radElemStride; // stride in Jacobian per dependence block
	const int offColumn = -nLeftRadElemDep * _radElemStride; // offset in Jacobian due to left element neighbours dependence
	const int offBlock = (2 - nLeftRadElemDep) * _radNNodes; // offset in block due to left element neighbours dependence
	int jac = offRow;

	for (unsigned int zNode = 0; zNode < _axNNodes; zNode++, jac += _axNodeStride - _radElemStride) // move the iterator to the next axial node, but account for row changes within the loop
	{
		for (unsigned int rNode = 0; rNode < _radNNodes; rNode++, jac += _strideBound) // move the iterator to the next radial node, but account for row changes within the loop
		{
			// inside the loop, the iterator is moved over all components
			for (unsigned int comp = 0; comp < _nComp; comp++, ++jac) // insert the entries for all components
			{
				for (unsigned int i = 0; i < _axNNodes; i++) // iterate over all axial node dependencies
				{
					for (unsigned int depBlock = 0; depBlock < depElem; depBlock++) // iterate over all radial block dependencies
					{
						for (unsigned int j = 0; j < _radNNodes; j++) // iterate over all radial node dependencies
						{
							const double entry = block(zNode * _radNNodes + rNode, offBlock + i * 5 * _radNNodes + depBlock * _radNNodes + j);
							if (std::abs(entry) > 1e-14)
								// row: at current node and component
								// col: add offset, go to current element, jump to current node
								addEntry(jac, offRow + offColumn + depBlock * depBlockStride + j * _radNodeStride + i * _axNodeStride + comp, entry);
						}
					}
				}
			}
		}
	}
}
/**
 * @brief adds an element block for all components to the system Jacobian
 * @detail the element block has axNNodes * radNNodes rows and arbitrary number of columns
 * @param [in] block to be added
 * @param [in] jacobian jacobian matrix
 * @param [in] offRow offset to first considered row
 * @param [in] nLeftRadElemDep number of left radial elements this block depends on
 * @param [in] depElem number of elements the element block element depends on
 */
void TwoDimensionalConvectionDispersionOperatorDG::addRadElemBlockToJac(const Eigen::MatrixXd& block, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, const int offRow, const int nLeftRadElemDep, const int depElem)
{
	auto addEntry = [&jacobian](int row, int col, double value) {
		jacobian.coeffRef(row, col) += value;
		};

	addRadElemBlockToJac(block, offRow, nLeftRadElemDep, depElem, addEntry);
}
/**
 * @brief adds an element block for all components to the system Jacobian
 * @detail the element block has axNNodes * radNNodes rows and arbitrary number of columns
 * @param [in] block to be added
 * @param [in] jacobian sparsity triplet list
 * @param [in] offRow offset to first considered row
 * @param [in] nLeftRadElemDep number of left radial elements this block depends on
 * @param [in] depElem number of elements the element block element depends on
 */
void TwoDimensionalConvectionDispersionOperatorDG::addRadElemBlockToJac(const Eigen::MatrixXd& block, std::vector<T>& tripletList, const int offRow, const int nLeftRadElemDep, const int depElem)
{
	auto addEntry = [&tripletList](int row, int col, double value) {
		//if (std::abs(value) > 1e-14)
			tripletList.push_back(T(row, col, 0.0));
		};

	addRadElemBlockToJac(block, offRow, nLeftRadElemDep, depElem, addEntry);
}
/**
 * @brief defines the convection dispersion Jacobian sparsity pattern in a tripletlist
 * @detail additionally sets the bulk reaction pattern. This way, the strides are always the same for the Bulk block
 * @param [in,out] tripletList triplet list with Jacobian entries and row, column positions
 * @param [in] bulkOffset offset to first bulk entry
 */
void TwoDimensionalConvectionDispersionOperatorDG::convDispJacPattern(std::vector<T>& tripletList, const int bulkOffset)
{
	// note that inlet DOFs are not included here as we define a dedicated inlet Jacobian in the unit operation

	const int Np = _elemNPoints; //<! number of points per 2D element
	const int uAxElem = std::min(5, static_cast<int>(_axNElem)); // number of unique axial Jacobian blocks (per radial element)

	for (int zElem = 0; zElem < _axNElem; zElem++)
	{
		for (int rElem = 0; rElem < _radNElem; rElem++)
		{
			const int offSetRow = bulkOffset + zElem * _axElemStride + rElem * _radElemStride;

			/* Axial dispersion pattern (which the convection pattern is a subset of) */
			{
				const int nRightAxElemDep = std::min(2, static_cast<int>(_axNElem) - 1 - zElem); //<! number of right axial elements, this element depends on
				const int nLeftAxElemDep = std::min(2, zElem); //<! number of left axial elements, this element depends on
				const int offSetColumnToRow = -nLeftAxElemDep * _axElemStride;
				int uAxBlockIdx = zElem; // unique block index
				if (zElem > 2)
				{
					if (zElem + 2 == _axNElem || _axNElem == 4) // next element is last element or special case of four elements
						uAxBlockIdx = 3;
					else if (zElem + 1 == _axNElem) // this element is the last element and we have at least five elements
						uAxBlockIdx = 4;
					else // this element has at least two neighbours in each axial direction
						uAxBlockIdx = 2;
				}

				addAxElemBlockToJac(-_jacAxDispersion[rElem * uAxElem + uAxBlockIdx].block(0, Np * (2 - nLeftAxElemDep), Np, Np * (nLeftAxElemDep + 1 + nRightAxElemDep)), tripletList, offSetRow, offSetColumnToRow, nLeftAxElemDep + 1 + nRightAxElemDep);
			}

			/* Radial dispersion pattern */
			const int nLeftRadElem = std::min(2, rElem);
			const int nRightRadElem = std::min(2, static_cast<int>(_radNElem) - 1 - rElem);
			addRadElemBlockToJac(-_jacRadDispersion[rElem], tripletList, offSetRow, nLeftRadElem, nLeftRadElem + 1 + nRightRadElem);
		}
	}

	/* Bulk reaction / binding pattern */
	// Also needs to be set no matter if we actually have reactions/bindings, so that we can assume the bulk pattern always be the same,
	// which we assume in the addAxElemBlockToJac and addRadElemBlockToJac functions.

	for (int bulk = 0; bulk < _axNElem * _axElemStride; bulk++)
	{
		for (int conc = 0; conc < _radNodeStride; conc++)
		{
			for (int concDep = 0; concDep < _radNodeStride; concDep++)
			{
				tripletList.push_back(
					T(bulkOffset + bulk * _radNodeStride + conc,
						bulkOffset + bulk * _radNodeStride + concDep,
						0.0)
				);
			}
		}
	}
}
/**
 * @brief Assembles the transport Jacobian
 * @param [in] jacobian Jacobian matrix of unit in sparse format
 * @param [in] jacInlet Inlet Jacobian matrix
 * @param [in] bulkOffset Offset to first entry of bulk phase in unit Jacobian, defaults to 0
 */
bool TwoDimensionalConvectionDispersionOperatorDG::assembleConvDispJacobian(Eigen::SparseMatrix<double, RowMajor>& jacobian, Eigen::MatrixXd& jacInlet, const int bulkOffset)
{
	const int Np = _elemNPoints; //<! number of points per 2D element
	const int uAxElem = std::min(5, static_cast<int>(_axNElem)); // number of unique axial Jacobian blocks (per radial element)

	for (int zElem = 0; zElem < _axNElem; zElem++)
	{
		// Note: Jacobian blocks *-1 for residual

		for (int rElem = 0; rElem < _radNElem; rElem++)
		{
			/* handle axial convection Jacobian */

			// inlet Jacobian
			if (zElem == 0)
			{
				for (int rNode = 0; rNode < _radNNodes; rNode++)
				{
					for (int zNode = 0; zNode < _axNNodes; zNode++)
					{
						for (int j = 0; j < _radNNodes; j++)
						{
							for (int comp = 0; comp < _nComp; comp++)
								jacInlet(rElem * _radElemStride + rNode * _radNodeStride + zNode * _axNodeStride + comp, rElem * _radElemStride + j * _radNodeStride + comp) += -_jacConvection[rElem].block(0, Np + Np - _radNNodes, Np, _radNNodes)(rNode + zNode * _radNNodes, j);
						}
					}
				}
			}
			// "inner" Jacobian
			const int nRightAxElemDep = std::min(2, static_cast<int>(_axNElem) - 1 - zElem); //<! number of right axial elements, this element depends on
			const int nLeftAxElemDep = std::min(2, zElem); //<! number of left axial elements, this element depends on
			const int offSetRow = bulkOffset + zElem * _axElemStride + rElem * _radElemStride;
			const int offSetColumnToRow = -nLeftAxElemDep * _axElemStride;

			addAxElemBlockToJac(-_jacConvection[rElem].block(0, Np * (2 - nLeftAxElemDep), Np, Np * (nLeftAxElemDep + 1 + nRightAxElemDep)), jacobian, offSetRow, offSetColumnToRow, nLeftAxElemDep + 1 + nRightAxElemDep);

			/* handle axial dispersion Jacobian */
			{
				int uAxBlockIdx = zElem; // unique block index
				if (zElem > 2)
				{
					if (zElem + 2 == _axNElem || _axNElem == 4) // next element is last element or special case of four elements
						uAxBlockIdx = 3;
					else if (zElem + 1 == _axNElem) // this element is the last element and we have at least five elements
						uAxBlockIdx = 4;
					else // this element has at least two neighbours in each axial direction
						uAxBlockIdx = 2;
				}

				addAxElemBlockToJac(-_jacAxDispersion[rElem * uAxElem + uAxBlockIdx].block(0, Np * (2 - nLeftAxElemDep), Np, Np * (nLeftAxElemDep + 1 + nRightAxElemDep)), jacobian, offSetRow, offSetColumnToRow, nLeftAxElemDep + 1 + nRightAxElemDep);
			}

			/* handle radial dispersion Jacobian */
			// @todo ? we could also think about handling axial and radial Jacobian at the same time, which requires a different insertion function
			const int nLeftRadElem = std::min(2, rElem);
			const int nRightRadElem = std::min(2, static_cast<int>(_radNElem) - 1 - rElem);
			addRadElemBlockToJac(-_jacRadDispersion[rElem], jacobian, offSetRow, nLeftRadElem, nLeftRadElem + 1 + nRightRadElem);

		}

	}

	return true;
}

void TwoDimensionalConvectionDispersionOperatorDG::setEquidistantRadialDisc()
{
	const active h = _colRadius / _radNElem;
	const double pi = 3.1415926535897932384626434;
	active subcellLeftEnd = 0.0;
	active subcellRightEnd = 0.0;

	std::fill(_radDelta.begin(), _radDelta.end(), h);

	_radialElemInterfaces[0] = 0.0;
	for (unsigned int r = 0; r < _radNElem; ++r)
	{
		// Set last edge to _colRadius for exact geometry
		if (r == _radNElem - 1)
			_radialElemInterfaces[r + 1] = _colRadius;
		else
			_radialElemInterfaces[r + 1] = h * (r + 1);

		subcellLeftEnd = dgtoolbox::mapRefToPhys<active>(_radDelta, r, -1.0);
		VectorXd radWeights = _radInvWeights.cwiseInverse();
		for (unsigned int node = 0; node < _radNNodes; ++node)
		{
			subcellRightEnd = subcellLeftEnd + radWeights[node] / 2.0 * _radDelta[r];
			_nodalCrossSections[r * _radNNodes + node] = pi * (pow(subcellRightEnd, 2.0) - pow(subcellLeftEnd, 2.0));
			subcellLeftEnd = subcellRightEnd;
		}
	}
}

void TwoDimensionalConvectionDispersionOperatorDG::setEquivolumeRadialDisc()
{
	const active volPerElement = _colRadius * _colRadius / _radNElem;
	const double pi = 3.1415926535897932384626434;
	active subcellLeftEnd = 0.0;
	active subcellRightEnd = 0.0;

	_radialElemInterfaces[0] = 0.0;
	for (unsigned int r = 0; r < _radNElem; ++r)
	{
		// Set last edge to _colRadius for exact geometry
		if (r == _radNElem - 1)
			_radialElemInterfaces[r+1] = _colRadius;
		else
			_radialElemInterfaces[r+1] = sqrt(volPerElement + _radialElemInterfaces[r] * _radialElemInterfaces[r]);

		_radDelta[r] = _radialElemInterfaces[r + 1] - _radialElemInterfaces[r];

		subcellLeftEnd = dgtoolbox::mapRefToPhys<active>(_radDelta, r, -1.0);
		VectorXd radWeights = _radInvWeights.cwiseInverse();
		for (unsigned int node = 0; node < _radNNodes; ++node)
		{
			subcellRightEnd = subcellLeftEnd + radWeights[node] / 2.0 * _radDelta[r];
			_nodalCrossSections[r * _radNNodes + node] = pi * (pow(subcellRightEnd, 2.0) - pow(subcellLeftEnd, 2.0));
			subcellLeftEnd = subcellRightEnd;
		}
	}
}

void TwoDimensionalConvectionDispersionOperatorDG::setUserdefinedRadialDisc()
{
	const double pi = 3.1415926535897932384626434;
	active subcellLeftEnd = 0.0;
	active subcellRightEnd = 0.0;

	for (unsigned int r = 0; r < _radNElem; ++r)
	{
		_radDelta[r] = _radialElemInterfaces[r + 1] - _radialElemInterfaces[r];

		subcellLeftEnd = dgtoolbox::mapRefToPhys<active>(_radDelta, r, -1.0);
		VectorXd radWeights = _radInvWeights.cwiseInverse();
		for (unsigned int node = 0; node < _radNNodes; ++node)
		{
			subcellRightEnd = subcellLeftEnd + radWeights[node] / 2.0 * _radDelta[r];
			_nodalCrossSections[r * _radNNodes + node] = pi * (pow(subcellRightEnd, 2.0) - pow(subcellLeftEnd, 2.0));
			subcellLeftEnd = subcellRightEnd;
		}
	}
}

void TwoDimensionalConvectionDispersionOperatorDG::updateRadialDisc()
{
	if (_radialDiscretizationMode == RadialDiscretizationMode::Equidistant)
		setEquidistantRadialDisc();
	else if (_radialDiscretizationMode == RadialDiscretizationMode::Equivolume)
		setEquivolumeRadialDisc();
	else if (_radialDiscretizationMode == RadialDiscretizationMode::UserDefined)
		setUserdefinedRadialDisc();
}

bool TwoDimensionalConvectionDispersionOperatorDG::setParameter(const ParameterId& pId, double value)
{
	if (_singlePorosity && (pId.name == hashString("COL_POROSITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep)
		&& (pId.reaction == ReactionIndep) && (pId.section == SectionIndep) && (pId.particleType == ParTypeIndep))
	{
		for (unsigned int i = 0; i < _radNPoints; ++i)
			_colPorosities[i].setValue(value);
		return true;
	}

	if (_singleVelocity && (pId.name == hashString("VELOCITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep) && (pId.reaction == ReactionIndep))
	{
		if (_velocity.size() > _radNPoints)
		{
			// Section dependent
			if (pId.section == SectionIndep)
				return false;

			for (unsigned int i = 0; i < _radNPoints; ++i)
				_velocity[pId.section * _radNPoints + i].setValue(value);
		}
		else
		{
			// Section independent
			if (pId.section != SectionIndep)
				return false;

			for (unsigned int i = 0; i < _radNPoints; ++i)
				_velocity[i].setValue(value);
		}
	}

	const bool ad = multiplexParameterValue(pId, hashString("COL_DISPERSION_AXIAL"), _axialDispersionMode, _axialDispersion, _nComp, _radNPoints, value, nullptr);
	if (ad)
		return true;

	const bool adr = multiplexParameterValue(pId, hashString("COL_DISPERSION_RADIAL"), _radialDispersionMode, _radialDispersion, _nComp, _radNPoints, value, nullptr);
	if (adr)
		return true;

	return false;
}

bool TwoDimensionalConvectionDispersionOperatorDG::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
{
	if (_singlePorosity && (pId.name == hashString("COL_POROSITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep)
		&& (pId.reaction == ReactionIndep) && (pId.section == SectionIndep) && (pId.particleType == ParTypeIndep))
	{
		if (contains(sensParams, &_colPorosities[0]))
		{
			for (unsigned int i = 0; i < _radNPoints; ++i)
				_colPorosities[i].setValue(value);
			return true;
		}
	}

	if (_singleVelocity && (pId.name == hashString("VELOCITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep) && (pId.reaction == ReactionIndep))
	{
		if (_velocity.size() > _radNPoints)
		{
			// Section dependent
			if ((pId.section == SectionIndep) || !contains(sensParams, &_velocity[pId.section * _radNPoints]))
				return false;

			for (unsigned int i = 0; i < _radNPoints; ++i)
				_velocity[pId.section * _radNPoints + i].setValue(value);
		}
		else
		{
			// Section independent
			if ((pId.section != SectionIndep) || !contains(sensParams, &_velocity[0]))
				return false;

			for (unsigned int i = 0; i < _radNPoints; ++i)
				_velocity[i].setValue(value);
		}
	}

	const bool ad = multiplexParameterValue(pId, hashString("COL_DISPERSION_AXIAL"), _axialDispersionMode, _axialDispersion, _nComp, _radNPoints, value, &sensParams);
	if (ad)
		return true;

	const bool adr = multiplexParameterValue(pId, hashString("COL_DISPERSION_RADIAL"), _radialDispersionMode, _radialDispersion, _nComp, _radNPoints, value, &sensParams);
	if (adr)
		return true;

	return false;
}

bool TwoDimensionalConvectionDispersionOperatorDG::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if (_singlePorosity && (pId.name == hashString("COL_POROSITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep)
		&& (pId.reaction == ReactionIndep) && (pId.section == SectionIndep) && (pId.particleType == ParTypeIndep))
	{
		sensParams.insert(&_colPorosities[0]);
		for (unsigned int i = 0; i < _radNPoints; ++i)
			_colPorosities[i].setADValue(adDirection, adValue);

		return true;
	}

	if (_singleVelocity && (pId.name == hashString("VELOCITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep) && (pId.reaction == ReactionIndep))
	{
		if (_velocity.size() > _radNPoints)
		{
			// Section dependent
			if (pId.section == SectionIndep)
				return false;

			sensParams.insert(&_velocity[pId.section * _radNPoints]);
			for (unsigned int i = 0; i < _radNPoints; ++i)
				_velocity[pId.section * _radNPoints + i].setADValue(adDirection, adValue);
		}
		else
		{
			// Section independent
			if (pId.section != SectionIndep)
				return false;

			sensParams.insert(&_velocity[0]);
			for (unsigned int i = 0; i < _radNPoints; ++i)
				_velocity[i].setADValue(adDirection, adValue);
		}
	}

	const bool ad = multiplexParameterAD(pId, hashString("COL_DISPERSION_AXIAL"), _axialDispersionMode, _axialDispersion, _nComp, _radNPoints, adDirection, adValue, sensParams);
	if (ad)
		return true;

	const bool adr = multiplexParameterAD(pId, hashString("COL_DISPERSION_RADIAL"), _radialDispersionMode, _radialDispersion, _nComp, _radNPoints, adDirection, adValue, sensParams);
	if (adr)
		return true;

	return false;
}

}  // namespace parts

}  // namespace model

}  // namespace cadet

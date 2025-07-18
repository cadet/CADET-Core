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

cadet::model::MultiplexMode readAndRegisterMultiplexParam(cadet::IParameterProvider& paramProvider, std::unordered_map<cadet::ParameterId, cadet::active*>& parameters, std::vector<cadet::active>& values, const std::string& name, unsigned int nComp, unsigned int radNElem, cadet::UnitOpIdx uoi)
{
	cadet::model::MultiplexMode mode = cadet::model::MultiplexMode::Independent;
	readParameterMatrix(values, paramProvider, name, nComp * radNElem, 1);
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
			if (values.size() != radNElem)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(radNElem) + ")");
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
			if (values.size() != nComp * radNElem)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nComp * radNElem) + ")");
		}
		else if (modeConfig == 4)
		{
			mode = cadet::model::MultiplexMode::Section;
			nSec = values.size();
		}
		else if (modeConfig == 5)
		{
			mode = cadet::model::MultiplexMode::RadialSection;
			if (values.size() % radNElem != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " is not a positive multiple of radNElem (" + std::to_string(radNElem) + ")");

			nSec = values.size() / radNElem;
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
			if (values.size() % (nComp * radNElem) != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " is not a positive multiple of NCOMP * radNElem (" + std::to_string(nComp * radNElem) + ")");

			nSec = values.size() / (nComp * radNElem);
		}
	}
	else
	{
		if (values.size() == 1)
			mode = cadet::model::MultiplexMode::Independent;
		else if (values.size() == nComp)
			mode = cadet::model::MultiplexMode::Component;
		else if (values.size() == radNElem)
			mode = cadet::model::MultiplexMode::Radial;
		else if (values.size() == radNElem * nComp)
			mode = cadet::model::MultiplexMode::ComponentRadial;
		else if (values.size() % nComp == 0)
		{
			mode = cadet::model::MultiplexMode::ComponentSection;
			nSec = values.size() / nComp;
		}
		else if (values.size() % radNElem == 0)
		{
			mode = cadet::model::MultiplexMode::RadialSection;
			nSec = values.size() / radNElem;
		}
		else if (values.size() % (radNElem * nComp) == 0)
		{
			mode = cadet::model::MultiplexMode::ComponentRadialSection;
			nSec = values.size() / (nComp * radNElem);
		}
		else
			throw cadet::InvalidParameterException("Could not infer multiplex mode of field " + name + ", set " + name + "_MULTIPLEX or change number of elements");

		// Do not infer cadet::model::MultiplexMode::Section in case of no matches (might hide specification errors)
	}

	const cadet::StringHash nameHash = cadet::hashStringRuntime(name);
	switch (mode)
	{
		case cadet::model::MultiplexMode::Independent: // nSec = 1
		case cadet::model::MultiplexMode::Section:
			{
				std::vector<cadet::active> p(nComp * radNElem * nSec);
				for (unsigned int s = 0; s < nSec; ++s)
					std::fill(p.begin() + s * radNElem * nComp, p.begin() + (s+1) * radNElem * nComp, values[s]);

				values = std::move(p);

				for (unsigned int s = 0; s < nSec; ++s)
					parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, (mode == cadet::model::MultiplexMode::Independent) ? cadet::SectionIndep : s)] = &values[s * radNElem * nComp];
			}
			break;
		case cadet::model::MultiplexMode::Component: // nSec = 1
		case cadet::model::MultiplexMode::ComponentSection:
			{
				std::vector<cadet::active> p(nComp * radNElem * nSec);
				for (unsigned int s = 0; s < nSec; ++s)
				{
					for (unsigned int r = 0; r < radNElem; ++r)
						std::copy(values.begin() + s * nComp, values.begin() + (s+1) * nComp, p.begin() + r * nComp + s * nComp * radNElem);
				}

				values = std::move(p);

				for (unsigned int s = 0; s < nSec; ++s)
				{
					for (unsigned int i = 0; i < nComp; ++i)
						parameters[cadet::makeParamId(nameHash, uoi, i, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, (mode == cadet::model::MultiplexMode::Component) ? cadet::SectionIndep : s)] = &values[s * radNElem * nComp + i];
				}
			}
			break;
		case cadet::model::MultiplexMode::Radial: // nSec = 1
		case cadet::model::MultiplexMode::RadialSection:
			{
				std::vector<cadet::active> p(nComp * radNElem * nSec);
				for (unsigned int idx = 0; idx < radNElem * nSec; ++idx)
					std::fill(p.begin() + idx * nComp, p.begin() + (idx+1) * nComp, values[idx]);

				values = std::move(p);

				for (unsigned int s = 0; s < nSec; ++s)
				{
					for (unsigned int r = 0; r < radNElem; ++r)
						parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, r, cadet::BoundStateIndep, cadet::ReactionIndep, (mode == cadet::model::MultiplexMode::Radial) ? cadet::SectionIndep : s)] = &values[s * radNElem * nComp + r * nComp];
				}
			}
			break;
		case cadet::model::MultiplexMode::ComponentRadial: // nSec = 1
		case cadet::model::MultiplexMode::ComponentRadialSection:
			cadet::registerParam3DArray(parameters, values, [=](bool multi, unsigned int sec, unsigned int compartment, unsigned int comp) { return cadet::makeParamId(nameHash, uoi, comp, compartment, cadet::BoundStateIndep, cadet::ReactionIndep, multi ? sec : cadet::SectionIndep); }, nComp, radNElem);
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

bool multiplexParameterValue(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nComp, unsigned int radNElem, double value, std::unordered_set<cadet::active*> const* sensParams)
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

				if (sensParams && !cadet::contains(*sensParams, &data[pId.section * nComp * radNElem]))
					return false;

				for (unsigned int i = 0; i < nComp * radNElem; ++i)
					data[i + pId.section * nComp * radNElem].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::Component:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component]))
					return false;

				for (unsigned int i = 0; i < radNElem; ++i)
					data[i * nComp + pId.component].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component + pId.section * nComp * radNElem]))
					return false;

				for (unsigned int i = 0; i < radNElem; ++i)
					data[i * nComp + pId.component + pId.section * nComp * radNElem].setValue(value);

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

				if (sensParams && !cadet::contains(*sensParams, &data[pId.particleType * nComp + pId.section * nComp * radNElem]))
					return false;

				for (unsigned int i = 0; i < nComp; ++i)
					data[i + pId.particleType * nComp + pId.section * nComp * radNElem].setValue(value);

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

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component + pId.particleType * nComp + pId.section * nComp * radNElem]))
					return false;

				data[pId.component + pId.particleType * nComp + pId.section * nComp * radNElem].setValue(value);

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

bool multiplexParameterAD(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nComp, unsigned int radNElem, unsigned int adDirection, double adValue, std::unordered_set<cadet::active*>& sensParams)
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

				sensParams.insert(&data[pId.section * nComp * radNElem]);

				for (unsigned int i = 0; i < nComp * radNElem; ++i)
					data[i + pId.section * nComp * radNElem].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::Component:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.component]);

				for (unsigned int i = 0; i < radNElem; ++i)
					data[i * nComp + pId.component].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.component + pId.section * nComp * radNElem]);

				for (unsigned int i = 0; i < radNElem; ++i)
					data[i * nComp + pId.component + pId.section * nComp * radNElem].setADValue(adDirection, adValue);

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

				sensParams.insert(&data[pId.particleType * nComp + pId.section * nComp * radNElem]);

				for (unsigned int i = 0; i < nComp; ++i)
					data[i + pId.particleType * nComp + pId.section * nComp * radNElem].setADValue(adDirection, adValue);

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

				sensParams.insert(&data[pId.component + pId.particleType * nComp + pId.section * nComp * radNElem]);

				data[pId.component + pId.particleType * nComp + pId.section * nComp * radNElem].setADValue(adDirection, adValue);

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
TwoDimensionalConvectionDispersionOperatorDG::TwoDimensionalConvectionDispersionOperatorDG() : _colPorosities(0), _dir(0), _dispersionDep(nullptr),
_radLiftMCyl(nullptr), _transMrCyl(nullptr), _invTransMrCyl(nullptr), _SrCyl(nullptr), _jacConvection(nullptr), _jacAxDispersion(nullptr), _jacRadDispersion(nullptr)
{
}

TwoDimensionalConvectionDispersionOperatorDG::~TwoDimensionalConvectionDispersionOperatorDG() CADET_NOEXCEPT
{
	delete _dispersionDep;
	delete[] _radLiftMCyl;
	delete[] _transMrCyl;
	delete[] _invTransMrCyl;
	delete[] _SrCyl;
	delete[] _jacConvection;
	delete[] _jacAxDispersion;
	delete[] _jacRadDispersion;
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
	dgtoolbox::writeDGCoordinates(coords, _radNElem, _radNNodes, static_cast<const double*>(&_radNodes[0]), static_cast<double>(_colRadius), leftElemBndries);
	delete[] leftElemBndries;
}

void TwoDimensionalConvectionDispersionOperatorDG::initializeDG()
{
	const bool firstConfig = _jacConvection == nullptr; // used to avoid multiply allocation

	_axInvWeights = VectorXd::Zero(_axNNodes);
	_axNodes = VectorXd::Zero(_axNNodes);
	_radInvWeights = VectorXd::Zero(_radNNodes);
	_radNodes = VectorXd::Zero(_radNNodes);
	// LGL nodes and weights for axial and radial reference element and radial quadrature
	dgtoolbox::lglNodesWeights(_axPolyDeg, _axNodes, _axInvWeights, true);
	dgtoolbox::lglNodesWeights(_radPolyDeg, _radNodes, _radInvWeights, true);

	// auxiliary variables
	_axAuxStateG.resize(_axNPoints * _radNPoints);
	_radAuxStateG.resize(_axNPoints * _radNPoints);

	// numerical fluxes
	_fStarAux1.resize(2 * _radNNodes);
	_fStarAux2.resize(_axNNodes * 2);
	_fStarConv.resize(2 * _radNNodes);
	_gZStarDisp.resize(2 * _radNNodes);
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

	const int uAxElem = std::min(5, static_cast<int>(_axNElem)); // number of unique axial Jacobian blocks (per radial element)

	if (firstConfig)
	{
		_transMrCyl = new MatrixXd[_radNElem];
		_invTransMrCyl = new MatrixXd[_radNElem];
		_SrCyl = new MatrixXd[_radNElem];
		_radLiftMCyl = new MatrixXd[_radNElem];
		// Jacobian blocks
		_jacConvection = new MatrixXd[_radNElem];
		_jacAxDispersion = new MatrixXd[_radNElem * uAxElem];
		_jacRadDispersion = new MatrixXd[_nComp * _radNElem];
	}

	for (int rElem = 0; rElem < _radNElem; rElem++)
	{
		_jacConvection[rElem] = MatrixXd::Zero(_elemNPoints, 5 * _elemNPoints);
		for (int i = 0; i < uAxElem; i++)
			_jacAxDispersion[rElem * uAxElem + i] = MatrixXd::Zero(_elemNPoints, 5 * _elemNPoints);

		for (int comp = 0; comp < _nComp; comp++)
			_jacRadDispersion[rElem * _nComp + comp] = MatrixXd::Zero(_elemNPoints, 5 * _elemNPoints);
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
	_bulkNPoints = _axNPoints * _radNPoints;

	paramProvider.popScope();

	_axNodeStride = radNodeStride * _radNPoints;
	_axElemStride = _axNNodes * _axNodeStride;
	_radNodeStride = radNodeStride;
	_radElemStride = _radNNodes * _radNodeStride;

	_radialElemInterfaces.resize(_radNElem + 1);
	_radDelta.resize(_radNElem);
	_elementCrossSections.resize(_radNElem);
	_nodalCrossSections.resize(_radNPoints);
	_curVelocity.resize(_radNElem);

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

	if ((_colPorosities.size() != 1) && (_colPorosities.size() != _radNElem))
		throw InvalidParameterException("Number of elements in field COL_POROSITY is neither 1 nor radNElem (" + std::to_string(_radNElem) + ")");

	_singlePorosity = (_colPorosities.size() == 1);
	if (_singlePorosity)
		_colPorosities = std::vector<active>(_radNElem, _colPorosities[0]);

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

			if (!_singleVelocity && (_velocity.size() % _radNElem != 0))
				throw InvalidParameterException("Number of elements in field VELOCITY is not a positive multiple of radNElem (" + std::to_string(_radNElem) + ")");
			if ((mode == 0) && (_velocity.size() != 1))
				throw InvalidParameterException("Number of elements in field VELOCITY inconsistent with VELOCITY_MULTIPLEX (should be 1)");
			if ((mode == 1) && (_velocity.size() != _radNElem))
				throw InvalidParameterException("Number of elements in field VELOCITY inconsistent with VELOCITY_MULTIPLEX (should be " + std::to_string(_radNElem) + ")");
		}
		else
		{
			// Infer radial dependence of VELOCITY:
			//   size not divisible by _radNElem -> radial independent
			_singleVelocity = ((_velocity.size() % _radNElem) != 0);
		}

		// Expand _velocity to make it component dependent
		if (_singleVelocity)
		{
			std::vector<active> expanded(_velocity.size() * _radNElem);
			for (std::size_t i = 0; i < _velocity.size(); ++i)
				std::fill(expanded.begin() + i * _radNElem, expanded.begin() + (i + 1) * _radNElem, _velocity[i]);

			_velocity = std::move(expanded);
		}
	}
	else
	{
		_singleVelocity = false;
		_velocity.resize(_radNElem, 1.0);
	}

	// Register VELOCITY
	if (_singleVelocity)
	{
		if (_velocity.size() > _radNElem)
		{
			// Register only the first item in each section
			for (std::size_t i = 0; i < _velocity.size() / _radNElem; ++i)
				parameters[makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, i)] = &_velocity[i * _radNElem];
		}
		else
		{
			// We have only one parameter
			parameters[makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_velocity[0];
		}
	}
	else
		registerParam2DArray(parameters, _velocity, [=](bool multi, unsigned int sec, unsigned int compartment) { return makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, compartment, BoundStateIndep, ReactionIndep, multi ? sec : SectionIndep); }, _radNElem);

	_dir = std::vector<int>(_radNElem, 1);

	_axialDispersionMode = readAndRegisterMultiplexParam(paramProvider, parameters, _axialDispersion, "COL_DISPERSION", _nComp, _radNElem, unitOpIdx);
	_radialDispersionMode = readAndRegisterMultiplexParam(paramProvider, parameters, _radialDispersion, "COL_DISPERSION_RADIAL", _nComp, _radNElem, unitOpIdx);

	// Add parameters to map
	parameters[makeParamId(hashString("COL_LENGTH"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colLength;
	parameters[makeParamId(hashString("COL_RADIUS"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colRadius;
	registerParam1DArray(parameters, _colPorosities, [=](bool multi, unsigned int i) { return makeParamId(hashString("COL_POROSITY"), unitOpIdx, CompIndep, multi ? i : ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep); });

	// compute standard DG operators
	initializeDG();

	// compute nodal cross section areas and radial element interfaces
	updateRadialDisc();

	// compute DG operators that depend on radial geometry and dispersion after updateRadialDisc
	const int comp = 0;
	for (unsigned int rElem = 0; rElem < _radNElem; rElem++)
	{
		// todo ? use active types for column radius sensitivity
		_radLiftMCyl[rElem] = MatrixXd::Zero(2, _radNNodes);
		_radLiftMCyl[rElem](0, 0) = -static_cast<double>(_radialElemInterfaces[rElem]);
		_radLiftMCyl[rElem](1, _radNNodes - 1) = static_cast<double>(_radialElemInterfaces[rElem + 1]);

		_transMrCyl[rElem] = (static_cast<double>(_radialElemInterfaces[rElem]) * dgtoolbox::mMatrix(_radPolyDeg, _radNodes, 0.0, 0.0) + static_cast<double>(_radDelta[rElem]) / 2.0 * dgtoolbox::mMatrix(_radPolyDeg, _radNodes, 0.0, 1.0)).transpose();
		_invTransMrCyl[rElem] = _transMrCyl[rElem].inverse();

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

	if (!_velocity.empty())
	{
		// _curVelocity has already been set to the network flow rate in setFlowRates()
		// the direction of the flow (i.e., sign of _curVelocity) is given by _velocity
		active const* const dirNew = getSectionDependentSlice(_velocity, _radNElem, secIdx);

		for (unsigned int i = 0; i < _radNElem; ++i)
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
	_curVelocity[compartment] = _dir[compartment] * in / (_elementCrossSections[compartment] * _colPorosities[compartment]);
}

void TwoDimensionalConvectionDispersionOperatorDG::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	for (unsigned int radElem = 0; radElem < _radNElem; ++radElem)
	{
		for (unsigned int node = 0; node < _radNNodes; ++node)
		{
			const unsigned int compartment = radElem * _radNNodes + node;

			active v = in[compartment] / (_nodalCrossSections[compartment] * _colPorosities[radElem]);

			if (node == 0)
				_curVelocity[radElem] = v;

			if (abs(v - _curVelocity[radElem]) > 1E-12) // flow rates are specified for every radial node, but velocities must be constant per radial zone, which the user has to guarantee
				throw cadet::InvalidParameterException("Inconsistent definition of flow rates: Must define constant velocity within a radial DG element. Please refer to the documentation for further information.");

		}
	}
}

double TwoDimensionalConvectionDispersionOperatorDG::inletFactor(unsigned int idxSec, int idxRad) const CADET_NOEXCEPT // todo is this function needed?
{
	const double h = static_cast<double>(_colLength) / static_cast<double>(_axNPoints);
	return -std::abs(static_cast<double>(_curVelocity[idxRad])) / h;
}

const active& TwoDimensionalConvectionDispersionOperatorDG::axialDispersion(unsigned int idxSec, int idxRad, int idxComp) const CADET_NOEXCEPT
{
	return *(getSectionDependentSlice(_axialDispersion, _radNElem * _nComp, idxSec) + idxRad * _nComp + idxComp);
}

const active& TwoDimensionalConvectionDispersionOperatorDG::radialDispersion(unsigned int idxSec, int idxRad, int idxComp) const CADET_NOEXCEPT
{
	return *(getSectionDependentSlice(_radialDispersion, _radNElem * _nComp, idxSec) + idxRad * _nComp + idxComp);
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

	const active* const curAxialDispersion = getSectionDependentSlice(_axialDispersion, _radNElem * _nComp, secIdx);
	const active* const curRadialDispersion = getSectionDependentSlice(_radialDispersion, _radNElem * _nComp, secIdx);

	for (unsigned int comp = 0; comp < _nComp; comp++)
	{
		/*	auxiliary equations	*/

		for (unsigned int zEidx = 0; zEidx < _axNElem; zEidx++)
		{
			for (unsigned int rEidx = 0; rEidx < _radNElem; rEidx++)
			{
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

				// auxiliary equation g^z
				_Gz = _axInvMM.template cast<StateType>() * (_axLiftM.template cast<StateType>() * _fAux1 - _axTransStiffM.template cast<StateType>() * _C);

				// auxiliary equation g^r
				_Gr = (_fAux2 * _radLiftM.template cast<StateType>() - _C * _radStiffM.template cast<StateType>()) * _radInvTransMM.template cast<StateType>();

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

		for (unsigned int zEidx = 0; zEidx < _axNElem; zEidx++)
		{
			for (unsigned int rEidx = 0; rEidx < _radNElem; rEidx++)
			{
				const int elemOffset = zEidx * _axElemStride + rEidx * _radElemStride;
				const int auxElemOffset = zEidx * auxAxElemStride + rEidx * _radNNodes;

				ConstMatrixMap<StateType> _C(y + offsetC + comp + elemOffset, _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(_radNodeStride, _axNodeStride));
				MatrixMap<ResidualType> _Res(res + offsetC + comp + elemOffset, _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(_radNodeStride, _axNodeStride));

				ConstMatrixMap<StateType> _Gz(reinterpret_cast<StateType*>(&_axAuxStateG[0]) + auxElemOffset, _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(1, _radNNodes * _radNElem));
				
				ConstMatrixMap<StateType> _Gr(reinterpret_cast<StateType*>(&_radAuxStateG[0]) + auxElemOffset, _axNNodes, _radNNodes, Stride<Dynamic, Dynamic>(1, _radNNodes * _radNElem));

				MatrixMap<ResidualType> _fStarConvZ(reinterpret_cast<ResidualType*>(&_fStarConv[0]), 2, _radNNodes, Stride<Dynamic, Dynamic>(1, _radNNodes));
				MatrixMap<StateType> _gStarDispZ(reinterpret_cast<StateType*>(&_gZStarDisp[0]), 2, _radNNodes, Stride<Dynamic, Dynamic>(1, _radNNodes));
				MatrixMap<ResidualType> _gStarDispR(reinterpret_cast<ResidualType*>(&_gRStarDisp[0]), _axNNodes, 2, Stride<Dynamic, Dynamic>(1, 2));
				
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
						_gStarDispR.col(1) = 0.5 * (_Gr.col(_radPolyDeg) + GRnextNode) * 0.5 * (static_cast<ParamType>(curRadialDispersion[rEidx * _nComp + comp]) + static_cast<ParamType>(curRadialDispersion[(rEidx + 1) * _nComp + comp])) * 0.5 * (static_cast<ParamType>(_colPorosities[rEidx]) + static_cast<ParamType>(_colPorosities[rEidx + 1]));
					}
				}

				// axial dispersion (central) flux and axial convection (upwind) flux
				if (zEidx == 0) // Danckwerts inlet BC
				{
					_gStarDispZ.row(0).setZero();
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> inlet(y + rEidx * _radElemStride + comp, _radNNodes, InnerStride<Dynamic>(_radNodeStride));

					_fStarConvZ.row(0) = -static_cast<ParamType>(_curVelocity[rEidx]) * inlet.template cast<StateType>();
				}
				else
				{
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> GZPrevNode(reinterpret_cast<StateType*>(&_axAuxStateG[0]) + auxElemOffset - auxAxNodeStride, _radNNodes, InnerStride<Dynamic>(1));
					_gStarDispZ.row(0) = 0.5 * (_Gz.row(0) + GZPrevNode.transpose());

					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> prevNodeZ(y + offsetC + comp + elemOffset - _axNodeStride, _radNNodes, InnerStride<Dynamic>(_radNodeStride));
					_fStarConvZ.row(0) = -static_cast<ParamType>(_curVelocity[rEidx]) * prevNodeZ;
				}

				_fStarConvZ.row(1) = -static_cast<ParamType>(_curVelocity[rEidx]) * _C.block(_axPolyDeg, 0, 1, _radNNodes);
				if (zEidx != _axNElem - 1)
				{
					Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> GZNextNode(reinterpret_cast<StateType*>(&_axAuxStateG[0]) + auxElemOffset + auxAxElemStride, _radNNodes, InnerStride<Dynamic>(1));
					_gStarDispZ.row(1) = 0.5 * (_Gz.row(_axPolyDeg) + GZNextNode.transpose());
				}
				else // Danckwert outflow boundary condition
					_gStarDispZ.row(1).setZero();

				// Axial convection
				_Res -= 2.0 / static_cast<ParamType>(_axDelta) * _axInvMM.template cast<ResidualType>() * (
					// axial surface integral
					_axLiftM.template cast<ResidualType>() * _fStarConvZ
					// axial volume integral
					- _axTransStiffM.template cast<ResidualType>() * (-static_cast<ParamType>(_curVelocity[rEidx]) * _C.template cast<ResidualType>())
					);

				// Axial dispersion
				_Res -= 2.0 / static_cast<ParamType>(_axDelta) * 2.0 / static_cast<ParamType>(_axDelta) * static_cast<ParamType>(curAxialDispersion[rEidx * _nComp + comp]) * _axInvMM.template cast<ResidualType>() * (
						// surface integral
					    _axLiftM.template cast<ResidualType>() * (
						_gStarDispZ.template cast<ResidualType>()
						)
					    // volume integral
						-_axTransStiffM.template cast<ResidualType>() * (
							_Gz.template cast<ResidualType>()
							)
						);

				const ParamType porosityFluxFactor = 1.0 / static_cast<ParamType>(_colPorosities[rEidx]);

				// Radial dispersion
				_Res -= 2.0 / static_cast<ParamType>(_radDelta[rEidx]) * 2.0 / static_cast<ParamType>(_radDelta[rEidx]) * (
					_gStarDispR * porosityFluxFactor * _radLiftMCyl[rEidx].template cast<ResidualType>()
					- static_cast<ParamType>(curRadialDispersion[rEidx * _nComp + comp]) * _Gr.template cast<ResidualType>() * _SrCyl[rEidx].template cast<ResidualType>()
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
		const double u = static_cast<double>(_curVelocity[rElem]);

		_jacConvection[rElem] = MzKronMrCylInv * 2.0 / deltaZ * u * (SzTKronMrCyl * cDer - BzKronMrCyl * cStarDer);

		// filter numerical noise
		_jacConvection[rElem] = _jacConvection[rElem].unaryExpr([](double val) {
			return (std::abs(val) < 1e-14) ? 0.0 : val;
			});

		/* axial dispersion block */
		const active* const curAxialDispersion = getSectionDependentSlice(_axialDispersion, _radNElem * _nComp, 0);

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

			_jacAxDispersion[rElem * uAxElem + i] = MzKronMrCylInv * 2.0 / deltaZ * (BzKronMrCyl * gStarZ);
			_jacAxDispersion[rElem * uAxElem + i].block(0, Np, Np, 3 * Np) -= MzKronMrCylInv * 2.0 / deltaZ * (SzTKronMrCyl * GzDer[auxIdx]);

			// filter numerical noise
			_jacAxDispersion[rElem * uAxElem + i] = _jacAxDispersion[rElem * uAxElem + i].unaryExpr([](double val) {
				return (std::abs(val) < 1e-14) ? 0.0 : val;
				});

			gStarZ.setZero();
		}

		/* radial dispersion block */
		const active* const curRadialDispersion = getSectionDependentSlice(_radialDispersion, _radNElem * _nComp, 0);
		const double curPorosity = static_cast<double>(_colPorosities[rElem]);

		MatrixXd gStarR = MatrixXd::Zero(Np, 5 * Np);

		for (int comp = 0; comp < _nComp; comp++)
		{
			const double radDisp = static_cast<double>(curRadialDispersion[rElem * _nComp + comp]);

			for (int i = 0; i < _axNNodes; i++)
			{
				if (rElem > 0) // else zero as per boundary condition
				{
					const int ll = rElem == 1 ? 0 : 1;
					const int lr = rElem == _radNElem - 1 ? 2 : 1;

					const double dispFactor = 0.5 * (static_cast<double>(curRadialDispersion[(rElem - 1) * _nComp + comp]) + radDisp);
					const double porosityFactor = 0.5 * (static_cast<double>(_colPorosities[rElem - 1]) + curPorosity) / curPorosity;

					gStarR.block(i * _radNNodes, i * 5 * _radNNodes, 1, 3 * _radNNodes) = porosityFactor * dispFactor * 0.5 * GrDer[ll].block((i + 1) * _radNNodes - 1, i * 3 * _radNNodes, 1, 3 * _radNNodes);
					gStarR.block(i * _radNNodes, i * 5 * _radNNodes + _radNNodes, 1, 3 * _radNNodes) += porosityFactor * dispFactor * 0.5 * GrDer[lr].block(i * _radNNodes, i * 3 * _radNNodes, 1, 3 * _radNNodes);
				}
				if (rElem < _radNElem - 1) // else zero as per boundary condition
				{
					const int rl = rElem == 0 ? 0 : 1;
					const int rr = rElem + 1 == _radNElem - 1 ? 2 : 1;

					const double dispFactor = 0.5 * (radDisp + static_cast<double>(curRadialDispersion[(rElem + 1) * _nComp + comp]));
					const double porosityFactor = 0.5 * (curPorosity + static_cast<double>(_colPorosities[rElem + 1])) / curPorosity;

					gStarR.block(i * _radNNodes + _radPolyDeg, i * 5 * _radNNodes + _radNNodes, 1, 3 * _radNNodes) = porosityFactor * dispFactor * 0.5 * GrDer[rl].block((i + 1) * _radNNodes - 1, i * 3 * _radNNodes, 1, 3 * _radNNodes);
					gStarR.block(i * _radNNodes + _radPolyDeg, i * 5 * _radNNodes + 2 * _radNNodes, 1, 3 * _radNNodes) += porosityFactor * dispFactor * 0.5 * GrDer[rr].block(i * _radNNodes, i * 3 * _radNNodes, 1, 3 * _radNNodes);
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

			_jacRadDispersion[rElem * _nComp + comp] = 2.0 / static_cast<double>(_radDelta[rElem]) * 2.0 / static_cast<double>(_radDelta[rElem]) * MzKronMrCylInv * (MzKronBrCyl * gStarR);
			for (int zNode = 0; zNode < _axNNodes; zNode++) // todo can this be done more elegantly?
				_jacRadDispersion[rElem * _nComp + comp].block(0, (5 * zNode + 1) * _radNNodes, Np, 3 * _radNNodes) -= 2.0 / static_cast<double>(_radDelta[rElem]) * 2.0 / static_cast<double>(_radDelta[rElem]) * radDisp * MzKronMrCylInv * (MzKronSrTCyl * GrDer[auxIdx].block(0, zNode * 3 * _radNNodes, Np, 3 * _radNNodes));

			// filter numerical noise
			_jacRadDispersion[rElem * _nComp + comp] = _jacRadDispersion[rElem * _nComp + comp].unaryExpr([](double val) {
				return (std::abs(val) < 1e-14) ? 0.0 : val;
				});

			gStarR.setZero();
		}
	}

	delete[] GrDer;
	delete[] GzDer;

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
void TwoDimensionalConvectionDispersionOperatorDG::addAxElemBlockToJac(const Eigen::MatrixXd& block, const int offRow, const int offColumn, const int depElem, Action addEntry, const active* const compFac)
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
							double entry = block(zNode * _radNNodes + rNode, depBlock * _elemNPoints + i * _radNNodes + j);
							if (compFac)
								entry *= static_cast<double>(compFac[comp]);

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
void TwoDimensionalConvectionDispersionOperatorDG::addAxElemBlockToJac(const Eigen::MatrixXd& block, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, const int offRow, const int offColumn, const int depElem, const active* const compFac)
{
	auto addEntry = [&jacobian](int row, int col, double value) {
		jacobian.coeffRef(row, col) += value;
		};

	addAxElemBlockToJac(block, offRow, offColumn, depElem, addEntry, compFac);
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
void TwoDimensionalConvectionDispersionOperatorDG::addRadElemBlockToJac(const Eigen::MatrixXd* block, const int offRow, const int nLeftRadElemDep, const int depElem, Action addEntry)
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
							const double entry = -block[comp](zNode * _radNNodes + rNode, offBlock + i * 5 * _radNNodes + depBlock * _radNNodes + j);

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
 * @param [in] jacobian jacobian matrix times -1.0
 * @param [in] offRow offset to first considered row
 * @param [in] nLeftRadElemDep number of left radial elements this block depends on
 * @param [in] depElem number of elements the element block element depends on
 */
void TwoDimensionalConvectionDispersionOperatorDG::addRadElemBlockToJac(const Eigen::MatrixXd* block, Eigen::SparseMatrix<double, Eigen::RowMajor>& jacobian, const int offRow, const int nLeftRadElemDep, const int depElem)
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
void TwoDimensionalConvectionDispersionOperatorDG::addRadElemBlockToJac(const Eigen::MatrixXd* block, std::vector<T>& tripletList, const int offRow, const int nLeftRadElemDep, const int depElem)
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
			/* Axial convection Pattern */

			const int nRightAxElemDep = std::min(2, static_cast<int>(_axNElem) - 1 - zElem); //<! number of right axial elements, this element depends on
			const int nLeftAxElemDep = std::min(2, zElem); //<! number of left axial elements, this element depends on
			const int offSetRow = bulkOffset + zElem * _axElemStride + rElem * _radElemStride;
			const int offSetColumnToRow = -nLeftAxElemDep * _axElemStride;

			addAxElemBlockToJac(-_jacConvection[rElem].block(0, Np * (2 - nLeftAxElemDep), Np, Np * (nLeftAxElemDep + 1 + nRightAxElemDep)), tripletList, offSetRow, offSetColumnToRow, nLeftAxElemDep + 1 + nRightAxElemDep);

			/* Axial dispersion Pattern */
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

				addAxElemBlockToJac(-_jacAxDispersion[rElem * uAxElem + uAxBlockIdx].block(0, Np * (2 - nLeftAxElemDep), Np, Np * (nLeftAxElemDep + 1 + nRightAxElemDep)), tripletList, offSetRow, offSetColumnToRow, nLeftAxElemDep + 1 + nRightAxElemDep);
			}

			/* Radial dispersion Pattern */
			const int nLeftRadElem = std::min(2, rElem);
			const int nRightRadElem = std::min(2, static_cast<int>(_radNElem) - 1 - rElem);
			addRadElemBlockToJac(_jacRadDispersion + rElem * _nComp, tripletList, offSetRow, nLeftRadElem, nLeftRadElem + 1 + nRightRadElem);
		}
	}

	/* Bulk reaction / binding pattern */
	// Also needs to be set no matter if we actually have reactions/bindings, so that we can assume the bulk pattern always be the same,
	// which we assume in the addAxElemBlockToJac and addRadElemBlockToJac functions.

	for (int bulk = 0; bulk < _bulkNPoints; bulk++)
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

	const active* const Dax = getSectionDependentSlice(_axialDispersion, _radNElem * _nComp, 0);
	const active* const Drad = getSectionDependentSlice(_radialDispersion, _radNElem * _nComp, 0);

	for (int zElem = 0; zElem < _axNElem; zElem++)
	{
		const active* curDax = Dax;
		const active* curDrad = Drad;

		// Note: Jacobian blocks *-1 for residual

		for (int rElem = 0; rElem < _radNElem; rElem++, curDax += _nComp, curDrad += _nComp)
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
								jacInlet(rElem * _radElemStride + rNode * _radNodeStride + zNode * _axNodeStride + comp, rElem * _radElemStride + j * _radNodeStride + comp) = -_jacConvection[rElem].block(0, Np + Np - _radNNodes, Np, _radNNodes)(rNode + zNode * _radNNodes, j);
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

				addAxElemBlockToJac(-_jacAxDispersion[rElem * uAxElem + uAxBlockIdx].block(0, Np * (2 - nLeftAxElemDep), Np, Np * (nLeftAxElemDep + 1 + nRightAxElemDep)), jacobian, offSetRow, offSetColumnToRow, nLeftAxElemDep + 1 + nRightAxElemDep, curDax);
			}

			/* handle radial dispersion Jacobian */
			const int nLeftRadElem = std::min(2, rElem);
			const int nRightRadElem = std::min(2, static_cast<int>(_radNElem) - 1 - rElem);
			addRadElemBlockToJac(_jacRadDispersion + rElem * _nComp, jacobian, offSetRow, nLeftRadElem, nLeftRadElem + 1 + nRightRadElem);
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

		_elementCrossSections[r] = pi * (pow(_radialElemInterfaces[r + 1], 2.0) - pow(_radialElemInterfaces[r], 2.0)); // equivalent to the sum of the nodal cross sections of that element
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

		_elementCrossSections[r] = pi * (pow(_radialElemInterfaces[r + 1], 2.0) - pow(_radialElemInterfaces[r], 2.0));
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

		_elementCrossSections[r] = pi * (pow(_radialElemInterfaces[r + 1], 2.0) - pow(_radialElemInterfaces[r], 2.0));
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
		for (unsigned int i = 0; i < _radNElem; ++i)
			_colPorosities[i].setValue(value);
		return true;
	}

	if (_singleVelocity && (pId.name == hashString("VELOCITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep) && (pId.reaction == ReactionIndep))
	{
		if (_velocity.size() > _radNElem)
		{
			// Section dependent
			if (pId.section == SectionIndep)
				return false;

			for (unsigned int i = 0; i < _radNElem; ++i)
				_velocity[pId.section * _radNElem + i].setValue(value);
		}
		else
		{
			// Section independent
			if (pId.section != SectionIndep)
				return false;

			for (unsigned int i = 0; i < _radNElem; ++i)
				_velocity[i].setValue(value);
		}
	}

	const bool ad = multiplexParameterValue(pId, hashString("COL_DISPERSION"), _axialDispersionMode, _axialDispersion, _nComp, _radNElem, value, nullptr);
	if (ad)
		return true;

	const bool adr = multiplexParameterValue(pId, hashString("COL_DISPERSION_RADIAL"), _radialDispersionMode, _radialDispersion, _nComp, _radNElem, value, nullptr);
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
			for (unsigned int i = 0; i < _radNElem; ++i)
				_colPorosities[i].setValue(value);
			return true;
		}
	}

	if (_singleVelocity && (pId.name == hashString("VELOCITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep) && (pId.reaction == ReactionIndep))
	{
		if (_velocity.size() > _radNElem)
		{
			// Section dependent
			if ((pId.section == SectionIndep) || !contains(sensParams, &_velocity[pId.section * _radNElem]))
				return false;

			for (unsigned int i = 0; i < _radNElem; ++i)
				_velocity[pId.section * _radNElem + i].setValue(value);
		}
		else
		{
			// Section independent
			if ((pId.section != SectionIndep) || !contains(sensParams, &_velocity[0]))
				return false;

			for (unsigned int i = 0; i < _radNElem; ++i)
				_velocity[i].setValue(value);
		}
	}

	const bool ad = multiplexParameterValue(pId, hashString("COL_DISPERSION"), _axialDispersionMode, _axialDispersion, _nComp, _radNElem, value, &sensParams);
	if (ad)
		return true;

	const bool adr = multiplexParameterValue(pId, hashString("COL_DISPERSION_RADIAL"), _radialDispersionMode, _radialDispersion, _nComp, _radNElem, value, &sensParams);
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
		for (unsigned int i = 0; i < _radNElem; ++i)
			_colPorosities[i].setADValue(adDirection, adValue);

		return true;
	}

	if (_singleVelocity && (pId.name == hashString("VELOCITY")) && (pId.component == CompIndep) && (pId.boundState == BoundStateIndep) && (pId.reaction == ReactionIndep))
	{
		if (_velocity.size() > _radNElem)
		{
			// Section dependent
			if (pId.section == SectionIndep)
				return false;

			sensParams.insert(&_velocity[pId.section * _radNElem]);
			for (unsigned int i = 0; i < _radNElem; ++i)
				_velocity[pId.section * _radNElem + i].setADValue(adDirection, adValue);
		}
		else
		{
			// Section independent
			if (pId.section != SectionIndep)
				return false;

			sensParams.insert(&_velocity[0]);
			for (unsigned int i = 0; i < _radNElem; ++i)
				_velocity[i].setADValue(adDirection, adValue);
		}
	}

	const bool ad = multiplexParameterAD(pId, hashString("COL_DISPERSION"), _axialDispersionMode, _axialDispersion, _nComp, _radNElem, adDirection, adValue, sensParams);
	if (ad)
		return true;

	const bool adr = multiplexParameterAD(pId, hashString("COL_DISPERSION_RADIAL"), _radialDispersionMode, _radialDispersion, _nComp, _radNElem, adDirection, adValue, sensParams);
	if (adr)
		return true;

	return false;
}

}  // namespace parts

}  // namespace model

}  // namespace cadet

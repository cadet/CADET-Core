// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//
//  Copyright © 2008-2019: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/GeneralRateModel2D.hpp"
#include "BindingModelFactory.hpp"
#include "ParamReaderHelper.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/ExternalFunction.hpp"
#include "cadet/SolutionRecorder.hpp"
#include "ConfigurationHelper.hpp"
#include "model/BindingModel.hpp"
#include "SimulationTypes.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "linalg/Norms.hpp"

#include "Stencil.hpp"
#include "Weno.hpp"
#include "AdUtils.hpp"
#include "SensParamUtil.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>
#include <iomanip>

#include "ParallelSupport.hpp"
#ifdef CADET_PARALLELIZE
	#include <tbb/tbb.h>
#endif


namespace
{

cadet::model::MultiplexMode readAndRegisterMultiplexParam(cadet::IParameterProvider& paramProvider, std::unordered_map<cadet::ParameterId, cadet::active*>& parameters, std::vector<cadet::active>& values, const std::string& name, unsigned int nAxial, unsigned int nRad, unsigned int nParType, cadet::UnitOpIdx uoi)
{
	cadet::model::MultiplexMode mode = cadet::model::MultiplexMode::Independent;
	readScalarParameterOrArray(values, paramProvider, name, 1);
	if (paramProvider.exists(name + "_MULTIPLEX"))
	{
		const int modeConfig = paramProvider.getInt(name + "_MULTIPLEX");
		if (modeConfig == 0)
		{
			mode = cadet::model::MultiplexMode::Independent;
			if (values.size() != nParType)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nParType) + ")");
		}
		else if (modeConfig == 1)
		{
			mode = cadet::model::MultiplexMode::Radial;
			if (values.size() != nRad * nParType)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nRad * nParType) + ")");
		}
		else if (modeConfig == 2)
		{
			mode = cadet::model::MultiplexMode::Axial;
			if (values.size() != nAxial * nParType)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nAxial * nParType) + ")");
		}
		else if (modeConfig == 3)
		{
			mode = cadet::model::MultiplexMode::AxialRadial;
			if (values.size() != nAxial * nRad * nParType)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nAxial * nRad * nParType) + ")");
		}
	}
	else
	{
		if (values.size() == nParType)
			mode = cadet::model::MultiplexMode::Independent;
		else if (values.size() == nRad * nParType)
			mode = cadet::model::MultiplexMode::Radial;
		else if (values.size() == nAxial * nParType)
			mode = cadet::model::MultiplexMode::Axial;
		else if (values.size() == nRad * nAxial * nParType)
			mode = cadet::model::MultiplexMode::AxialRadial;
		else
			throw cadet::InvalidParameterException("Could not infer multiplex mode of field " + name + ", set " + name + "_MULTIPLEX or change number of elements");
	}

	const cadet::StringHash nameHash = cadet::hashStringRuntime(name);
	switch (mode)
	{
		case cadet::model::MultiplexMode::Independent:
			{
				std::vector<cadet::active> p(nAxial * nRad * nParType);
				for (unsigned int s = 0; s < nAxial * nRad; ++s)
					std::copy(values.begin(), values.end(), p.begin() + s * nParType);

				values = std::move(p);

				for (unsigned int s = 0; s < nParType; ++s)
					parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, s, cadet::BoundStateIndep, cadet::ReactionIndep, cadet::SectionIndep)] = &values[s];
			}
			break;
		case cadet::model::MultiplexMode::Radial:
			{
				std::vector<cadet::active> p(nAxial * nRad * nParType);
				for (unsigned int s = 0; s < nAxial; ++s)
					std::copy(values.begin(), values.end(), p.begin() + s * nParType * nRad);

				values = std::move(p);

				for (unsigned int s = 0; s < nRad; ++s)
				{
					for (unsigned int i = 0; i < nParType; ++i)
						parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, i, cadet::BoundStateIndep, s, cadet::SectionIndep)] = &values[s * nParType + i];
				}
			}
			break;
		case cadet::model::MultiplexMode::Axial:
			{
				std::vector<cadet::active> p(nAxial * nRad * nParType);
				for (unsigned int i = 0; i < nAxial; ++i)
				{
					for (unsigned int j = 0; j < nRad; ++j)
						std::copy(values.begin() + i * nParType, values.begin() + (i+1) * nParType, p.begin() + i * nRad * nParType + j * nParType);
				}

				values = std::move(p);

				for (unsigned int s = 0; s < nAxial; ++s)
				{
					for (unsigned int i = 0; i < nParType; ++i)
						parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, i, cadet::BoundStateIndep, cadet::ReactionIndep, s)] = &values[s * nParType * nRad + i];
				}
			}
			break;
		case cadet::model::MultiplexMode::AxialRadial:
			cadet::registerParam3DArray(parameters, values, [=](bool multi, unsigned int ax, unsigned int rad, unsigned int pt) { return cadet::makeParamId(nameHash, uoi, cadet::CompIndep, pt, cadet::BoundStateIndep, rad, ax); }, nParType, nRad);
			break;
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Component:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::ComponentSection:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			cadet_assert(false);
			break;
	}

	return mode;
}

bool multiplexParameterValue(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nAxial, unsigned int nRad, unsigned int nParType, double value, std::unordered_set<cadet::active*> const* sensParams)
{
	if (pId.name != nameHash)
		return false;

	switch (mode)
	{
		case cadet::model::MultiplexMode::Independent:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.particleType]))
					return false;

				for (unsigned int i = 0; i < nAxial * nRad; ++i)
					data[i * nParType + pId.particleType].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::Radial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction == cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.reaction * nParType + pId.particleType]))
					return false;

				for (unsigned int i = 0; i < nAxial; ++i)
					data[i * nRad * nParType + pId.reaction * nParType + pId.particleType].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::Axial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.section * nParType * nRad + pId.particleType]))
					return false;

				for (unsigned int i = 0; i < nRad; ++i)
					data[pId.section * nParType * nRad + i * nParType + pId.particleType].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::AxialRadial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction == cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.section * nParType * nRad + pId.reaction * nParType + pId.particleType]))
					return false;

				data[pId.section * nParType * nRad + pId.reaction * nParType + pId.particleType].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Component:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::ComponentSection:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			cadet_assert(false);
			break;
	}

	return false;
}

bool multiplexParameterAD(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nAxial, unsigned int nRad, unsigned int nParType, unsigned int adDirection, double adValue, std::unordered_set<cadet::active*>& sensParams)
{
	if (pId.name != nameHash)
		return false;

	switch (mode)
	{
		case cadet::model::MultiplexMode::Independent:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.particleType]);

				for (unsigned int i = 0; i < nAxial * nRad; ++i)
					data[i * nParType + pId.particleType].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::Radial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction == cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.reaction * nParType + pId.particleType]);

				for (unsigned int i = 0; i < nAxial; ++i)
					data[i * nRad * nParType + pId.reaction * nParType + pId.particleType].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::Axial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.section * nParType * nRad + pId.particleType]);

				for (unsigned int i = 0; i < nRad; ++i)
					data[pId.section * nParType * nRad + i * nParType + pId.particleType].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::AxialRadial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction == cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.section * nParType * nRad + pId.reaction * nParType + pId.particleType]);
				data[pId.section * nParType * nRad + pId.reaction * nParType + pId.particleType].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Component:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::ComponentSection:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			cadet_assert(false);
			break;
	}

	return false;
}



cadet::model::MultiplexMode readAndRegisterMultiplexParam(cadet::IParameterProvider& paramProvider, std::unordered_map<cadet::ParameterId, cadet::active*>& parameters, std::vector<cadet::active>& values, const std::string& name, unsigned int nParType, unsigned int nComp, cadet::UnitOpIdx uoi)
{
	cadet::model::MultiplexMode mode = cadet::model::MultiplexMode::Independent;
	readScalarParameterOrArray(values, paramProvider, name, 1);
	if (paramProvider.exists(name + "_MULTIPLEX"))
	{
		const int modeConfig = paramProvider.getInt(name + "_MULTIPLEX");
		if (modeConfig == 0)
		{
			mode = cadet::model::MultiplexMode::Component;
			if (values.size() != nComp)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nComp) + ")");
		}
		else if (modeConfig == 1)
		{
			mode = cadet::model::MultiplexMode::ComponentSection;
			if ((values.size() % nComp) != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be positive multiple of " + std::to_string(nComp) + ")");
		}
		else if (modeConfig == 2)
		{
			mode = cadet::model::MultiplexMode::ComponentType;
			if (values.size() != nComp * nParType)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nComp * nParType) + ")");
		}
		else if (modeConfig == 3)
		{
			mode = cadet::model::MultiplexMode::ComponentSectionType;
			if ((values.size() % (nComp * nParType)) != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be positive multiple of " + std::to_string(nComp * nParType) + ")");
		}
	}
	else
	{
		if (values.size() == nComp)
			mode = cadet::model::MultiplexMode::Component;
		else if (values.size() == nComp * nParType)
			mode = cadet::model::MultiplexMode::ComponentType;
		else if ((values.size() % (nComp * nParType)) == 0)
			mode = cadet::model::MultiplexMode::ComponentSectionType;
		else if ((values.size() % nComp) == 0)
			mode = cadet::model::MultiplexMode::ComponentSection;
		else
			throw cadet::InvalidParameterException("Could not infer multiplex mode of field " + name + ", set " + name + "_MULTIPLEX or change number of elements");
	}

	const cadet::StringHash nameHash = cadet::hashStringRuntime(name);
	switch (mode)
	{
		case cadet::model::MultiplexMode::Component:
			{
				std::vector<cadet::active> p(nComp * nParType);
				for (unsigned int s = 0; s < nParType; ++s)
					std::copy(values.begin(), values.end(), p.begin() + s * nComp);

				values = std::move(p);

				for (unsigned int s = 0; s < nComp; ++s)
					parameters[cadet::makeParamId(nameHash, uoi, s, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, cadet::SectionIndep)] = &values[s];
			}
			break;
		case cadet::model::MultiplexMode::ComponentSection:
			{
				std::vector<cadet::active> p(values.size() * nParType);
				for (unsigned int s = 0; s < values.size() / nComp; ++s)
				{
					for (unsigned int pt = 0; pt < nParType; ++pt)
					{
						std::copy(values.begin() + s * nComp, values.begin() + (s+1) * nComp, p.begin() + s * nParType * nComp + pt * nComp);
					}
				}

				values = std::move(p);

				for (unsigned int s = 0; s < values.size() / nComp; ++s)
				{
					for (unsigned int i = 0; i < nComp; ++i)
						parameters[cadet::makeParamId(nameHash, uoi, i, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, s)] = &values[s * nParType * nComp + i];
				}
			}
			break;
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			cadet::registerParam3DArray(parameters, values, [=](bool multi, unsigned int sec, unsigned int pt, unsigned int comp) { return cadet::makeParamId(nameHash, uoi, comp, pt, cadet::BoundStateIndep, cadet::ReactionIndep, multi ? sec : cadet::SectionIndep); }, nComp, nParType);
			break;
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Independent:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::Axial:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::Radial:
		case cadet::model::MultiplexMode::AxialRadial:
			cadet_assert(false);
			break;
	}

	return mode;
}

bool multiplexParameterValue(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nParType, unsigned int nComp, double value, std::unordered_set<cadet::active*> const* sensParams)
{
	if (pId.name != nameHash)
		return false;

	switch (mode)
	{
		case cadet::model::MultiplexMode::Component:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.component]))
					return false;

				for (unsigned int i = 0; i < nParType; ++i)
					data[i * nComp + pId.component].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.section * nComp * nParType + pId.component]))
					return false;

				for (unsigned int i = 0; i < nParType; ++i)
					data[i * nComp + pId.section * nComp * nParType + pId.component].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Independent:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::Axial:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::Radial:
		case cadet::model::MultiplexMode::AxialRadial:
			cadet_assert(false);
			break;
	}

	return false;
}

bool multiplexParameterAD(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nParType, unsigned int nComp, unsigned int adDirection, double adValue, std::unordered_set<cadet::active*>& sensParams)
{
	if (pId.name != nameHash)
		return false;

	switch (mode)
	{
		case cadet::model::MultiplexMode::Component:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.component]);

				for (unsigned int i = 0; i < nParType; ++i)
					data[i * nComp + pId.component].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.section * nComp * nParType + pId.component]);

				for (unsigned int i = 0; i < nParType; ++i)
					data[i * nComp + pId.section * nComp * nParType + pId.component].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Independent:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::Axial:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::Radial:
		case cadet::model::MultiplexMode::AxialRadial:
			cadet_assert(false);
			break;
	}

	return false;
}



cadet::model::MultiplexMode readAndRegisterMultiplexParam(cadet::IParameterProvider& paramProvider, std::unordered_map<cadet::ParameterId, cadet::active*>& parameters, std::vector<cadet::active>& values, const std::string& name, unsigned int nParType, unsigned int nComp, unsigned int const* strideBound, unsigned int const* nBound, cadet::UnitOpIdx uoi)
{
	const unsigned int nTotalBound = strideBound[nParType];
	if (nTotalBound == 0)
		return cadet::model::MultiplexMode::Component;

	cadet::model::MultiplexMode mode = cadet::model::MultiplexMode::Independent;
	readScalarParameterOrArray(values, paramProvider, name, 1);
	if (paramProvider.exists(name + "_MULTIPLEX"))
	{
		const int modeConfig = paramProvider.getInt(name + "_MULTIPLEX");
		if (modeConfig == 0)
		{
			mode = cadet::model::MultiplexMode::Component;
			if (values.size() != strideBound[0])
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(strideBound[0]) + ")");
		}
		else if (modeConfig == 1)
		{
			mode = cadet::model::MultiplexMode::ComponentSection;
			if ((values.size() % strideBound[0]) != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be positive multiple of " + std::to_string(strideBound[0]) + ")");
		}
		else if (modeConfig == 2)
		{
			mode = cadet::model::MultiplexMode::ComponentType;
			if (values.size() != nTotalBound)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nTotalBound) + ")");
		}
		else if (modeConfig == 3)
		{
			mode = cadet::model::MultiplexMode::ComponentSectionType;
			if ((values.size() % nTotalBound) != 0)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be positive multiple of " + std::to_string(nTotalBound) + ")");
		}
	}
	else
	{
		if (values.size() == strideBound[0])
			mode = cadet::model::MultiplexMode::Component;
		else if (values.size() == nTotalBound)
			mode = cadet::model::MultiplexMode::ComponentType;
		else if ((values.size() % nTotalBound) == 0)
			mode = cadet::model::MultiplexMode::ComponentSectionType;
		else if ((values.size() % strideBound[0]) == 0)
			mode = cadet::model::MultiplexMode::ComponentSection;
		else
			throw cadet::InvalidParameterException("Could not infer multiplex mode of field " + name + ", set " + name + "_MULTIPLEX or change number of elements");
	}

	const cadet::StringHash nameHash = cadet::hashStringRuntime(name);
	switch (mode)
	{
		case cadet::model::MultiplexMode::Component:
			{
				std::vector<cadet::active> p(nTotalBound);
				for (unsigned int s = 0; s < nParType; ++s)
					std::copy(values.begin(), values.end(), p.begin() + s * strideBound[0]);

				values = std::move(p);

				unsigned int idx = 0;
				for (unsigned int c = 0; c < nComp; ++c)
				{
					for (unsigned int s = 0; s < nBound[c]; ++s, ++idx)
						parameters[cadet::makeParamId(nameHash, uoi, c, cadet::ParTypeIndep, s, cadet::ReactionIndep, cadet::SectionIndep)] = &values[idx];
				}
			}
			break;
		case cadet::model::MultiplexMode::ComponentSection:
			{
				std::vector<cadet::active> p(values.size() / strideBound[0] * nTotalBound);
				for (unsigned int s = 0; s < values.size() / strideBound[0]; ++s)
				{
					for (unsigned int pt = 0; pt < nParType; ++pt)
					{
						std::copy(values.begin() + s * strideBound[0], values.begin() + (s+1) * strideBound[0], p.begin() + s * nTotalBound + pt * strideBound[0]);
					}
				}

				values = std::move(p);

				for (unsigned int sec = 0; sec < values.size() / strideBound[0]; ++sec)
				{
					unsigned int idx = nTotalBound * sec;
					for (unsigned int c = 0; c < nComp; ++c)
					{
						for (unsigned int s = 0; s < nBound[c]; ++s, ++idx)
							parameters[cadet::makeParamId(nameHash, uoi, c, cadet::ParTypeIndep, s, cadet::ReactionIndep, sec)] = &values[idx];
					}
				}
			}
			break;
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			{
				unsigned int idx = 0;
				const bool multiSec = (values.size() != nTotalBound);
				for (unsigned int sec = 0; sec < values.size() / nTotalBound; ++sec)
				{
					for (unsigned int type = 0; type < nParType; ++type)
					{
						for (unsigned int c = 0; c < nComp; ++c)
						{
							for (unsigned int s = 0; s < nBound[type * nComp + c]; ++s, ++idx)
								parameters[cadet::makeParamId(nameHash, uoi, c, type, s, cadet::ReactionIndep, multiSec ? sec : cadet::SectionIndep)] = &values[idx];
						}
					}
				}
			}
			break;
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Independent:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::Axial:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::Radial:
		case cadet::model::MultiplexMode::AxialRadial:
			cadet_assert(false);
			break;
	}

	return mode;
}

bool multiplexParameterValue(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data,
	unsigned int nParType, unsigned int nComp, unsigned int const* strideBound, unsigned int const* nBound, unsigned int const* boundOffset, double value, std::unordered_set<cadet::active*> const* sensParams)
{
	if (pId.name != nameHash)
		return false;

	if (strideBound[nParType] == 0)
		return true;

	switch (mode)
	{
		case cadet::model::MultiplexMode::Component:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState == cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[boundOffset[pId.component] + pId.boundState]))
					return false;

				for (unsigned int i = 0; i < nParType; ++i)
					data[boundOffset[pId.component] + pId.boundState + i * strideBound[0]].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState == cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.section * strideBound[nParType] + boundOffset[pId.component] + pId.boundState]))
					return false;

				for (unsigned int i = 0; i < nParType; ++i)
					data[pId.section * strideBound[nParType] + boundOffset[pId.component] + pId.boundState + i * strideBound[0]].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Independent:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::Axial:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::Radial:
		case cadet::model::MultiplexMode::AxialRadial:
			cadet_assert(false);
			break;
	}

	return false;
}

bool multiplexParameterAD(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data,
	unsigned int nParType, unsigned int nComp, unsigned int const* strideBound, unsigned int const* nBound, unsigned int const* boundOffset, unsigned int adDirection, double adValue, std::unordered_set<cadet::active*>& sensParams)
{
	if (pId.name != nameHash)
		return false;

	if (strideBound[nParType] == 0)
		return true;

	switch (mode)
	{
		case cadet::model::MultiplexMode::Component:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState == cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[boundOffset[pId.component] + pId.boundState]);

				for (unsigned int i = 0; i < nParType; ++i)
					data[boundOffset[pId.component] + pId.boundState + i * strideBound[0]].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentSection:
			{
				if ((pId.component == cadet::CompIndep) || (pId.particleType != cadet::ParTypeIndep) || (pId.boundState == cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.section * strideBound[nParType] + boundOffset[pId.component] + pId.boundState]);

				for (unsigned int i = 0; i < nParType; ++i)
					data[pId.section * strideBound[nParType] + boundOffset[pId.component] + pId.boundState + i * strideBound[0]].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Independent:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::Axial:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::Radial:
		case cadet::model::MultiplexMode::AxialRadial:
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

int schurComplementMultiplierGRM2D(void* userData, double const* x, double* z)
{
	GeneralRateModel2D* const grm = static_cast<GeneralRateModel2D*>(userData);
	return grm->schurComplementMatrixVector(x, z);
}


GeneralRateModel2D::GeneralRateModel2D(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
	_jacP(nullptr), _jacPdisc(nullptr), _jacPF(nullptr), _jacFP(nullptr), _jacInlet(),
	_analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr),
	_initC(0), _singleRadiusInitC(true), _initCp(0), _singleRadiusInitCp(true), _initQ(0), _singleRadiusInitQ(true), _initState(0), _initStateDot(0)
{
}

GeneralRateModel2D::~GeneralRateModel2D() CADET_NOEXCEPT
{
	delete[] _tempState;

	delete[] _jacPF;
	delete[] _jacFP;

	delete[] _jacP;
	delete[] _jacPdisc;

	delete[] _bindingWorkspaceOffset;

	delete[] _disc.nParCell;
	delete[] _disc.parTypeOffset;
	delete[] _disc.nParCellsBeforeType;
	delete[] _disc.nBound;
	delete[] _disc.boundOffset;
	delete[] _disc.strideBound;
	delete[] _disc.nBoundBeforeType;
}

unsigned int GeneralRateModel2D::numDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nCol * nComp * nRad
	// Particle DOFs: nCol * nRad * nParType particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	//                in each shell; there are nParCell shells for each particle type
	// Flux DOFs: nCol * nComp * nRad * nParType (as many as column bulk DOFs)
	// Inlet DOFs: nComp * nRad
	return _disc.nCol * _disc.nRad * (_disc.nComp * (1 + _disc.nParType)) + _disc.parTypeOffset[_disc.nParType] + _disc.nComp * _disc.nRad;
}

unsigned int GeneralRateModel2D::numPureDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nCol * nComp * nRad
	// Particle DOFs: nCol * nRad * nParType particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	//                in each shell; there are nParCell shells for each particle type
	// Flux DOFs: nCol * nComp * nRad * nParType (as many as column bulk DOFs)
	return _disc.nCol * _disc.nRad * (_disc.nComp * (1 + _disc.nParType)) + _disc.parTypeOffset[_disc.nParType];
}


bool GeneralRateModel2D::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

bool GeneralRateModel2D::configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper)
{
	// ==== Read discretization
	_disc.nComp = paramProvider.getInt("NCOMP");

	paramProvider.pushScope("discretization");

	_disc.nCol = paramProvider.getInt("NCOL");
	_disc.nRad = paramProvider.getInt("NRAD");

	const std::vector<int> nParCell = paramProvider.getIntArray("NPAR");

	const std::vector<int> nBound = paramProvider.getIntArray("NBOUND");
	if (nBound.size() < _disc.nComp)
		throw InvalidParameterException("Field NBOUND contains too few elements (NCOMP = " + std::to_string(_disc.nComp) + " required)");

	if (paramProvider.exists("NPARTYPE"))
		_disc.nParType = paramProvider.getInt("NPARTYPE");
	else
	{
		// Infer number of particle types
		_disc.nParType = std::max(nBound.size() / _disc.nComp, nParCell.size());
	}

	if ((nParCell.size() > 1) && (nParCell.size() < _disc.nParType))
		throw InvalidParameterException("Field NPAR must have 1 or NPARTYPE (" + std::to_string(_disc.nParType) + ") entries");

	_disc.nParCell = new unsigned int[_disc.nParType];
	if (nParCell.size() < _disc.nParType)
	{
		// Multiplex number of particle shells to all particle types
		for (unsigned int i = 0; i < _disc.nParType; ++i)
			std::fill(_disc.nParCell, _disc.nParCell + _disc.nParType, nParCell[0]);
	}
	else
		std::copy_n(nParCell.begin(), _disc.nParType, _disc.nParCell);

	if ((nBound.size() > _disc.nComp) && (nBound.size() < _disc.nComp * _disc.nParType))
		throw InvalidParameterException("Field NBOUND must have NCOMP (" + std::to_string(_disc.nComp) + ") or NCOMP * NPARTYPE (" + std::to_string(_disc.nComp * _disc.nParType) + ") entries");

	_disc.nBound = new unsigned int[_disc.nComp * _disc.nParType];
	if (nBound.size() < _disc.nComp * _disc.nParType)
	{
		// Multiplex number of bound states to all particle types
		for (unsigned int i = 0; i < _disc.nParType; ++i)
			std::copy_n(nBound.begin(), _disc.nComp, _disc.nBound + i * _disc.nComp);
	}
	else
		std::copy_n(nBound.begin(), _disc.nComp * _disc.nParType, _disc.nBound);

	const unsigned int nTotalBound = std::accumulate(_disc.nBound, _disc.nBound + _disc.nComp * _disc.nParType, 0u);

	// Precompute offsets and total number of bound states (DOFs in solid phase)
	_disc.boundOffset = new unsigned int[_disc.nComp * _disc.nParType];
	_disc.strideBound = new unsigned int[_disc.nParType + 1];
	_disc.nBoundBeforeType = new unsigned int[_disc.nParType];
	_disc.strideBound[_disc.nParType] = nTotalBound;
	_disc.nBoundBeforeType[0] = 0;
	for (unsigned int j = 0; j < _disc.nParType; ++j)
	{
		unsigned int* const ptrOffset = _disc.boundOffset + j * _disc.nComp;
		unsigned int* const ptrBound = _disc.nBound + j * _disc.nComp;
		
		ptrOffset[0] = 0;
		for (unsigned int i = 1; i < _disc.nComp; ++i)
		{
			ptrOffset[i] = ptrOffset[i - 1] + ptrBound[i - 1];
		}
		_disc.strideBound[j] = ptrOffset[_disc.nComp - 1] + ptrBound[_disc.nComp - 1];

		if (j != _disc.nParType - 1)
			_disc.nBoundBeforeType[j + 1] = _disc.nBoundBeforeType[j] + _disc.strideBound[j];
	}

	// Precompute offsets of particle type DOFs
	_disc.parTypeOffset = new unsigned int[_disc.nParType + 1];
	_disc.nParCellsBeforeType = new unsigned int[_disc.nParType + 1];
	_disc.parTypeOffset[0] = 0;
	_disc.nParCellsBeforeType[0] = 0;
	unsigned int nTotalParCells = 0;
	for (unsigned int j = 1; j < _disc.nParType + 1; ++j)
	{
		_disc.parTypeOffset[j] = _disc.parTypeOffset[j-1] + (_disc.nComp + _disc.strideBound[j-1]) * _disc.nParCell[j-1] * _disc.nCol * _disc.nRad;
		_disc.nParCellsBeforeType[j] = _disc.nParCellsBeforeType[j-1] + _disc.nParCell[j-1];
		nTotalParCells += _disc.nParCell[j-1];
	}
	_disc.nParCellsBeforeType[_disc.nParType] = nTotalParCells;

	// Configure particle discretization
	_parCellSize.resize(nTotalParCells);
	_parCenterRadius.resize(nTotalParCells);
	_parOuterSurfAreaPerVolume.resize(nTotalParCells);
	_parInnerSurfAreaPerVolume.resize(nTotalParCells);

	// Read particle discretization mode and default to "EQUIDISTANT_PAR"
	_parDiscType = std::vector<ParticleDiscretizationMode>(_disc.nParType, ParticleDiscretizationMode::Equidistant);
	std::vector<std::string> pdt = paramProvider.getStringArray("PAR_DISC_TYPE");
	if ((pdt.size() == 1) && (_disc.nParType > 1))
	{
		// Multiplex using first value
		pdt.resize(_disc.nParType, pdt[0]);
	}
	else if (pdt.size() < _disc.nParType)
		throw InvalidParameterException("Field PAR_DISC_TYPE contains too few elements (" + std::to_string(_disc.nParType) + " required)");

	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		if (pdt[i] == "EQUIVOLUME_PAR")
			_parDiscType[i] = ParticleDiscretizationMode::Equivolume;
		else if (pdt[i] == "USER_DEFINED_PAR")
			_parDiscType[i] = ParticleDiscretizationMode::UserDefined;
	}

	if (paramProvider.exists("PAR_DISC_VECTOR"))
	{
		_parDiscVector = paramProvider.getDoubleArray("PAR_DISC_VECTOR");
		if (_parDiscVector.size() < nTotalParCells + _disc.nParType)
			throw InvalidParameterException("Field PAR_DISC_VECTOR contains too few elements (Sum [NPAR + 1] = " + std::to_string(nTotalParCells + _disc.nParType) + " required)");
	}

	// Determine whether analytic Jacobian should be used but don't set it right now.
	// We need to setup Jacobian matrices first.
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	const bool analyticJac = paramProvider.getBool("USE_ANALYTIC_JACOBIAN");
#else
	const bool analyticJac = false;
#endif

	// Initialize and configure GMRES for solving the Schur-complement
	_gmres.initialize(_disc.nCol * _disc.nRad * _disc.nComp * _disc.nParType, paramProvider.getInt("MAX_KRYLOV"), linalg::toOrthogonalization(paramProvider.getInt("GS_TYPE")), paramProvider.getInt("MAX_RESTARTS"));
	_gmres.matrixVectorMultiplier(&schurComplementMultiplierGRM2D, this);
	_schurSafety = paramProvider.getDouble("SCHUR_SAFETY");

	// Allocate space for initial conditions
	_initC.resize(_disc.nComp * _disc.nRad);
	_initCp.resize(_disc.nComp * _disc.nRad * _disc.nParType);
	_initQ.resize(nTotalBound * _disc.nRad);

	paramProvider.popScope();

	const bool transportSuccess = _convDispOp.configureModelDiscretization(paramProvider, _disc.nComp, _disc.nCol, _disc.nRad);

	// Allocate memory
	Indexer idxr(_disc);

	_jacInlet.resize(_disc.nComp * _disc.nRad);

	_jacP = new linalg::BandMatrix[_disc.nCol * _disc.nRad * _disc.nParType];
	_jacPdisc = new linalg::FactorizableBandMatrix[_disc.nCol * _disc.nRad * _disc.nParType];
	for (unsigned int j = 0; j < _disc.nParType; ++j)
	{
		linalg::BandMatrix* const ptrJac = _jacP + _disc.nCol * _disc.nRad * j;
		linalg::FactorizableBandMatrix* const ptrJacDisc = _jacPdisc + _disc.nCol * _disc.nRad * j;
		for (unsigned int i = 0; i < _disc.nCol * _disc.nRad; ++i)
		{
			ptrJacDisc[i].resize(_disc.nParCell[j] * (_disc.nComp + _disc.strideBound[j]), _disc.nComp + _disc.strideBound[j], _disc.nComp + 2 * _disc.strideBound[j]);
			ptrJac[i].resize(_disc.nParCell[j] * (_disc.nComp + _disc.strideBound[j]), _disc.nComp + _disc.strideBound[j], _disc.nComp + 2 * _disc.strideBound[j]);
		}
	}

	_jacPF = new linalg::DoubleSparseMatrix[_disc.nCol * _disc.nRad * _disc.nParType];
	_jacFP = new linalg::DoubleSparseMatrix[_disc.nCol * _disc.nRad * _disc.nParType];
	for (unsigned int i = 0; i < _disc.nCol * _disc.nRad * _disc.nParType; ++i)
	{
		_jacPF[i].resize(_disc.nComp);
		_jacFP[i].resize(_disc.nComp);
	}

	_jacCF.resize(_disc.nComp * _disc.nCol * _disc.nRad * _disc.nParType);
	_jacFC.resize(_disc.nComp * _disc.nCol * _disc.nRad * _disc.nParType);

	_discParFlux.resize(sizeof(active) * _disc.nComp);

	// Set whether analytic Jacobian is used
	useAnalyticJacobian(analyticJac);

	// ==== Construct and configure binding model
	clearBindingModels();
	_binding = std::vector<IBindingModel*>(_disc.nParType, nullptr);

	const std::vector<std::string> bindModelNames = paramProvider.getStringArray("ADSORPTION_MODEL");

	if (paramProvider.exists("ADSORPTION_MODEL_MULTIPLEX"))
		_singleBinding = (paramProvider.getInt("ADSORPTION_MODEL_MULTIPLEX") == 1);
	else
	{
		// Infer multiplex mode
		_singleBinding = (bindModelNames.size() == 1);
	}

	if (!_singleBinding && (bindModelNames.size() < _disc.nParType))
		throw InvalidParameterException("Field ADSORPTION_MODEL contains too few elements (" + std::to_string(_disc.nParType) + " required)");
	else if (_singleBinding && (bindModelNames.size() != 1))
		throw InvalidParameterException("Field ADSORPTION_MODEL requires (only) 1 element");

	bool bindingConfSuccess = true;
	unsigned int bindingRequiredMem = 0;
	_bindingWorkspaceOffset = new unsigned int[_disc.nParType];
	_bindingWorkspaceOffset[0] = 0;
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		if (_singleBinding && (i > 0))
		{
			// Reuse first binding model
			_binding[i] = _binding[0];
		}
		else
		{
			_binding[i] = helper.createBindingModel(bindModelNames[i]);
			if (!_binding[i])
				throw InvalidParameterException("Unknown binding model " + bindModelNames[i]);

			bindingConfSuccess = _binding[i]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound + i * _disc.nComp, _disc.boundOffset + i * _disc.nComp) && bindingConfSuccess;
		}

		// In case of a single binding model, we still require the additional memory per type because of parallelization
		if (_binding[i]->requiresWorkspace())
		{
			// Required memory (number of doubles) for nonlinear solvers
			const unsigned int requiredMem = (_binding[i]->workspaceSize(_disc.nComp, _disc.strideBound[i], _disc.nBound + i * _disc.nComp) + sizeof(double) - 1) / sizeof(double) * _disc.nCol * _disc.nRad;
			bindingRequiredMem += requiredMem;

			if (i != _disc.nParType - 1)
				_bindingWorkspaceOffset[i+1] = _bindingWorkspaceOffset[i] + requiredMem;
		}
		else
		{
			if (i != _disc.nParType - 1)
				_bindingWorkspaceOffset[i+1] = _bindingWorkspaceOffset[i];
		}
	}

	// setup the memory for tempState based on state vector or memory needed for consistent initialization of isotherms, whichever is larger
	const unsigned int size = std::max(numDofs(), bindingRequiredMem);
	_tempState = new double[size];

	return transportSuccess && bindingConfSuccess;
}

bool GeneralRateModel2D::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	const bool transportSuccess = _convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

	// Read geometry parameters
	_singleParRadius = readScalarParameterOrArray(_parRadius, paramProvider, "PAR_RADIUS", _disc.nParType);
	_singleParPorosity = readScalarParameterOrArray(_parPorosity, paramProvider, "PAR_POROSITY", _disc.nParType);

	// Let PAR_CORERADIUS default to 0.0 for backwards compatibility
	if (paramProvider.exists("PAR_CORERADIUS"))
		_singleParCoreRadius = readScalarParameterOrArray(_parCoreRadius, paramProvider, "PAR_CORERADIUS", _disc.nParType);
	else
	{
		_singleParCoreRadius = true;
		_parCoreRadius = std::vector<active>(_disc.nParType, 0.0);
	}

	// Check whether PAR_TYPE_VOLFRAC is required or not
	if ((_disc.nParType > 1) && !paramProvider.exists("PAR_TYPE_VOLFRAC"))
		throw InvalidParameterException("The required parameter \"PAR_TYPE_VOLFRAC\" was not found");

	// Let PAR_TYPE_VOLFRAC default to 1.0 for backwards compatibility
	if (paramProvider.exists("PAR_TYPE_VOLFRAC"))
		_parTypeVolFracMode = readAndRegisterMultiplexParam(paramProvider, _parameters, _parTypeVolFrac, "PAR_TYPE_VOLFRAC", _disc.nCol, _disc.nRad, _disc.nParType, _unitOpIdx);
	else
	{
		// Only one particle type present
		_parTypeVolFrac.resize(_disc.nCol * _disc.nRad, 1.0);
		_parTypeVolFracMode = MultiplexMode::Independent;
	}

	// Check whether all sizes are matched
	if (_disc.nParType != _parRadius.size())
		throw InvalidParameterException("Number of elements in field PAR_RADIUS does not match number of particle types");
	if (_disc.nParType * _disc.nCol * _disc.nRad != _parTypeVolFrac.size())
		throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types times number of bulk cells");
	if (_disc.nParType != _parPorosity.size())
		throw InvalidParameterException("Number of elements in field PAR_POROSITY does not match number of particle types");
	if (_disc.nParType != _parCoreRadius.size())
		throw InvalidParameterException("Number of elements in field PAR_CORERADIUS does not match number of particle types");

	// Check that particle volume fractions sum to 1.0
	for (unsigned int i = 0; i < _disc.nCol * _disc.nRad; ++i)
	{
		const double volFracSum = std::accumulate(_parTypeVolFrac.begin() + i * _disc.nParType, _parTypeVolFrac.begin() + (i+1) * _disc.nParType, 0.0, 
			[](double a, const active& b) -> double { return a + static_cast<double>(b); });
		if (std::abs(1.0 - volFracSum) > 1e-10)
			throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ") in axial cell " + std::to_string(i / _disc.nRad) + " radial cell " + std::to_string(i % _disc.nRad));
	}

	// Read vectorial parameters (which may also be section dependent; transport)
	_filmDiffusionMode = readAndRegisterMultiplexParam(paramProvider, _parameters, _filmDiffusion, "FILM_DIFFUSION", _disc.nParType, _disc.nComp, _unitOpIdx);
	_parDiffusionMode = readAndRegisterMultiplexParam(paramProvider, _parameters, _parDiffusion, "PAR_DIFFUSION", _disc.nParType, _disc.nComp, _unitOpIdx);
	
	if (paramProvider.exists("PAR_SURFDIFFUSION"))
		_parSurfDiffusionMode = readAndRegisterMultiplexParam(paramProvider, _parameters, _parSurfDiffusion, "PAR_SURFDIFFUSION", _disc.nParType, _disc.nComp, _disc.strideBound, _disc.nBound, _unitOpIdx);
	else
	{
		_parSurfDiffusionMode = MultiplexMode::Component;
		_parSurfDiffusion.resize(_disc.strideBound[_disc.nParType], 0.0);
	}

	if ((_filmDiffusion.size() < _disc.nComp * _disc.nParType) || (_filmDiffusion.size() % (_disc.nComp * _disc.nParType) != 0))
		throw InvalidParameterException("Number of elements in field FILM_DIFFUSION is not a positive multiple of NCOMP * NPARTYPE (" + std::to_string(_disc.nComp * _disc.nParType) + ")");
	if ((_parDiffusion.size() < _disc.nComp * _disc.nParType) || (_parDiffusion.size() % (_disc.nComp * _disc.nParType) != 0))
		throw InvalidParameterException("Number of elements in field PAR_DIFFUSION is not a positive multiple of NCOMP * NPARTYPE (" + std::to_string(_disc.nComp * _disc.nParType) + ")");
	if ((_parSurfDiffusion.size() < _disc.strideBound[_disc.nParType]) || ((_disc.strideBound[_disc.nParType] > 0) && (_parSurfDiffusion.size() % _disc.strideBound[_disc.nParType] != 0)))
		throw InvalidParameterException("Number of elements in field PAR_SURFDIFFUSION is not a positive multiple of NTOTALBND (" + std::to_string(_disc.strideBound[_disc.nParType]) + ")");

	if (paramProvider.exists("PORE_ACCESSIBILITY"))
		_poreAccessFactorMode = readAndRegisterMultiplexParam(paramProvider, _parameters, _poreAccessFactor, "PORE_ACCESSIBILITY", _disc.nParType, _disc.nComp, _unitOpIdx);
	else
	{
		_poreAccessFactorMode = MultiplexMode::ComponentType;
		_poreAccessFactor = std::vector<cadet::active>(_disc.nComp * _disc.nParType, 1.0);
	}

	if (_disc.nComp * _disc.nParType != _poreAccessFactor.size())
		throw InvalidParameterException("Number of elements in field PORE_ACCESSIBILITY differs from NCOMP * NPARTYPE (" + std::to_string(_disc.nComp * _disc.nParType) + ")");

	// Add parameters to map
	if (_singleParRadius)
		_parameters[makeParamId(hashString("PAR_RADIUS"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parRadius[0];
	else
		registerParam1DArray(_parameters, _parRadius, [=](bool multi, unsigned int type) { return makeParamId(hashString("PAR_RADIUS"), _unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, SectionIndep); });

	if (_singleParCoreRadius)
		_parameters[makeParamId(hashString("PAR_CORERADIUS"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parCoreRadius[0];
	else
		registerParam1DArray(_parameters, _parCoreRadius, [=](bool multi, unsigned int type) { return makeParamId(hashString("PAR_CORERADIUS"), _unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, SectionIndep); });

	if (_singleParPorosity)
		_parameters[makeParamId(hashString("PAR_POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parPorosity[0];
	else
		registerParam1DArray(_parameters, _parPorosity, [=](bool multi, unsigned int type) { return makeParamId(hashString("PAR_POROSITY"), _unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, SectionIndep); });

	// Calculate the particle radial discretization variables (_parCellSize, _parCenterRadius, etc.)
	updateRadialDisc();

	// Register initial conditions parameters
	registerParam1DArray(_parameters, _initC, [=](bool multi, unsigned int comp) { return makeParamId(hashString("INIT_C"), _unitOpIdx, comp, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep); });
	if (_disc.nRad > 1)
		registerParam2DArray(_parameters, _initC, [=](bool multi, unsigned int rad, unsigned int comp) { return makeParamId(hashString("INIT_C"), _unitOpIdx, comp, ParTypeIndep, BoundStateIndep, rad, SectionIndep); }, _disc.nComp);

	if (_singleBinding)
	{
		for (unsigned int c = 0; c < _disc.nComp; ++c)
			_parameters[makeParamId(hashString("INIT_CP"), _unitOpIdx, c, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_initCp[c];

		if (_disc.nRad > 1)
		{
			for (unsigned int r = 0; r < _disc.nRad; ++r)
			{
				for (unsigned int c = 0; c < _disc.nComp; ++c)
					_parameters[makeParamId(hashString("INIT_CP"), _unitOpIdx, c, ParTypeIndep, BoundStateIndep, r, SectionIndep)] = &_initCp[r * _disc.nComp * _disc.nParType + c];
			}
		}
	}
	else
	{
		registerParam2DArray(_parameters, _initCp, [=](bool multi, unsigned int type, unsigned int comp) { return makeParamId(hashString("INIT_CP"), _unitOpIdx, comp, type, BoundStateIndep, ReactionIndep, SectionIndep); }, _disc.nComp);
		if (_disc.nRad > 1)
			registerParam3DArray(_parameters, _initCp, [=](bool multi, unsigned int rad, unsigned int type, unsigned int comp) { return makeParamId(hashString("INIT_CP"), _unitOpIdx, comp, type, BoundStateIndep, rad, SectionIndep); }, _disc.nComp, _disc.nParType);
	}

	if (!_binding.empty())
	{
		const unsigned int maxBoundStates = *std::max_element(_disc.strideBound, _disc.strideBound + _disc.nParType);
		std::vector<ParameterId> initParams(maxBoundStates);

		// Register radially independent
		if (_singleBinding)
		{
			_binding[0]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, ParTypeIndep);

			active* const iq = _initQ.data() + _disc.nBoundBeforeType[0];
			for (unsigned int i = 0; i < _disc.strideBound[0]; ++i)
				_parameters[initParams[i]] = iq + i;
		}
		else
		{
			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
				_binding[type]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, type);

				active* const iq = _initQ.data() + _disc.nBoundBeforeType[type];
				for (unsigned int i = 0; i < _disc.strideBound[type]; ++i)
					_parameters[initParams[i]] = iq + i;
			}
		}

		// Register radially dependent
		if (_disc.nRad > 1)
		{
			if (_singleBinding)
			{
				_binding[0]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, ParTypeIndep);

				for (unsigned int r = 0; r < _disc.nRad; ++r)
				{
					for (ParameterId& pId : initParams)
						pId.reaction = r;

					active* const iq = _initQ.data() + _disc.nBoundBeforeType[0] + r * _disc.strideBound[_disc.nParType];
					for (unsigned int i = 0; i < _disc.strideBound[0]; ++i)
						_parameters[initParams[i]] = iq + i;
				}
			}
			else
			{
				for (unsigned int r = 0; r < _disc.nRad; ++r)
				{
					for (unsigned int type = 0; type < _disc.nParType; ++type)
					{
						_binding[type]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, type);
						for (ParameterId& pId : initParams)
							pId.reaction = r;

						active* const iq = _initQ.data() + _disc.nBoundBeforeType[type] + r * _disc.strideBound[_disc.nParType];
						for (unsigned int i = 0; i < _disc.strideBound[type]; ++i)
							_parameters[initParams[i]] = iq + i;
					}
				}
			}
		}
	}

	// Reconfigure binding model
	if (!_binding.empty())
	{
		bool bindingConfSuccess = true;
		if (_singleBinding)
		{
			if (_binding[0] && _binding[0]->requiresConfiguration())
			{
				if (paramProvider.exists("adsorption"))
					paramProvider.pushScope("adsorption");
				else if (paramProvider.exists("adsorption_000"))
					paramProvider.pushScope("adsorption_000");
				else
					throw InvalidParameterException("Group \"adsorption\" or \"adsorption_000\" required");

				bindingConfSuccess = _binding[0]->configure(paramProvider, _unitOpIdx, ParTypeIndep);
				paramProvider.popScope();
			}
		}
		else
		{
			std::ostringstream oss;
			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
	 			if (!_binding[type] || !_binding[type]->requiresConfiguration())
	 				continue;

				oss.str("");
				oss << "adsorption_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << type;

	 			// If there is just one type, allow legacy "adsorption" scope
				if (!paramProvider.exists(oss.str()) && (_disc.nParType == 1))
					oss.str("adsorption");

				if (!paramProvider.exists(oss.str()))
					continue;

				paramProvider.pushScope(oss.str());
				bindingConfSuccess = _binding[type]->configure(paramProvider, _unitOpIdx, type) && bindingConfSuccess;
				paramProvider.popScope();
			}
		}

		return transportSuccess && bindingConfSuccess;
	}

	return transportSuccess;
}

unsigned int GeneralRateModel2D::numAdDirsForJacobian() const CADET_NOEXCEPT
{
	// We need as many directions as the highest bandwidth of the diagonal blocks:
	// The bandwidth of the column block depends on the size of the WENO stencil, whereas
	// the bandwidth of the particle blocks are given by the number of components and bound states.

	// Get maximum stride of particle type blocks
	unsigned int maxStride = 0;
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		maxStride = std::max(maxStride, _jacP[type * _disc.nCol * _disc.nRad].stride());
	}

	return maxStride;
}

void GeneralRateModel2D::useAnalyticJacobian(const bool analyticJac)
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	_analyticJac = analyticJac;
	if (!_analyticJac)
		_jacobianAdDirs = numAdDirsForJacobian();
	else
		_jacobianAdDirs = 0;
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always enable AD for comparison and use it in simulation
	_analyticJac = false;
	_jacobianAdDirs = numAdDirsForJacobian();
#endif
}

void GeneralRateModel2D::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const AdJacobianParams& adJac)
{
	// Setup flux Jacobian blocks at the beginning of the simulation or in case of
	// section dependent film or particle diffusion coefficients
	if ((secIdx == 0) || (_filmDiffusion.size() > _disc.nComp * _disc.nParType) || (_parDiffusion.size() > _disc.nComp * _disc.nParType))
		assembleOffdiagJac(t, secIdx);

	Indexer idxr(_disc);

	// ConvectionDispersionOperator tells us whether flow direction has changed
	if (!_convDispOp.notifyDiscontinuousSectionTransition(t, secIdx))
		return;

	// Setup the matrix connecting inlet DOFs to first column cells
	_jacInlet.clear();

	for (unsigned int rad = 0; rad < _disc.nRad; ++rad)
	{
		const double f = _convDispOp.inletFactor(secIdx, rad);
		if (_convDispOp.currentVelocity(rad) >= 0.0)
		{
			// Forwards flow

			// Place entries for inlet DOF to first axial cell conversion
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				_jacInlet.addElement(comp * idxr.strideColComp() + rad * idxr.strideColRadialCell(), comp, f);
		}
		else
		{
			// Backwards flow

			// Place entries for inlet DOF to last column cell conversion
			const unsigned int offset = (_disc.nCol - 1) * idxr.strideColAxialCell();
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				_jacInlet.addElement(offset + comp * idxr.strideColComp() + rad * idxr.strideColRadialCell(), comp, f);
		}
	}
}

void GeneralRateModel2D::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	_convDispOp.setFlowRates(in, out);
}

void GeneralRateModel2D::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, *this, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void GeneralRateModel2D::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, *this, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}


unsigned int GeneralRateModel2D::requiredADdirs() const CADET_NOEXCEPT
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	return _jacobianAdDirs;
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always need the AD directions for the Jacobian
	return numAdDirsForJacobian();
#endif
}

void GeneralRateModel2D::prepareADvectors(const AdJacobianParams& adJac) const
{
	// Early out if AD is disabled
	if (!adJac.adY)
		return;

	Indexer idxr(_disc);

	// Particle blocks
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		const unsigned int lowerParBandwidth = _jacP[type * _disc.nCol * _disc.nRad].lowerBandwidth();
		const unsigned int upperParBandwidth = _jacP[type * _disc.nCol * _disc.nRad].upperBandwidth();

		for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nRad; ++pblk)
		{
			ad::prepareAdVectorSeedsForBandMatrix(adJac.adY + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{pblk}), adJac.adDirOffset, idxr.strideParBlock(type), lowerParBandwidth, upperParBandwidth, lowerParBandwidth);
		}
	}
}

/**
 * @brief Extracts the system Jacobian from band compressed AD seed vectors
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void GeneralRateModel2D::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	Indexer idxr(_disc);

	// Particles
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nRad; ++pblk)
		{
			linalg::BandMatrix& jacMat = _jacP[_disc.nCol * _disc.nRad * type + pblk];
			ad::extractBandedJacobianFromAd(adRes + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{pblk}), adDirOffset, jacMat.lowerBandwidth(), jacMat);
		}
	}
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the corresponding band matrices.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void GeneralRateModel2D::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	Indexer idxr(_disc);

	LOG(Debug) << "AD dir offset: " << adDirOffset << " DiagDirCol: " << _convDispOp.jacobian().lowerBandwidth() << " DiagDirPar: " << _jacP[0].lowerBandwidth();

	// Particles
	double maxDiffPar = 0.0;
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nRad; ++pblk)
		{
			linalg::BandMatrix& jacMat = _jacP[_disc.nCol * _disc.nRad * type + pblk];
			const double localDiff = ad::compareBandedJacobianWithAd(adRes + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{pblk}), adDirOffset, jacMat.lowerBandwidth(), jacMat);
			LOG(Debug) << "-> Par type " << type << " block " << pblk << " diff: " << localDiff;
			maxDiffPar = std::max(maxDiffPar, localDiff);
		}
	}
}

#endif

int GeneralRateModel2D::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual do not compute Jacobian or parameter sensitivities
	return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simTime.timeFactor, simState.vecStateY, simState.vecStateYdot, res);
}

int GeneralRateModel2D::residualWithJacobian(const ActiveSimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual, use AD for Jacobian if required but do not evaluate parameter derivatives
	return residual(simTime, simState, res, adJac, true, false);
}

int GeneralRateModel2D::residual(const ActiveSimulationTime& simTime, const ConstSimulationState& simState, double* const res, 
	const AdJacobianParams& adJac, bool updateJacobian, bool paramSensitivity)
{
	if (updateJacobian)
	{
		_factorizeJacobian = true;

#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
		if (_analyticJac)
		{
			if (paramSensitivity)
			{
				const int retCode = residualImpl<double, active, active, true>(simTime.t, simTime.secIdx, simTime.timeFactor, simState.vecStateY, simState.vecStateYdot, adJac.adRes);

				// Copy AD residuals to original residuals vector
				if (res)
					ad::copyFromAd(adJac.adRes, res, numDofs());

				return retCode;
			}
			else
				return residualImpl<double, double, double, true>(static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor), simState.vecStateY, simState.vecStateYdot, res);
		}
		else
		{
			// Compute Jacobian via AD

			// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
			// and initalize residuals with zero (also resetting directional values)
			ad::copyToAd(simState.vecStateY, adJac.adY, numDofs());
			// @todo Check if this is necessary
			ad::resetAd(adJac.adRes, numDofs());

			// Evaluate with AD enabled
			int retCode = 0;
			if (paramSensitivity)
				retCode = residualImpl<active, active, active, false>(simTime.t, simTime.secIdx, simTime.timeFactor, adJac.adY, simState.vecStateYdot, adJac.adRes);
			else
				retCode = residualImpl<active, active, double, false>(static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor), adJac.adY, simState.vecStateYdot, adJac.adRes);

			// Copy AD residuals to original residuals vector
			if (res)
				ad::copyFromAd(adJac.adRes, res, numDofs());

			// Extract Jacobian
			extractJacobianFromAD(adJac.adRes, adJac.adDirOffset);

			return retCode;
		}
#else
		// Compute Jacobian via AD

		// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
		// and initalize residuals with zero (also resetting directional values)
		ad::copyToAd(simState.vecStateY, adJac.adY, numDofs());
		// @todo Check if this is necessary
		ad::resetAd(adJac.adRes, numDofs());

		// Evaluate with AD enabled
		int retCode = 0;
		if (paramSensitivity)
			retCode = residualImpl<active, active, active, false>(simTime.t, simTime.secIdx, simTime.timeFactor, adJac.adY, simState.vecStateYdot, adJac.adRes);
		else
			retCode = residualImpl<active, active, double, false>(static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor), adJac.adY, simState.vecStateYdot, adJac.adRes);

		// Only do comparison if we have a residuals vector (which is not always the case)
		if (res)
		{
			// Evaluate with analytical Jacobian which is stored in the band matrices
			retCode = residualImpl<double, double, double, true>(static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor), simState.vecStateY, simState.vecStateYdot, res);

			// Compare AD with anaytic Jacobian
			checkAnalyticJacobianAgainstAd(adJac.adRes, adJac.adDirOffset);
		}

		// Extract Jacobian
		extractJacobianFromAD(adJac.adRes, adJac.adDirOffset);

		return retCode;
#endif
	}
	else
	{
		if (paramSensitivity)
		{
			// Initalize residuals with zero
			// @todo Check if this is necessary
			ad::resetAd(adJac.adRes, numDofs());

			const int retCode = residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simTime.timeFactor, simState.vecStateY, simState.vecStateYdot, adJac.adRes);

			// Copy AD residuals to original residuals vector
			if (res)
				ad::copyFromAd(adJac.adRes, res, numDofs());

			return retCode;
		}
		else
			return residualImpl<double, double, double, false>(static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor), simState.vecStateY, simState.vecStateYdot, res);
	}
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int GeneralRateModel2D::residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* const y, double const* const yDot, ResidualType* const res)
{
	BENCH_START(_timerResidualPar);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), size_t(_disc.nCol * _disc.nRad * _disc.nParType + 1), [&](size_t pblk)
#else
	for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nRad * _disc.nParType + 1; ++pblk)
#endif
	{
		if (cadet_unlikely(pblk == 0))
			_convDispOp.residual(t, secIdx, timeFactor, y, yDot, res, wantJac);
		else
		{
			const unsigned int type = (pblk - 1) / (_disc.nCol * _disc.nRad);
			const unsigned int par = (pblk - 1) % (_disc.nCol * _disc.nRad);
			residualParticle<StateType, ResidualType, ParamType, wantJac>(t, type, par, secIdx, timeFactor, y, yDot, res);
		}
	} CADET_PARFOR_END;

	BENCH_STOP(_timerResidualPar);

	residualFlux<StateType, ResidualType, ParamType>(t, secIdx, y, yDot, res);

	// Handle inlet DOFs, which are simply copied to res
	for (unsigned int i = 0; i < _disc.nComp * _disc.nRad; ++i)
	{
		res[i] = y[i];
	}

	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int GeneralRateModel2D::residualParticle(const ParamType& t, unsigned int parType, unsigned int colCell, unsigned int secIdx, const ParamType& timeFactor, StateType const* yBase, double const* yDotBase, ResidualType* resBase)
{
	Indexer idxr(_disc);

	// Go to the particle block of the given column cell
	StateType const* y = yBase + idxr.offsetCp(ParticleTypeIndex{parType}, ParticleIndex{colCell});
	double const* yDot = yDotBase + idxr.offsetCp(ParticleTypeIndex{parType}, ParticleIndex{colCell});
	ResidualType* res = resBase + idxr.offsetCp(ParticleTypeIndex{parType}, ParticleIndex{colCell});

	const unsigned int requiredMem = (_binding[parType]->workspaceSize(_disc.nComp, _disc.strideBound[parType], _disc.nBound + parType * _disc.nComp) + sizeof(double) - 1) / sizeof(double);
	double* const buffer = _tempState + _bindingWorkspaceOffset[parType] + requiredMem * colCell;

	// Prepare parameters
	active const* const parDiff = getSectionDependentSlice(_parDiffusion, _disc.nComp * _disc.nParType, secIdx) + parType * _disc.nComp;

	// Ordering of particle surface diffusion:
	// bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2
	active const* const parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _disc.strideBound[_disc.nParType], secIdx) + _disc.nBoundBeforeType[parType];

	// Midpoint of current column cell (z, rho coordinate) - needed in externally dependent adsorption kinetic
	const unsigned int axialCell = colCell / _disc.nRad;
	const unsigned int radialCell = colCell % _disc.nRad;
	const double r = static_cast<double>(_convDispOp.radialCenters()[radialCell]) / static_cast<double>(_convDispOp.columnRadius());
	const double z = (0.5 + static_cast<double>(axialCell)) / static_cast<double>(_disc.nCol);

	// Reset Jacobian
	if (wantJac)
		_jacP[_disc.nCol * _disc.nRad * parType + colCell].setAll(0.0);

	// The RowIterator is always centered on the main diagonal.
	// This means that jac[0] is the main diagonal, jac[-1] is the first lower diagonal,
	// and jac[1] is the first upper diagonal. We can also access the rows from left to
	// right beginning with the last lower diagonal moving towards the main diagonal and
	// continuing to the last upper diagonal by using the native() method.
	linalg::BandMatrix::RowIterator jac = _jacP[_disc.nCol * _disc.nRad * parType + colCell].row(0);

	active const* const outerSurfPerVol = _parOuterSurfAreaPerVolume.data() + _disc.nParCellsBeforeType[parType];
	active const* const innerSurfPerVol = _parInnerSurfAreaPerVolume.data() + _disc.nParCellsBeforeType[parType];
	active const* const parCenterRadius = _parCenterRadius.data() + _disc.nParCellsBeforeType[parType];

	// Loop over particle cells
	for (unsigned int par = 0; par < _disc.nParCell[parType]; ++par)
	{
		// Geometry
		const ParamType outerAreaPerVolume = static_cast<ParamType>(outerSurfPerVol[par]);
		const ParamType innerAreaPerVolume = static_cast<ParamType>(innerSurfPerVol[par]);

		// Mobile phase
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++res, ++y, ++yDot, ++jac)
		{
			*res = 0.0;
			const unsigned int nBound = _disc.nBound[_disc.nComp * parType + comp];
			const ParamType invBetaP = (1.0 - static_cast<ParamType>(_parPorosity[parType])) / (static_cast<ParamType>(_poreAccessFactor[_disc.nComp * parType + comp]) * static_cast<ParamType>(_parPorosity[parType]));

			// Add time derivatives
			if (yDotBase)
			{
				// Ultimately, we need dc_{p,comp} / dt + 1 / beta_p * [ sum_i  dq_comp^i / dt ]
				// Compute the sum in the brackets first, then divide by beta_p and add dc_p / dt

				// Sum dq_comp^1 / dt + dq_comp^2 / dt + ... + dq_comp^{N_comp} / dt
				for (unsigned int i = 0; i < nBound; ++i)
					// Index explanation:
					//   -comp -> go back to beginning of liquid phase
					//   + strideParLiquid() skip to solid phase
					//   + offsetBoundComp() jump to component (skips all bound states of previous components)
					//   + i go to current bound state
					// Remember this, you'll see it quite a lot ...
					*res += yDot[idxr.strideParLiquid() - comp + idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i];

				// Divide by beta_p and add dcp_i / dt
				*res = timeFactor * (yDot[0] + invBetaP * res[0]);
			}

			const ParamType dp = static_cast<ParamType>(parDiff[comp]);

			// Add flow through outer surface
			// Note that inflow boundary conditions are handled in residualFlux().
			if (cadet_likely(par != 0))
			{
				// Difference between two cell-centers
				const ParamType dr = static_cast<ParamType>(parCenterRadius[par - 1]) - static_cast<ParamType>(parCenterRadius[par]);

				// Molecular diffusion contribution
				const ResidualType gradCp = (y[-idxr.strideParShell(parType)] - y[0]) / dr;
				*res -= outerAreaPerVolume * dp * gradCp;

				// Surface diffusion contribution
				for (unsigned int i = 0; i < nBound; ++i)
				{
					// See above for explanation of curIdx value
					const int curIdx = idxr.strideParLiquid() - comp + idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i;
					const ResidualType gradQ = (y[-idxr.strideParShell(parType) + curIdx] - y[curIdx]) / dr;
					*res -= outerAreaPerVolume * static_cast<ParamType>(parSurfDiff[idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i]) * invBetaP * gradQ;
				}

				if (wantJac)
				{
					const double localInvBetaP = static_cast<double>(invBetaP);
					const double ouApV = static_cast<double>(outerAreaPerVolume);
					const double ldr = static_cast<double>(dr);

					// Liquid phase
					jac[0] += ouApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j)
					jac[-idxr.strideParShell(parType)] = -ouApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j-1)

					// Solid phase
					for (unsigned int i = 0; i < nBound; ++i)
					{
						// See above for explanation of curIdx value
						const int curIdx = idxr.strideParLiquid() - comp + idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i;
						jac[curIdx] += ouApV * localInvBetaP * static_cast<double>(parSurfDiff[idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i]) / ldr; // dres / dq_i^(p,j)
						jac[-idxr.strideParShell(parType) + curIdx] = -ouApV * localInvBetaP * static_cast<double>(parSurfDiff[idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i]) / ldr; // dres / dq_i^(p,j-1)
					}
				}
			}

			// Add flow through inner surface
			// Note that this term vanishes for the most inner shell due to boundary conditions
			if (cadet_likely(par != _disc.nParCell[parType] - 1))
			{
				// Difference between two cell-centers
				const ParamType dr = static_cast<ParamType>(parCenterRadius[par]) - static_cast<ParamType>(parCenterRadius[par + 1]);

				// Molecular diffusion contribution
				const ResidualType gradCp = (y[0] - y[idxr.strideParShell(parType)]) / dr;
				*res += innerAreaPerVolume * dp * gradCp;

				// Surface diffusion contribution
				for (unsigned int i = 0; i < nBound; ++i)
				{
					// See above for explanation of curIdx value
					const unsigned int curIdx = idxr.strideParLiquid() - comp + idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i;
					const ResidualType gradQ = (y[curIdx] - y[idxr.strideParShell(parType) + curIdx]) / dr;
					*res += innerAreaPerVolume * static_cast<ParamType>(parSurfDiff[idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i]) * invBetaP * gradQ;
				}

				if (wantJac)
				{
					const double localInvBetaP = static_cast<double>(invBetaP);
					const double inApV = static_cast<double>(innerAreaPerVolume);
					const double ldr = static_cast<double>(dr);

					// Liquid phase
					jac[0] += inApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j)
					jac[idxr.strideParShell(parType)] = -inApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j+1)

					// Solid phase
					for (unsigned int i = 0; i < nBound; ++i)
					{
						// See above for explanation of curIdx value
						const int curIdx = idxr.strideParLiquid() - comp + idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i;
						jac[curIdx] += inApV * localInvBetaP * static_cast<double>(parSurfDiff[idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i]) / ldr; // dres / dq_i^(p,j)
						jac[idxr.strideParShell(parType) + curIdx] = -inApV * localInvBetaP * static_cast<double>(parSurfDiff[idxr.offsetBoundComp(ParticleTypeIndex{parType}, ComponentIndex{comp}) + i]) / ldr; // dres / dq_i^(p,j-1)
					}
				}
			}
		}

		// Bound phases
		if (!yDotBase)
			yDot = nullptr;

		_binding[parType]->residual(t, secIdx, timeFactor, ColumnPosition{z, r, static_cast<double>(parCenterRadius[par]) / static_cast<double>(_parRadius[parType])}, y, y - _disc.nComp, yDot, res, buffer);
		if (wantJac)
		{
			// static_cast should be sufficient here, but this statement is also analyzed when wantJac = false
			_binding[parType]->analyticJacobian(static_cast<double>(t), secIdx, ColumnPosition{z, r, static_cast<double>(parCenterRadius[par]) / static_cast<double>(_parRadius[parType])}, reinterpret_cast<double const*>(y), _disc.nComp, jac, buffer);
		}

		// Advance pointers over all bound states
		y += idxr.strideParBound(parType);
		yDot += idxr.strideParBound(parType);
		res += idxr.strideParBound(parType);
		jac += idxr.strideParBound(parType);
	}
	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType>
int GeneralRateModel2D::residualFlux(const ParamType& t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase)
{
	Indexer idxr(_disc);

	// Get offsets
	ResidualType* const resCol = resBase + idxr.offsetC();
	ResidualType* const resFlux = resBase + idxr.offsetJf();

	StateType const* const yCol = yBase + idxr.offsetC();
	StateType const* const yFlux = yBase + idxr.offsetJf();

	// J_f block (identity matrix), adds flux state to flux equation
	for (unsigned int i = 0; i < _disc.nComp * _disc.nCol * _disc.nRad * _disc.nParType; ++i)
		resFlux[i] = yFlux[i];

	// Discretized film diffusion kf for finite volumes
	ParamType* const kf_FV = _discParFlux.create<ParamType>(_disc.nComp);

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		ResidualType* const resParType = resBase + idxr.offsetCp(ParticleTypeIndex{type});
		ResidualType* const resFluxType = resBase + idxr.offsetJf(ParticleTypeIndex{type});

		StateType const* const yParType = yBase + idxr.offsetCp(ParticleTypeIndex{type});
		StateType const* const yFluxType = yBase + idxr.offsetJf(ParticleTypeIndex{type});

		const ParamType epsP = static_cast<ParamType>(_parPorosity[type]);

		// Ordering of diffusion:
		// sec0type0comp0, sec0type0comp1, sec0type0comp2, sec0type1comp0, sec0type1comp1, sec0type1comp2,
		// sec1type0comp0, sec1type0comp1, sec1type0comp2, sec1type1comp0, sec1type1comp1, sec1type1comp2, ...
		active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp * _disc.nParType, secIdx) + type * _disc.nComp;
		active const* const parDiff = getSectionDependentSlice(_parDiffusion, _disc.nComp * _disc.nParType, secIdx) + type * _disc.nComp;

		const ParamType surfaceToVolumeRatio = 3.0 / static_cast<ParamType>(_parRadius[type]);
		const ParamType outerAreaPerVolume = static_cast<ParamType>(_parOuterSurfAreaPerVolume[_disc.nParCellsBeforeType[type]]);

		const ParamType jacPF_val = -outerAreaPerVolume / epsP;

		// Discretized film diffusion kf for finite volumes
		const ParamType absOuterShellHalfRadius = 0.5 * static_cast<ParamType>(_parCellSize[_disc.nParCellsBeforeType[type]]);
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			kf_FV[comp] = 1.0 / (absOuterShellHalfRadius / epsP / static_cast<ParamType>(_poreAccessFactor[type * _disc.nComp + comp]) / static_cast<ParamType>(parDiff[comp]) + 1.0 / static_cast<ParamType>(filmDiff[comp]));
		}

		// J_{0,f} block, adds flux to column void / bulk volume equations
		unsigned int idx = 0;
		for (unsigned int i = 0; i < _disc.nCol; ++i)
		{
			for (unsigned int j = 0; j < _disc.nRad; ++j)
			{
				const ParamType invBetaC = 1.0 / static_cast<ParamType>(_convDispOp.columnPorosity(j)) - 1.0;
				const ParamType jacCF_val = invBetaC * surfaceToVolumeRatio;
				for (unsigned int k = 0; k < _disc.nComp; ++k, ++idx)
				{
					resCol[idx] += jacCF_val * static_cast<ParamType>(_parTypeVolFrac[type + i * _disc.nParType * _disc.nRad + j * _disc.nParType]) * yFluxType[idx];
				}
			}
		}

		// J_{f,0} block, adds bulk volume state c_i to flux equation
		for (unsigned int bnd = 0; bnd < _disc.nCol * _disc.nRad; ++bnd)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int eq = bnd * idxr.strideColRadialCell() + comp * idxr.strideColComp();
				resFluxType[eq] -= kf_FV[comp] * yCol[eq];
			}
		}

		// J_{p,f} block, implements bead boundary condition in outer bead shell equation
		for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nRad; ++pblk)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int eq = pblk * idxr.strideColRadialCell() + comp * idxr.strideColComp();
				resParType[pblk * idxr.strideParBlock(type) + comp] += jacPF_val / static_cast<ParamType>(_poreAccessFactor[type * _disc.nComp + comp]) * yFluxType[eq];
			}
		}

		// J_{f,p} block, adds outer bead shell state c_{p,i} to flux equation
		for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nRad; ++pblk)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int eq = pblk * idxr.strideColRadialCell() + comp * idxr.strideColComp();
				resFluxType[eq] += kf_FV[comp] * yParType[comp + pblk * idxr.strideParBlock(type)];
			}
		}
	}

	_discParFlux.destroy<ParamType>();
	return 0;
}

/**
 * @brief Assembles off diagonal Jacobian blocks
 * @details Assembles the fixed blocks @f$ J_{0,f}, \dots, J_{N_p,f} @f$ and @f$ J_{f,0}, \dots, J_{f, N_p}. @f$
 *          The blocks are fixed for each section.
 * @param [in] t Current time
 * @param [in] secIdx Index of the current section
 */
void GeneralRateModel2D::assembleOffdiagJac(double t, unsigned int secIdx)
{
	// Clear matrices for new assembly
	_jacCF.clear();
	_jacFC.clear();
	for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nRad * _disc.nParType; ++pblk)
	{
		_jacPF[pblk].clear();
		_jacFP[pblk].clear();
	}

	// Note that the J_f block, which is the identity matrix, is treated in the linear solver

	Indexer idxr(_disc);

	// Discretized film diffusion kf for finite volumes
	double* const kf_FV = _discParFlux.create<double>(_disc.nComp);

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		const unsigned int typeOffset = type * _disc.nCol * _disc.nComp * _disc.nRad;
		const double epsP = static_cast<double>(_parPorosity[type]);

		// Ordering of diffusion:
		// sec0type0comp0, sec0type0comp1, sec0type0comp2, sec0type1comp0, sec0type1comp1, sec0type1comp2,
		// sec1type0comp0, sec1type0comp1, sec1type0comp2, sec1type1comp0, sec1type1comp1, sec1type1comp2, ...
		active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp * _disc.nParType, secIdx) + type * _disc.nComp;
		active const* const parDiff = getSectionDependentSlice(_parDiffusion, _disc.nComp * _disc.nParType, secIdx) + type * _disc.nComp;

		const double surfaceToVolumeRatio = 3.0 / static_cast<double>(_parRadius[type]);
		const double outerAreaPerVolume = static_cast<double>(_parOuterSurfAreaPerVolume[_disc.nParCellsBeforeType[type]]);

		const double jacPF_val = -outerAreaPerVolume / epsP;
		const double absOuterShellHalfRadius = 0.5 * static_cast<double>(_parCellSize[_disc.nParCellsBeforeType[type]]);

		// Discretized film diffusion kf for finite volumes
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			kf_FV[comp] = 1.0 / (absOuterShellHalfRadius / epsP / static_cast<double>(_poreAccessFactor[type * _disc.nComp + comp]) / static_cast<double>(parDiff[comp]) + 1.0 / static_cast<double>(filmDiff[comp]));
		}

		// J_{0,f} block, adds flux to column void / bulk volume equations
		unsigned int idx = 0;
		for (unsigned int col = 0; col < _disc.nCol; ++col)
		{
			for (unsigned int rad = 0; rad < _disc.nRad; ++rad)
			{
				const double invBetaC = 1.0 / static_cast<double>(_convDispOp.columnPorosity(rad)) - 1.0;
				const double jacCF_val = invBetaC * surfaceToVolumeRatio;

				for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++idx)
				{
					// Main diagonal corresponds to j_{f,i} (flux) state variable
					_jacCF.addElement(idx, idx + typeOffset, jacCF_val * static_cast<double>(_parTypeVolFrac[type + col * _disc.nRad * _disc.nParType + rad * _disc.nParType]));
				}
			}
		}

		// J_{f,0} block, adds bulk volume state c_i to flux equation
		for (unsigned int col = 0; col < _disc.nCol * _disc.nRad; ++col)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				// Main diagonal corresponds to c_i state variable in each column cell
				const unsigned int eq = col * idxr.strideColRadialCell() + comp * idxr.strideColComp();
				_jacFC.addElement(eq + typeOffset, eq, -kf_FV[comp]);
			}
		}

		// J_{p,f} block, implements bead boundary condition in outer bead shell equation
		linalg::DoubleSparseMatrix* const jacPFtype = _jacPF + type * _disc.nCol * _disc.nRad;
		for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nRad; ++pblk)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int eq = typeOffset + pblk * idxr.strideColRadialCell() + comp * idxr.strideColComp();
				jacPFtype[pblk].addElement(comp, eq, jacPF_val / static_cast<double>(_poreAccessFactor[comp]));
			}
		}

		// J_{f,p} block, adds outer bead shell state c_{p,i} to flux equation
		linalg::DoubleSparseMatrix* const jacFPtype = _jacFP + type * _disc.nCol * _disc.nRad;
		for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nRad; ++pblk)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int eq = typeOffset + pblk * idxr.strideColRadialCell() + comp * idxr.strideColComp();
				jacFPtype[pblk].addElement(eq, comp, kf_FV[comp]);
			}
		}
	}

	_discParFlux.destroy<double>();
}

int GeneralRateModel2D::residualSensFwdWithJacobian(const ActiveSimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode and at the same time update the 
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(simTime, simState, nullptr, adJac, true, true);
}

int GeneralRateModel2D::residualSensFwdAdOnly(const ActiveSimulationTime& simTime, const ConstSimulationState& simState, active* const adRes)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simTime.timeFactor, simState.vecStateY, simState.vecStateYdot, adRes); 
}

int GeneralRateModel2D::residualSensFwdCombine(const ActiveSimulationTime& simTime, const ConstSimulationState& simState, 
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes, 
	double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_SCOPE(_timerResidualSens);

	// tmp1 stores result of (dF / dy) * s
	// tmp2 stores result of (dF / dyDot) * sDot

	const SimulationTime cst{static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor)};
	const ConstSimulationState css{nullptr, nullptr};
	for (unsigned int param = 0; param < yS.size(); ++param)
	{

		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(cst, css, yS[param], 1.0, 0.0, tmp1);

		// Directional derivative (dF / dyDot) * sDot
		multiplyWithDerivativeJacobian(cst, css, ySdot[param], tmp2);

		double* const ptrResS = resS[param];

		BENCH_START(_timerResidualSensPar);

		// Complete sens residual is the sum:
		// TODO: Chunk TBB loop
#ifdef CADET_PARALLELIZE
		tbb::parallel_for(size_t(0), size_t(numDofs()), [&](size_t i)
#else
		for (unsigned int i = 0; i < numDofs(); ++i)
#endif
		{
			ptrResS[i] = tmp1[i] + tmp2[i] + adRes[i].getADValue(param);
		} CADET_PARFOR_END;

		BENCH_STOP(_timerResidualSensPar);
	}

	return 0;
}

/**
 * @brief Multiplies the given vector with the system Jacobian (i.e., @f$ \frac{\partial F}{\partial y}\left(t, y, \dot{y}\right) @f$)
 * @details Actually, the operation @f$ z = \alpha \frac{\partial F}{\partial y} x + \beta z @f$ is performed.
 * 
 *          Note that residual() or one of its cousins has to be called with the requested point @f$ (t, y, \dot{y}) @f$ once
 *          before calling multiplyWithJacobian() as this implementation ignores the given @f$ (t, y, \dot{y}) @f$.
 * @param [in] simTime Current simulation time point
 * @param [in] simState Simulation state vectors
 * @param [in] yS Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial y} @f$
 * @param [in] alpha Factor @f$ \alpha @f$ in front of @f$ \frac{\partial F}{\partial y} @f$
 * @param [in] beta Factor @f$ \beta @f$ in front of @f$ z @f$
 * @param [in,out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void GeneralRateModel2D::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Handle identity matrix of inlet DOFs
	for (unsigned int i = 0; i < _disc.nComp * _disc.nRad; ++i)
	{
		ret[i] = alpha * yS[i] + beta * ret[i];
	}

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), size_t(_disc.nCol * _disc.nRad * _disc.nParType + 1), [&](size_t idx)
#else
	for (unsigned int idx = 0; idx < _disc.nCol * _disc.nRad * _disc.nParType + 1; ++idx)
#endif
	{
		if (cadet_unlikely(idx == 0))
		{
			_convDispOp.jacobian().multiplyVector(yS + idxr.offsetC(), alpha, beta, ret + idxr.offsetC());
			_jacCF.multiplyVector(yS + idxr.offsetJf(), alpha, 1.0, ret + idxr.offsetC());
		}
		else
		{
			const unsigned int pblk = idx - 1;
			const unsigned int type = pblk / (_disc.nCol * _disc.nRad);
			const unsigned int par = pblk % (_disc.nCol * _disc.nRad);

			const int localOffset = idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{par});
			_jacP[pblk].multiplyVector(yS + localOffset, alpha, beta, ret + localOffset);
			_jacPF[pblk].multiplyVector(yS + idxr.offsetJf(), alpha, 1.0, ret + localOffset);
		}
	} CADET_PARFOR_END;

	// Handle flux equation

	// Set fluxes(ret) = fluxes(yS)
	// This applies the identity matrix in the bottom right corner of the Jaocbian (flux equation)
	for (unsigned int i = idxr.offsetJf(); i < numDofs(); ++i)
		ret[i] = alpha * yS[i] + beta * ret[i];

	double* const retJf = ret + idxr.offsetJf();
	_jacFC.multiplyVector(yS + idxr.offsetC(), alpha, 1.0, retJf);

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		for (unsigned int par = 0; par < _disc.nCol * _disc.nRad; ++par)
		{
			_jacFP[type * _disc.nCol * _disc.nRad + par].multiplyVector(yS + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{par}), alpha, 1.0, retJf);
		}
	}

	// Map inlet DOFs to the column inlet (first bulk cells)
	_jacInlet.multiplyAdd(yS, ret + idxr.offsetC(), alpha);
}

/**
 * @brief Multiplies the time derivative Jacobian @f$ \frac{\partial F}{\partial \dot{y}}\left(t, y, \dot{y}\right) @f$ with a given vector
 * @details The operation @f$ z = \frac{\partial F}{\partial \dot{y}} x @f$ is performed.
 *          The matrix-vector multiplication is transformed matrix-free (i.e., no matrix is explicitly formed).
 * @param [in] simTime Current simulation time point
 * @param [in] simState Simulation state vectors
 * @param [in] sDot Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void GeneralRateModel2D::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
{
	Indexer idxr(_disc);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), size_t(_disc.nCol * _disc.nRad * _disc.nParType + 1), [&](size_t idx)
#else
	for (unsigned int idx = 0; idx < _disc.nCol * _disc.nRad * _disc.nParType + 1; ++idx)
#endif
	{
		if (cadet_unlikely(idx == 0))
		{
			_convDispOp.multiplyWithDerivativeJacobian(simTime, sDot, ret);
		}
		else
		{
			const unsigned int idxParLoop = idx - 1;
			const unsigned int pblk = idxParLoop % (_disc.nCol * _disc.nRad);
			const unsigned int type = idxParLoop / (_disc.nCol * _disc.nRad);

			const double invBetaP = (1.0 / static_cast<double>(_parPorosity[type]) - 1.0) * simTime.timeFactor;
			unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;
			unsigned int const* const boundOffset = _disc.boundOffset + type * _disc.nComp;

			// Particle shells
			for (unsigned int shell = 0; shell < _disc.nParCell[type]; ++shell)
			{
				double const* const localSdot = sDot + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{pblk}) + shell * idxr.strideParShell(type);
				double* const localRet = ret + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{pblk}) + shell * idxr.strideParShell(type);

				// Mobile phase
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					// Add derivative with respect to dc_p / dt to Jacobian
					localRet[comp] = simTime.timeFactor * localSdot[comp];

					// Add derivative with respect to dq / dt to Jacobian (normal equations)
					for (unsigned int i = 0; i < nBound[comp]; ++i)
					{
						// Index explanation:
						//   nComp -> skip mobile phase
						//   + boundOffset[comp] skip bound states of all previous components
						//   + i go to current bound state
						localRet[comp] += invBetaP * localSdot[_disc.nComp + boundOffset[comp] + i];
					}
				}

				// Solid phase
				_binding[type]->multiplyWithDerivativeJacobian(localSdot + _disc.nComp, localRet + _disc.nComp, simTime.timeFactor);
			}
		}
	} CADET_PARFOR_END;

	// Handle fluxes (all algebraic)
	double* const dFdyDot = ret + idxr.offsetJf();
	std::fill(dFdyDot, dFdyDot + _disc.nCol * _disc.nComp * _disc.nRad * _disc.nParType, 0.0);

	// Handle inlet DOFs (all algebraic)
	std::fill_n(ret, _disc.nComp * _disc.nRad, 0.0);
}

void GeneralRateModel2D::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	for (IBindingModel* bm : _binding)
	{
		if (bm)
			bm->setExternalFunctions(extFuns, size);
	}
}

unsigned int GeneralRateModel2D::localOutletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	// Inlets are duplicated so need to be accounted for
	if (static_cast<double>(_convDispOp.currentVelocity(port)) >= 0.0)
		// Forward Flow: outlet is last cell
		return _disc.nComp * _disc.nRad + (_disc.nCol - 1) * _disc.nComp * _disc.nRad + port * _disc.nComp;
	else
		// Backward flow: Outlet is first cell
		return _disc.nComp * _disc.nRad + _disc.nComp * port;
}

unsigned int GeneralRateModel2D::localInletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	// Always 0 due to dedicated inlet DOFs
	return _disc.nComp * port;
}

unsigned int GeneralRateModel2D::localOutletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

unsigned int GeneralRateModel2D::localInletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

void GeneralRateModel2D::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// @todo Write this function
}

/**
 * @brief Computes equidistant radial nodes in the beads
 */
void GeneralRateModel2D::setEquidistantRadialDisc(unsigned int parType)
{
	const active radius = _parRadius[parType] - _parCoreRadius[parType];
	const active dr = radius / static_cast<double>(_disc.nParCell[parType]);
	std::fill(_parCellSize.data() + _disc.nParCellsBeforeType[parType], _parCellSize.data() + _disc.nParCellsBeforeType[parType] + _disc.nParCell[parType], dr);

	active* const ptrCenterRadius = _parCenterRadius.data() + _disc.nParCellsBeforeType[parType];
	active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data() + _disc.nParCellsBeforeType[parType];
	active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data() + _disc.nParCellsBeforeType[parType];

	for (unsigned int cell = 0; cell < _disc.nParCell[parType]; ++cell)
	{
		const active r_out = radius - static_cast<active>(cell) * dr;
		const active r_in = radius - static_cast<active>(cell + 1) * dr;

		ptrCenterRadius[cell] = radius - (0.5 + static_cast<active>(cell)) * dr;

		// Compute denominator -> corresponding to cell volume
		const active vol = pow(r_out, 3.0) - pow(r_in, 3.0);

		ptrOuterSurfAreaPerVolume[cell] = 3.0 * sqr(r_out) / vol;
		ptrInnerSurfAreaPerVolume[cell] = 3.0 * sqr(r_in) / vol;
	}
}

/**
 * @brief Computes the radial nodes in the beads in such a way that all shells have the same volume
 */
void GeneralRateModel2D::setEquivolumeRadialDisc(unsigned int parType)
{
	active r_out = _parRadius[parType];
	active r_in = _parCoreRadius[parType];
	const active volumePerShell = (pow(_parRadius[parType], 3.0) - pow(_parCoreRadius[parType], 3.0)) / static_cast<double>(_disc.nParCell[parType]);

	active* const ptrCellSize = _parCellSize.data() + _disc.nParCellsBeforeType[parType];
	active* const ptrCenterRadius = _parCenterRadius.data() + _disc.nParCellsBeforeType[parType];
	active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data() + _disc.nParCellsBeforeType[parType];
	active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data() + _disc.nParCellsBeforeType[parType];

	for (unsigned int cell = 0; cell < _disc.nParCell[parType]; ++cell)
	{
		if (cell != (_disc.nParCell[parType] - 1))
			r_in = pow(pow(r_out, 3.0) - volumePerShell, (1.0 / 3.0));
		else
			r_in = _parCoreRadius[parType];

		ptrCellSize[cell] = r_out - r_in;
		ptrCenterRadius[cell] = (r_out + r_in) * 0.5;

		ptrOuterSurfAreaPerVolume[cell] = 3.0 * sqr(r_out) / volumePerShell;
		ptrInnerSurfAreaPerVolume[cell] = 3.0 * sqr(r_in) / volumePerShell;

		// For the next cell: r_out == r_in of the current cell
		r_out = r_in;
	}
}

/**
 * @brief Computes all helper quantities for radial bead discretization from given radial cell boundaries
 * @details Calculates surface areas per volume for every shell and the radial shell centers.
 */
void GeneralRateModel2D::setUserdefinedRadialDisc(unsigned int parType)
{
	active* const ptrCellSize = _parCellSize.data() + _disc.nParCellsBeforeType[parType];
	active* const ptrCenterRadius = _parCenterRadius.data() + _disc.nParCellsBeforeType[parType];
	active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data() + _disc.nParCellsBeforeType[parType];
	active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data() + _disc.nParCellsBeforeType[parType];

	// Care for the right ordering and include 0.0 / 1.0 if not already in the vector.
	std::vector<double> orderedInterfaces = std::vector<double>(_parDiscVector.begin() + _disc.nParCellsBeforeType[parType] + parType, 
		_parDiscVector.begin() + _disc.nParCellsBeforeType[parType] + parType + _disc.nParCell[parType] + 1);

	// Sort in descending order
	std::sort(orderedInterfaces.begin(), orderedInterfaces.end(), std::greater<double>());

	// Force first and last element to be 1.0 and 0.0, respectively
	orderedInterfaces[0] = 1.0;
	orderedInterfaces.back() = 0.0;

	// Map [0, 1] -> [core radius, particle radius] via linear interpolation
	for (unsigned int cell = 0; cell < _disc.nParCell[parType]; ++cell)
		orderedInterfaces[cell] = orderedInterfaces[cell] * (static_cast<double>(_parRadius[parType]) - static_cast<double>(_parCoreRadius[parType])) + static_cast<double>(_parCoreRadius[parType]);

	for (unsigned int cell = 0; cell < _disc.nParCell[parType]; ++cell)
	{
		ptrCellSize[cell] = orderedInterfaces[cell] - orderedInterfaces[cell + 1];
		ptrCenterRadius[cell] = (orderedInterfaces[cell] + orderedInterfaces[cell + 1]) * 0.5;

		// Compute denominator -> corresponding to cell volume
		const active vol = std::pow(orderedInterfaces[cell], 3.0) - std::pow(orderedInterfaces[cell + 1], 3.0);

		ptrOuterSurfAreaPerVolume[cell] = 3.0 * sqr(orderedInterfaces[cell]) / vol;
		ptrInnerSurfAreaPerVolume[cell] = 3.0 * sqr(orderedInterfaces[cell + 1]) / vol;
	}
}

void GeneralRateModel2D::updateRadialDisc()
{
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		if (_parDiscType[i] == ParticleDiscretizationMode::Equidistant)
			setEquidistantRadialDisc(i);
		else if (_parDiscType[i] == ParticleDiscretizationMode::Equivolume)
			setEquivolumeRadialDisc(i);
		else if (_parDiscType[i] == ParticleDiscretizationMode::UserDefined)
			setUserdefinedRadialDisc(i);
	}
}

bool GeneralRateModel2D::setParameter(const ParameterId& pId, double value)
{
	if (multiplexParameterValue(pId, hashString("PAR_TYPE_VOLFRAC"), _parTypeVolFracMode, _parTypeVolFrac, _disc.nCol, _disc.nRad, _disc.nParType, value, nullptr))
		return true;
	if (multiplexParameterValue(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _disc.nParType, _disc.nComp, value, nullptr))
		return true;
	if (multiplexParameterValue(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _disc.nParType, _disc.nComp, value, nullptr))
		return true;
	if (multiplexParameterValue(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _disc.nParType, _disc.nComp, value, nullptr))
		return true;
	if (multiplexParameterValue(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _disc.nParType, _disc.nComp, _disc.strideBound, _disc.nBound, _disc.boundOffset, value, nullptr))
		return true;
	const int mpIc = multiplexInitialConditions(pId, value, false);
	if (mpIc > 0)
		return true;
	else if (mpIc < 0)
		return false;

	if (_singleParRadius && (pId.name == hashString("PAR_RADIUS")))
	{
		if ((pId.unitOperation != _unitOpIdx) || (pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
			return false;

		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parRadius[i].setValue(value);

		return true;
	}

	if (_singleParCoreRadius && (pId.name == hashString("PAR_CORERADIUS")))
	{
		if ((pId.unitOperation != _unitOpIdx) || (pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
			return false;

		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parCoreRadius[i].setValue(value);

		return true;
	}

	if (_singleParPorosity && (pId.name == hashString("PAR_POROSITY")))
	{
		if ((pId.unitOperation != _unitOpIdx) || (pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
			return false;

		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parPorosity[i].setValue(value);

		return true;
	}

	if (_convDispOp.setParameter(pId, value))
		return true;

	const bool result = UnitOperationBase::setParameter(pId, value);

	// Check whether particle radius or core radius has changed and update radial discretization if necessary
	if (result && ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS"))))
		updateRadialDisc();

	return result;
}

void GeneralRateModel2D::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if (multiplexParameterValue(pId, hashString("PAR_TYPE_VOLFRAC"), _parTypeVolFracMode, _parTypeVolFrac, _disc.nCol, _disc.nRad, _disc.nParType, value, &_sensParams))
		return;
	if (multiplexParameterValue(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _disc.nParType, _disc.nComp, value, &_sensParams))
		return;
	if (multiplexParameterValue(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _disc.nParType, _disc.nComp, value, &_sensParams))
		return;
	if (multiplexParameterValue(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _disc.nParType, _disc.nComp, value, &_sensParams))
		return;
	if (multiplexParameterValue(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _disc.nParType, _disc.nComp, _disc.strideBound, _disc.nBound, _disc.boundOffset, value, &_sensParams))
		return;
	if (multiplexInitialConditions(pId, value, true) != 0)
		return;

	if (_singleParRadius && (pId.name == hashString("PAR_RADIUS")))
	{
		if ((pId.unitOperation != _unitOpIdx) || (pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
			return;
		if (!contains(_sensParams, &_parRadius[0]))
			return;

		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parRadius[i].setValue(value);

		return;
	}

	if (_singleParCoreRadius && (pId.name == hashString("PAR_CORERADIUS")))
	{
		if ((pId.unitOperation != _unitOpIdx) || (pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
			return;
		if (!contains(_sensParams, &_parCoreRadius[0]))
			return;

		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parCoreRadius[i].setValue(value);

		return;
	}

	if (_singleParPorosity && (pId.name == hashString("PAR_POROSITY")))
	{
		if ((pId.unitOperation != _unitOpIdx) || (pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
			return;
		if (!contains(_sensParams, &_parPorosity[0]))
			return;

		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parPorosity[i].setValue(value);

		return;
	}

	if (_convDispOp.setSensitiveParameterValue(_sensParams, pId, value))
		return;

	UnitOperationBase::setSensitiveParameterValue(pId, value);

	// Check whether particle radius or core radius has changed and update radial discretization if necessary
	if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
		updateRadialDisc();
}

bool GeneralRateModel2D::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if (multiplexParameterAD(pId, hashString("PAR_TYPE_VOLFRAC"), _parTypeVolFracMode, _parTypeVolFrac, _disc.nCol, _disc.nRad, _disc.nParType, adDirection, adValue, _sensParams))
	{
		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
		return true;
	}

	if (multiplexParameterAD(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _disc.nParType, _disc.nComp, adDirection, adValue, _sensParams))
	{
		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
		return true;
	}

	if (multiplexParameterAD(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _disc.nParType, _disc.nComp, adDirection, adValue, _sensParams))
	{
		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
		return true;
	}

	if (multiplexParameterAD(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _disc.nParType, _disc.nComp, adDirection, adValue, _sensParams))
	{
		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
		return true;
	}

	if (multiplexParameterAD(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _disc.nParType, _disc.nComp, _disc.strideBound, _disc.nBound, _disc.boundOffset, adDirection, adValue, _sensParams))
	{
		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
		return true;
	}

	const int mpIc = multiplexInitialConditions(pId, adDirection, adValue);
	if (mpIc > 0)
	{
		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
		return true;
	}
	else if (mpIc < 0)
		return false;

	if (_singleParRadius && (pId.name == hashString("PAR_RADIUS")))
	{
		if ((pId.unitOperation != _unitOpIdx) || (pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
			return false;

		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;

		// Register parameter and set AD seed / direction
		_sensParams.insert(&_parRadius[0]);
		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parRadius[i].setADValue(adDirection, adValue);

		return true;
	}

	if (_singleParCoreRadius && (pId.name == hashString("PAR_CORERADIUS")))
	{
		if ((pId.unitOperation != _unitOpIdx) || (pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
			return false;

		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;

		// Register parameter and set AD seed / direction
		_sensParams.insert(&_parCoreRadius[0]);
		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parCoreRadius[i].setADValue(adDirection, adValue);

		return true;
	}

	if (_singleParPorosity && (pId.name == hashString("PAR_POROSITY")))
	{
		if ((pId.unitOperation != _unitOpIdx) || (pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
			return false;

		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;

		// Register parameter and set AD seed / direction
		_sensParams.insert(&_parPorosity[0]);
		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parPorosity[i].setADValue(adDirection, adValue);

		return true;
	}

	if (_convDispOp.setSensitiveParameter(_sensParams, pId, adDirection, adValue))
	{
		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
		return true;
	}

	const bool result = UnitOperationBase::setSensitiveParameter(pId, adDirection, adValue);

	// Check whether particle radius or core radius has been set active and update radial discretization if necessary
	// Note that we need to recompute the radial discretization variables (_parCellSize, _parCenterRadius, _parOuterSurfAreaPerVolume, _parInnerSurfAreaPerVolume)
	// because their gradient has changed (although their nominal value has not changed).
	if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
		updateRadialDisc();

	return result;
}

void registerGeneralRateModel2D(std::unordered_map<std::string, std::function<IUnitOperation*(UnitOpIdx)>>& models)
{
	models[GeneralRateModel2D::identifier()] = [](UnitOpIdx uoId) { return new GeneralRateModel2D(uoId); };
	models["GRM2D"] = [](UnitOpIdx uoId) { return new GeneralRateModel2D(uoId); };
}

}  // namespace model

}  // namespace cadet

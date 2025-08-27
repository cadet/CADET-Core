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

#include "model/ParameterMultiplexing.hpp"
#include "cadet/ParameterProvider.hpp"
#include "ParamIdUtil.hpp"
#include "SensParamUtil.hpp"
#include "AutoDiff.hpp"
#include "ParamReaderHelper.hpp"

namespace cadet
{

namespace model
{

bool setTypeParameterValue(const ParameterId& pId, StringHash nameHash, bool typeDep, active& parameter, unsigned int parTypeIdx, double value, std::unordered_set<active*> const* sensParams)
{
	if (!typeDep || (pId.name != nameHash) || (pId.particleType != ParTypeIndep && pId.particleType != parTypeIdx))
		return false;

	if ((pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
		return false;
	if (sensParams && !contains(*sensParams, &parameter))
		return false;

	parameter.setValue(value);

	return true;
}

bool setTypeParameterAD(const ParameterId& pId, StringHash nameHash, bool typeDep, active& parameter, unsigned int parTypeIdx, unsigned int adDirection, double adValue, std::unordered_set<active*>& sensParams)
{
	if (!typeDep || (pId.name != nameHash) || (pId.particleType != ParTypeIndep && pId.particleType != parTypeIdx))
		return false;

	if ((pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
		return false;

	sensParams.insert(&parameter);
	parameter.setADValue(adDirection, adValue);

	return true;
}


MultiplexMode readAndRegisterMultiplexCompSecParam(IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, std::vector<active>& values, const std::string& name, unsigned int nComp, unsigned int parTypeIdx, const bool parTypeDep, UnitOpIdx uoi, const int valTypeOffset)
{
	MultiplexMode mode = MultiplexMode::Independent;
	std::vector<active> tmpVals;
	readScalarParameterOrArray(tmpVals, paramProvider, name, 1);
	if (paramProvider.exists(name + "_MULTIPLEX"))
	{
		const int modeConfig = paramProvider.getInt(name + "_MULTIPLEX");
		if (modeConfig == 0)
		{
			mode = parTypeDep ? MultiplexMode::ComponentType : MultiplexMode::Component;
			if (tmpVals.size() != nComp)
				throw InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nComp) + ")");
		}
		else if (modeConfig == 1)
		{
			mode = parTypeDep ? MultiplexMode::ComponentSectionType : MultiplexMode::ComponentSection;
			if ((tmpVals.size() % nComp) != 0)
				throw InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be positive multiple of " + std::to_string(nComp) + ")");
		}
		else
			throw InvalidParameterException("Unsupported multiplex mode of field " + name + "_MULTIPLEX, must be 0 or 1");
	}
	else
	{
		if (tmpVals.size() == nComp)
			mode = parTypeDep ? MultiplexMode::ComponentType : MultiplexMode::Component;
		else if ((tmpVals.size() % nComp) == 0)
			mode = parTypeDep ? MultiplexMode::ComponentSectionType : MultiplexMode::ComponentSection;
		else
			throw InvalidParameterException("Could not infer multiplex mode of field " + name + ", set " + name + "_MULTIPLEX or change number of elements");
	}

	const StringHash nameHash = hashStringRuntime(name);
	switch (mode)
	{
	case MultiplexMode::Component:
	case MultiplexMode::ComponentType:
	{
		values.resize(valTypeOffset + tmpVals.size());
		std::copy_n(tmpVals.begin(), tmpVals.size(), values.begin() + valTypeOffset);
		if (parTypeDep || parTypeIdx == 0) // only register once for parTypeIndep
		{
			for (unsigned int s = 0; s < nComp; ++s)
				parameters[makeParamId(nameHash, uoi, s, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &values[valTypeOffset + s];
		}
	}
	break;
	case MultiplexMode::ComponentSection:
	case MultiplexMode::ComponentSectionType:
	{
		values.resize(valTypeOffset + tmpVals.size());
		std::move(tmpVals.begin(), tmpVals.end(), values.begin() + valTypeOffset);

		if (parTypeDep || parTypeIdx == 0) // only register once for parTypeIndep
		{
			for (std::size_t s = 0; s < tmpVals.size() / nComp; ++s)
			{
				for (unsigned int i = 0; i < nComp; ++i)
					parameters[makeParamId(nameHash, uoi, i, ParTypeIndep, BoundStateIndep, ReactionIndep, s)] = &values[valTypeOffset + s * nComp + i];
			}
		}
	}
	break;
	case MultiplexMode::RadialSection:
	case MultiplexMode::Independent:
	case MultiplexMode::ComponentRadial:
	case MultiplexMode::ComponentRadialSection:
	case MultiplexMode::Axial:
	case MultiplexMode::Section:
	case MultiplexMode::Type:
	case MultiplexMode::Radial:
	case MultiplexMode::AxialRadial:
		cadet_assert(false);
		break;
	}

	return mode;
}

MultiplexMode readAndRegisterMultiplexBndCompSecParam(IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, std::vector<active>& values, const std::string& name,
	unsigned int nComp, unsigned int strideBound, unsigned int const* nBound, const unsigned int parTypeIdx, const bool parTypeDep, UnitOpIdx uoi, const int valTypeOffset)
{
	if (strideBound == 0)
		return MultiplexMode::Component;

	MultiplexMode mode = MultiplexMode::Independent;
	std::vector<active> tmpVals;
	readScalarParameterOrArray(tmpVals, paramProvider, name, 1);
	if (paramProvider.exists(name + "_MULTIPLEX"))
	{
		const int modeConfig = paramProvider.getInt(name + "_MULTIPLEX");
		if (modeConfig == 0)
		{
			mode = parTypeDep ? MultiplexMode::ComponentType : MultiplexMode::Component;
			if (tmpVals.size() != strideBound)
				throw InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(strideBound) + ")");
		}
		else if (modeConfig == 1)
		{
			mode = parTypeDep ? MultiplexMode::ComponentSectionType : MultiplexMode::ComponentSection;
			if ((tmpVals.size() % strideBound) != 0)
				throw InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be positive multiple of " + std::to_string(strideBound) + ")");
		}
		else
			throw InvalidParameterException("Unsupported multiplex mode of field " + name + "_MULTIPLEX, must be 0 or 1");
	}
	else
	{
		if (tmpVals.size() == strideBound)
			mode = parTypeDep ? MultiplexMode::ComponentType : MultiplexMode::Component;
		else if ((tmpVals.size() % strideBound) == 0)
			mode = parTypeDep ? MultiplexMode::ComponentSectionType : MultiplexMode::ComponentSection;
		else
			throw InvalidParameterException("Could not infer multiplex mode of field " + name + ", set " + name + "_MULTIPLEX or change number of elements");
	}

	const StringHash nameHash = hashStringRuntime(name);
	switch (mode)
	{
	case MultiplexMode::Component:
	case MultiplexMode::ComponentType:
	{
		values.resize(valTypeOffset + tmpVals.size());
		std::move(tmpVals.begin(), tmpVals.end(), values.begin() + valTypeOffset);

		if (parTypeDep || parTypeIdx == 0) // only register once for parTypeIndep
		{
			unsigned int idx = 0;
			for (unsigned int c = 0; c < nComp; ++c)
			{
				for (unsigned int s = 0; s < nBound[c]; ++s, ++idx)
					parameters[makeParamId(nameHash, uoi, c, ParTypeIndep, s, ReactionIndep, SectionIndep)] = &values[valTypeOffset + idx];
			}
		}
	}
	break;
	case MultiplexMode::ComponentSection:
	case MultiplexMode::ComponentSectionType:
	{
		values.resize(valTypeOffset + tmpVals.size());
		std::move(tmpVals.begin(), tmpVals.end(), values.begin() + valTypeOffset);

		if (parTypeDep || parTypeIdx == 0) // only register once for parTypeIndep
		{
			for (std::size_t sec = 0; sec < values.size() / strideBound; ++sec)
			{
				unsigned int idx = strideBound * sec;
				for (unsigned int c = 0; c < nComp; ++c)
				{
					for (unsigned int s = 0; s < nBound[c]; ++s, ++idx)
						parameters[makeParamId(nameHash, uoi, c, ParTypeIndep, s, ReactionIndep, sec)] = &values[valTypeOffset + idx];
				}
			}
		}
	}
	break;
	case MultiplexMode::RadialSection:
	case MultiplexMode::Independent:
	case MultiplexMode::ComponentRadial:
	case MultiplexMode::ComponentRadialSection:
	case MultiplexMode::Axial:
	case MultiplexMode::Section:
	case MultiplexMode::Type:
	case MultiplexMode::Radial:
	case MultiplexMode::AxialRadial:
		cadet_assert(false);
		break;
	}

	return mode;
}

bool singleTypeMultiplexTypeParameterValue(const ParameterId& pId, StringHash nameHash, bool mode, active& data, unsigned int parTypeIdx, double value, std::unordered_set<active*> const* sensParams)
{
	if (!mode || (pId.name != nameHash) || (pId.particleType != ParTypeIndep && pId.particleType != parTypeIdx))
		return false;

	if ((pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
		return false;
	if (sensParams && !contains(*sensParams, &data))
		return false;

	data.setValue(value);

	return true;
}

bool singleTypeMultiplexTypeParameterAD(const ParameterId& pId, StringHash nameHash, bool mode, active& data, unsigned int parTypeIdx, unsigned int adDirection, double adValue, std::unordered_set<active*>& sensParams)
{
	if (!mode || (pId.name != nameHash) || (pId.particleType != ParTypeIndep && pId.particleType != parTypeIdx))
		return false;

	if ((pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep) || (pId.particleType != ParTypeIndep))
		return false;

	sensParams.insert(&data);
	data.setADValue(adDirection, adValue);

	return true;
}

bool singleTypeMultiplexCompTypeSecParameterValue(const ParameterId& pId, StringHash nameHash, MultiplexMode mode, std::vector<active>& data,
	unsigned int nComp, unsigned int parTypeIdx, double value, std::unordered_set<active*> const* sensParams, const int valTypeOffset)
{
	if (pId.name != nameHash || (pId.particleType != ParTypeIndep && pId.particleType != parTypeIdx))
		return false;

	switch (mode)
	{
	case MultiplexMode::Component:
	{
		if ((pId.component == CompIndep) || (pId.particleType != ParTypeIndep) || (pId.boundState != BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section != SectionIndep))
			return false;

		if (sensParams && !contains(*sensParams, &data[pId.component]))
			if (parTypeIdx == 0)
				return false; // sensParams is only expected to contain this address for the first particle type since the parameter is parTypeIndep

		data[pId.component].setValue(value);

		return true;
	}
	case MultiplexMode::ComponentSection:
	{
		if ((pId.component == CompIndep) || (pId.particleType != ParTypeIndep) || (pId.boundState != BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section == SectionIndep))
			return false;

		if (sensParams && !contains(*sensParams, &data[pId.section * nComp + pId.component]))
			if (parTypeIdx == 0)
				return false; // sensParams is only expected to contain this address for the first particle type since the parameter is parTypeIndep

		data[pId.section * nComp + pId.component].setValue(value);

		return true;
	}
	case MultiplexMode::ComponentType:
	{
		if ((pId.component == CompIndep) || (pId.particleType == ParTypeIndep) || (pId.boundState != BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section != SectionIndep))
			return false;

		if (sensParams && !contains(*sensParams, &data[pId.component]))
			return false;

		for (unsigned int s = 0; s < data.size() / nComp; ++s)
			data[s * nComp + pId.component].setValue(value);

		return true;
	}
	case MultiplexMode::ComponentSectionType:
	{
		if ((pId.component == CompIndep) || (pId.particleType == ParTypeIndep) || (pId.boundState != BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section == SectionIndep))
			return false;

		if (sensParams && !contains(*sensParams, &data[pId.section * nComp + pId.component]))
			if (parTypeIdx == 0)
				return false; // sensParams is only expected to contain this address for the first particle type since the parameter is parTypeIndep

		data[pId.section * nComp + pId.component].setValue(value);

		return true;
	}
	case MultiplexMode::RadialSection:
	case MultiplexMode::Independent:
	case MultiplexMode::ComponentRadial:
	case MultiplexMode::ComponentRadialSection:
	case MultiplexMode::Axial:
	case MultiplexMode::Section:
	case MultiplexMode::Type:
	case MultiplexMode::Radial:
	case MultiplexMode::AxialRadial:
		cadet_assert(false);
		break;
	}

	return false;
}

bool singleTypeMultiplexCompTypeSecParameterAD(const ParameterId& pId, StringHash nameHash, MultiplexMode mode, std::vector<active>& data,
	unsigned int nComp, unsigned int parTypeIdx, unsigned int adDirection, double adValue, std::unordered_set<active*>& sensParams, const int valTypeOffset)
{
	if (pId.name != nameHash || (pId.particleType != ParTypeIndep && pId.particleType != parTypeIdx))
		return false;

	switch (mode)
	{
	case MultiplexMode::Component:
	{
		if ((pId.component == CompIndep) || (pId.particleType != ParTypeIndep) || (pId.boundState != BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section != SectionIndep))
			return false;

		if (parTypeIdx == 0)
			sensParams.insert(&data[valTypeOffset + pId.component]);

		data[valTypeOffset + pId.component].setADValue(adDirection, adValue);

		return true;
	}
	case MultiplexMode::ComponentSection:
	{
		if ((pId.component == CompIndep) || (pId.particleType != ParTypeIndep) || (pId.boundState != BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section == SectionIndep))
			return false;

		if (parTypeIdx == 0)
			sensParams.insert(&data[valTypeOffset + pId.section * nComp + pId.component]);

		data[valTypeOffset + pId.section * nComp + pId.component].setADValue(adDirection, adValue);

		return true;
	}
	case MultiplexMode::ComponentType:
	{
		if ((pId.component == CompIndep) || (pId.particleType == ParTypeIndep) || (pId.boundState != BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section != SectionIndep))
			return false;

		sensParams.insert(&data[valTypeOffset + pId.component]);

		const int nSec = valTypeOffset > 0 ? valTypeOffset / parTypeIdx / nComp : data.size() / nComp;

		for (unsigned int s = 0; s < nSec; ++s)
			data[valTypeOffset + s * nComp + pId.component].setADValue(adDirection, adValue);

		return true;
	}
	case MultiplexMode::ComponentSectionType:
	{
		if ((pId.component == CompIndep) || (pId.particleType == ParTypeIndep) || (pId.boundState != BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section == SectionIndep))
			return false;

		sensParams.insert(&data[valTypeOffset + pId.section * nComp + pId.component]);

		data[valTypeOffset + pId.section * nComp + pId.component].setADValue(adDirection, adValue);

		return true;
	}
	case MultiplexMode::RadialSection:
	case MultiplexMode::Independent:
	case MultiplexMode::ComponentRadial:
	case MultiplexMode::ComponentRadialSection:
	case MultiplexMode::Axial:
	case MultiplexMode::Section:
	case MultiplexMode::Type:
	case MultiplexMode::Radial:
	case MultiplexMode::AxialRadial:
		cadet_assert(false);
		break;
	}

	return false;
}

bool singleTypeMultiplexBndCompTypeSecParameterValue(const ParameterId& pId, StringHash nameHash, MultiplexMode mode, std::vector<active>& data,
	unsigned int nComp, unsigned int strideBound, unsigned int const* boundOffset, unsigned int parTypeIdx, double value, std::unordered_set<active*> const* sensParams, const int valTypeOffset)
{
	if (pId.name != nameHash || (pId.particleType != ParTypeIndep && pId.particleType != parTypeIdx))
		return false;

	if (strideBound == 0)
		return true;

	switch (mode)
	{
	case MultiplexMode::Component:
	{
		if ((pId.component == CompIndep) || (pId.particleType != ParTypeIndep) || (pId.boundState == BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section != SectionIndep))
			return false;

		if (sensParams && !contains(*sensParams, &data[valTypeOffset + boundOffset[pId.component] + pId.boundState]))
			if (parTypeIdx == 0)
				return false;

		data[valTypeOffset + boundOffset[pId.component] + pId.boundState].setValue(value);

		return true;
	}
	case MultiplexMode::ComponentSection:
	{
		if ((pId.component == CompIndep) || (pId.particleType != ParTypeIndep) || (pId.boundState == BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section == SectionIndep))
			return false;

		if (sensParams && !contains(*sensParams, &data[valTypeOffset + pId.section * strideBound + boundOffset[pId.component] + pId.boundState]))
			if (parTypeIdx == 0)
				return false;

		data[valTypeOffset + pId.section * strideBound + boundOffset[pId.component] + pId.boundState].setValue(value);

		return true;
	}
	case MultiplexMode::ComponentType:
	{
		if ((pId.component == CompIndep) || (pId.particleType != ParTypeIndep) || (pId.boundState == BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section != SectionIndep))
			return false;

		if (sensParams && !contains(*sensParams, &data[valTypeOffset + boundOffset[pId.component] + pId.boundState]))
			if (parTypeIdx == 0)
				return false;

		const int nSec = valTypeOffset > 0 ? valTypeOffset / parTypeIdx / nComp : data.size() / nComp;

		for (unsigned int s = 0; s < nSec; ++s)
			data[valTypeOffset + s * strideBound + boundOffset[pId.component] + pId.boundState].setValue(value);

		return true;
	}
	case MultiplexMode::ComponentSectionType:
	{
		if ((pId.component == CompIndep) || (pId.particleType != ParTypeIndep) || (pId.boundState == BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section == SectionIndep))
			return false;

		if (sensParams && !contains(*sensParams, &data[valTypeOffset + pId.section * strideBound + boundOffset[pId.component] + pId.boundState]))
			if (parTypeIdx == 0)
				return false;

		data[valTypeOffset + pId.section * strideBound + boundOffset[pId.component] + pId.boundState].setValue(value);

		return true;
	}
	case MultiplexMode::RadialSection:
	case MultiplexMode::Independent:
	case MultiplexMode::ComponentRadial:
	case MultiplexMode::ComponentRadialSection:
	case MultiplexMode::Axial:
	case MultiplexMode::Section:
	case MultiplexMode::Type:
	case MultiplexMode::Radial:
	case MultiplexMode::AxialRadial:
		cadet_assert(false);
		break;
	}

	return false;
}

bool singleTypeMultiplexBndCompTypeSecParameterAD(const ParameterId& pId, StringHash nameHash, MultiplexMode mode, std::vector<active>& data,
	unsigned int nComp, unsigned int strideBound, unsigned int const* boundOffset, unsigned int parTypeIdx, unsigned int adDirection, double adValue, std::unordered_set<active*>& sensParams, const int valTypeOffset)
{
	if (pId.name != nameHash || (pId.particleType != ParTypeIndep && pId.particleType != parTypeIdx))
		return false;

	if (strideBound == 0)
		return true;

	switch (mode)
	{
	case MultiplexMode::Component:
	{
		if ((pId.component == CompIndep) || (pId.particleType != ParTypeIndep) || (pId.boundState == BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section != SectionIndep))
			return false;

		if (parTypeIdx == 0)
			sensParams.insert(&data[valTypeOffset + boundOffset[pId.component] + pId.boundState]);

		data[valTypeOffset + boundOffset[pId.component] + pId.boundState].setADValue(adDirection, adValue);

		return true;
	}
	case MultiplexMode::ComponentSection:
	{
		if ((pId.component == CompIndep) || (pId.particleType != ParTypeIndep) || (pId.boundState == BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section == SectionIndep))
			return false;

		if (parTypeIdx == 0)
			sensParams.insert(&data[valTypeOffset + pId.section * strideBound + boundOffset[pId.component] + pId.boundState]);

		data[valTypeOffset + pId.section * strideBound + boundOffset[pId.component] + pId.boundState].setADValue(adDirection, adValue);

		return true;
	}
	case MultiplexMode::ComponentType:
	{
		if ((pId.component == CompIndep) || (pId.particleType == ParTypeIndep) || (pId.boundState == BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section != SectionIndep))
			return false;

		sensParams.insert(&data[valTypeOffset + boundOffset[pId.component] + pId.boundState]);

		const int nSec = valTypeOffset > 0 ? valTypeOffset / parTypeIdx / nComp : data.size() / nComp;

		for (unsigned int s = 0; s < nSec; ++s)
			data[valTypeOffset + s * strideBound + boundOffset[pId.component] + pId.boundState].setADValue(adDirection, adValue);

		return true;
	}
	case MultiplexMode::ComponentSectionType:
	{
		if ((pId.component == CompIndep) || (pId.particleType == ParTypeIndep) || (pId.boundState == BoundStateIndep)
			|| (pId.reaction != ReactionIndep) || (pId.section == SectionIndep))
			return false;

		sensParams.insert(&data[valTypeOffset + pId.section * strideBound + boundOffset[pId.component] + pId.boundState]);

		data[valTypeOffset + pId.section * strideBound + boundOffset[pId.component] + pId.boundState].setADValue(adDirection, adValue);

		return true;
	}
	case MultiplexMode::RadialSection:
	case MultiplexMode::Independent:
	case MultiplexMode::ComponentRadial:
	case MultiplexMode::ComponentRadialSection:
	case MultiplexMode::Axial:
	case MultiplexMode::Section:
	case MultiplexMode::Type:
	case MultiplexMode::Radial:
	case MultiplexMode::AxialRadial:
		cadet_assert(false);
		break;
	}

	return false;
}

}  // namespace model

}  // namespace cadet

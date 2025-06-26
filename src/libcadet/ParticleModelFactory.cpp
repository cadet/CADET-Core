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

#include "ParticleModelFactory.hpp"
#include "cadet/Exceptions.hpp"

namespace cadet
{
	namespace model
	{
			void registerGeneralRateParticleModel(std::unordered_map<std::string, std::function<model::IParticleModel*()>>& particles);
			void registerHomogeneousParticleModel(std::unordered_map<std::string, std::function<model::IParticleModel*()>>& particles);
	}

	ParticleModelFactory::ParticleModelFactory()
	{
		// Register all ParticleModels here
		model::registerGeneralRateParticleModel(_particleModels);
		model::registerHomogeneousParticleModel(_particleModels);
	}

	ParticleModelFactory::~ParticleModelFactory() { }

	template <class ParticleModel_t>
	void ParticleModelFactory::registerModel(const std::string& name)
	{
		_particleModels[name] = []() { return new ParticleModel_t(); };
	}

	template <class ParticleModel_t>
	void ParticleModelFactory::registerModel()
	{
		registerModel<ParticleModel_t>(ParticleModel_t::identifier());
	}

	model::IParticleModel* ParticleModelFactory::create(const std::string& name) const
	{
		const auto it = _particleModels.find(name);
		if (it == _particleModels.end())
		{
			// ParticleModel was not found
			return nullptr;
		}

		// Call factory function (thanks to type erasure of std::function we can store 
		// all factory functions in one container)
		return it->second();
	}

	void ParticleModelFactory::registerModel(const std::string& name, std::function<model::IParticleModel* ()> factory)
	{
		if (_particleModels.find(name) == _particleModels.end())
			_particleModels[name] = std::move(factory);
		else
			throw InvalidParameterException("IParticleModel implementation with the name " + name + " is already registered and cannot be overwritten");
	}

	bool ParticleModelFactory::exists(const std::string& name) const
	{
		return _particleModels.find(name) != _particleModels.end();
	}
} // namespace cadet

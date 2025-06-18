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
 * Defines the ParticleModelFactory
 */

#ifndef LIBCADET_PARTICLEMODELFACTORY_HPP_
#define LIBCADET_PARTICLEMODELFACTORY_HPP_

#include <string>
#include <unordered_map>
#include <functional>

namespace cadet
{

	namespace model
	{
		class IParticleModel;
	}

	/**
	 * @brief Creates particle models
	 */
	class ParticleModelFactory
	{
	public:
		/**
		 * @brief Construct the ParticleModelFactory
		 * @details All internal particle models are registered here.
		 */
		ParticleModelFactory();

		~ParticleModelFactory();

		/**
		 * @brief Creates particle models with the given @p name
		 * @param [in] name Name of the particle model
		 * @return The particle model or @c NULL if a particle model with this name does not exist
		 */
		model::IParticleModel* create(const std::string& name) const;

		/**
		 * @brief Registers the given particle model implementation
		 * @param [in] name Name of the IParticleModel implementation
		 * @param [in] factory Function that creates an object of the IParticleModel class
		 */
		void registerModel(const std::string& name, std::function<model::IParticleModel*()> factory);

		/**
		 * @brief Returns whether a particle model of the given name @p name exists
		 * @param [in] name Name of the particle model
		 * @return @c true if a particle model of this name exists, otherwise @c false
		 */
		bool exists(const std::string& name) const;
	protected:

		/**
		 * @brief Registers an IParticleModel
		 * @param [in] name Name of the particle model
		 * @tparam ParticleModel_t Type of the particle model
		 */
		template <class ParticleModel_t>
		void registerModel(const std::string& name);

		/**
		 * @brief Registers an IParticleModel
		 * @details The name of the particle model is inferred from the static function IParticleModel::identifier().
		 * @tparam ParticleModel_t Type of the particle model
		 */
		template <class ParticleModel_t>
		void registerModel();

		std::unordered_map<std::string, std::function<model::IParticleModel* ()>> _particleModels; //!< Map with factory functions
	};

} // namespace cadet

#endif  // LIBCADET_ParticleModelFactory_HPP_

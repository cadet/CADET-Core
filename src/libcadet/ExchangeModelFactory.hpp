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
 * Defines the ExchangeModelFactory
 */

#ifndef LIBCADET_EXCHANGEMODELFACTORY_HPP_
#define LIBCADET_EXCHANGEMODELFACTORY_HPP_

#include <string>
#include <unordered_map>
#include <functional>

namespace cadet
{

	namespace model
	{
		class IExchangeModel;
	}

	/**
	 * @brief Creates Exchange models
	 */
	class ExchangeModelFactory
	{
	public:
		/**
		 * @brief Construct the ExchangeModelFactory
		 * @details All internal exchange models are registered here.
		 */
		ExchangeModelFactory();

		~ExchangeModelFactory();

		/**
		 * @brief Creates exchange models with the given @p name
		 * @param [in] name Name of the exchange model
		 * @return The exchange model or @c NULL if a exchange model with this name does not exist
		 */
		model::IExchangeModel* create(const std::string& name) const;

		/**
		 * @brief Registers the given exchange model implementation
		 * @param [in] name Name of the IExchangeModel implementation
		 * @param [in] factory Function that creates an object of the IExchangeModel class
		 */
		void registerModel(const std::string& name, std::function<model::IExchangeModel*()> factory);

		/**
		 * @brief Returns whether a exchange model of the given name @p name exists
		 * @param [in] name Name of the exchange model
		 * @return @c true if a exchange model of this name exists, otherwise @c false
		 */
		bool exists(const std::string& name) const;
	protected:

		/**
		 * @brief Registers an IExchangeModel
		 * @param [in] name Name of the exchange model
		 * @tparam ExchangeModel_t Type of the exchange model
		 */
		template <class ExchangeModel_t>
		void registerModel(const std::string& name);

		/**
		 * @brief Registers an IExchangeModel
		 * @details The name of the exchange model is inferred from the static function IExchangeModel::identifier().
		 * @tparam ExchangeModel_t Type of the exchange model
		 */
		template <class ExchangeModel_t>
		void registerModel();

		std::unordered_map<std::string, std::function<model::IExchangeModel* ()>> _exchangeModels; //!< Map with factory functions
	};

} // namespace cadet

#endif  // LIBCADET_EXCHANGEMODELFACTORY_HPP_

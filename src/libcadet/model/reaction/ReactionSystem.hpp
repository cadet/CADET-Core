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

#pragma once

#include "model/reaction/ReactionModelBase.hpp"
#include "model/ReactionModel.hpp"
#include "cadet/Exceptions.hpp"
#include "ConfigurationHelper.hpp"
#include "SensParamUtil.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>
#include <map>
#include <vector>
#include <string>
#include <cstdio>
#include <limits>
#include <cmath>
#include <unordered_set>

namespace cadet
{

namespace model
{

/**
 * @brief Manages reaction models across different phases in a unit operation
 * 
 * The ReactionSystem organizes and manages dynamic reaction models for different
 * phases (cross_phase, solid, bulk, liquid) within a unit operation.
 * It handles configuration, discretization, and memory management of reaction models.
 */
struct ReactionSystem
	{   
        private:
        /**
         * @brief Data structure holding reaction information for a specific phase
         */
        struct PhaseData
        {
            std::vector<IDynamicReactionModel*> dynReactions; //!< Dynamic reactions in the phase

            /**
             * @brief Default constructor initializing phase data with safe defaults
             */
            PhaseData()
            {
                dynReactions = std::vector<IDynamicReactionModel*>{};
            }
        }; 

        //!< Maps phase type strings to their corresponding phase data
        std::map<std::string, PhaseData> _phaseMap = {
            {"cross_phase", PhaseData{}}, 
            {"solid", PhaseData{}},
            {"liquid", PhaseData{}}, // unit liquid phase
        };

        /**
         * @brief Retrieves phase data for a given phase type
         * @param phaseType The type of phase to retrieve data for
         * @return Reference to the phase data
         * @throws InvalidParameterException if phase type is unknown
         */
        PhaseData& getPhaseData(const std::string& phaseType)
        {
            auto phase = _phaseMap.find(phaseType);
            if (phase == _phaseMap.end())
                throw InvalidParameterException("Unknown phase type: " + phaseType);

            return phase->second;
        }

        /**
         * @brief Const version of getPhaseData
         * @param phaseType The type of phase to retrieve data for
         * @return Const reference to the phase data
         * @throws InvalidParameterException if phase type is unknown
         */
        const PhaseData& getPhaseData(const std::string& phaseType) const
        {
            auto phase = _phaseMap.find(phaseType);
            if (phase == _phaseMap.end())
                throw InvalidParameterException("Unknown phase type: " + phaseType);

            return phase->second;
        }

        public:

        /**
         * @brief Checks if any reactions are present in the system
         * @return True if at least one reaction exists, false otherwise
         */
        bool hasReactions()
        {
            for(const auto& phasePair : _phaseMap)
            {
                const auto& dynReactions = phasePair.second.dynReactions;
                if (!dynReactions.empty())
                {
                    for (const auto* reaction : dynReactions)
                    {
                        if (reaction != nullptr)
                        {
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        /**
         * @brief Adds a dynamic reaction model to a specific phase
         * @param phaseType The phase type to add the reaction to
         * @param reactionModel Pointer to the dynamic reaction model to add
         */
        void addReactionModel(const std::string& phaseType, IDynamicReactionModel* reactionModel)
        {
            auto& dynReactions = getDynReactionVector(phaseType);
            dynReactions.push_back(reactionModel);
        }
        
        /**
         * @brief Gets const reference to dynamic reaction vector for a phase
         * @param phase_type The phase type to retrieve reactions for
         * @return Const reference to the vector of dynamic reaction models
         */
        const std::vector<IDynamicReactionModel*>& getDynReactionVector(const std::string& phase_type) const
        {
            return getPhaseData(phase_type).dynReactions;
        }

        /**
         * @brief Gets mutable reference to dynamic reaction vector for a phase
         * @param phase_type The phase type to retrieve reactions for
         * @return Mutable reference to the vector of dynamic reaction models
         */
        std::vector<IDynamicReactionModel*>& getDynReactionVector(const std::string& phase_type) 
        {
            return getPhaseData(phase_type).dynReactions;
        }
                
        /**
         * @brief Configures the dimensions of the dynamic reaction vector for a phase
         * @param phaseType The phase type to configure
         * @param numReac Number of reactions to allocate space for
         */
        void configureDimensions(const std::string& phaseType, unsigned int numReac)
        {
            auto& reactionVector = getPhaseData(phaseType).dynReactions;
            reactionVector.resize(numReac, nullptr);
        }


        /**
         * @brief Configures discretization for reaction models in a specific phase
         * @param phaseType The phase type to configure
         * @param nReactions Number of reactions to configure
         * @param nComp Number of components
         * @param nBound Array of bound states per component
         * @param boundOffset Array of bound state offsets per component
         * @param paramProvider Parameter provider for configuration
         * @param helper Configuration helper for creating reaction models
         * @return True if configuration succeeded, false otherwise
         */
        bool configureDiscretization(std::string phaseType, unsigned int nReactions, unsigned int nComp, unsigned int* nBound, unsigned int* boundOffset,  IParameterProvider& paramProvider, const IConfigHelper& helper)
	    {
            auto& dynReaction = getDynReactionVector(phaseType);
            configureDimensions(phaseType, nReactions);

            bool reactionConfSuccess = true;
            // Determine the number of reactions to configure
            const unsigned int nReac = nReactions;

            // Configure each reaction model
            for (unsigned int i = 0; i < nReac; ++i)
            {                    
                char reactionKey[32];
                snprintf(reactionKey, sizeof(reactionKey), "%s_reaction_%03d", phaseType.c_str(), i);


                paramProvider.pushScope(reactionKey); //scope reaction_xxx

                // Check if reaction type is specified
                if (!paramProvider.exists("TYPE"))
                {
                    paramProvider.popScope(); //scope reaction_xxx
                    throw InvalidParameterException("Missing 'type' parameter for " + std::string(reactionKey));
                }

                // Create reaction model based on type
                std::string reactionType = paramProvider.getString("TYPE");
                dynReaction[i] = helper.createDynamicReactionModel(reactionType);

                // Validate reaction model creation
                if (!dynReaction[i]) 
                {
                    paramProvider.popScope(); //scope reaction_xxx
                    throw InvalidParameterException("Unknown dynamic reaction model " + reactionType +
                        " for " + reactionKey);
                }

                // Configure the reaction model discretization
                reactionConfSuccess = dynReaction[i]->configureModelDiscretization(paramProvider, nComp, nBound, boundOffset) && reactionConfSuccess;

                // Handle configuration failure
                if (!reactionConfSuccess) 
                {
                    paramProvider.popScope(); //scope reaction_xxx
                    throw InvalidParameterException("Failed to configure reaction model " + reactionType +
                        " for " + reactionKey);
                }

                paramProvider.popScope(); //scope reaction_xxx
            }
            
            return reactionConfSuccess;
	    }

        /**
         * @brief Configures reaction models for a specific phase and particle type
         * @param phaseType The phase type to configure
         * @param parType The particle type index
         * @param unitOpIdx Unit operation index
         * @param paramProvider Parameter provider for configuration
         * @return True if configuration succeeded, false otherwise
         */
        bool configure(std::string phaseType, unsigned int parType, unsigned int unitOpIdx, IParameterProvider& paramProvider)
        {
            auto& dynReaction = getDynReactionVector(phaseType);
            bool dynReactionConfSuccess = true;
             
            // Determine the number of reactions to configure
            const unsigned int nReac = dynReaction.size();
            // Configure each reaction
            for (unsigned int reac = 0; reac < nReac; ++reac)
            {
                // Skip null reactions or those that don't require configuration
                if (!dynReaction[reac] || !dynReaction[ reac]->requiresConfiguration())
                    continue;

                // Create reaction scope key
                char reactionKey[32];
                snprintf(reactionKey, sizeof(reactionKey), "%s_reaction_%03d", phaseType.c_str(), reac);
                paramProvider.pushScope(reactionKey); //scope reaction_xxx

                // Configure the reaction model
                dynReactionConfSuccess = dynReaction[reac]->configure(paramProvider, unitOpIdx, parType) && dynReactionConfSuccess;

                paramProvider.popScope();//scope reaction_xxx
            }
        
            return dynReactionConfSuccess;
        }

        /**
         * @brief Sets workspace requirements for all reaction models across all phases
         * @param phaseType The phase type to configure
         * @param nComp Number of components
         * @param strideBound Stride of bound states
         * @param lms Linear memory sizer to configure workspace requirements
         */
        void setWorkspaceRequirements(std::string phaseType, unsigned int nComp, unsigned int const strideBound, LinearMemorySizer& lms) const
        {
            auto& dynReactionVector = getDynReactionVector(phaseType);
            for (unsigned int i = 0; i < dynReactionVector.size(); i++)
            {
                if (dynReactionVector[i] && dynReactionVector[i]->requiresWorkspace())
                {
                    lms.fitBlock(dynReactionVector[i]->workspaceSize(nComp, strideBound, nullptr));
                }
            }
        }

        /**
         * @brief Clears and deletes all dynamic reaction models in all phases
         * Properly deallocates memory and resets vectors to safe state
         */
		void clearDynamicReactionModels()
		{
			for (auto& phase: _phaseMap)
			{
				auto& dynReactionVector = getDynReactionVector(phase.first);
                // Delete all reaction model pointers
				for (auto* reac: dynReactionVector) 
				{
					delete reac;
				}
                // Clear vector and reset to safe initial state
				dynReactionVector.clear();
            }
		}

        /**
         * @brief Empties all reaction vectors without deleting the models
         * Used for cleanup without memory deallocation
         */
        void empty()
        {
            for (auto& phase : _phaseMap)
            {
                auto& dynReacVec = getPhaseData(phase.first).dynReactions;
                dynReacVec.clear();
                dynReacVec.resize(1, nullptr);
            }
        }

        void getAllParameterValues(std::unordered_map<ParameterId, double>& data)
        {
            for (const auto& phasePair : _phaseMap)
            {
                const std::vector<IDynamicReactionModel*>& dynReactions = phasePair.second.dynReactions;
                for (const auto* reaction : dynReactions)
                {
                    if (reaction)
                    {
                        const std::unordered_map<ParameterId, double> localData = reaction->getAllParameterValues();
                        for (const auto& pair : localData)
                            data[pair.first] = pair.second;
                    }
                }
            }
        }

        double getParameterDouble(const ParameterId& pId) const
        {
            for (const auto& phasePair : _phaseMap)
            {
                const std::vector<IDynamicReactionModel*>& dynReactions = phasePair.second.dynReactions;
                for (const auto* reaction : dynReactions)
                {
                    if (reaction)
                    {
                        active const* const val = const_cast<IDynamicReactionModel*>(reaction)->getParameter(pId);
                        if (val)
                        {
                            return static_cast<double>(*val);
                        }
                    }
                }
            }
            
            // Not found - return NaN
            return std::numeric_limits<double>::quiet_NaN();
        }

        bool hasParameter(const ParameterId& pId) const
        {
            for (const auto& phasePair : _phaseMap)
            {
                const std::vector<IDynamicReactionModel*>& dynReactions = phasePair.second.dynReactions;
                for (const auto* reaction : dynReactions)
                {
                    if (reaction && const_cast<IDynamicReactionModel*>(reaction)->hasParameter(pId))
                        return true;
                }
            }
            return false;
        }

        template <typename param_t>
        bool setParameter(const ParameterId& pId, param_t value)
        {
            for (auto& phasePair : _phaseMap)
            {
                std::vector<IDynamicReactionModel*>& dynReactions = phasePair.second.dynReactions;
                for (auto* reaction : dynReactions)
                {
                    if (reaction && reaction->setParameter(pId, value))
                        return true;
                }
            }
            return false;
        }

        bool setSensitiveParameterValue(const ParameterId& pId, double value, const std::unordered_set<active*>& sensParams)
        {
            for (auto& phasePair : _phaseMap)
            {
                std::vector<IDynamicReactionModel*>& dynReactions = phasePair.second.dynReactions;
                for (auto* reaction : dynReactions)
                {
                    if (!reaction)
                        continue;

                    active* const val = reaction->getParameter(pId);
                    if (val && contains(sensParams, val))
                    {
                        val->setValue(value);
                        return true;
                    }
                }
            }
            return false;
        }

        bool setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue, std::unordered_set<active*>& sensParams)
        {
            for (auto& phasePair : _phaseMap)
            {
                std::vector<IDynamicReactionModel*>& dynReactions = phasePair.second.dynReactions;
                for (auto* reaction : dynReactions)
                {
                    if (!reaction)
                        continue;

                    active* const paramReaction = reaction->getParameter(pId);
                    if (paramReaction)
                    {
                        // Register parameter and set AD seed / direction
                        sensParams.insert(paramReaction);
                        paramReaction->setADValue(adDirection, adValue);
                        return true;
                    }
                }
            }
            return false;
        }

        /**
         * @brief Default constructor - initializes phase map with default values
         */
        ReactionSystem(){ }

        /**
         * @brief Creates and configures a new ReactionSystem for each particle type
         * @param reacParticle Vector to hold the created ReactionSystem pointers
         */
        static void create(std::vector<ReactionSystem*>& reacParticle)
        {
            for (unsigned int par = 0; par < reacParticle.size(); par++)
            {
                reacParticle[par] = new ReactionSystem();
            }
        }
        /**
         * @brief Destructor - cleans up dynamic reaction models
         */
        ~ReactionSystem()
        {
            clearDynamicReactionModels();
        }

	};



} // namespace model
} // namespace cadet
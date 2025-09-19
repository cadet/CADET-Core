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

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>
#include <map>
#include <vector>
#include <string>
#include <cstdio>
namespace cadet
{

namespace model
{

/**
 * @brief Manages reaction models across different phases in a unit operation
 * 
 * The ReactionSystem organizes and manages dynamic reaction models for different
 * phases (cross_phase, pore, solid, bulk, liquid) within a unit operation.
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
            std::vector<int> offsets; //!< Offsets for reactions in the phase
            std::vector<int> nReacParType; //!< Number of reactions per particle type
            unsigned int totalReactions; //!< Total reactions in the phase

            /**
             * @brief Default constructor initializing phase data with safe defaults
             */
            PhaseData()
            {
                dynReactions = { nullptr };
                offsets = { 0 };
                nReacParType = { 0 };
                totalReactions = 0;

            }
        }; 

        //!< Maps phase type strings to their corresponding phase data
        std::map<std::string, PhaseData> _phaseMap = {
            {"cross_phase", PhaseData{}}, 
            {"solid", PhaseData{}},
            {"liquid", PhaseData{}}, // unit liquid phase
            {"pore", PhaseData{}} // particle liquid phase
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

        unsigned int _parTypes; //!< Internal storage for the number of particles in the unit

        /**
         * @brief Gets the total number of reactions for a specific phase
         * @param phaseType The phase type to query
         * @return Total number of reactions in the phase
         */
        unsigned int getTotalReactions(const std::string& phaseType) const
        {
            return getPhaseData(phaseType).totalReactions;
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
         * @brief Updates and returns the total number of reactions for a phase
         * @param phaseType The phase type to update
         * @param nReactions Number of reactions to add to the total
         * @return Updated total number of reactions
         */
        unsigned int updateTotalReactions(const std::string& phaseType, unsigned int nReactions)
        {
            PhaseData& phaseData = getPhaseData(phaseType);
            phaseData.totalReactions += nReactions;
            return phaseData.totalReactions;
        }
        
        /**
         * @brief Configures dimensions for offset and reaction-per-particle-type arrays
         * @param nParTypes Number of particle types to allocate for
         */
        void configureDimOfSetAndReacParType(unsigned int nParTypes)
        {
            unsigned int dim = std::max(1u, nParTypes);
            
            for (auto& [phaseType, phaseData] : _phaseMap)
            {
                auto& offset = phaseData.offsets;
                auto& nReacParTyps = phaseData.nReacParType;
                offset.resize(dim, 0);
                nReacParTyps.resize(dim, 0);

            }
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
         * @brief Computes offsets and stores number of reactions per particle type
         * @param phaseType The phase type to compute offsets for
         * @param numOfReaOfParType Number of reactions for this particle type
         * @param parType The particle type index
         */
        void computeOffsetsAndReaOfParType(const std::string& phaseType, unsigned int numOfReaOfParType, unsigned int parType)
        {
            auto& offsets = getPhaseData(phaseType).offsets;
            auto& nReacParTyps = getPhaseData(phaseType).nReacParType;

            if (parType == 0)
            {
                offsets[parType] = 0;
                nReacParTyps[parType] = numOfReaOfParType;
            }
            else
            {
                offsets[parType] = offsets[parType - 1] + nReacParTyps[parType - 1];
                nReacParTyps[parType] = numOfReaOfParType;
            }
            
        }

        /**
         * @brief Gets the offset for a specific particle type in a phase
         * @param phaseType The phase type to query
         * @param parType The particle type index
         * @return Offset value for the given particle type
         * @throws InvalidParameterException if parType is out of bounds
         */
        int getOffsetForPhase(const std::string& phaseType, unsigned int parType) const
        {
            const auto& offsets = getPhaseData(phaseType).offsets;
            if (parType >= offsets.size()) {
                throw InvalidParameterException("ParType index " + std::to_string(parType) + 
                    " out of bounds for phase " + phaseType);
            }
            return offsets[parType];
        }

        /**
         * @brief Gets the number of reactions for a specific particle type in a phase
         * @param phaseType The phase type to query
         * @param parType The particle type index
         * @return Number of reactions for the given particle type
         * @throws InvalidParameterException if parType is out of bounds
         */
        int getnReactionOfParType(const std::string& phaseType, unsigned int parType) const
        {
            const auto& nReacType = getPhaseData(phaseType).nReacParType;
            if (parType >= nReacType.size()) {
                throw InvalidParameterException("ParType index " + std::to_string(parType) + 
                    " out of bounds for phase " + phaseType);
            }
            return nReacType[parType];
        }

        /**
         * @brief Configures discretization for reaction models in a specific phase
         * @param phaseType The phase type to configure
         * @param parType The particle type index
         * @param nReactions Number of reactions to configure
         * @param nComp Number of components
         * @param nBound Array of bound states per component
         * @param boundOffset Array of bound state offsets per component
         * @param paramProvider Parameter provider for configuration
         * @param helper Configuration helper for creating reaction models
         * @return True if configuration succeeded, false otherwise
         */
        bool configureDiscretization(std::string phaseType, unsigned int parType, unsigned int nReactions, unsigned int nComp, unsigned int* nBound, unsigned int* boundOffset,  IParameterProvider& paramProvider, const IConfigHelper& helper)
	    {
            auto& dynReaction = getDynReactionVector(phaseType);
            int offSet = getOffsetForPhase(phaseType, parType);
            configureDimensions(phaseType, nReactions);

            bool reactionConfSuccess = true;
            // Determine the number of reactions to configure
            const unsigned int nReac =
                (getPhaseData(phaseType).nReacParType[parType] == 0)
                ? static_cast<unsigned int>(getPhaseData(phaseType).dynReactions.size())
                : getPhaseData(phaseType).nReacParType[parType];

            // Configure each reaction model
            for (unsigned int i = 0; i < nReac; ++i)
				{
                    // for finite volumen units where the particles are not handling their own reactions.  
                    // neet to map between the unser interphase and the phase types 
                    std::string interphase_type = phaseType;
                    if (interphase_type == "pore")
                        interphase_type = "liquid";
                    
                    char reactionKey[32];
                    snprintf(reactionKey, sizeof(reactionKey), "%s_reaction_%03d", interphase_type.c_str(), i);

                    paramProvider.pushScope(reactionKey);

                    // Check if reaction type is specified
					if (!paramProvider.exists("TYPE")) 
					{
						paramProvider.popScope();
						throw InvalidParameterException("Missing 'type' parameter for " + std::string(reactionKey));
					}

                    // Create reaction model based on type
					std::string reactionType = paramProvider.getString("TYPE");
					dynReaction[offSet + i] = helper.createDynamicReactionModel(reactionType);

                    // Validate reaction model creation
					if (!dynReaction[offSet + i]) 
					{
						paramProvider.popScope();
						throw InvalidParameterException("Unknown dynamic reaction model " + reactionType +
							" for " + reactionKey);
					}

                    // Configure the reaction model discretization
					reactionConfSuccess = dynReaction[offSet + i]->configureModelDiscretization(paramProvider, nComp, nBound + parType * nComp, boundOffset + parType * nComp) && reactionConfSuccess;

                    // Handle configuration failure
					if (!reactionConfSuccess) 
					{
						if (dynReaction[offSet + i]->usesParamProviderInDiscretizationConfig())
							paramProvider.popScope();
						paramProvider.popScope();
						throw InvalidParameterException("Failed to configure reaction model " + reactionType +
							" for " + reactionKey);
					}

                    // Pop scope if reaction model used parameter provider
					if (dynReaction[offSet + i]->usesParamProviderInDiscretizationConfig())
						paramProvider.popScope();
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
            int offSet = getOffsetForPhase(phaseType, parType);
            unsigned int nTotalReactions = dynReaction.size();

            bool dynReactionConfSuccess = true;
             
            // Determine the number of reactions to configure
            const unsigned int nReac =
                (getPhaseData(phaseType).nReacParType[parType] == 0)
                ? static_cast<unsigned int>(getPhaseData(phaseType).dynReactions.size())
                : getPhaseData(phaseType).nReacParType[parType];
            
            // Configure each reaction
            for (unsigned int reac = 0; reac < nReac; ++reac)
            {
                // Skip null reactions or those that don't require configuration
                if (!dynReaction[offSet + reac] || !dynReaction[offSet + reac]->requiresConfiguration())
                    continue;

                // Create reaction scope key

                std::string interphase_type = phaseType;
                if (interphase_type == "pore")
                    interphase_type = "liquid";

                char reactionKey[32];
                snprintf(reactionKey, sizeof(reactionKey), "%s_reaction_%03d", interphase_type.c_str(), reac);
                paramProvider.pushScope(reactionKey); //scope reaction_xxx

                // Configure the reaction model
                dynReactionConfSuccess = dynReaction[offSet + reac]->configure(paramProvider, unitOpIdx, parType) && dynReactionConfSuccess;

                paramProvider.popScope();//scope reaction_xxx
            }
        
            return dynReactionConfSuccess;
        }

        /**
         * @brief Sets workspace requirements for all reaction models across all phases
         * @param lms Linear memory sizer to configure workspace requirements
         * @param nComp Number of components in the system
         */
        void setWorkspaceRequirements(std::string phaseType, unsigned int nParType, unsigned int nComp, unsigned int const* strideBound, LinearMemorySizer& lms) const
        {
            auto& dynReactionVector = getDynReactionVector(phaseType);
                
            for (auto par = 0; par < nParType; par++)
            {   
                int numReacOfPartical = getnReactionOfParType(phaseType, par);
                int offSet = getOffsetForPhase(phaseType,par);

                for (auto i = 0; i < numReacOfPartical; i++)
                {
                    if (dynReactionVector[offSet + i] && dynReactionVector[offSet + i]->requiresWorkspace())
                    {
                        lms.fitBlock(dynReactionVector[offSet+ i]->workspaceSize(nComp, strideBound[i], nullptr));
                    }
                }
            }
        }

        void setWorkspaceRequirements(std::string phaseType, unsigned int nComp, unsigned int const strideBound, LinearMemorySizer& lms) const
        {
            auto& dynReactionVector = getDynReactionVector(phaseType);
            for (auto i = 0; i < dynReactionVector.size(); i++)
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
				dynReactionVector.resize(1, nullptr);
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

        /**
         * @brief Default constructor - initializes phase map with default values
         */
        ReactionSystem(){ }

	};

} // namespace model
} // namespace cadet
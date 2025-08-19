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

struct ReactionSystem
	{   
        private:
        struct PhaseData
        {
            std::vector<IDynamicReactionModel*> dynReactions; //!< Dynamic reactions in the phase
            std::vector<int> offsets; //!< Offsets for reactions in the phase
            unsigned int totalReactions = 0; //!< Total reactions in the phase
        }; 

        std::map<std::string, PhaseData> phaseMap = {
            {"cross_phase", PhaseData{}},
            {"pore", PhaseData{}},
            {"solid", PhaseData{}}
        };

        PhaseData& getPhaseData(const std::string& phaseType)
        {
            auto phase = phaseMap.find(phaseType);
            if (phase == phaseMap.end())
                throw InvalidParameterException("Unknown phase type: " + phaseType);

            return phase->second;
        }

        // Const version of getPhaseData
        const PhaseData& getPhaseData(const std::string& phaseType) const
        {
            auto phase = phaseMap.find(phaseType);
            if (phase == phaseMap.end())
                throw InvalidParameterException("Unknown phase type: " + phaseType);

            return phase->second;
        }

        public:

        bool oldReactionInterface; //!< Flag to distinguish between old and new reaction interface
        unsigned int numberOfParticles; //!< Internal storage for the number of particles in the unit

        unsigned int getTotalReactions(const std::string& phaseType) const
        {
            return getPhaseData(phaseType).totalReactions;
        }

        std::vector<IDynamicReactionModel*>& getDynReactionVector(const std::string& phase_type)
        {
            return getPhaseData(phase_type).dynReactions;
        }

        unsigned int updateTotalReactions(const std::string& phaseType, unsigned int nReactions)
        {
            PhaseData& phaseData = getPhaseData(phaseType);
            phaseData.totalReactions += nReactions;
            return phaseData.totalReactions;
        }

        void configureDimensions(const std::string& phaseType, unsigned int numReac)
        {
            auto& reactionVector = getPhaseData(phaseType).dynReactions;
            reactionVector.resize(numReac, nullptr);
        }

        void configureDimensionsOffSet(unsigned int nParTypes)
        {

            for (auto phase : phaseMap)
            {
                auto& offset = getPhaseData(phase.first).offsets;
                offset.resize(nParTypes, 0);
            }           
        }

        void computeOffsets(const std::string& phaseType, unsigned int numOfReaOfParType, unsigned int parType)
        {
            auto& offsets = getPhaseData(phaseType).offsets;
            for (size_t i = 0; i <= parType; ++i)
            {
                if (i == 0)
                    offsets[i] = 0;
                else
                    offsets[i] = offsets[i - 1] + numOfReaOfParType;
            }
        }

        int getOffsetForPhase(const std::string& phaseType, unsigned int parType) const
        {
            const auto& offsets = getPhaseData(phaseType).offsets;
            return (parType < offsets.size()) ? offsets[parType] : 0;
        }

        bool configureDiscretization(std::string phaseType, unsigned int parType, unsigned int nReactions, unsigned int nComp, unsigned int* nBound, unsigned int* boundOffset,  IParameterProvider& paramProvider, const IConfigHelper& helper)
	    {
            auto& dynReaction = getDynReactionVector(phaseType);
            int offSet = getOffsetForPhase(phaseType, parType);
            configureDimensions(phaseType, nReactions);

            bool reactionConfSuccess = true;
            for (unsigned int i = 0; i < nReactions; ++i) 
				{
					char reactionKey[32];
                    snprintf(reactionKey, sizeof(reactionKey), "%s_reaction_%03d", phaseType.c_str(), i);

                    paramProvider.pushScope(reactionKey);

					if (!paramProvider.exists("TYPE")) 
					{
						paramProvider.popScope();
						throw InvalidParameterException("Missing 'type' parameter for " + std::string(reactionKey));
					}

					std::string reactionType = paramProvider.getString("TYPE");
					dynReaction[offSet + i] = helper.createDynamicReactionModel(reactionType);

					if (!dynReaction[offSet + i]) 
					{
						paramProvider.popScope();
						throw InvalidParameterException("Unknown dynamic reaction model " + reactionType +
							" for " + reactionKey);
					}

					reactionConfSuccess = dynReaction[offSet + i]->configureModelDiscretization(paramProvider, nComp, nBound + parType * nComp, boundOffset + parType * nComp) && reactionConfSuccess;

					if (!reactionConfSuccess) 
					{
						if (dynReaction[offSet + i]->usesParamProviderInDiscretizationConfig())
							paramProvider.popScope();
						paramProvider.popScope();
						throw InvalidParameterException("Failed to configure reaction model " + reactionType +
							" for " + reactionKey);
					}

					if (dynReaction[offSet + i]->usesParamProviderInDiscretizationConfig())
						paramProvider.popScope();
				}

            return reactionConfSuccess;
	    }

        bool configure(std::string phaseType, unsigned int parType, unsigned int unitOpIdx, unsigned int nReactions, IParameterProvider& paramProvider)
        {
            auto& dynReaction = getDynReactionVector(phaseType);
            int offSet = getOffsetForPhase(phaseType, parType);
            
            bool dynReactionConfSuccess = true;

            for (unsigned int reac = 0; reac < nReactions; ++reac)
            {
                if (!dynReaction[offSet + reac] || !dynReaction[offSet + reac]->requiresConfiguration())
                    continue;

                char reactionKey[32];
                snprintf(reactionKey, sizeof(reactionKey), "reaction_model_%03d", reac);
                paramProvider.pushScope(reactionKey);

                dynReactionConfSuccess = dynReaction[offSet + reac]->configure(paramProvider, unitOpIdx, parType) && dynReactionConfSuccess;

                paramProvider.popScope();
            }

            paramProvider.popScope();
        
            return dynReactionConfSuccess;
        }

        void setWorkspaceRequirements(LinearMemorySizer lms, unsigned int nComp)
        {
            for (auto phase: phaseMap)
            {   
                auto& dynReactionVector = getDynReactionVector(phase.first);
                for (auto i = 0; i < dynReactionVector.size(); i++)
                {
                        if (dynReactionVector[i] && dynReactionVector[i]->requiresWorkspace())
                            lms.fitBlock(dynReactionVector[i]->workspaceSize(nComp, 0, nullptr));
                }
            }
        }

		void clearDynamicReactionModels()
		{
			for (auto phase: phaseMap)
			{
				auto& dynReactionVector = getDynReactionVector(phase.first);
				for (auto* reac : dynReactionVector) delete reac;
			}
		}

        void empty()
        {
            for (auto phase : phaseMap)
            {
                auto& dynReacVec = getPhaseData(phase.first).dynReactions;
                dynReacVec.resize(1);
            }
        }

	};

} // namespace model
} // namespace cadet
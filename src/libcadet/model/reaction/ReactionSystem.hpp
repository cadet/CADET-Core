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

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>

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

        PhaseData& getPhaseData(std::string phaseType)
        {
            auto phase = phaseMap.find(phaseType);
            if (phase == phaseMap.end())
                throw InvalidParameterException("Unknown phase type: " + phaseType);

            return phase->second;
        }

        // Const version of getPhaseData
        const PhaseData& getPhaseData(std::string phaseType) const
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
            return getPhaseData(phaseType).totalReactions; // ← Kein const_cast mehr nötig!
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

        void configureDimensions(unsigned int numberOfParticles)
        {
            _reactions.numberOfParticles = numberOfParticles;
            _crossPhaseOffsets.resize(numberOfParticles,0);
            _poreOffsets.resize(numberOfParticles,0);
            _solidOffsets.resize(numberOfParticles,0);
        }

        int setCrossPhaseOffset(unsigned int numOfReaOfParType,  unsigned int parType) const
		{
			computeOffsets(numOfReaOfParType, parType, crossPhaseOffsets);
			return (parType < crossPhaseOffsets.size()) ? crossPhaseOffsets[parType] : 0;
		}

        int setPoreOffset(unsigned int numOfReaOfParType, unsigned int parType) const
		{
			computeOffsets(numOfReaOfParType, parType, poreOffsets);
			return (parType < poreOffsets.size()) ? poreOffsets[parType] : 0;
		}

		int setSolidOffset(unsigned int numOfReaOfParType, unsigned int parType) const
		{
			computeOffsets(numOfReaOfParType, parType, solidOffsets);
			return (parType < solidOffsets.size()) ? solidOffsets[parType] : 0;
		}

		static void computeOffsets(const std::string& phaseType, unsigned int numOfReaOfParType, unsigned int parType)
		{
            auto& offsets = getPhaseData(phaseType).offsets;
            for (size_t i = 0; i < parType; ++i)
			{
				if (i == 0)
					offsets[i] = numOfReaOfParType;
				else
					offsets[i] = offsets[i - 1] + numOfReaOfParType;
			}
		}

        std::vector<IDynamicReactionModel*>& getDynReactionVector(const std::string& phase_type)
        {
            if (phase_type == "CROSSPHASE" || phase_type == "crossphase")
                return dynReactionCrossPhase;
            else if (phase_type == "PORE" || phase_type == "pore")
                return dynReactionPore;
            else if (phase_type == "SOLID" || phase_type == "solid")
                return dynReactionSolid;
            else
                throw InvalidParameterException("Unknown phase type: " + phase_type);
        }

        int getOffsetForPhase(const std::string& phase_type, unsigned int parType) const
        {
            const auto& offsets = const_cast<ReactionSystem*>(this)->getPhaseData(phase_type).offsets;
            return (parType < offsets.size()) ? offsets[parType] : 0;
        }

        bool configureDiscretization(std::string phase_type, unsigned int parType, unsigned int nReactions, IParameterProvider& paramProvider, const IConfigHelper& helper)
	    {
            auto& dynReaction = getDynReactionVector(phase_type);
            int offSet = getOffsetForPhase(phase_type, parType);

            for (unsigned int i = 0; i < nReactions; ++i) 
				{
					char reactionKey[32];
					snprintf(reactionKey, sizeof(reactionKey), phase_type +"_%03d", i);

					if (!paramProvider.exists(reactionKey)) 
					{
						paramProvider.popScope();
						throw InvalidParameterException("Missing reaction model definition for " + std::string(reactionKey));
					}

					paramProvider.pushScope(reactionKey);

					if (!paramProvider.exists("TYPE")) 
					{
						paramProvider.popScope();
						throw InvalidParameterException("Missing 'type' parameter for " + std::string(reactionKey));
					}

					std::string reactionType = paramProvider.getString("TYPE");
					paramProvider.popScope();
					dynReaction[offSet + i] = helper.createDynamicReactionModel(reactionType);

					if (!dynReaction[offSet + i]) 
					{
						paramProvider.popScope();
						throw InvalidParameterException("Unknown dynamic reaction model " + reactionType +
							" for " + reactionKey);
					}

					if (dynReaction[offSet + i]->usesParamProviderInDiscretizationConfig())
						paramProvider.pushScope(reactionKey);

					reactionConfSuccess = dynReaction[offSet + i]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound + parType * _disc.nComp, _disc.boundOffset + parType * _disc.nComp) && reactionConfSuccess;

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
	    }

        bool configure(std::string reactionType, unsigned int parType, unsigned int unitOpIdx, IParameterProvider& paramProvider)
        {
            auto& dynReaction = getDynReactionVector(phase_type);
            int offSet = getOffsetForPhase(phase_type, parType);
            
            bool dynReactionConfSuccess = true;

            for (int reac = 0; reac < nReactions; ++reac)
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

        void setWorkspaceRequirements(LinearMemorySizer lms)
        {
            for (auto phase: phaseMap)
            {
                auto& dynReactionVector = _reactions.getDynReactionVector(phase);
                    for (auto i = 0; i < dynReactionVector.size(); i++)
                    {
                            if (dynReactionVector[i] && dynReactionVector[i]->requiresWorkspace())
                                lms.fitBlock(dynReactionVector[i]->workspaceSize(_disc.nComp, 0, nullptr));
                    }
            }
        }

		void clearDynamicReactionModels()
		{
			for (auto phase: phaseMap)
			{
				auto& dynReactionVector = _reactions.getDynReactionVector(phase);
				for (auto* reac : dynReactionVector) delete reac;
			}
		}

		~ReactionSystem()
		{
			for (auto* reac : _dynReactionBulk) delete reac;
			for (auto* reac : _dynReactionPoreParticle) delete reac;
			for (auto* reac : _dynReactionSolidParticle) delete reac;
		}
	};

} // namespace model
} // namespace cadet
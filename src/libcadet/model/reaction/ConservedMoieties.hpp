#pragma once

#include "model/ReactionModel.hpp"
#include "cadet/Exceptions.hpp"
#include "ConfigurationHelper.hpp"
#include "SensParamUtil.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <Eigen/Core>
#include <Eigen/SVD>

#include <utility>
#include <vector>

namespace cadet
{

namespace model
{



class ConservedMoieties
    {
        
        private:
        bool _enabled = false;

        unsigned int _nStates = 0;
        unsigned int _nReactionTotal = 0; // 1. only one reaction model later: for all reaction models
        unsigned int _nReactionEquilibirum = 0;
        
        std::vector<unsigned int> _reactionColumnOffset;
        std::vector<bool> _eqReactionMask;

        Eigen::MatrixXd _S; // N dim: nStates x nReaction
        Eigen::MatrixXd _eqS; // N_eq dim: nStates x nEqreaction
        
        Eigen::MatrixXd _L; // L dim: nMoities x nStates
        unsigned int _rank = 0;

        double TOL = 1e-15; 

        public:
        const Eigen::MatrixXd& getConservedMoietiesMatrix() const { return _L; }
        const Eigen::MatrixXd& conservedMoietyMatrix() const { return _L; }
        const Eigen::MatrixXd& stoichiometryMatrix() const { return _S; }
        const Eigen::MatrixXd& equilibriumStoichiometryMatrix() const { return _eqS; }


        bool isEnabled() const {return _enabled; }
        unsigned int rank() const {return _rank; }
        unsigned int numMoieties() const  {return static_cast<unsigned int>(_L.rows()); }
        unsigned int numEquilibriumReactions() const {return static_cast<unsigned int>(_eqS.cols()); }

        const std::vector<bool>& equilibriumReactionFlags() const { return _eqReactionFlags;}

        bool configure(unsigned int states, std::vector<unsigned int>&& reactionColumnOffset,
            std::vector<bool>&& eqReactionFlags, Eigen::MatrixXd&& stoichiometry, double rankTol)
        {
            _enabled = false;
            _nStates = states;
            _nReactionTotal = static_cast<unsigned int>(eqReactionFlags.size());
            _reactionColumnOffset = std::move(reactionColumnOffset);
            _eqReactionMask = std::move(eqReactionFlags);
            _S = std::move(stoichiometry);
            TOL = rankTol;

            extractEquilibriumStoichiometry();
            computeLeftNullspace();

            _enabled = true;
            return true;
        }

        void extractEquilibriumStoichiometry()
        {
            unsigned int nEq = 0;
            for (bool isEq : _eqReactionMask)
            {
                if (isEq) ++nEq;
            }

            _nReactionEquilibirum = nEq;
            _eqS.resize(_S.rows(), nEq);
            unsigned int col = 0;

            for (unsigned int r = 0; r < _eqReactionMask.size(); ++r)
            {
                if (!_eqReactionMask[r])
                    continue;
                _eqS.col(col) = _S.col(r);
                col++;
            }

        }

        void computeLeftNullspace()
        {
            if (_eqS.cols() == 0)
            { 
                _rank = 0;
                _L = Eigen::MatrixXd::Identity(0, _eqS.rows());
                return;
            }

            Eigen::MatrixXd A = _eqS.transpose();
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
            
            const auto& singularValues = svd.singularValues();

            _rank = 0;
            for (int i = 0; i < singularValues.size(); ++i)
            {
                if (singularValues(i) > TOL)
                    ++_rank;
            }

            if (_nReactionEquilibirum != _rank)
	            throw InvalidParameterException("Conserved Moieties: Redundant equilibrium reactions are not supported");

            const unsigned int nullity = static_cast<unsigned int>(A.cols()) - _rank;
            
            if (nullity == 0)
            {
                _L = Eigen::MatrixXd::Zero(0, _eqS.rows());
                return;
            }
            
            const Eigen::MatrixXd V = svd.matrixV();
            
            // Columns rank ... end span null(A)
            const Eigen::MatrixXd nullspace = V.rightCols(nullity);
            // L such that L * N = 0
            _L = nullspace.transpose();
        }

    };


} //model

} //cadet

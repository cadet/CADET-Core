// SPDX-License-Identifier: AGPL-3.0-or-later
// =================================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Affero General Public
//  License v3.0 (or, at your option, any later version).
// =================================================================================

#pragma once

#include "model/ReactionModel.hpp"
#include "cadet/Exceptions.hpp"
#include "ConfigurationHelper.hpp"
#include "SensParamUtil.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>
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
        unsigned int _nReactionTotal = 0;
        unsigned int _nReactionEquilibirum = 0;
        
        std::vector<unsigned int> _reactionColumnOffset;
        std::vector<bool> _eqReactionMask;

        Eigen::MatrixXd _S; // N dim: nStates x nReaction
        Eigen::MatrixXd _eqS; // N_eq dim: nStates x nEqreaction
        
        Eigen::MatrixXd _L; // L dim: nMoities x nStates
        unsigned int _rank = 0;

        double _TOL = 1e-15;

        public:
        const Eigen::MatrixXd& getConservedMoietiesMatrix() const { return _L; }
        const Eigen::MatrixXd& conservedMoietyMatrix() const { return _L; }
        const Eigen::MatrixXd& stoichiometryMatrix() const { return _S; }
        const Eigen::MatrixXd& equilibriumStoichiometryMatrix() const { return _eqS; }


        bool isEnabled() const {return _enabled; }
        unsigned int rank() const {return _rank; }
        unsigned int numMoieties() const  {return static_cast<unsigned int>(_L.rows()); }
        unsigned int numEquilibriumReactions() const {return static_cast<unsigned int>(_eqS.cols()); }

        const std::vector<bool>& equilibriumReactionFlags() const { return _eqReactionMask;}

        void clear()
        {
            _enabled = false;
            _nStates = 0;
            _nReactionTotal = 0;
            _nReactionEquilibirum = 0;
            _reactionColumnOffset.clear();
            _eqReactionMask.clear();
            _S.resize(0, 0);
            _eqS.resize(0, 0);
            _L.resize(0, 0);
            _rank = 0;
        }

        bool configure(unsigned int states, std::vector<unsigned int>&& reactionColumnOffset,
            std::vector<bool>&& eqReactionFlags, Eigen::MatrixXd&& stoichiometry, double rankTol)
        {
            _enabled = false;
            _nStates = states;
            _nReactionTotal = static_cast<unsigned int>(eqReactionFlags.size());
            _reactionColumnOffset = std::move(reactionColumnOffset);
            _eqReactionMask = std::move(eqReactionFlags);
            _S = std::move(stoichiometry);
            _TOL = rankTol;

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
                _eqS.col(col++) = _S.col(r);
            }

        }

        void computeLeftNullspace()
        {
            if (_eqS.cols() == 0)
            { 
                _rank = 0;
                _L = Eigen::MatrixXd::Identity(_eqS.rows(), _eqS.rows());
                return;
            }

            Eigen::MatrixXd A = _eqS.transpose();
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
            
            const auto& singularValues = svd.singularValues();

            _rank = 0;
            for (int i = 0; i < singularValues.size(); ++i)
            {
                if (singularValues(i) > _TOL)
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

        /**
         * @brief Multiplies the conserved-moiety matrix with a state-sized vector
         * @param [out] targetVec Result vector with numMoieties() entries
         * @param [in] sourceVec Source vector with @p size entries, which must not alias @p targetVec
         * @param [in] size Number of source entries, equal to the number of states
         */
        template <typename TargetType, typename SourceType>
        void multiplyToVector(TargetType* const targetVec, SourceType const* const sourceVec, unsigned int size) const
        {
            cadet_assert(size == static_cast<unsigned int>(_L.cols()));

            for (unsigned int moiety = 0; moiety < numMoieties(); ++moiety)
            {
                targetVec[moiety] = TargetType{0.0};
                for (unsigned int state = 0; state < size; ++state)
                    targetVec[moiety] += static_cast<double>(_L(moiety, state)) * sourceVec[state];
            }
        }

        /**
         * @brief Replaces a sparse row block by its conserved-moiety transformation
         * @details All @p numStates rows in the target block are cleared. The first numMoieties()
         *          rows are then set to L times the corresponding source rows. The source and
         *          target matrices must not alias. The target pattern must contain the union
         *          of all source row patterns.
         */
        void multiplyToMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& targetMat,
            Eigen::SparseMatrix<double, Eigen::RowMajor> const& sourceMat, unsigned int numStates, unsigned int rowOffset) const
        {
            multiplyToMatrix(targetMat, sourceMat, numStates, rowOffset, 0, targetMat.cols());
        }

        void multiplyToMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& targetMat,
            Eigen::SparseMatrix<double, Eigen::RowMajor> const& sourceMat, unsigned int numStates,
            unsigned int rowOffset, Eigen::Index firstColumn, Eigen::Index lastColumn) const
        {
            cadet_assert(numStates == static_cast<unsigned int>(_L.cols()));
            cadet_assert(rowOffset + numStates <= static_cast<unsigned int>(targetMat.rows()));
            cadet_assert(rowOffset + numStates <= static_cast<unsigned int>(sourceMat.rows()));
            cadet_assert((firstColumn >= 0) && (firstColumn <= lastColumn));
            cadet_assert(lastColumn <= targetMat.cols());
            cadet_assert(lastColumn <= sourceMat.cols());
            cadet_assert(&targetMat != &sourceMat);

            // Clear the transformed row block once.
            for (unsigned int row = 0; row < numStates; ++row)
            {
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(targetMat, rowOffset + row); it; ++it)
                {
                    if ((it.col() >= firstColumn) && (it.col() < lastColumn))
                        it.valueRef() = 0.0;
                }
            }

            // Apply L to complete sparse rows.
            for (unsigned int moiety = 0; moiety < numMoieties(); ++moiety)
            {
                const Eigen::Index targetRow = rowOffset + moiety;
                for (unsigned int state = 0; state < numStates; ++state)
                {
                    const double factor = _L(moiety, state);
                    if (factor == 0.0)
                        continue;

                    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(sourceMat, rowOffset + state); it; ++it)
                    {
                        if ((it.col() >= firstColumn) && (it.col() < lastColumn))
                            targetMat.coeffRef(targetRow, it.col()) += factor * it.value();
                    }
                }
            }
        }

        struct EigenSparseMatrixEntry
        {
            unsigned int sourceRow;
            Eigen::Index column;
            double value;
        };

        /**
         * @brief Replaces a sparse row block in place by its conserved-moiety transformation
         * @details The affected source entries are preserved in @p scratch before their rows are
         *          cleared. The scratch capacity is retained between calls. The matrix pattern must
         *          contain the union of all source row patterns in every conserved-moiety row.
         */
        void multiplyToMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix, unsigned int numStates, unsigned int rowOffset,
            Eigen::Index firstColumn, Eigen::Index lastColumn, std::vector<EigenSparseMatrixEntry>& orig) const
        {
            cadet_assert(numStates == static_cast<unsigned int>(_L.cols()));
            cadet_assert(rowOffset + numStates <= static_cast<unsigned int>(matrix.rows()));
            cadet_assert((firstColumn >= 0) && (firstColumn <= lastColumn));
            cadet_assert(lastColumn <= matrix.cols());

            orig.clear();

            // Save and clear source rows
            for (unsigned int state = 0; state < numStates; ++state)
            {
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(matrix, rowOffset + state); it; ++it)
                {
                    if ((it.col() < firstColumn) || (it.col() >= lastColumn))
                        continue;

                    orig.push_back({state, it.col(), it.value()});
                    it.valueRef() = 0.0;
                }
            }

            // Rebuild conserved rows from the preserved entries
            for (const EigenSparseMatrixEntry& entry : orig)
            {
                for (unsigned int m = 0; m < numMoieties(); ++m)
                {
                    const double factor = _L(m, entry.sourceRow);
                    if (factor == 0.0)
                        continue;

                    matrix.coeffRef(rowOffset + m, entry.column) += factor * entry.value;
                }
            }
        }

        void multiplyToMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix, unsigned int numStates, unsigned int rowOffset,
            std::vector<EigenSparseMatrixEntry>& scratch) const
        {
            multiplyToMatrix(
                matrix, numStates, rowOffset, 0, matrix.cols(), scratch);
        }

    };


} //model

} //cadet

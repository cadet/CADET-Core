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

#include "cadet/Exceptions.hpp"
#include "common/CompilerSpecific.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SVD>

#include <algorithm>
#include <limits>
#include <vector>

namespace cadet
{

namespace model
{

class ConservedMoieties
{
public:
	struct EigenSparseMatrixEntry
	{
		unsigned int sourceRow;
		Eigen::Index column;
		double value;
	};

	const Eigen::MatrixXd& conservedMoietyMatrix() const CADET_NOEXCEPT { return _L; }

	bool isEnabled() const CADET_NOEXCEPT { return _enabled; }
	unsigned int numMoieties() const CADET_NOEXCEPT { return static_cast<unsigned int>(_L.rows()); }
	unsigned int numEquilibriumReactions() const CADET_NOEXCEPT
	{
		return static_cast<unsigned int>(_L.cols() - _L.rows());
	}

	void clear()
	{
		_enabled = false;
		_L.resize(0, 0);
	}

	bool configure(unsigned int numStates, const std::vector<bool>& equilibriumReactionFlags,
		const Eigen::MatrixXd& stoichiometry, double rankTolerance)
	{
		cadet_assert(numStates == static_cast<unsigned int>(stoichiometry.rows()));
		cadet_assert(equilibriumReactionFlags.size() == static_cast<std::size_t>(stoichiometry.cols()));

		_enabled = false;
		computeConservedMoieties(equilibriumReactionFlags, stoichiometry, rankTolerance);
		_enabled = true;
		return true;
	}

	/**
	 * @brief Multiplies the conserved-moiety matrix with a state-sized vector
	 * @param [out] targetVec Result vector with numMoieties() entries
	 * @param [in] sourceVec Source vector with @p size entries, which must not alias @p targetVec
	 * @param [in] size Number of source entries, equal to the number of states
	 */
	template <typename TargetType, typename SourceType>
	void applyToVector(TargetType* const targetVec, SourceType const* const sourceVec, unsigned int size) const
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
	 * @brief Replaces a vector prefix by its conserved-moiety transformation
	 * @details The transformed values are first written to @p scratch so that @p vector can
	 *          also serve as the source. The scratch buffer must provide numMoieties() entries.
	 * @param [in,out] vector State-sized source vector whose first numMoieties() entries are replaced
	 * @param [in] size Number of source entries, equal to the number of states
	 * @param [out] scratch Temporary buffer with numMoieties() entries
	 */
	template <typename ValueType>
	void applyToVector(ValueType* const vector, unsigned int size, ValueType* const scratch) const
	{
		applyToVector(scratch, vector, size);
		std::copy_n(scratch, numMoieties(), vector);
	}

	/**
	 * @brief Replaces a vector prefix by its conserved-moiety transformation
	 * @details The scratch capacity is retained between calls.
	 */
	template <typename ValueType>
	void applyToVector(ValueType* const vector, unsigned int size, std::vector<ValueType>& scratch) const
	{
		scratch.resize(numMoieties());
		applyToVector(vector, size, scratch.data());
	}

	/**
	 * @brief Applies the conserved-moiety transformation to a time-derivative result
	 * @details Conserved-moiety rows are transformed and equilibrium-reaction rows are
	 *          set to zero because equilibrium equations are algebraic.
	 */
	template <typename TargetType, typename SourceType>
	void applyToDerivativeVector(TargetType* const targetVec, SourceType const* const sourceVec, unsigned int size) const
	{
		cadet_assert(numMoieties() + numEquilibriumReactions() == size);
		applyToVector(targetVec, sourceVec, size);
		std::fill_n(targetVec + numMoieties(), numEquilibriumReactions(), TargetType{0.0});
	}

	template <typename ValueType>
	void applyToDerivativeVector(ValueType* const vector, unsigned int size, std::vector<ValueType>& scratch) const
	{
		cadet_assert(numMoieties() + numEquilibriumReactions() == size);
		applyToVector(vector, size, scratch);
		std::fill_n(vector + numMoieties(), numEquilibriumReactions(), ValueType{0.0});
	}

	/**
	 * @brief Replaces rows of a dense matrix by their conserved-moiety transformation
	 * @details Each column is transformed in place. The scratch buffer must provide
	 *          numMoieties() entries, independent of the matrix size.
	 */
	template <typename MatrixType>
	void applyToMatrix(MatrixType& matrix, unsigned int numStates, unsigned int rowOffset,
		unsigned int firstColumn, unsigned int lastColumn, double* const scratch) const
	{
		cadet_assert(numStates == static_cast<unsigned int>(_L.cols()));

		for (unsigned int column = firstColumn; column < lastColumn; ++column)
		{
			for (unsigned int moiety = 0; moiety < numMoieties(); ++moiety)
			{
				scratch[moiety] = 0.0;
				for (unsigned int state = 0; state < numStates; ++state)
					scratch[moiety] += _L(moiety, state) * matrix.native(rowOffset + state, column);
			}

			for (unsigned int moiety = 0; moiety < numMoieties(); ++moiety)
				matrix.native(rowOffset + moiety, column) = scratch[moiety];
		}
	}

	std::size_t matrixScratchSize(Eigen::SparseMatrix<double, Eigen::RowMajor> const& matrix,
		unsigned int numStates, unsigned int rowOffset) const
	{
		return matrixScratchSize(matrix, numStates, rowOffset, 0, matrix.cols());
	}

	/**
	 * @brief Replaces a sparse row block in place by its conserved-moiety transformation
	 * @details The affected source entries are preserved in @p scratch before their rows are
	 *          cleared. The scratch capacity is retained between calls. The matrix pattern must
	 *          contain the union of all source row patterns in every conserved-moiety row.
	 */
	void applyToMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix, unsigned int numStates,
		unsigned int rowOffset, Eigen::Index firstColumn, Eigen::Index lastColumn,
		std::vector<EigenSparseMatrixEntry>& scratch) const
	{
		scratch.resize(matrixScratchSize(matrix, numStates, rowOffset, firstColumn, lastColumn));
		applyToMatrix(matrix, numStates, rowOffset, firstColumn, lastColumn, scratch.data(), scratch.size());
	}

	/**
	 * @brief Replaces a sparse row block in place using caller-provided scratch memory
	 * @details The scratch array must contain at least matrixScratchSize() entries.
	 */
	void applyToMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix, unsigned int numStates,
		unsigned int rowOffset, Eigen::Index firstColumn, Eigen::Index lastColumn,
		EigenSparseMatrixEntry* const scratch, std::size_t scratchSize) const
	{
		cadet_assert(numStates == static_cast<unsigned int>(_L.cols()));
		cadet_assert(rowOffset + numStates <= static_cast<unsigned int>(matrix.rows()));
		cadet_assert((firstColumn >= 0) && (firstColumn <= lastColumn));
		cadet_assert(lastColumn <= matrix.cols());

		std::size_t numEntries = 0;
		for (unsigned int state = 0; state < numStates; ++state)
		{
			for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(matrix, rowOffset + state); it; ++it)
			{
				if ((it.col() < firstColumn) || (it.col() >= lastColumn))
					continue;

				cadet_assert(numEntries < scratchSize);
				scratch[numEntries++] = {state, it.col(), it.value()};
				it.valueRef() = 0.0;
			}
		}
		cadet_assert(numEntries == scratchSize);

		for (std::size_t entryIndex = 0; entryIndex < numEntries; ++entryIndex)
		{
			const EigenSparseMatrixEntry& entry = scratch[entryIndex];
			for (unsigned int moiety = 0; moiety < numMoieties(); ++moiety)
			{
				const double factor = _L(moiety, entry.sourceRow);
				if (factor != 0.0)
					matrix.coeffRef(rowOffset + moiety, entry.column) += factor * entry.value;
			}
		}
	}

	void applyToMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix, unsigned int numStates,
		unsigned int rowOffset, std::vector<EigenSparseMatrixEntry>& scratch) const
	{
		applyToMatrix(matrix, numStates, rowOffset, 0, matrix.cols(), scratch);
	}

	/**
	 * @brief Adds the row-pattern required by a conserved-moiety transformation
	 * @details Every conserved-moiety row receives the union of the patterns of all
	 *          source rows. Only entries present when the method is called are used as
	 *          source entries, so generated entries are never transformed recursively.
	 */
	void applyToPattern(std::vector<Eigen::Triplet<double>>& entries, unsigned int numStates,
		unsigned int rowOffset) const
	{
		applyToPattern(entries, numStates, rowOffset,
			std::numeric_limits<Eigen::Index>::lowest(), std::numeric_limits<Eigen::Index>::max());
	}

	void applyToPattern(std::vector<Eigen::Triplet<double>>& entries, unsigned int numStates,
		unsigned int rowOffset, Eigen::Index firstColumn, Eigen::Index lastColumn) const
	{
		cadet_assert(numStates == static_cast<unsigned int>(_L.cols()));
		cadet_assert(firstColumn <= lastColumn);

		const std::size_t originalSize = entries.size();
		std::size_t sourceEntries = 0;
		for (std::size_t i = 0; i < originalSize; ++i)
		{
			const Eigen::Index row = entries[i].row();
			const Eigen::Index column = entries[i].col();
			if ((row >= rowOffset) && (row < rowOffset + numStates)
				&& (column >= firstColumn) && (column < lastColumn))
				++sourceEntries;
		}

		entries.reserve(entries.size() + static_cast<std::size_t>(numMoieties()) * sourceEntries);
		for (std::size_t i = 0; i < originalSize; ++i)
		{
			const Eigen::Index sourceRow = entries[i].row();
			const Eigen::Index sourceColumn = entries[i].col();
			if ((sourceRow < rowOffset) || (sourceRow >= rowOffset + numStates)
				|| (sourceColumn < firstColumn) || (sourceColumn >= lastColumn))
				continue;

			for (unsigned int moiety = 0; moiety < numMoieties(); ++moiety)
				entries.emplace_back(rowOffset + moiety, sourceColumn, 0.0);
		}
	}

private:
	std::size_t matrixScratchSize(Eigen::SparseMatrix<double, Eigen::RowMajor> const& matrix,
		unsigned int numStates, unsigned int rowOffset, Eigen::Index firstColumn,
		Eigen::Index lastColumn) const
	{
		cadet_assert(numStates == static_cast<unsigned int>(_L.cols()));
		cadet_assert(rowOffset + numStates <= static_cast<unsigned int>(matrix.rows()));

		std::size_t numEntries = 0;
		for (unsigned int state = 0; state < numStates; ++state)
		{
			for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(matrix, rowOffset + state); it; ++it)
			{
				if ((it.col() >= firstColumn) && (it.col() < lastColumn))
					++numEntries;
			}
		}
		return numEntries;
	}

	void computeConservedMoieties(const std::vector<bool>& equilibriumReactionFlags,
		const Eigen::MatrixXd& stoichiometry, double rankTolerance)
	{
		const unsigned int numEquilibriumReactions = static_cast<unsigned int>(
			std::count(equilibriumReactionFlags.begin(), equilibriumReactionFlags.end(), true));

		Eigen::MatrixXd equilibriumStoichiometry(stoichiometry.rows(), numEquilibriumReactions);
		unsigned int equilibriumColumn = 0;
		for (unsigned int reaction = 0; reaction < equilibriumReactionFlags.size(); ++reaction)
		{
			if (equilibriumReactionFlags[reaction])
				equilibriumStoichiometry.col(equilibriumColumn++) = stoichiometry.col(reaction);
		}

		if (numEquilibriumReactions == 0)
		{
			_L = Eigen::MatrixXd::Identity(stoichiometry.rows(), stoichiometry.rows());
			return;
		}

		const Eigen::MatrixXd transposedStoichiometry = equilibriumStoichiometry.transpose();
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(transposedStoichiometry, Eigen::ComputeFullV);

		unsigned int rank = 0;
		for (int singularValue = 0; singularValue < svd.singularValues().size(); ++singularValue)
		{
			if (svd.singularValues()(singularValue) > rankTolerance)
				++rank;
		}

		if (numEquilibriumReactions != rank)
			throw InvalidParameterException("Conserved Moieties: Redundant equilibrium reactions are not supported");

		const unsigned int nullity = static_cast<unsigned int>(transposedStoichiometry.cols()) - rank;
		_L = svd.matrixV().rightCols(nullity).transpose();
	}

	bool _enabled = false;
	Eigen::MatrixXd _L;
};

} // namespace model

} // namespace cadet

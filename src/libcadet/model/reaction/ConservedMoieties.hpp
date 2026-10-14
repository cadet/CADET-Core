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
#include "cadet/cadetCompilerInfo.hpp"
#include "common/CompilerSpecific.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SVD>

#include <algorithm>
#include <cstddef>
#include <vector>

namespace cadet
{

namespace model
{


class ConservedMoieties
{
public:
	const Eigen::MatrixXd& conservedMoietyMatrix() const CADET_NOEXCEPT { return _L; }

	bool isEnabled() const CADET_NOEXCEPT { return _enabled; }
	unsigned int rank() const CADET_NOEXCEPT { return numEquilibriumReactions(); }
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

	bool configure(unsigned int numStates, const std::vector<bool>& equilibriumReactionFlags, const Eigen::MatrixXd& stoichiometry, double rankTolerance)
	{
		cadet_assert(numStates == static_cast<unsigned int>(stoichiometry.rows()));
		cadet_assert(equilibriumReactionFlags.size() == static_cast<std::size_t>(stoichiometry.cols()));

		_enabled = false;
		computeLeftNullspace(equilibriumReactionFlags, stoichiometry, rankTolerance);
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

	std::size_t matrixBufferSize(Eigen::SparseMatrix<double, Eigen::RowMajor> const& matrix, unsigned int numStates,
		unsigned int rowOffset, Eigen::Index firstColumn, Eigen::Index lastColumn) const
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
				if ((it.col() >= firstColumn) && (it.col() < lastColumn))
					++numEntries;
			}
		}
		return numEntries;
	}

	/**
	 * @brief Replaces a sparse row block in place using caller-provided buffer memory
	 * @details The affected source entries are preserved in @p buffer before their rows are
	 *          cleared. The buffer array must contain at least matrixBufferSize() entries. The
	 *          matrix pattern must contain the union of all source row patterns in every
	 *          conserved-moiety row.
	 */
	void applyToMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix, unsigned int numStates, unsigned int rowOffset, Eigen::Index firstColumn, Eigen::Index lastColumn, Eigen::Triplet<double>* const buffer, std::size_t bufferSize) const
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

				cadet_assert(numEntries < bufferSize);
				buffer[numEntries++] = Eigen::Triplet<double>(state, it.col(), it.value());
				it.valueRef() = 0.0;
			}
		}
		cadet_assert(numEntries <= bufferSize);

		for (std::size_t entryIndex = 0; entryIndex < numEntries; ++entryIndex)
		{
			const Eigen::Triplet<double>& entry = buffer[entryIndex];
			for (unsigned int moiety = 0; moiety < numMoieties(); ++moiety)
			{
				const double factor = _L(moiety, entry.row());
				if (factor != 0.0)
					matrix.coeffRef(rowOffset + moiety, entry.col()) += factor * entry.value();
			}
		}
	}

	/**
	 * @brief Multiplies the conserved-moiety transformation to a time-derivative result
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

	/**
	 * @brief Adds the row-pattern for one or more conserved-moiety row blocks
	 * @details Every conserved-moiety row receives the union of the patterns of all source rows
	 *          in its block and within the selected column range. All blocks are processed in one
	 *          pass over @p entries. Only entries present when the method is called are used as
	 *          source entries, so generated entries are never transformed recursively.
	 */
	void addPatternToBlocks(std::vector<Eigen::Triplet<double>>& entries, unsigned int numStates, unsigned int firstRowOffset, unsigned int numBlocks, unsigned int rowStride, Eigen::Index firstColumn, Eigen::Index lastColumn) const
	{
		cadet_assert(numStates == static_cast<unsigned int>(_L.cols()));
		cadet_assert(rowStride >= numStates);
		cadet_assert(firstColumn <= lastColumn);
		if ((numBlocks == 0) || (numStates == 0))
			return;

		const std::size_t originalSize = entries.size();
		std::size_t sourceEntries = 0;
		for (std::size_t i = 0; i < originalSize; ++i)
		{
			const Eigen::Index column = entries[i].col();
			if ((column < firstColumn) || (column >= lastColumn))
				continue;

			const Eigen::Index relativeRow = entries[i].row() - static_cast<Eigen::Index>(firstRowOffset);
			if (relativeRow < 0)
				continue;

			const unsigned int block = static_cast<unsigned int>(relativeRow) / rowStride;
			const unsigned int localRow = static_cast<unsigned int>(relativeRow) % rowStride;
			if ((block < numBlocks) && (localRow < numStates))
				++sourceEntries;
		}

		entries.reserve(entries.size() + static_cast<std::size_t>(numMoieties()) * sourceEntries);
		for (std::size_t i = 0; i < originalSize; ++i)
		{
			const Eigen::Index sourceColumn = entries[i].col();
			if ((sourceColumn < firstColumn) || (sourceColumn >= lastColumn))
				continue;

			const Eigen::Index relativeRow = entries[i].row() - static_cast<Eigen::Index>(firstRowOffset);
			if (relativeRow < 0)
				continue;

			const unsigned int block = static_cast<unsigned int>(relativeRow) / rowStride;
			const unsigned int localRow = static_cast<unsigned int>(relativeRow) % rowStride;
			if ((block >= numBlocks) || (localRow >= numStates))
				continue;

			const Eigen::Index targetRowOffset = firstRowOffset + block * rowStride;
			for (unsigned int moiety = 0; moiety < numMoieties(); ++moiety)
				entries.emplace_back(targetRowOffset + moiety, sourceColumn, 0.0);
		}
	}

private:
	bool _enabled = false;
	Eigen::MatrixXd _L; // L dim: nMoities x nStates

	void computeLeftNullspace(const std::vector<bool>& equilibriumReactionFlags, const Eigen::MatrixXd& stoichiometry, double rankTolerance)
	{
		const unsigned int numEquilibriumReactions = static_cast<unsigned int>(std::count(equilibriumReactionFlags.begin(), equilibriumReactionFlags.end(), true));
		if (numEquilibriumReactions == 0)
		{
			_L = Eigen::MatrixXd::Identity(stoichiometry.rows(), stoichiometry.rows());
			return;
		}

		Eigen::MatrixXd transposedStoichiometry(numEquilibriumReactions, stoichiometry.rows());
		unsigned int equilibriumRow = 0;
		for (unsigned int reaction = 0; reaction < equilibriumReactionFlags.size(); ++reaction)
		{
			if (equilibriumReactionFlags[reaction])
				transposedStoichiometry.row(equilibriumRow++) = stoichiometry.col(reaction).transpose();
		}

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

};


} // namespace model

} // namespace cadet

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

#ifndef LIBCADET_EigenSolverWrapper_HPP_
#define LIBCADET_EigenSolverWrapper_HPP_

#include <Eigen/Sparse>

#include <algorithm>
#include <memory>
#include <vector>

namespace cadet
{

namespace linalg
{

/**
 * @brief Storage order the models assemble their Jacobians in
 */
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMatrixRM;

/**
 * @brief Square diagonal block of a row-major Jacobian, in the storage order a solver needs
 * @details The models assemble their Jacobian row-major, while Eigen's SparseLU requires a
 *          column-major matrix. Letting Eigen convert at the call site re-sorts and
 *          re-allocates the whole block at every Newton iteration, although only the values
 *          change: the sparsity pattern is fixed for the whole simulation.
 *
 *          This class does the symbolic work once in prepare() and records, for every entry
 *          of the block, where it sits in the source value array and where it sits in the
 *          target value array. Every later update() is a scatter of doubles through that
 *          permutation, without allocation and without sorting. Should the pattern or the
 *          requested block change after all, update() notices and redoes the symbolic work.
 *
 * @tparam TargetOptions Storage order the solver needs, Eigen::ColMajor or Eigen::RowMajor
 */
template <int TargetOptions>
class DiagonalBlockView
{
public:
	typedef Eigen::SparseMatrix<double, TargetOptions> TargetMatrix;

	/**
	 * @brief Lays out the block @p src[offset, offset+size) and builds the permutation
	 * @param [in] src Row-major source matrix
	 * @param [in] offset Index of the first row and column of the block
	 * @param [in] size Number of rows and columns of the block
	 */
	void prepare(const SparseMatrixRM& src, int offset, int size)
	{
		_offset = offset;
		_size = size;
		_srcNonZeros = static_cast<int>(src.nonZeros());
		_srcCompressed = src.isCompressed();

		if (!_srcCompressed)
		{
			// The permutation needs the compressed arrays, so leave this case to Eigen
			_target = src.block(offset, offset, size, size);
			_target.makeCompressed();
			_srcIdx.clear();
			_srcForDst.clear();
			return;
		}

		int const* const outer = src.outerIndexPtr();
		int const* const inner = src.innerIndexPtr();

		_srcIdx.clear();
		_rows.clear();
		_cols.clear();

		std::vector<Eigen::Triplet<double>> triplets;
		for (int row = offset; row < offset + size; ++row)
		{
			for (int k = outer[row]; k < outer[row + 1]; ++k)
			{
				const int col = inner[k];
				if ((col < offset) || (col >= offset + size))
					continue;

				_srcIdx.push_back(k);
				_rows.push_back(row - offset);
				_cols.push_back(col - offset);
				triplets.emplace_back(row - offset, col - offset, 0.0);
			}
		}

		_target.resize(size, size);
		_target.setFromTriplets(triplets.begin(), triplets.end());
		_target.makeCompressed();

		// Record, for every entry of the target in its own storage order, where to read it
		// from in the source. Indexing by the destination keeps the writes of update()
		// sequential: a scatter over the target instead would miss cache on every entry
		// once the block outgrows the last level cache, which costs more than the
		// conversion this class exists to avoid.
		int const* const tOuter = _target.outerIndexPtr();
		int const* const tInner = _target.innerIndexPtr();

		_srcForDst.assign(static_cast<std::size_t>(_target.nonZeros()), 0);
		for (std::size_t e = 0; e < _srcIdx.size(); ++e)
		{
			const int o = (TargetOptions & Eigen::RowMajor) ? _rows[e] : _cols[e];
			const int i = (TargetOptions & Eigen::RowMajor) ? _cols[e] : _rows[e];

			int const* const begin = tInner + tOuter[o];
			int const* const end = tInner + tOuter[o + 1];
			const std::size_t dst = static_cast<std::size_t>(std::lower_bound(begin, end, i) - tInner);
			_srcForDst[dst] = _srcIdx[e];
		}
	}

	/**
	 * @brief Copies the current values of the block into the target matrix
	 */
	void update(const SparseMatrixRM& src, int offset, int size)
	{
		if ((offset != _offset) || (size != _size) || (static_cast<int>(src.nonZeros()) != _srcNonZeros)
			|| (src.isCompressed() != _srcCompressed) || !src.isCompressed())
		{
			prepare(src, offset, size);
			if (!_srcCompressed)
				return; // prepare() has copied the values already
		}

		double const* const sv = src.valuePtr();
		double* const tv = _target.valuePtr();

		const std::size_t n = _srcForDst.size();
		for (std::size_t d = 0; d < n; ++d)
			tv[d] = sv[_srcForDst[d]];
	}

	const TargetMatrix& matrix() const { return _target; }

private:
	TargetMatrix _target;
	std::vector<int> _srcIdx; //!< Index into src.valuePtr() of each entry of the block
	std::vector<int> _srcForDst; //!< For each entry of _target, its index in src.valuePtr()
	std::vector<int> _rows; //!< Block local row of each entry
	std::vector<int> _cols; //!< Block local column of each entry

	int _offset = 0;
	int _size = 0;
	int _srcNonZeros = -1;
	bool _srcCompressed = false;
};

class EigenSolverBase {
public:
    virtual ~EigenSolverBase() = default;

    /**
     * @brief Analyzes the pattern of the square diagonal block @p mat[offset, offset+size)
     */
    virtual void analyzePattern(const SparseMatrixRM& mat, int offset, int size) = 0;

    /**
     * @brief Factorizes the square diagonal block @p mat[offset, offset+size)
     */
    virtual void factorize(const SparseMatrixRM& mat, int offset, int size) = 0;

    /**
     * @brief Solves the linear system in place
     * @param [in,out] rhs Right hand side on entry, solution on exit
     * @param [in] size Number of unknowns
     */
    virtual void solve(double* rhs, int size) = 0;

    virtual Eigen::ComputationInfo info() const = 0;
};

template <typename OrderingType>
class SparseLU : public EigenSolverBase {
public:
    void analyzePattern(const SparseMatrixRM& mat, int offset, int size) override {
        _view.prepare(mat, offset, size);
        solver.analyzePattern(_view.matrix());
    }

    void factorize(const SparseMatrixRM& mat, int offset, int size) override {
        _view.update(mat, offset, size);
        solver.factorize(_view.matrix());
    }

    void solve(double* rhs, int size) override {
        Eigen::Map<Eigen::VectorXd> b(rhs, size);
        _x = solver.solve(b);
        b = _x;
    }

    Eigen::ComputationInfo info() const override {
        return solver.info();
    }

private:
    Eigen::SparseLU<Eigen::SparseMatrix<double>, OrderingType> solver;
    DiagonalBlockView<Eigen::ColMajor> _view;
    Eigen::VectorXd _x;
};

template <typename OrderingType>
class SparseQR : public EigenSolverBase {
public:
    void analyzePattern(const SparseMatrixRM& mat, int offset, int size) override {
        _view.prepare(mat, offset, size);
        //solver.analyzePattern(_view.matrix());
    }

    void factorize(const SparseMatrixRM& mat, int offset, int size) override {
        _view.update(mat, offset, size);
        solver.factorize(_view.matrix());
    }

    void solve(double* rhs, int size) override {
        Eigen::Map<Eigen::VectorXd> b(rhs, size);
        _x = solver.solve(b);
        b = _x;
    }

    Eigen::ComputationInfo info() const override {
        return solver.info();
    }

private:
    Eigen::SparseQR<Eigen::SparseMatrix<double, Eigen::RowMajor>, OrderingType> solver;
    DiagonalBlockView<Eigen::RowMajor> _view;
    Eigen::VectorXd _x;
};

template <typename PreConditioner>
class BiCGSTAB : public EigenSolverBase {
public:
    void analyzePattern(const SparseMatrixRM& mat, int offset, int size) override {
        _view.prepare(mat, offset, size);
        solver.analyzePattern(_view.matrix());
    }

    void factorize(const SparseMatrixRM& mat, int offset, int size) override {

        _view.update(mat, offset, size);
        solver.compute(_view.matrix());

        if (solver.info() != Eigen::Success)
        {
            throw std::runtime_error("BiCGSTAB decomposition failed");
        }
    }

    void solve(double* rhs, int size) override {
        Eigen::Map<Eigen::VectorXd> b(rhs, size);
        _x = solver.solve(b);
        b = _x;
    }

    Eigen::ComputationInfo info() const override {
        return solver.info();
    }

private:
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, PreConditioner> solver;
    DiagonalBlockView<Eigen::RowMajor> _view;
    Eigen::VectorXd _x;
};

template <typename PreConditioner>
class LeastSquaresConjugateGradient : public EigenSolverBase {
public:
    void analyzePattern(const SparseMatrixRM& mat, int offset, int size) override {
        _view.prepare(mat, offset, size);
        //solver.analyzePattern(_view.matrix());
    }

    void factorize(const SparseMatrixRM& mat, int offset, int size) override {

        _view.update(mat, offset, size);
        solver.compute(_view.matrix());

        if (solver.info() != Eigen::Success)
        {
            throw std::runtime_error("LeastSquaresConjugateGradient decomposition failed");
        }
    }

    void solve(double* rhs, int size) override {
        Eigen::Map<Eigen::VectorXd> b(rhs, size);
        _x = solver.solve(b);
        b = _x;
    }

    Eigen::ComputationInfo info() const override {
        return solver.info();
    }

private:
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>, PreConditioner> solver;
    DiagonalBlockView<Eigen::RowMajor> _view;
    Eigen::VectorXd _x;
};

cadet::linalg::EigenSolverBase* setLinearSolver(const std::string solverName);

} // namespace linalg

} // namespace cadet

#endif  // LIBCADET_EigenSolverWrapper_HPP_

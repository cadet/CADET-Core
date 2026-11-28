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
 * Controls the floating point mode of the current thread.
 */

#ifndef LIBCADET_FLOATINGPOINTMODE_HPP_
#define LIBCADET_FLOATINGPOINTMODE_HPP_

#include "cadet/cadetCompilerInfo.hpp"

#if defined(_M_X64) || defined(_M_IX86) || defined(__x86_64__) || defined(__i386__)
	#define CADET_HAS_SSE_CSR 1
	#include <xmmintrin.h>
	#include <pmmintrin.h>
#endif

namespace cadet
{

/**
 * @brief Flushes subnormal floating point values to zero while in scope
 * @details Entries far below the scale of the system carry no information for a linear
 *          solve, but every operation that takes a subnormal operand or produces a
 *          subnormal result costs a microcode assist orders of magnitude more expensive
 *          than the operation itself. A factorization performs a fixed amount of work for
 *          given matrix dimensions, so those assists are pure overhead and can make up the
 *          larger part of its run time.
 *
 *          The mode is set per thread. The previous floating point state is restored on
 *          destruction, so it never leaves the scope this object guards.
 */
class ScopedFlushSubnormals
{
public:
#ifdef CADET_HAS_SSE_CSR
	ScopedFlushSubnormals() CADET_NOEXCEPT : _csr(_mm_getcsr()) { _mm_setcsr(_csr | (_MM_FLUSH_ZERO_ON | _MM_DENORMALS_ZERO_ON)); }
	~ScopedFlushSubnormals() CADET_NOEXCEPT { _mm_setcsr(_csr); }

	ScopedFlushSubnormals(const ScopedFlushSubnormals&) = delete;
	ScopedFlushSubnormals& operator=(const ScopedFlushSubnormals&) = delete;

private:
	unsigned int _csr;
#endif
};

} // namespace cadet

#endif  // LIBCADET_FLOATINGPOINTMODE_HPP_

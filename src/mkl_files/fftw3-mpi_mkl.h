/*******************************************************************************
* Copyright 2010-2019 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
 *
 * Definitions for MPI FFTW3 wrappers to Intel(R) Math Kernel Library (Intel(R) MKL).
 *
 ******************************************************************************
 */

#ifndef FFTW3_MPI_MKL_H
#define FFTW3_MPI_MKL_H

#include "fftw3-mpi.h"

#if defined(MKL_SINGLE)
typedef float real_t;
typedef fftwf_complex complex_t;
#define MPI_PREC MPI_FLOAT
#define MKL_PREC DFTI_SINGLE
#define FFTW_MPI_MANGLE(name) FFTW_MPI_MANGLE_FLOAT(name)
#define FFTW_MANGLE(name) FFTW_MANGLE_FLOAT(name)
#else
typedef double real_t;
typedef fftw_complex complex_t;
#define MPI_PREC MPI_DOUBLE
#define MKL_PREC DFTI_DOUBLE
#define FFTW_MPI_MANGLE(name) FFTW_MPI_MANGLE_DOUBLE(name)
#define FFTW_MANGLE(name) FFTW_MANGLE_DOUBLE(name)
#endif

#include "fftw3_mkl.h"
#include "mkl_cdft.h"

#define WANT_FAST_INPLACE_CLUSTER_FFT 1
/* if WANT_FAST_INPLACE_CLUSTER_FFT set to 1, FFTW3 MPI wrappers internally
 * allocate additional memory(workspace) needed for fast inplace Intel(R) MKL CDFT
 * otherwise, no additional memory is used, though the perfomance would be
 * worse, because of many MPI communications */

#endif /* FFTW3_MPI_MKL_H */


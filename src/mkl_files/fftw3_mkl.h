/*******************************************************************************
* Copyright 2005-2019 Intel Corporation.
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
 * Definitions for FFTW3 wrappers to Intel(R) Math Kernel Library (Intel(R) MKL).
 *
 ******************************************************************************
 */

#ifndef FFTW3_MKL_H
#define FFTW3_MKL_H

#include <stdlib.h>
#include "fftw3.h"
#include "mkl_dfti.h"
#include "mkl_trig_transforms.h"
#include "mkl_service.h"

typedef struct fftw_mkl_plan_s *fftw_mkl_plan;
typedef struct fftw3_mkl_s      fftw3_mkl_s;

/* Plan holder for the wrappers */
struct fftw_mkl_plan_s
{
    DFTI_DESCRIPTOR_HANDLE desc;
    void *io[4];
    MKL_INT *ipar;
    double *dpar;
    float *spar;
    void (*execute) (fftw_mkl_plan p);
    void (*destroy) (fftw_mkl_plan p);
    void *mpi_plan; /* placeholder for FFTW3 MPI Intel(R) MKL plan */
};

/* Global helper structure */
struct fftw3_mkl_s
{
    int verbose;
    int nthreads;
    double timelimit;
    int number_of_user_threads; /* Will be deprecated in nearest future */
    fftw_mkl_plan (*new_plan) (void);
    int default_alignment;
};

FFTW_EXTERN fftw3_mkl_s fftw3_mkl;

#define MKL_MAXRANK 7
#define MKL_ONE     1
#define MKL_RODFT00 413

#define BAD(status) ((status) && !DftiErrorClass((status),DFTI_NO_ERROR))

#ifndef UNUSED
#define UNUSED(p) (void)p
#endif

#endif /* FFTW3_MKL_H */

// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2014: Eric von Lieres¹, Joel Andersson,
//                         Andreas Puettmann¹, Sebastian Schnittert¹,
//                         Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <stdio.h>
#include <stdlib.h>

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* Modified version of NETLIB's dgbmv */
void mydgbmv(char *trans, int *m, int *n, int *kl, int *ku, double *alpha,
        double *a, int *lda, double *x, int *incx, double *beta, double *y,
        int *incy)
{

    /* Local variables */
    double temp;
    int lenx, leny, i, j, k;
    int ix, iy, jx, jy, kx, ky;
    int kup1;

#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set

     up the start points in  X  and  Y. */

    if (0)
    {
        lenx = *n;
        leny = *m;
    }
    else
    {
        lenx = *m;
        leny = *n;
    }
    if (*incx > 0)
    {
        kx = 1;
    }
    else
    {
        kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0)
    {
        ky = 1;
    }
    else
    {
        ky = 1 - (leny - 1) * *incy;
    }

    /*     Sta~rt the operations. In this version the elements of A are
     accessed sequentially with one pass through the band part of A.

     First form  y := beta*y. */

    if (*beta != 1.)
    {
        if (*incy == 1)
        {
            if (*beta == 0.)
            {
                for (i = 1; i <= leny; ++i)
                {
                    Y(i) = 0.;
                }
            }
            else
            {
                for (i = 1; i <= leny; ++i)
                {
                    Y(i) = *beta * Y(i);
                }
            }
        }
        else
        {
            iy = ky;
            if (*beta == 0.)
            {
                for (i = 1; i <= leny; ++i)
                {
                    Y(iy) = 0.;
                    iy += *incy;
                }
            }
            else
            {
                for (i = 1; i <= leny; ++i)
                {
                    Y(iy) = *beta * Y(iy);
                    iy += *incy;
                }
            }
        }
    }
    if (*alpha == 0.)
    {
        return;
    }
    kup1 = *ku + 1;
    if (0)
    {

        /*        Form  y := alpha*A*x + y. */

        jx = kx;
        if (*incy == 1)
        {
            for (j = 1; j <= *n; ++j)
            {
                if (X(jx) != 0.)
                {
                    temp = *alpha * X(jx);
                    k = kup1 - j;
                    /* Computing MAX */
                    /* Computing MIN */
                    for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i)
                    {
                        Y(i) += temp * A(k+i,j);
                    }
                }
                jx += *incx;
            }
        }
        else
        {
            for (j = 1; j <= *n; ++j)
            {
                if (X(jx) != 0.)
                {
                    temp = *alpha * X(jx);
                    iy = ky;
                    k = kup1 - j;
                    /* Computing MAX */
                    /* Computing MIN */
                    for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i)
                    {
                        Y(iy) += temp * A(k+i,j);
                        iy += *incy;
                    }
                }
                jx += *incx;
                if (j > *ku)
                {
                    ky += *incy;
                }
            }
        }
    }
    else
    {

        /*        Form  y := alpha*A'*x + y. */

        jy = ky;
        if (*incx == 1)
        {
            for (j = 1; j <= *n; ++j)
            {
                temp = 0.;
                k = kup1 - j;
                /* Computing MAX */
                /* Computing MIN */
                for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i)
                {
                    temp += A(k+i,j) * X(i);
                }
                Y(jy) += *alpha * temp;
                jy += *incy;
            }
        }
        else
        {
            for (j = 1; j <= *n; ++j)
            {
                temp = 0.;
                ix = kx;
                k = kup1 - j;
                /* Computing MAX */
                /* Computing MIN */
                for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i)
                {
                    temp += A(k+i,j) * X(ix);
                    ix += *incx;
                }
                Y(jy) += *alpha * temp;
                jy += *incy;
                if (j > *ku)
                {
                    kx += *incx;
                }
            }
        }
    }

    /*     End of DGBMV . */

} /* dgbmv_ */

/* Modified version of NETLIB's dgbtrs */
void mydgbtrs(char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab,
        int *ldab, int *piv, double *b, int *ldb, int *info)
{
    int kul = *ku + *kl;

    int i_1, i_2, i_3;
    double d_1;
    double *a = ab;
    double *x = b;

    /*  Solve U'*X = B, overwriting B with X. */
    a += kul;
    for (i_2 = 0; i_2 < kul; ++i_2)
    {
        d_1 = *x;
        for (i_1 = -i_2; i_1 < 0; ++i_1)
        {
            d_1 -= a[i_1] * x[i_1];
        }
        d_1 /= *a;
        *x = d_1;
        a += *ldab;
        x++;
    }
    for (; i_2 < *n; ++i_2)
    {
        d_1 = *x;
        for (i_1 = -kul; i_1 < 0; ++i_1)
        {
            d_1 -= a[i_1] * x[i_1];
        }

        d_1 /= *a;
        *x = d_1;
        a += *ldab;
        x++;
    }

    /*  Solve L'*X = B, overwriting B with X. */
    ab += kul + 1 + (*n - 2) * *ldab;

    for (i_2 = *n - 2; i_2 >= 0; --i_2)
    {
        i_3 = *n - 1 - i_2;
        if (i_3 > *kl)
            i_3 = *kl;
        a = b + i_2 + 1;
        d_1 = b[i_2];
        for (i_1 = 0; i_1 < i_3; ++i_1)
        {
            d_1 -= *a++ * *ab++;
        }
        ab -= *ldab + i_3;
        i_3 = piv[i_2] - 1;
        if (i_3 == i_2)
        {
            b[i_2] = d_1;
        }
        else
        {
            b[i_2] = b[i_3];
            b[i_3] = d_1;
        }
    }

    *info = 0;
}

/* "LAPACK" function for factorization of band matrices */
void mydgbtrf(int *m, int *n, int *kl, int *ku, double *ab, int *ldab,
        int *ipiv, int *info)
{

    /* System generated locals */
    int i__1, i__2, i__3, i__4;
    double d__1, d__2;
    /* Local variables */
    int i__, i, j, k;
    int km, jp, ju, kv;
#define ab_ref(a_1,a_2) ab[(a_2)**ldab + a_1]

    /* Function Body */
    kv = *ku + *kl;

    *info = 0;
    /*     Gaussian eliMINation with partial pivoting

     Set fill-in elements in columns KU+2 to KV to zero. */

    i__1 = MIN(kv,*n);
    for (j = *ku + 1; j < i__1; ++j)
    {
        i__2 = *kl;
        for (i__ = kv - j; i__ < i__2; ++i__)
        {
            ab_ref(i__, j) = 0.;
        }
    }

    /*     JU is the index of the last column affected by the current stage
     of the factorization. */

    ju = 0;

    i__1 = MIN(*m,*n);
    for (j = 0; j < i__1; ++j)
    {

        /*        Set fill-in elements in column J+KV to zero. */

        if (j + kv < *n)
        {
            i__2 = *kl;
            for (i__ = 0; i__ < i__2; ++i__)
            {
                ab_ref(i__,kv) = 0.;
            }
        }

        /*        Find pivot and test for singularity. KM is the number of
         subdiagonal elements in the current column.

         Computing MIN */
        i__2 = *kl, i__3 = *m - j - 1;
        km = MIN(i__2,i__3);

        d__2 = -1.;
        jp = -1;
        for (i = 0; i <= km; i++)
        {
            d__1 = ab[kv + i];
            d__1 *= d__1;
            if (d__1 > d__2)
            {
                d__2 = d__1;
                jp = i;
            }

        }

        *ipiv = jp + 1 + j;
        if (ab[kv + jp] != 0.)
        {
            /* Computing MAX
             Computing MIN */
            i__2 = *n - 1;
            i__3 = j + *ku + jp;
            if (i__3 < i__2)
                i__2 = i__3;
            if (i__2 > ju)
                ju = i__2;

            /*           Apply interchange to columns J to JU. */

            if (jp != 0)
            {
                i__2 = ju - j;
                i__3 = *ldab - 1;
                i__4 = *ldab - 1;
                for (i = 0; i <= i__2; i++)
                {
                    d__1 = ab[kv + jp + i * i__3];
                    ab[kv + jp + i * i__3] = ab[kv + i * i__4];
                    ab[kv + i * i__4] = d__1;
                }
            }

            if (km > 0)
            {

                /*              Compute multipliers. */

                d__1 = 1. / ab[kv];
                i__2 = km + kv;

                for (i = 1 + kv; i <= i__2; i++)
                    ab[i] *= d__1;

                /*              Update trailing submatrix within the band. */

                if (ju > j)
                {
                    i__2 = ju - j;
                    i__3 = km + kv;
                    i__4 = 1 + kv;
                    for (k = 1; k <= i__2; k++)
                    {
                        d__1 = ab_ref(kv-k, k);
                        for (i = i__4; i <= i__3; i++)
                            ab_ref(i-k, k) -= ab[i] * d__1;
                    }

                }
            }
        }
        else
        {

            /*           If pivot is zero, set INFO to the index of the pivot
             unless a zero pivot has already been found. */

            if (*info == 0)
            {
                *info = j + 1;
            }
        }
        ab += *ldab;
        ipiv++;
    }
    /*     return 0; */

}

#undef ab_ref



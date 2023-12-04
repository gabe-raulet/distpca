#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "kiss.h"
#include "cblas.h"
#include "lapacke.h"
#include "svd_utils.h"
#include "mmio_dense.h"
#include "utils.h"

double l2dist(const double *x, const double *y, int n)
{
    double v = 0;

    for (int i = 0; i < n; ++i)
    {
        v += ((x[i] - y[i]) * (x[i] - y[i]));
    }

    return sqrt(v);
}

int svds_naive(double *A, double *Up, double *Sp, double *Vpt, int m, int n, int p)
{
    int r = m < n? m : n;

    assert(A != NULL && Up != NULL && Sp != NULL && Vpt != NULL && r >= p && p >= 1);

    double *S, *U, *Vt, *work;

    work = dalloc(5*r, 0);
    S = dalloc(r, 0);
    U = dalloc(m*r, 0);
    Vt = dalloc(r*n, 0);

    LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, A, m, S, U, m, Vt, r, work);

    memcpy(Up, U, p*m*sizeof(double));
    memcpy(Sp, S, p*sizeof(double));

    double *Vt_ptr = Vt;
    double *Vpt_ptr = Vpt;

    for (int j = 0; j < r; ++j)
    {
        memcpy(Vpt_ptr, Vt_ptr, p*sizeof(double));

        Vpt_ptr += p;
        Vt_ptr += r;
    }

    free(S);
    free(U);
    free(Vt);

    return 0;
}

double* combine_routine(double *Ak_2i_0, double *Vtk_2i_0, double *Ak_2i_1, double *Vtk_2i_1, double *USki, int m, int n, int k, int q, int p)
{
    int b, s, d;
    double *Aki, *Vhtki, *Ski, *Vtki, *W, *tau;

    b = 1 << q;
    s = n / b;
    d = (1 << (k-1)) * s;

    Aki = dalloc(m*(2*p), 0);

    memcpy(&Aki[0],   Ak_2i_0, m*p*sizeof(double));
    memcpy(&Aki[m*p], Ak_2i_1, m*p*sizeof(double));

    Vhtki = dalloc((2*p)*(2*d), 1);

    for (int j = 0; j < d; ++j)
    {
        memcpy(&Vhtki[j*(2*p)], &Vtk_2i_0[j*p], p*sizeof(double));
        memcpy(&Vhtki[(j+d)*(2*p)+p], &Vtk_2i_1[j*p], p*sizeof(double));
    }

    Ski = dalloc(p, 0);
    Vtki = dalloc(p*(2*p), 0);

    svds_naive(Aki, USki, Ski, Vtki, m, 2*p, p);

    for (int j = 0; j < p; ++j)
        for (int i = 0; i < m; ++i)
            USki[i + j*m] *= Ski[j];

    free(Ski);

    W = dalloc((2*d)*p, 0);

    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, 2*d, p, 2*p, 1.0, Vhtki, 2*p, Vtki, p, 0.0, W, 2*d);

    free(Vtki);
    free(Vhtki);

    assert(2*d >= p);
    tau = dalloc(p, 0);

    LAPACKE_dgeqrf(LAPACK_COL_MAJOR, 2*d, p, W, 2*d, tau);
    LAPACKE_dtrtri(LAPACK_COL_MAJOR, 'U', 'N', p, W, 2*d);

    cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, p, 1.0, W, 2*d, USki, m);

    LAPACKE_dorgqr(LAPACK_COL_MAJOR, 2*d, p, p, W, 2*d, tau);

    free(tau);
    free(Aki);

    return W;
}

int combine_node(double *Ak_2i_0, double *Vtk_2i_0, double *Ak_2i_1, double *Vtk_2i_1, double *Ak1_lj, double *Vtk1_lj, int m, int n, int k, int q, int p)
{
    int b, s, d;
    double *W;

    b = 1 << q;
    s = n / b;
    d = (1 << (k-1)) * s;

    W = combine_routine(Ak_2i_0, Vtk_2i_0, Ak_2i_1, Vtk_2i_1, Ak1_lj, m, n, k, q, p);

    for (int j = 0; j < p; ++j)
        for (int i = 0; i < 2*d; ++i)
            Vtk1_lj[j + i*p] = W[i + j*2*d]; /* Vtk1_lj[j,i] = W[i,j]; W is 2d-by-p */

    free(W);

    return 0;
}

int seed_node(double *Ai, double *A1i, double *Vt1i, int m, int n, int q, int p)
{
    int b = 1 << q;
    int s = n / b;

    double *Sp = dalloc(p, 0);

    svds_naive(Ai, A1i, Sp, Vt1i, m, s, p);

    for (int j = 0; j < p; ++j)
        for (int i = 0; i < m; ++i)
            A1i[i + m*j] *= Sp[j];

    free(Sp);
    return 0;
}

int extract_node(double *Aq1_11, double *Vtq1_11, double *Aq1_12, double *Vtq1_12, double *U, double *S, double *Vt, int m, int n, int q, int p)
{
    double *USq, *Qq, *Vtp;

    USq = dalloc(m*p, 0);
    Vtp = dalloc(p*n, 0);

    Qq = combine_routine(Aq1_11, Vtq1_11, Aq1_12, Vtq1_12, USq, m, n, q, q, p);

    svds_naive(USq, U, S, Vtp, m, p, p);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, p, n, p, 1.0, Vtp, p, Qq, n, 0.0, Vt, p );

    free(USq);
    free(Qq);
    free(Vtp);

    return 0;
}



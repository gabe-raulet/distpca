#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <fmi.h>
#include <string>
#include "kiss.h"
#include "cblas.h"
#include "lapacke.h"
#include "mmio_dense.h"
#include "svd_utils.h"
#include "utils.h"
#include "fmi_wrapper.h"
#include "svd_algs_fmi.h"

int compute_errors
(
    const double *A,
    const double *U, const double *Up,
    const double *S, const double *Sp,
    const double *Vt, const double *Vtp,
    int m, int n, int p,
    double errs[4]
);

int main(int argc, char *argv[])
{
    int myrank = atoi(argv[1]);
    int nprocs = atoi(argv[2]);

    if (nprocs == 1 || (nprocs&(nprocs-1)))
    {
        if (!myrank) fprintf(stderr, "[error::main][nprocs=%d] nprocs must be a power of 2 greater than 1\n", nprocs);
        return 1;
    }

    /*
     * NOTE: MUST ASSUME THAT p <= r AND p <= s
     */

    int m; /* (rank[0..nprocs-1]) rows of A */
    int n; /* (rank[0..nprocs-1]) columns of A */
    int r; /* (rank[0..nprocs-1]) = min(m, n) */
    int p; /* (rank[0..nprocs-1]) SVD truncation parameter */
    int s; /* (rank[0..nprocs-1]) columns of Aloc; equals n/nprocs */
    double *A; /* (rank[0]) m-by-n A matrix */
    double *S; /* (rank[0]) r-by-r diagonal S matrix (ground truth) */
    double *U; /* (rank[0]) m-by-r U matrix (ground truth) */
    double *Vt; /* (rank[0]) r-by-n Vt matrix (ground truth) */
    double *Sp; /* (rank[0]) p-by-p diagonal Sp matrix (computed) */
    double *Up; /* (rank[0]) m-by-p Up matrix (computed) */
    double *Vtp; /* (rank[0]) p-by-n Vtp matrix (computed) */
    double *Aloc; /* (rank[myrank]) m-by-s) matrix A[:,myrank*s:(myrank+1)*s] */

    if (argc != 10)
    {
        if (!myrank) fprintf(stderr, "Usage: %s <myrank> <nprocs> <m> <n> <A.mtx> <S.diag> <U.mtx> <Vt.mtx> <p>\n", argv[0]);
        return 1;
    }

    m = atoi(argv[3]);
    n = atoi(argv[4]);
    p = atoi(argv[9]);

    std::string comm_name = std::to_string(std::time(nullptr));
    std::string config_path = "fmi/config/fmi_test.json";
    auto comm = FMI::Communicator(myrank, nprocs, config_path, comm_name);
    comm.hint(FMI::Utils::Hint::fast);
    comm.barrier();

    if (!myrank)
    {
        /*
         * Read A, S, U, and Vt
         */
        char *Afname, *Sfname, *Ufname, *Vtfname;
        int tmp[2];

        Afname = argv[5];
        A = mmread(Afname, tmp, tmp+1);
        assert(tmp[0] == m && tmp[1] == n);
        r = m < n? m : n;

        Sfname = argv[6];
        S = dalloc(r, 0);
        FILE *f = fopen(Sfname, "r");
        for (int i = 0; i < r; ++i)
            fscanf(f, "%lg\n", S+i);
        fclose(f);

        Ufname = argv[7];
        U = mmread(Ufname, tmp, tmp+1);
        assert(m==tmp[0] && r==tmp[1]);

        Vtfname = argv[8];
        Vt = mmread(Vtfname, tmp, tmp+1);
        assert(r==tmp[0] && n==tmp[1]);

        fprintf(stderr, "[sanity::main][m=%d,n=%d,p=%d,nprocs=%d]\n", m, n, p, nprocs);
    }

    s = n / nprocs;
    r = m < n? m : n; /* now everyone has it */

    if (!(m >= 1 && n >= 1) || (m&(m-1)) || (n&(n-1)))
    {
        if (!myrank) fprintf(stderr, "[error::main][m=%d,n=%d] must have m,n >= 1 with n and m both being powers of 2\n", m, n);
        return 1;
    }

    if (p > r || n % nprocs != 0 || p > s)
    {
        if (!myrank) fprintf(stderr, "[error::main][p=%d,min(m,n)=%d,n=%d,nprocs=%d] must have p <= min(m,n), n %% nprocs == 0, and p <= n/nprocs\n", p, r, n, nprocs);
        return 1;
    }

    /*
     * Distribute A to Aloc
     */

    if (!myrank)
    {
        Up = dalloc(m*p, 0);
        Sp = dalloc(p, 0);
        Vtp = dalloc(p*n, 0);
    }
    else
    {
        Up = Sp = Vtp = NULL;
    }

    Aloc = dalloc(m*s, 0);
    fmi_scatter((void*)A, m*n*sizeof(double), (void*)Aloc, m*s*sizeof(double), 0, comm);

    /*
     * Run distributed SVD
     */

    if (svd_fmi(Aloc, Up, Sp, Vtp, m, n, p, 0, myrank, nprocs, comm) != 0)
    {
        return 1;
    }

    /*
     * Compute and report errors
     */

    if (!myrank)
    {
        double errs[4];
        compute_errors(A, U, Up, S, Sp, Vt, Vtp, m, n, p, errs);

        fprintf(stderr, "[main:compute_errors] Aerr=%.18e [normF(A - Up@Sp@Vtp)]\n", errs[0]);
        fprintf(stderr, "[main:compute_errors] Serr=%.18e [normF(S[:p,:p] - Sp)]\n", errs[1]);
        fprintf(stderr, "[main:compute_errors] Uerr=%.18e [normF(U[:,:p]@U[:,:p].T - Up@Up.T)]\n", errs[2]);
        fprintf(stderr, "[main:compute_errors] Verr=%.18e [normF(Vt[:p,:].T@Vt[:p,:] - Vtp.T@Vtp)]\n", errs[3]);
    }

    if (!myrank)
    {
        free(A);
        free(S);
        free(U);
        free(Vt);
    }

    free(Aloc);
    (void)Sp;
    (void)Up;
    (void)Vtp;

    return 0;
}

int compute_errors
(
    const double *A,
    const double *U, const double *Up,
    const double *S, const double *Sp,
    const double *Vt, const double *Vtp,
    int m, int n, int p,
    double errs[4]
)
{
    int r = m < n? m : n;
    int d = m > n? m : n;

    double Aerr, Serr, Uerr, Verr;
    double *mem = dalloc(d*d, 0);

    /*
     * Compute Aerr = normF(A - Up@Sp@Vtp).
     */

    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; ++i)
        {
            double acc = 0;

            for (int k = 0; k < p; ++k)
            {
                acc += Up[i + k*m]*Sp[k]*Vtp[k + j*p];
            }

            mem[i + j*m] = acc - A[i + j*m];
        }

    Aerr = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', m, n, mem, m);

    /*
     * Compute Serr = normF(S[:p,:p] - Sp).
     */

    Serr = 0;
    for (int i = 0; i < p; ++i)
    {
        double delta = S[i] - Sp[i];
        Serr += delta*delta;
    }
    Serr = sqrt(Serr);

    /*
     * Compute Uerr = normF(U[:,:p]@U[:,:p].T - Up@Up.T).
     *
     * U is m-by-r, Up is m-by-p, p <= r.
     */

    for (int j = 0; j < m; ++j)
        for (int i = 0; i < m; ++i)
        {
            double acc = 0;

            for (int k = 0; k < p; ++k)
            {
                acc += U[i + k*m]*U[j + k*m];
                acc -= Up[i + k*m]*Up[j + k*m];
            }

            mem[i + j*m] = acc;
        }

    Uerr = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', m, m, mem, m);

    /*M = Vt.T@Vt - Vtp.T@Vtp n-by-n */

    /*
     * Compute Verr = normF(Vt[:p,:].T@Vt[:p,:] - Vtp.T@Vtp).
     *
     * Vt is r-by-n, Vtp is p-by-n, p <= r.
     */

    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
        {
            double acc = 0;

            for (int k = 0; k < p; ++k)
            {
                acc += Vt[k + i*r]*Vt[k + j*r];
                acc -= Vtp[k + i*p]*Vtp[k + j*p];
            }

            mem[i + j*n] = acc;
        }

    Verr = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', n, n, mem, n);

    errs[0] = Aerr, errs[1] = Serr, errs[2] = Uerr, errs[3] = Verr;

    free(mem);
    return 0;

}

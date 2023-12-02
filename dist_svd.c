#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "distpca.h"

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
    int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (nprocs == 1 || (nprocs&(nprocs-1)))
    {
        if (!myrank) fprintf(stderr, "[error::main][nprocs=%d] nprocs must be a power of 2 greater than 1\n", nprocs);
        MPI_Finalize();
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

    if (argc != 6)
    {
        if (!myrank) fprintf(stderr, "Usage: %s <A.mtx> <S.diag> <U.mtx> <Vt.mtx> <p>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    int intparams[3];

    mpi_timer_t timer;
    double maxtime, proctime;

    mpi_timer_init(&timer, MPI_COMM_WORLD);
    mpi_timer_start(&timer);

    if (!myrank)
    {
        /*
         * Read A, S, U, and Vt
         */
        char *Afname, *Sfname, *Ufname, *Vtfname;
        int tmp[2];

        Afname = argv[1];
        A = mmread(Afname, &m, &n);
        r = m < n? m : n;

        Sfname = argv[2];
        S = malloc(r*sizeof(double));
        FILE *f = fopen(Sfname, "r");
        for (int i = 0; i < r; ++i)
            fscanf(f, "%lg\n", S+i);
        fclose(f);

        Ufname = argv[3];
        U = mmread(Ufname, tmp, tmp+1);
        assert(m==tmp[0] && r==tmp[1]);

        Vtfname = argv[4];
        Vt = mmread(Vtfname, tmp, tmp+1);
        assert(r==tmp[0] && n==tmp[1]);

        p = atoi(argv[5]);

        intparams[0] = m, intparams[1] = n, intparams[2] = p;

        fprintf(stderr, "[sanity::main][m=%d,n=%d,p=%d,nprocs=%d]\n", m, n, p, nprocs);
    }

    MPI_Bcast(intparams, 3, MPI_INT, 0, MPI_COMM_WORLD);

    m = intparams[0];
    n = intparams[1];
    p = intparams[2];
    s = n / nprocs;
    r = m < n? m : n; /* now everyone has it */

    if (!(m >= 1 && n >= 1) || (m&(m-1)) || (n&(n-1)))
    {
        if (!myrank) fprintf(stderr, "[error::main][m=%d,n=%d] must have m,n >= 1 with n and m both being powers of 2\n", m, n);
        MPI_Finalize();
        return 1;
    }

    if (p > r || n % nprocs != 0 || p > s)
    {
        if (!myrank) fprintf(stderr, "[error::main][p=%d,min(m,n)=%d,n=%d,nprocs=%d] must have p <= min(m,n), n %% nprocs == 0, and p <= n/nprocs\n", p, r, n, nprocs);
        MPI_Finalize();
        return 1;
    }

    mpi_timer_stop(&timer);
    mpi_timer_query(&timer, &maxtime, &proctime);

    if (!myrank) fprintf(stderr, "[read_input::main::maxtime=%.5f(s)::proctime=%.5f(s)] finished\n", maxtime, proctime);

    mpi_timer_start(&timer);

    /*
     * Distribute A to Aloc
     */

    if (!myrank)
    {
        Up = malloc(m*p*sizeof(double));
        Sp = malloc(p*sizeof(double));
        Vtp = malloc(p*n*sizeof(double));
    }
    else
    {
        Up = Sp = Vtp = NULL;
    }

    Aloc = malloc(m*s*sizeof(double));
    MPI_Scatter(A, m*s, MPI_DOUBLE, Aloc, m*s, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    mpi_timer_stop(&timer);
    mpi_timer_query(&timer, &maxtime, &proctime);

    if (!myrank) fprintf(stderr, "[distribute_Amat::main::maxtime=%.5f(s)::proctime=%.5f(s)] finished\n", maxtime, proctime);

    mpi_timer_start(&timer);

    /*
     * Run distributed SVD
     */

    if (svd_dist(Aloc, Up, Sp, Vtp, m, n, p, 0, MPI_COMM_WORLD) != 0)
    {
        MPI_Finalize();
        return 1;
    }

    mpi_timer_stop(&timer);
    mpi_timer_query(&timer, &maxtime, &proctime);

    if (!myrank) fprintf(stderr, "[svd_dist::main::maxtime=%.5f(s)::proctime=%.5f(s)] finished\n", maxtime, proctime);

    mpi_timer_start(&timer);

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

    mpi_timer_stop(&timer);
    mpi_timer_query(&timer, &maxtime, &proctime);

    if (!myrank) fprintf(stderr, "[compute_errors::main::maxtime=%.5f(s)::proctime=%.5f(s)] finished\n", maxtime, proctime);

    /*
     * Clean up
     */

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

    MPI_Finalize();
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
    double *mem = malloc(d*d*sizeof(double));

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

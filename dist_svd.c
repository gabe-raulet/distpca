#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "distpca.h"

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

    int m; /* (rank[0..nprocs-1]) rows of A */
    int n; /* (rank[0..nprocs-1]) columns of A */
    int p; /* (rank[0..nprocs-1]) SVD truncation parameter */
    int s; /* (rank[0..nprocs-1]) columns of Aloc; equals n/nprocs */
    double *A; /* (rank[0]) m-by-n A matrix */
    double *S; /* (rank[0]) n-by-n diagonal S matrix (ground truth) */
    double *U; /* (rank[0]) m-by-n U matrix (ground truth) */
    double *Vt; /* (rank[0]) n-by-n Vt matrix (ground truth) */
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

    if (!myrank)
    {
        /*
         * Read A, S, U, and Vt
         */
        char *Afname, *Sfname, *Ufname, *Vtfname;
        int tmp[2];

        Afname = argv[1];
        A = mmread(Afname, &m, &n);

        Sfname = argv[2];
        S = malloc(n*sizeof(double));
        FILE *f = fopen(Sfname, "r");
        for (int i = 0; i < n; ++i)
            fscanf(f, "%lg\n", S+i);
        fclose(f);

        Ufname = argv[3];
        U = mmread(Ufname, tmp, tmp+1);
        assert(m==tmp[0] && n==tmp[1]);

        Vtfname = argv[4];
        Vt = mmread(Vtfname, tmp, tmp+1);
        assert(n==tmp[0] && n==tmp[1]);

        p = atoi(argv[5]);

        intparams[0] = m, intparams[1] = n, intparams[2] = p;

        fprintf(stderr, "[sanity::main][m=%d,n=%d,p=%d,nprocs=%d]\n", m, n, p, nprocs);
    }

    MPI_Bcast(intparams, 3, MPI_INT, 0, MPI_COMM_WORLD);

    m = intparams[0];
    n = intparams[1];
    p = intparams[2];
    s = n / nprocs;

    if (!(1 <= n && n <= m) || (m&(m-1)) || (n&(n-1)))
    {
        if (!myrank) fprintf(stderr, "[error::main][m=%d,n=%d] must have 1 <= n <= m with n and m both being powers of 2\n", m, n);
        MPI_Finalize();
        return 1;
    }

    if (p > n || n % nprocs != 0 || p > s)
    {
        if (!myrank) fprintf(stderr, "[error::main][p=%d,n=%d,nprocs=%d] must have p <= n, n %% nprocs == 0, and p <= n/nprocs\n", p, n, nprocs);
        MPI_Finalize();
        return 1;
    }

    /*
     * Distribute A to Aloc
     */

    if (!myrank)
    {
        Up = malloc(m*p*sizeof(double));
        Sp = malloc(p*sizeof(double));
        Vtp = malloc(p*n*sizeof(double));
    }

    Aloc = malloc(m*s*sizeof(double));
    MPI_Scatter(A, m*s, MPI_DOUBLE, Aloc, m*s, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*
     * Run distributed SVD
     */

    /*
     * Compute and report errors
     */

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


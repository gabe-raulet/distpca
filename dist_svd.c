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

    int m; /* (rank[0..nprocs-1]) rows of A */
    int n; /* (rank[0..nprocs-1]) columns of A */
    int p; /* (rank[0..nprocs-1]) SVD truncation parameter */
    double *A; /* (rank[0]) m-by-n A matrix */
    double *S; /* (rank[0]) n-by-n diagonal S matrix (ground truth) */
    double *U; /* (rank[0]) m-by-n U matrix (ground truth) */
    double *Vt; /* (rank[0]) n-by-n Vt matrix (ground truth) */
    double *Sp; /* (rank[0]) p-by-p diagonal Sp matrix (computed) */
    double *Up; /* (rank[0]) m-by-p Up matrix (computed) */
    double *Vtp; /* (rank[0]) p-by-n Vtp matrix (computed) */
    double *Aloc; /* (rank[myrank]) m-by-(n/nprocs) matrix A[:,myrank*(n/nprocs):(myrank+1)*(n/nprocs)] */

    /*
     * Read inputs.
     */

    if (!myrank)
    {
        /*
         * Read A, S, U, and Vt
         */
    }

    /*
     * Distribute A to Aloc
     */

    /*
     * Run distributed SVD
     */

    /*
     * Compute and report errors
     */

    /*
     * Clean up
     */

    MPI_Finalize();
    return 0;
}


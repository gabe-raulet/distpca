#ifndef SVD_ALGS_H_
#define SVD_ALGS_H_

#include <mpi.h>

int svd_serial
(
    double const *A, /* input m-by-n matrx */
    double *Up, /* output m-by-p matrix */
    double *Sp, /* output p-by-p diagonal matrix */
    double *Vtp, /* output p-by-n matrix */
    int m, /* rows of A */
    int n, /* columns of A */
    int p, /* rank approximation */
    int b /* number of seed nodes in binary topology */
);

int svd_dist
(
    double const *Aloc, /* (rank[myrank]) input m-by-(n/nprocs) matrix */
    double *Up, /* (rank[root]) output m-by-p matrix */
    double *Sp, /* (rank[root]) output p-by-p diagonal matrix */
    double *Vtp, /* (rank[root]) output p-by-n matrix */
    int m, /* rows of A */
    int n, /* columns of A */
    int p, /* rank approximation */
    int root, /* root rank */
    MPI_Comm comm /* communicator */
);

#endif

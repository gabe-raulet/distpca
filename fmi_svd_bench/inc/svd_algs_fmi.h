#ifndef SVD_ALGS_FMI_H_
#define SVD_ALGS_FMI_H_

#include <fmi.h>

int svd_fmi
(
    double *Aloc, /* (rank[myrank]) input m-by-(n/nprocs) matrix */
    double *Up, /* (rank[root]) output m-by-p matrix */
    double *Sp, /* (rank[root]) output p-by-p diagonal matrix */
    double *Vtp, /* (rank[root]) output p-by-n matrix */
    int m, /* rows of A */
    int n, /* columns of A */
    int p, /* rank approximation */
    int root, /* root rank */
    int myrank,
    int nprocs,
    FMI::Communicator& comm /* FMI communicator */
);

#endif

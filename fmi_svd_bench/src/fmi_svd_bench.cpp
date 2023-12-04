#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "kiss.h"
#include "cblas.h"
#include "lapacke.h"
#include "mmio_dense.h"
#include "svd_utils.h"
#include "utils.h"
#include "fmi_wrapper.h"
#include "svd_algs_fmi.h"

int main(int argc, char *argv[])
{
    int myrank = atoi(argv[1]);
    int nprocs = atoi(argv[2]);

    int m, n, p, seed;
    double *Aloc, *Up, *Sp, *Vtp;

    if (argc != 7)
    {
        if (!myrank) fprintf(stderr, "usage: %s <myrank> <nprocs> <m> <n> <p> <seed>\n", argv[0]);
        return 1;
    }

    std::string comm_name = std::to_string(std::time(nullptr));
    std::string config_path = "fmi/config/fmi_test.json";
    auto comm = FMI::Communicator(myrank, nprocs, config_path, comm_name);
    comm.hint(FMI::Utils::Hint::fast);
    comm.barrier();

    m = atoi(argv[3]);
    n = atoi(argv[4]);
    p = atoi(argv[5]);
    seed = atoi(argv[6]);

    double maxtime, proctime, elapsed;

    kiss_seed(myrank*seed);

    int s = n / nprocs;

    Aloc = (double*)malloc(m*s*sizeof(double));

    for (int i = 0; i < m*s; ++i)
    {
        Aloc[i] = 10*kiss_unirandf();
    }

    if (!myrank)
    {
        Up = (double*)malloc(m*p*sizeof(double));
        Sp = (double*)malloc(p*sizeof(double));
        Vtp = (double*)malloc(p*n*sizeof(double));
    }
    else
    {
        Up = Sp = Vtp = NULL;
    }

    if (svd_fmi(Aloc, Up, Sp, Vtp, m, n, p, 0, myrank, nprocs, comm) != 0)
    {
        return 1;
    }

    if (!myrank)
    {
        free(Up);
        free(Sp);
        free(Vtp);
    }

    free(Aloc);

    return 0;
}

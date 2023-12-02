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
#include "mmio.h"
#include "mmio_dense.h"
#include "svd_utils.h"
#include "utils.h"
#include "fmi_wrapper.h"
#include "svd_algs_fmi.h"

int p = 10;
char *Afname, *PCfname, *expvarfname;

int usage(char *argv[]);
double* mmread_centered(char const *fname, int *n, int *d);

int main(int argc, char *argv[])
{
    int myrank = atoi(argv[1]);
    int nprocs = atoi(argv[2]);

    if (nprocs == 1 || (nprocs&(nprocs-1)))
    {
        if (!myrank) fprintf(stderr, "[error::main][nprocs=%d] nprocs must be a power of 2 greater than 1\n", nprocs);
        return 1;
    }

    if (argc != 9)
    {
        if (!myrank) fprintf(stderr, "Usage: %s <myrank> <nprocs> <n> <d> <A.mtx> <PC.mtx> <expvar.txt> <p>\n", argv[0]);
        return 1;
    }

    int n, d;
    double *A, *Aloc;

    n = atoi(argv[3]);
    d = atoi(argv[4]);

    Afname = argv[5];
    PCfname = argv[6];
    expvarfname = argv[7];
    p = atoi(argv[8]);

    std::string comm_name = std::to_string(std::time(nullptr));
    std::string config_path = "fmi/config/fmi_test.json";
    auto comm = FMI::Communicator(myrank, nprocs, config_path, comm_name);
    comm.hint(FMI::Utils::Hint::fast);
    comm.barrier();

    //mpi_timer_t timer;
    //double maxtime, proctime;

    //mpi_timer_init(&timer, MPI_COMM_WORLD);
    //mpi_timer_start(&timer);

    if (!myrank)
    {
        int tmp[2];
        A = mmread_centered(Afname, tmp, tmp+1);
        assert(tmp[0] == n && tmp[1] == d);
    }

    int s = d / nprocs;

    //if (!myrank)
    //{
        //if ((n&(n-1)) || (d&(d-1)))
        //{
            //if (!myrank) fprintf(stderr, "[error::main][n=%d,d=%d] must have d and n both being powers of 2\n", n, d);
            //return 1;
        //}

        //if (!(1 <= d && d <= n) || (n&(n-1)) || (d&(d-1)))
        //{
            //if (!myrank) fprintf(stderr, "[error::main][n=%d,d=%d] must have 1 <= d <= n with d and n both being powers of 2\n", n, d);
            //return 1;
        //}

        //if (p > d || d % nprocs != 0 || p > s)
        //{
            //if (!myrank) fprintf(stderr, "[error::main][p=%d,d=%d,nprocs=%d] must have p <= d, d %% nprocs == 0, and p <= d/nprocs\n", p, d, nprocs);
            //return 1;
        //}
    //}

    //mpi_timer_stop(&timer);
    //mpi_timer_query(&timer, &maxtime, &proctime);

    //if (!myrank) fprintf(stderr, "[read_input::main::maxtime=%.5f(s)::proctime=%.5f(s)] finished\n", maxtime, proctime);

    //mpi_timer_start(&timer);

    double *Up, *Sp, *Vtp;

    if (!myrank)
    {
        Up = dalloc(n*p, 0);
        Sp = dalloc(p, 0);
        Vtp = dalloc(p*d, 0);
    }
    else
    {
        Up = Sp = Vtp = NULL;
    }

    Aloc = dalloc(n*s, 0);
    /*MPI_Scatter(A, n*s, MPI_DOUBLE, Aloc, n*s, MPI_DOUBLE, 0, MPI_COMM_WORLD);*/
    fmi_scatter((void*)A, n*d*sizeof(double), (void*)Aloc, n*s*sizeof(double), 0, comm);

    /*mpi_timer_stop(&timer);*/
    /*mpi_timer_query(&timer, &maxtime, &proctime);*/

    /*if (!myrank) fprintf(stderr, "[distribute_Amat::main::maxtime=%.5f(s)::proctime=%.5f(s)] finished\n", maxtime, proctime);*/

    /*mpi_timer_start(&timer);*/

    if (svd_fmi(Aloc, Up, Sp, Vtp, n, d, p, 0, myrank, nprocs, comm) != 0)
    {
        return 1;
    }

    /*mpi_timer_stop(&timer);*/
    /*mpi_timer_query(&timer, &maxtime, &proctime);*/

    /*if (!myrank) fprintf(stderr, "[svd_dist::main::maxtime=%.5f(s)::proctime=%.5f(s)] finished\n", maxtime, proctime);*/

    /*mpi_timer_start(&timer);*/

    if (!myrank)
    {
        for (int i = 0; i < p; ++i)
        {
            Sp[i] = (Sp[i]*Sp[i])/(n-1.0);
        }

        FILE *f = fopen(expvarfname, "w");
        for (int i = 0; i < p; ++i)
            fprintf(f, "%.18e\n", Sp[i]);
        fclose(f);

        mmwrite(PCfname, Vtp, p, d);
    }

    /*mpi_timer_stop(&timer);*/
    /*mpi_timer_query(&timer, &maxtime, &proctime);*/

    /*if (!myrank) fprintf(stderr, "[write_files::main::maxtime=%.5f(s)::proctime=%.5f(s)] finished\n", maxtime, proctime);*/

    if (!myrank)
    {
        free(A);
        free(Up);
        free(Sp);
        free(Vtp);
    }
    free(Aloc);

    return 0;
}

double* mmread_centered(char const *fname, int *n, int *d)
{
    MM_typecode matcode;
    int N, D, i, j;
    double *A;
    FILE *f;

    f = fopen(fname, "r");
    mm_read_banner(f, &matcode);
    assert(mm_is_dense(matcode));
    mm_read_mtx_array_size(f, &N, &D);
    A = dalloc(N*D, 0);
    assert(A != NULL);

    for (j = 0; j < D; ++j)
    {
        double *col = &A[j*N];
        double mean = 0.;

        for (i = 0; i < N; ++i)
        {
            fscanf(f, "%lg", &col[i]);
            mean += col[i];
        }

        mean /= N;

        for (i = 0; i < N; ++i)
            col[i] -= mean;
    }

    fclose(f);

    *n = N;
    *d = D;

    return A;
}

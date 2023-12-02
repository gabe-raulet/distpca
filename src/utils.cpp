#include "utils.h"
#include "kiss.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

int log2i(int v)
{
    int x = 0;
    while (v >>= 1) ++x;
    return x;
}

double* dalloc(size_t n, int clear)
{
    return (double*)(clear? calloc(n, sizeof(double)) : malloc(n*sizeof(double)));
}

static int iseed_prv[4];
static int initialized = 0;

static int iseed_init_prv()
{
    if (initialized) return 1;
    for (int i = 0; i < 4; ++i) iseed_prv[i] = (kiss_rand() % 4096);
    iseed_prv[3] &= (iseed_prv[3]^1); /* iseed_prv[3] must be odd */
    initialized = 1;
    return 0;
}


int iseed_init_dev()
{
    kiss_init();
    iseed_init_prv();
    return 0;
}

int iseed_init_usr(int seed)
{
    if (seed <= 0) kiss_init();
    else kiss_seed((uint32_t)seed);
    iseed_init_prv();
    return 0;
}

int iseed_get(int iseed[4])
{
    memcpy(iseed, iseed_prv, 4*sizeof(int));
    return 0;
}

//typedef struct
//{
//    int isroot;
//    double telapsed;
//    MPI_Comm comm;
//} mpi_timer_t;

//int mpi_timer_init(mpi_timer_t *timer, MPI_Comm comm)
//{
    //int myrank;

    //MPI_Comm_rank(comm, &myrank);
    //timer->isroot = (myrank == 0);
    //timer->telapsed = 0;
    //timer->comm = comm;

    //return 0;
//}

//int mpi_timer_start(mpi_timer_t *timer)
//{
    //MPI_Barrier(timer->comm);
    //timer->telapsed = -MPI_Wtime();
    //return 0;
//}

//int mpi_timer_stop(mpi_timer_t *timer)
//{
    //timer->telapsed += MPI_Wtime();
    //return 0;
//}

//int mpi_timer_query(mpi_timer_t *timer, double *maxtime, double *proctime)
//{
    //MPI_Reduce(&timer->telapsed, maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, timer->comm);
    //MPI_Reduce(&timer->telapsed, proctime, 1, MPI_DOUBLE, MPI_SUM, 0, timer->comm);
    //return 0;
//}

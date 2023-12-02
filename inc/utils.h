#ifndef UTILS_H_
#define UTILS_H_

#include <stddef.h>

int log2i(int v);
int iseed_init_dev();
int iseed_init_usr(int seed);
int iseed_get(int iseed[4]);

double* dalloc(size_t n, int clear);

#include <mpi.h>

typedef struct
{
    int isroot;
    double telapsed;
    MPI_Comm comm;
} mpi_timer_t;

int mpi_timer_init(mpi_timer_t *timer, MPI_Comm comm);
int mpi_timer_start(mpi_timer_t *timer);
int mpi_timer_stop(mpi_timer_t *timer);
int mpi_timer_query(mpi_timer_t *timer, double *maxtime, double *proctime);

#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "distpca.h"
#include "mmio.h"

int p = 10;
char *fname = NULL;

int usage(char *argv[]);
double* mmread_centered(char const *fname, int *n, int *d);

int main(int argc, char *argv[])
{
    int c;

    while ((c = getopt(argc, argv, "p:h")) >= 0)
    {
        if      (c == 'h') return usage(argv);
        else if (c == 'p') p = atoi(optarg);
    }

    if (optind >= argc)
    {
        fprintf(stderr, "error: missing <input.mtx>\n");
        return usage(argv);
    }

    int n, d;
    double *A;

    A = mmread_centered(argv[optind], &n, &d);

    mmwrite("Acentered.mtx", A, n, d);

    free(A);
    return 0;
}

int usage(char *argv[])
{
    fprintf(stderr, "Usage: %s [options] <input.mtx>\n", argv[0]);
    fprintf(stderr, "Options: -p INT   number of principal components [%d]\n", p);
    fprintf(stderr, "         -h       help message\n");
    return 1;
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
    A = malloc(N*D*sizeof(double));
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

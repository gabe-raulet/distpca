#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "kiss.h"
#include "cblas.h"
#include "lapacke.h"
#include "mmio_dense.h"
#include "svd_algs.h"
#include "utils.h"

typedef struct
{
    int m; /* number of rows */
    int n; /* number of columns */
    int mode; /* DLATMS mode (1..5) or < 0 for cond, cond/dmax, cond/dmax^2, ... */
    int seed; /* random seed */
    double cond; /* DLATMS cond (mode 1..5) or condition number (mode < 0) */
    double dmax; /* DLATMS dmax (mode 1..5) or damping factor (mode < 0) */
    char *label; /* generate files {A,U,Vt}_{label}.mtx and S_{label}.diag */
} params_t;

int params_init(params_t *ps);
int usage(char *argv[], const params_t *ps);
int parse_params(int argc, char *argv[], params_t *ps);
int param_check(const params_t *ps);

int main(int argc, char *argv[])
{
    params_t ps;

    if (parse_params(argc, argv, &ps) != 0)
        return 1;

    iseed_init_usr(ps.seed);

    if (param_check(&ps) != 0)
        return 1;

    return 0;
}

int params_init(params_t *ps)
{
    ps->m = 256;
    ps->n = 128;
    ps->mode = -1;
    ps->cond = 100.f;
    ps->dmax = 2.f;
    ps->seed = -1;
    ps->label = NULL;
    return 0;
}

int usage(char *argv[], const params_t *ps)
{
    fprintf(stderr, "Usage: %s [options]\n", argv[0]);
    fprintf(stderr, "Options: -m INT    rows of test matrix [%d]\n", ps->m);
    fprintf(stderr, "         -n INT    columns of test matrix [%d]\n", ps->n);
    fprintf(stderr, "         -u INT    matrix generation mode [%d]\n", ps->mode);
    fprintf(stderr, "         -c FLOAT  DLATMS cond or condition number [%.3e]\n", ps->cond);
    fprintf(stderr, "         -d FLOAT  DLATMS dmax or damping factor [%.3e]\n", ps->dmax);
    fprintf(stderr, "         -s INT    RNG seed [%d]\n", ps->seed);
    fprintf(stderr, "         -o STR    output label (required)\n");
    fprintf(stderr, "         -h        help message\n");
    return 1;
}

int parse_params(int argc, char *argv[], params_t *ps)
{
    int c;
    params_init(ps);

    while ((c = getopt(argc, argv, "m:n:u:c:d:s:o:h")) >= 0)
    {
        if      (c == 'm') ps->m = atoi(optarg);
        else if (c == 'n') ps->n = atoi(optarg);
        else if (c == 'u') ps->mode = atoi(optarg);
        else if (c == 'c') ps->cond = atof(optarg);
        else if (c == 'd') ps->dmax = atof(optarg);
        else if (c == 's') ps->seed = atoi(optarg);
        else if (c == 'o') ps->label = optarg;
        else if (c == 'h')
        {
            params_init(ps);
            return usage(argv, ps);
        }
    }

    if (!ps->label)
    {
        fprintf(stderr, "error: missing required -o parameter\n");
        return usage(argv, ps);
    }

    return 0;
}

int param_check(const params_t *ps)
{
    int m = ps->m, n = ps->n, mode = ps->mode;
    double cond = ps->cond, dmax = ps->dmax;

    if (!(1 <= n && n <= m) || (m&(m-1)) || (n&(n-1)))
    {
        fprintf(stderr, "[error::param_check][m=%d,n=%d] must have 1 <= n <= m with n and m both being powers of 2\n", m, n);
        return 1;
    }

    if (mode == 0 || mode > 5)
    {
        fprintf(stderr, "[error::param_check][mode=%d] mode must be an integer between 1 and 5 inclusive or it must be negative\n", mode);
        return 1;
    }

    if (mode != 0 && (cond < 1 || dmax <= 0))
    {
        fprintf(stderr, "[error::param_check][mode=%d,cond=%.5e,dmax=%.5e] must have cond >= 1 and dmax > 0 when mode=1,2,..,5\n", mode, cond, dmax);
        return 1;
    }

    if (mode < 0 && (cond <= 0 || dmax < 1))
    {
        fprintf(stderr, "[error::param_check][mode=%d,cond=%.5e,damp=%.5e] must have cond > 0 and damp >= 1 when mode < 0\n", mode, cond, dmax);
        return 1;
    }

    return 0;
}

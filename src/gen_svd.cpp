#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "distpca.h"

typedef struct
{
    int m; /* number of rows */
    int n; /* number of columns */
    int mode; /* DLATMS mode (1..5) or < 0 for cond, cond/dmax, cond/dmax^2, ... */
    int seed; /* random seed */
    double cond; /* DLATMS cond (mode 1..5) or condition number (mode < 0) */
    double dmax; /* DLATMS dmax (mode 1..5) or damping factor (mode < 0) */
    char *label; /* generate files {label}_{A,U,Vt}.mtx and {label}_S.diag */
} params_t;

int params_init(params_t *ps);
int usage(char *argv[], const params_t *ps);
int parse_params(int argc, char *argv[], params_t *ps);
int param_check(const params_t *ps);

int gen_test_mat(double *A, double *S, int m, int n, int mode, double cond, double dmax, int iseed[4]);
int gen_uv_mats(const double *A, double *S, double *U, double *Vt, int m, int n);

int main(int argc, char *argv[])
{
    params_t ps;

    if (parse_params(argc, argv, &ps) != 0)
        return 1;

    iseed_init_usr(ps.seed);

    if (param_check(&ps) != 0)
        return 1;

    int iseed[4];
    int m=ps.m, n=ps.n, mode=ps.mode;
    double cond=ps.cond, dmax=ps.dmax;
    char *label=ps.label;

    assert(label != NULL);

    int r = m < n? m : n;

    double *A = dalloc(m*n, 0);
    double *S = dalloc(r, 0);

    iseed_get(iseed);
    gen_test_mat(A, S, m, n, mode, cond, dmax, iseed);

    double *Scheck = dalloc(r, 0);
    double *U = dalloc(m*r, 0);
    double *Vt = dalloc(r*n, 0);

    gen_uv_mats(A, Scheck, U, Vt, m, n);

    int show = r < 5? r : 5;
    fprintf(stderr, "[gen_svd:main] diag(S)[0..%d] = ", show-1);
    for (int i = 0; i < show; ++i) fprintf(stderr, "%.3e,", S[i]);
    fprintf(stderr, "%s\n", r < 5? "" : "...");

    fprintf(stderr, "[gen_svd:main] err=%.18e (DLATMS-vs-DGESVD) [S :: singular values]\n", l2dist(S, Scheck, r));
    free(Scheck);

    char fname[1024];

    snprintf(fname, 1024, "%s_A.mtx", label);
    mmwrite(fname, A, m, n);

    snprintf(fname, 1024, "%s_U.mtx", label);
    mmwrite(fname, U, m, r);

    snprintf(fname, 1024, "%s_Vt.mtx", label);
    mmwrite(fname, Vt, r, n);

    snprintf(fname, 1024, "%s_S.diag", label);

    FILE *f = fopen(fname, "w");
    for (int i = 0; i < r; ++i)
        fprintf(f, "%.18e\n", S[i]);
    fclose(f);

    free(A);
    free(S);
    free(U);
    free(Vt);

    return 0;
}

int gen_uv_mats(const double *A, double *S, double *U, double *Vt, int m, int n)
{
    int r = m < n? m : n;

    double *Al = dalloc(m*n, 0);
    double *work = dalloc(5*r, 0);

    memcpy(Al, A, m*n*sizeof(double));

    LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, Al, m, S, U, m, Vt, r, work);

    free(Al);
    free(work);

    return 0;
}

int gen_test_mat(double *A, double *S, int m, int n, int mode, double cond, double dmax, int iseed[4])
{
    int r = m < n? m : n;

    if (mode < 0)
    {
        S[0] = cond;

        for (int i = 1; i < r; ++i)
            S[i] = S[i-1] / dmax;

        mode = 0;
    }
    else
    {
        assert(1 <= mode && mode <= 5);
    }

    LAPACKE_dlatms(LAPACK_COL_MAJOR, m, n, 'U', iseed, 'N', S, mode, cond, dmax, m, n, 'N', A, m);

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

    if (!(m >= 1 && n >= 1) || (m&(m-1)) || (n&(n-1)))
    {
        fprintf(stderr, "[error::param_check][m=%d,n=%d] must have m,n >= 1 with n and m both being powers of 2\n", m, n);
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

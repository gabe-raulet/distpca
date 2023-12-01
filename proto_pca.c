#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "distpca.h"

int p = 10;
char *fname = NULL;

int usage(char *argv[]);

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

    return 0;
}

int usage(char *argv[])
{
    fprintf(stderr, "Usage: %s [options] <input.mtx>\n", argv[0]);
    fprintf(stderr, "Options: -p INT   number of principal components [%d]\n", p);
    fprintf(stderr, "         -h       help message\n");
    return 1;
}

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


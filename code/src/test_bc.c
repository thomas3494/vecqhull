#include <stdio.h>
#include <stdlib.h>

#include "blockcyclic.h"

int main(int argc, char **argv)
{
    if (argc != 4) {
        printf("Usage: %s <n> <block> <nthreads>\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n = atol(argv[1]);
    size_t block = atol(argv[2]);
    unsigned int nthreads = atoi(argv[3]);

    double *x = (double *)malloc(n * sizeof(double));    

    for (unsigned int t = 0; t < nthreads; t++) {
        BC_index i = {t * block, 0};
        while (i.k + i.j < n) {
            x[i.k + i.j] = t;
            BC_Add(&i, 1, block, nthreads);
        }
    }

    printf("Block cyclic distribution with n = %zu, b = %zu, nt = %x.\n",
                n, block, nthreads);
    for (size_t i = 0; i < n; i++) {
        printf("%d", (int)x[i]);
    }
    printf("\n");

    for (unsigned int t = 0; t < nthreads; t++) {
        BC_index i = BC_Upper(t, n / 2, block, nthreads);
        printf("BC_Upper(%zu), thread %d = %zu\n", 
                    n / 2, t, i.k + i.j);
    }

    printf("Copy thread 0 to thread 1\n");
    BC_index low0  = BC_Upper(0, 0, block, nthreads);
    BC_index low1  = BC_Upper(1, 0, block, nthreads);
    BC_index high1 = BC_Upper(1, n, block, nthreads);
    
    printf("Copy from %zu + %zu to %zu + %zu\n", low0.k, low0.j, low1.k, low1.j);

    size_t count = BC_Dist(low1, high1, nthreads);
    BC_Copy(x, low1, x, low0, count, block, nthreads);
    for (size_t i = 0; i < n; i++) {
        printf("%d", (int)x[i]);
    }
    printf("\n");

    free(x);

    return EXIT_SUCCESS;
}

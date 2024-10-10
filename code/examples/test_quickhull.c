#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vqhull.h>

#ifdef RAPL_ENERGY
#include <rapl_energy.h>
#endif

#include "util.h"

int main(int argc, char **argv)
{
    if (argc == 1 || argc > 5) {
        printf("Usage: Provide points on stdin in PBBS format.\n"
               "s will print summary of results to stderr\n"
               "b accepts binary input\n"
               "m will print runtime in seconds to stdout\n"
               "p will print points to stdout");
        return EXIT_FAILURE;
    }

    bool summary = false;
    bool print = false;
    bool binary = false;
    bool measure = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "s") == 0) {
            summary = true;
        } else if (strcmp(argv[i], "p") == 0) {
            print = true;
        } else if (strcmp(argv[i], "b") == 0) {
            binary = true;
        } else if (strcmp(argv[i], "m") == 0) {
            measure = true;
        } else {
            printf("%s is not a valid option\n", argv[0]);
            return EXIT_FAILURE;
        }
    }

    size_t n;
    Points P = (binary) ? input_b(&n) : input(&n);

    #ifdef RAPL_ENERGY
    struct RaplEnergy *rapl;
    rapl = rapl_start();
    #endif

    double start = wtime();
    size_t count = VecQuickhull(n, P.x, P.y);
    double stop = wtime();

    #ifdef RAPL_ENERGY
    struct RaplElapsed *elapsed;
    elapsed = rapl_elapsed(rapl);
    #endif

    double duration = stop - start;

    if (summary) {
        fprintf(stderr, "We found %zu elements on the hull\n", count);
        fprintf(stderr, "Runtime in ms: %lf\n", duration * 1e3);
        fprintf(stderr, "Bandwidth in MB/s: %lf\n",
                    n * 2 * sizeof(double) / duration / 1e6);
    }

    if (measure) {
        printf("%lf", duration);
        #ifdef RAPL_ENERGY
        double total_energy = 0;
        for (uintptr_t i = 0; i < elapsed->len; i++) {
            total_energy += elapsed->energy[i];
        }
        printf(" %lf", total_energy);
        #endif
        printf("\n");
    }

    #ifdef RAPL_ENERGY
    elapsed_free(elapsed);
    rapl_free(rapl);
    #endif

    if (print) {
        PrintPoints(count, P);
    }

    free(P.x);
    free(P.y);

    return EXIT_SUCCESS;
}

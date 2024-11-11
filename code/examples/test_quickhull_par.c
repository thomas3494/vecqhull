#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vqhull.h>

#ifdef ENERGY
#include <energy_bench.h>
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

    #ifdef ENERGY
    EnergyInfo *start = start_energy_measure();
    #else
    double start = wtime();
    #endif

    size_t count = VecQuickhullP(n, P.x, P.y);

    #ifdef ENERGY
    EnergyResult *stop = stop_energy_measure(start);
    double duration = energy_duration(&stop);
    #else
    double stop = wtime();
    double duration = stop - start;
    #endif

    if (summary) {
        fprintf(stderr, "We found %zu elements on the hull\n", count);
        fprintf(stderr, "Runtime in ms: %lf\n", duration * 1e3);
        fprintf(stderr, "Bandwidth in MB/s: %lf\n",
                    n * 2 * sizeof(double) / duration / 1e6);
    }

    if (measure) {
        #ifdef ENERGY
        print_energy_results(&stop);
        #else
        printf("%lf\n", duration);
        #endif
    }

    if (print) {
        PrintPoints(count, P);
    }

    free(P.x);
    free(P.y);
    #ifdef ENERGY
    free_energy_results(stop);
    #endif

    return EXIT_SUCCESS;
}

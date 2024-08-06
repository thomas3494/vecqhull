#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "quickhull.h"

int main(int argc, char **argv)
{
    if (argc == 1 || argc > 4) {
        printf("Usage: Provide points on stdin in PBBS format.\n"
               "s will print summary of results to stderr\n"
               "b will print runtime in seconds to stdout\n"
               "p will print points to stdout");
        return EXIT_FAILURE;
    }

    bool summary = false;
    bool print = false;
    bool bench = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "s") == 0) {
            summary = true;
        } else if (strcmp(argv[i], "p") == 0) {
            print = true;
        } else if (strcmp(argv[i], "b") == 0) {
            bench = true;
        } else {
            printf("%s is not a valid option\n", argv[0]);
            return EXIT_FAILURE;
        }
    }

    size_t n;
    Points P = (bench) ? input_b(&n) : input(&n);

    double start = wtime();
    size_t count = QuickhullP(n, P);
    double stop = wtime();

    double duration = stop - start;

    if (summary) {
        fprintf(stderr, "We found %zu elements on the hull\n", count);
        fprintf(stderr, "Runtime in ms: %lf\n", duration * 1e3);
        fprintf(stderr, "Bandwidth in MB/s: %lf\n",
                    n * 2 * sizeof(double) / duration / 1e6);
    }

    if (bench) {
        printf("%lf\n", duration);
    }

    if (print) {
        PrintPoints(count, P);
    }

    free(P.x);
    free(P.y);

    return EXIT_SUCCESS;
}

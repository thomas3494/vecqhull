/**
 * Generates random points on a disk
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char **argv)
{
    if (argc != 2) {
        printf("Usage: %s <n>\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n = atol(argv[1]);

    srand(314159);

    fprintf(stderr, "Creating %ld points\n", n);

    for (size_t i = 0; i < n; i++) {
        double theta = M_2_PI * ((double)rand() / RAND_MAX);
        double r     = (double)rand() / RAND_MAX;
        double x     = r * cos(theta);
        double y     = r * sin(theta);
        fwrite(&x, sizeof(double), 1, stdout);
        fwrite(&y, sizeof(double), 1, stdout);
    }

    return EXIT_SUCCESS;
}

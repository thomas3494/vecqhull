/**
 * Generates equidistanced points on a circle.
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

void equi_circle(long n)
{
    for (long i = 0; i < n; i++) {
        double theta = M_PI * i / (n / 2);
        double x = cos(theta);
        double y = sin(theta);
        fwrite(&x, sizeof(double), 1, stdout);
        fwrite(&y, sizeof(double), 1, stdout);
    }
}

long *permute(long n)
{
    long *indices = (long *)malloc(n * sizeof(long));
    for (long i = 0; i < n; i++) {
        indices[i] = i;
    }
    for (long i = 0; i <= n - 2; i++) {
        /* Pretty bad PRNG, but not that important */
        long j = i + rand() % (n - i);
        long temp = indices[i];
        indices[i] = indices[j];
        indices[j] = temp;
    }

    return indices;
}

void random_equi_circle(long n)
{
    long *indices = permute(n);

    for (long i = 0; i < n; i++) {
        double theta = M_PI * indices[i] / (n / 2);
        double x = cos(theta);
        double y = sin(theta);
        fwrite(&x, sizeof(double), 1, stdout);
        fwrite(&y, sizeof(double), 1, stdout);
    }

    free(indices);
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        printf("Usage: %s <power>\n", argv[0]);
        return EXIT_FAILURE;
    }

    /* Power of two so that 2 * pi / n is rounded correctly. */
    long n = 1;
    for (int i = 0; i < atol(argv[1]); i++) {
        n *= 2;
    }

    fprintf(stderr, "Creating %ld points\n", n);
    random_equi_circle(n);

    return EXIT_SUCCESS;
}

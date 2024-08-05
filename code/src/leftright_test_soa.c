#include "common.h"
#include <cstdio>
#include <cstdlib>

int main(int argc, char **argv)
{
    if (argc != 1) {
        printf("Usage: %s < <input in binary format>\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n;
    Points P = input_b_soa(&n);

    size_t left, right;
    Point p, q;

    double time1 = wtime();

    FindLeftRightV_soa(n, P, &left, &right);
    p = {P.x[left],  P.y[left]};
    q = {P.x[right], P.y[right]};

    double time2 = wtime();
    
    printf("We found left-point (%e, %e) and right point (%e, %e)\n",
            p.x, p.y, q.x, q.y);
    /* Bandwidth computed from the entire input, in reality we will not
     * read in most y-coordinates. */
    printf("This took %lf ms, or %lf GB/s\n", (time2 - time1) * 1e3,
                (double)n * sizeof(Point) / (time2 - time1) / 1e9);

    free(P.x);
    free(P.y);

    return EXIT_SUCCESS;
}

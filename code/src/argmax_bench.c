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
    Point *P = input_b(&n);

    Point p, q;

    FindLeftRight(n, P, &p, &q);

    double time1 = wtime();

    Point argmax1, argmax2;
    MinMaxV(n, P, p, q, p, &argmax1, &argmax2);

    double time2 = wtime();
    
    printf("Point furthest to the left is (%e, %e)\n"
           "Point furthest to the right is (%e, %e)\n",
            argmax1.x, argmax1.y, argmax2.x, argmax2.y);
    double duration = time2 - time1;
    printf("This took %lf ms, or %lf gflops or %lf GB/s\n", duration * 1e3,
                14.0 * n / duration / 1e9,
                (double)n * sizeof(Point) / duration / 1e9);

    free(P);

    return EXIT_SUCCESS;
}

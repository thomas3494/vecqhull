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

    double time1 = wtime();

    FindLeftRight(n, P, &p, &q);

    double time2 = wtime();
    
    printf("We found left-point (%e, %e) and right point (%e, %e)\n",
            p.x, p.y, q.x, q.y);
    printf("This took %lf ms, or %lf GB/s\n", (time2 - time1) * 1e3,
                (double)n * sizeof(Point) / (time2 - time1) / 1e9);

    double time3 = wtime();

    FindLeftRightV(n, P, &p, &q);

    double time4 = wtime();
    
    printf("We found left-point (%e, %e) and right point (%e, %e)\n",
            p.x, p.y, q.x, q.y);
    printf("This took %lf ms, or %lf GB/s\n", (time4 - time3) * 1e3,
                (double)n * sizeof(Point) / (time4 - time3) / 1e9);

    free(P);

    return EXIT_SUCCESS;
}

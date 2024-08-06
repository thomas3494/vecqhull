#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstring>

#include "common.h"

void TestPartition(size_t n, Points P, Point p, Point r, Point q, 
                   Point *r1_out, Point *r2_out, size_t *total1_out,
                   size_t *total2_out)
{
    Point r1, r2;
    size_t total1, total2;

    double time3 = wtime();
    TriPartitionV(n, P, p, r, q, &r1, &r2, &total1, &total2);
    double time4 = wtime();

    double duration = time4 - time3;
    printf("This took %lf ms, or %lf gflops or %lf GB/s\n", duration * 1e3,
            14.0 * n / 1e9 / duration,
            2.0 * sizeof(Point) * n / 1e9 / duration);
    printf("r1 = (%e, %e), r2 = (%e, %e)\n", r1.x, r1.y, r2.x, r2.y);
    printf("%zu total, %zu points to the left, %zu points to the right\n", 
            n, total1, total2);

    size_t wrong = 0;
    for (size_t i = 0; i < total1; i++) {
        Point u = {P.x[i], P.y[i]};
        if (!(orient(p, u, r) > 0)) {
            wrong++;
            printf("Error at %zu\n", i);
        }
    }
    for (size_t i = total1; i < total1 + total2; i++) {
        Point u = {P.x[i], P.y[i]};
        if (!(orient(r, u, q) > 0)) {
            wrong++;
            printf("Error at %zu\n", i);
        }
    }

    printf("%zu elements are in the wrong position\n", wrong);

    *r1_out = r1;
    *r2_out = r2;
    *total1_out = total1;
    *total2_out = total2;
}

int main(int argc, char **argv)
{
    if (argc != 1) {
        printf("Usage: %s < <points in PBBS format>\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n;
    Points P = input_b(&n);

    double time1 = wtime();

    size_t left, right;
    FindLeftRightV(n, P, &left, &right);
    Point p = {P.x[left],  P.y[left]};
    Point q = {P.x[right], P.y[right]};

    double time2 = wtime();
    printf("Finding left and right took %lf ms\n", (time2 - time1) * 1e3);

    Point r1, r2;
    size_t total1, total2;

    printf("First TriPartition\n");

    TestPartition(n, P, p, q, p, &r1, &r2, &total1, &total2);

    printf("Second TriPartition (left)\n");
    Point r1l, r2l;
    size_t total1l, total2l;
    Points Q = {P.x, P.y};
    TestPartition(total1, Q, p, r1, q, &r1l, &r2l, &total1l, &total2l);

    free(P.x);
    free(P.y);

    return EXIT_SUCCESS;
}

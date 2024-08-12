#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstring>
#include <omp.h>

#include "common.h"

void TestPartition(size_t n, Points P, Points P2, Point p, Point r, Point q, 
                   Point *r1_out, Point *r2_out, size_t *total1_out,
                   size_t *total2_out)
{
    Point r1, r2;
    size_t total1, total2;

    double time3 = wtime();
    TriPartitionV(n, P, p, r, q, &r1, &r2, &total1, &total2);
    double time4 = wtime();

    double duration = time4 - time3;
    printf("Sequential took %lf ms, or %lf gflops or %lf GB/s\n", 
            duration * 1e3,
            14.0 * n / 1e9 / duration,
            2.0 * sizeof(Point) * n / 1e9 / duration);
    printf("r1 = (%e, %e), r2 = (%e, %e)\n", r1.x, r1.y, r2.x, r2.y);
    printf("%zu total, %zu points to the left, %zu points to the right\n", 
            n, total1, total2);

    unsigned int nthreads;
    #pragma omp parallel master
    {
        nthreads = omp_get_num_threads();
    }
    time3 = wtime();
    TriPartitionP(n, P2, p, r, q, &r1, &r2, &total1, &total2, nthreads);
    time4 = wtime();

    duration = time4 - time3;
    printf("Parallel took %lf ms, or %lf gflops or %lf GB/s\n", 
            duration * 1e3,
            14.0 * n / 1e9 / duration,
            2.0 * sizeof(Point) * n / 1e9 / duration);
    printf("r1 = (%e, %e), r2 = (%e, %e)\n", r1.x, r1.y, r2.x, r2.y);
    printf("%zu total, %zu points to the left, %zu points to the right\n", 
            n, total1, total2);

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
    Points P2;
    P2.x = (double *)aligned_alloc(64, n * sizeof(double));
    P2.y = (double *)aligned_alloc(64, n * sizeof(double));
    memcpy(P2.x, P.x, n * sizeof(double));
    memcpy(P2.y, P.y, n * sizeof(double));

    double time1 = wtime();

    size_t left, right;
    FindLeftRightV(n, P, &left, &right);
    Point p = {P.x[left],  P.y[left]};
    Point q = {P.x[right], P.y[right]};

    double time2 = wtime();
    printf("Finding left and right took %lf ms\n", (time2 - time1) * 1e3);

    Point r1, r2;
    size_t total1, total2;

    TestPartition(n, P, P2, p, q, p, &r1, &r2, &total1, &total2);

    free(P.x);
    free(P.y);
    free(P2.x);
    free(P2.y);

    return EXIT_SUCCESS;
}

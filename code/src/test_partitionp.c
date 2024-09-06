#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstring>
#include <omp.h>

#include "common.h"

void TestPartition(size_t n, Points P, Points P2, Point p, Point r, Point q,
                   Point *r1_out, Point *r2_out, size_t *c1_out,
                   size_t *c2_out)
{
    Point r1_seq, r2_seq;
    Point r1_par, r2_par;
    size_t c1_seq, c2_seq;
    size_t c1_par, c2_par;

    printf("\n========== Sequential ==========\n");

    double time3 = wtime();
    TriPartitionV(n, P, p, r, q, &r1_seq, &r2_seq, &c1_seq, &c2_seq);
    double time4 = wtime();

    //Points S2 = {P.x + c2, P.y + c2};
    //PrintPoints(n - c2, S2);

    double duration = time4 - time3;
    printf("Sequential took %lf ms, or %lf gflops or %lf GB/s\n",
            duration * 1e3,
            14.0 * n / 1e9 / duration,
            2.0 * sizeof(Point) * n / 1e9 / duration);
    printf("r1 = (%e, %e), r2 = (%e, %e)\n", r1_seq.x, r1_seq.y, r2_seq.x, r2_seq.y);
    printf("c1 = %zu, c2 = %zu\n", c1_seq, c2_seq);

    printf("\n========== Parallel ==========\n");

    unsigned int nthreads;
    #pragma omp parallel master
    {
        nthreads = omp_get_num_threads();
    }
    time3 = wtime();
    TriPartitionP(n, P2, p, r, q, &r1_par, &r2_par, &c1_par, &c2_par, nthreads);
    time4 = wtime();

//    Points S2 = {P2.x + c2, P2.y + c2};
//    PrintPoints(n - c2, S2);

    duration = time4 - time3;
    printf("Parallel took %lf ms, or %lf gflops or %lf GB/s\n",
            duration * 1e3,
            14.0 * n / 1e9 / duration,
            2.0 * sizeof(Point) * n / 1e9 / duration);
    printf("r1 = (%e, %e), r2 = (%e, %e)\n", r1_par.x, r1_par.y, r2_par.x, r2_par.y);
    printf("c1 = %zu, c2 = %zu\n", c1_par, c2_par);

    bool success = true;
    for (size_t i = 0; i < c1_par; i++) {
        Point u = {P2.x[i], P2.y[i]};
        if (!(orient(p, u, r) > 0)) {
            printf("Error at S1 (%zu not in S1)\n", i);
            success = false;
        }
    }

    for (size_t i = c2_par; i < n; i++) {
        Point u = {P2.x[i], P2.y[i]};
        if (!(orient(r, u, q) > 0)) {
            printf("Error at S2 (%zu not in S2)\n", i);
            success = false;
        }

        if (orient(p, u, r) > 0) {
            printf("Error at S2 (%zu already in S1)\n", i);
            success = false;
        }
    }

    if (r1_seq.x != r1_par.x || r1_seq.y != r1_par.y) {
        printf("R1 not equal\n");
        success = false;
    }

    if (r2_seq.x != r2_par.x || r2_seq.y != r2_par.y) {
        printf("R1 not equal\n");
        success = false;
    }

    if (c1_seq != c1_par) {
        printf("C1 not equal\n");
        success = false;
    }

    if (c2_seq != c2_par) {
        printf("C2 not equal\n");
        success = false;
    }

    printf("Partition was %ssuccesfull\n", (success) ? "" : "un");

    *r1_out = r1_par;
    *r2_out = r2_par;
    *c1_out = c1_par;
    *c2_out = c2_par;
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
    P2.x = (double *)aligned_alloc(64, (n * sizeof(double) + 63) / 64 * 64);
    P2.y = (double *)aligned_alloc(64, (n * sizeof(double) + 63) / 64 * 64);
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
    size_t c1, c2;

    TestPartition(n, P, P2, p, q, p, &r1, &r2, &c1, &c2);

//    printf("Partition left\n");
//    size_t c1l, c2l;
//    Point r1l, r2l;
//    TestPartition(c1, P, P2, p, r1, q, &r1l, &r2l, &c1l, &c2l);
//
//    printf("Partition right\n");
//    Points S2 = {P.x + c2, P.y + c2};
//    Points S22 = {P2.x + c2, P2.y + c2};
//    size_t c1r, c2r;
//    Point r1r, r2r;
//    TestPartition(n - c2, S2, S22, q, r2, p, &r1r, &r2r, &c1r, &c2r);

    free(P.x);
    free(P.y);
    free(P2.x);
    free(P2.y);

    return EXIT_SUCCESS;
}

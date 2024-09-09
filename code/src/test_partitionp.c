#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstring>
#include <omp.h>

#include "common.h"

void TestPartition(size_t n, Points P1, Points P2, Point p, Point r, Point q,
                   Point *r1_out, Point *r2_out, size_t *c1_out,
                   size_t *c2_out)
{
    assert(n > 1);

    Point r1_seq, r2_seq, r1_par, r2_par;
    size_t c1_seq, c2_seq, c1_par, c2_par;

    printf("\n========== Sequential ==========\n");

    double time1 = wtime();
    TriPartitionV(n, P1, p, r, q, &r1_seq, &r2_seq, &c1_seq, &c2_seq);
    double time2 = wtime();

    //Points S2 = {P1.x + c2, P1.y + c2};
    //PrintPoints(n - c2, S2);

    double duration = time2 - time1;
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

    time1 = wtime();
    TriPartitionP(n, P2, p, r, q, &r1_par, &r2_par, &c1_par, &c2_par, nthreads);
    time2 = wtime();

    //Points S2 = {P2.x + c2, P2.y + c2};
    //PrintPoints(n - c2, S2);

    duration = time2 - time1;
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

    bool s1_is_empty = (c1_seq == 0);
    bool s2_is_empty = (c2_seq == n - 2);

    if (!s1_is_empty && (r1_seq.x != r1_par.x || r1_seq.y != r1_par.y)) {
        printf("R1 not equal\n");
        success = false;
    }

    if (!s2_is_empty && (r2_seq.x != r2_par.x || r2_seq.y != r2_par.y)) {
        printf("R2 not equal\n");
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

    printf("Partition was %ssuccesfull\n\n", (success) ? "" : "un");

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

    size_t n, left, right;
    Points P1 = input_b(&n);

    double time1 = wtime();
    FindLeftRightV(n, P1, &left, &right);
    double time2 = wtime();

    Point p = {P1.x[left],  P1.y[left]};
    Point q = {P1.x[right], P1.y[right]};

    printf("Finding left and right took %lf ms\n", (time2 - time1) * 1e3);
    printf("p = (%e, %e), q = (%e, %e)\n", p.x, p.y, q.x, q.y);

    if (right != 0) {
        swap(P1, 0, left);
        swap(P1, n - 1, right);
    } else {
        swap(P1, n - 1, right);
        swap(P1, 0, left);
    }

    Points P2;
    P2.x = (double *)aligned_alloc(64, (n * sizeof(double) + 63) / 64 * 64);
    P2.y = (double *)aligned_alloc(64, (n * sizeof(double) + 63) / 64 * 64);
    memcpy(P2.x, P1.x, n * sizeof(double));
    memcpy(P2.y, P1.y, n * sizeof(double));

    printf("Partition:\n");
    Point r1, r2;
    size_t c1, c2;
    Points S1_1 = {P1.x + 1, P1.y + 1};
    Points S1_2 = {P2.x + 1, P2.y + 1};
    TestPartition(n - 2, S1_1, S1_2, p, q, p, &r1, &r2, &c1, &c2);

    size_t total1 = c1;
    size_t total2 = n - 2 - c2;

    printf("Partition left:\n");
    Point r1l, r2l;
    size_t c1l, c2l;
    TestPartition(total1, S1_1, S1_2, p, r1, q, &r1l, &r2l, &c1l, &c2l);

    printf("Partition right:\n");
    Point r1r, r2r;
    size_t c1r, c2r;
    Points S2_1 = {S1_1.x + c2, S1_1.y + c2};
    Points S2_2 = {S1_2.x + c2, S1_2.y + c2};
    TestPartition(total2, S2_1, S2_2, q, r2, p, &r1r, &r2r, &c1r, &c2r);

    free(P1.x);
    free(P1.y);
    free(P2.x);
    free(P2.y);

    return EXIT_SUCCESS;
}

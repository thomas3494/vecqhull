#include <omp.h>
#include <string.h>
#include <math.h>

#include "quickhull.h"
#include "common.h"

//#define MEASURE_BW

#ifdef MEASURE_BW
size_t read_bw     = 0;
size_t write_bw    = 0;
size_t condense_bw = 0;
#endif

/* Forward declarations */
size_t FindHull(size_t n, Points P, Point p, Point r, Point q);
size_t FindHullP(size_t n, Points P, Point p, Point r, Point q,
                 unsigned int nthreads);

size_t Quickhull(size_t n, Points P)
{
    /* Find the points with left-most and right-most x-coordinate.
     * (In case of ties, with bottom and top y-coordinate.)
     * These are guaranteed to be on the convex hull, and will be our
     * first bisection. */
#ifdef MEASURE_BW
    read_bw += n * sizeof(Point);
#endif
    size_t left, right;
    FindLeftRightV(n, P, &left, &right);
    Point p = {P.x[left], P.y[left]};
    Point q = {P.x[right], P.y[right]};

    if (right != 0) {
        swap(P, 0, left);
        swap(P, n - 1, right);
    } else {
        swap(P, n - 1, right);
        swap(P, 0, left);
    }

    Point r1, r2;
    size_t c1, c2;
    Points S1 = {P.x + 1, P.y + 1};
    TriPartitionV(n - 2, S1, p, q, p, &r1, &r2, &c1, &c2);
    size_t total1 = c1;
    size_t total2 = n - 2 - c2;

#ifdef MEASURE_BW
    read_bw  += (n - 2) * sizeof(Point);
    write_bw += (total1 + total2) * sizeof(Point);
#endif

    Points S2 = {S1.x + c2, S1.y + c2};
    size_t lcount = FindHull(total1, S1, p, r1, q);
    size_t rcount = FindHull(total2, S2, q, r2, p);

    /* Put <p> <left hull> <q> <right hull> together */
    memmove(S1.x + lcount + 1, S2.x, rcount * sizeof(double));
    memmove(S1.y + lcount + 1, S2.y, rcount * sizeof(double));

#ifdef MEASURE_BW
    condense_bw  += 2 * rcount * sizeof(Point);
#endif

    P.x[lcount + 1] = q.x;
    P.y[lcount + 1] = q.y;

#ifdef MEASURE_BW
    printf("Total bw: %zu, read bw: %zu bytes, write bw: %zu, condese bw: %zu\n",
            read_bw + write_bw + condense_bw, read_bw, write_bw, condense_bw);
#endif

    return 2 + lcount + rcount;
}

size_t FindHull(size_t n, Points P, Point p, Point r, Point q)
{
    if (n <= 1) return n;

    Point r1, r2;
    size_t c1, c2;
    TriPartitionV(n, P, p, r, q, &r1, &r2, &c1, &c2);
    size_t total1 = c1;
    size_t total2 = n - c2;

#ifdef MEASURE_BW
    read_bw  += n * sizeof(Point);
    write_bw += (total1 + total2) * sizeof(Point);
#endif

    Points S1 = P;
    size_t lcount = FindHull(total1, S1, p, r1, r);

    Points S2 = {P.x + c2, P.y + c2};
    size_t rcount = FindHull(total2, S2, r, r2, q);

    /* Condense left and right hull into contiguous memory */
    /* TODO: we have to move everything because the order of right hull
     * matters. We could also have cut the last part of S2 before,
     * the recursive call, potentially moving less data. */
    memmove(S1.x + lcount + 1, S2.x, rcount * sizeof(double));
    memmove(S1.y + lcount + 1, S2.y, rcount * sizeof(double));

#ifdef MEASURE_BW
    condense_bw  += 2 * rcount * sizeof(Point);
#endif

    P.x[lcount] = r.x;
    P.y[lcount] = r.y;

    return 1 + lcount + rcount;
}

size_t QuickhullP(size_t n, Points P)
{
    unsigned int nthreads;
    #pragma omp parallel master
    {
        nthreads = omp_get_num_threads();
    }

    omp_set_max_active_levels(nthreads);

    /* Find the points with left-most and right-most x-coordinate.
     * (In case of ties, with bottom and top y-coordinate.)
     * These are guaranteed to be on the convex hull, and will be our
     * first bisection. */
    size_t left, right;
    FindLeftRightVP(n, P, &left, &right);
    Point p = {P.x[left], P.y[left]};
    Point q = {P.x[right], P.y[right]};

    if (right != 0) {
        swap(P, 0, left);
        swap(P, n - 1, right);
    } else {
        swap(P, n - 1, right);
        swap(P, 0, left);
    }

    Point r1, r2;
    size_t c1, c2;
    Points S1 = {P.x + 1, P.y + 1};
    TriPartitionP(n - 2, S1, p, q, p, &r1, &r2, &c1, &c2, nthreads);
    size_t total1 = c1;
    size_t total2 = n - 2 - c2;
    Points S2 = {S1.x + c2, S1.y + c2};

    unsigned int threads1 = roundf((double)total1 /
                                        (total1 + total2) * nthreads);
    unsigned int threads2 = nthreads - threads1;

    size_t lcount, rcount;
    if (threads1 == 0 || threads2 == 0) {
        lcount = FindHull(total1, S1, p, r1, q);
        rcount = FindHull(total2, S2, q, r2, p);
    } else {
        #pragma omp parallel shared(threads1, threads2) num_threads(2)
        {
            #pragma omp single nowait
            {
                #pragma omp task
                lcount = FindHullP(total1, S1, p, r1, q, threads1);
            }

            #pragma omp single nowait
            {
                #pragma omp task
                rcount = FindHullP(total2, S2, q, r2, p, threads2);
            }
        }
    }

    /* Put <p> <left hull> <q> <right hull> together */
    memmove(S1.x + lcount + 1, S2.x, rcount * sizeof(double));
    memmove(S1.y + lcount + 1, S2.y, rcount * sizeof(double));
    P.x[lcount + 1] = q.x;
    P.y[lcount + 1] = q.y;

    return 2 + lcount + rcount;
}

size_t FindHullP(size_t n, Points P, Point p, Point r, Point q,
                 unsigned int nthreads)
{
    if (n <= 1) return n;

    Point r1, r2;
    size_t c1, c2;
    TriPartitionP(n, P, p, r, q, &r1, &r2, &c1, &c2, nthreads);
    size_t total1 = c1;
    size_t total2 = n - c2;

    Points S1 = P;
    Points S2 = {P.x + c2, P.y + c2};

    unsigned int threads1 = (total1 + total2 == 0) ?
                        0 :
                        roundf((double)total1 / (total1 + total2) * nthreads);
    unsigned int threads2 = nthreads - threads1;

    size_t lcount, rcount;
    if (threads1 == 0 || threads2 == 0) {
        lcount = FindHull(total1, S1, p, r1, r);
        rcount = FindHull(total2, S2, r, r2, q);
    } else {
        #pragma omp parallel shared(threads1, threads2) num_threads(2)
        {
            #pragma omp single nowait
            {
                #pragma omp task
                {
                    lcount = (threads1 > 1) ?
                                 FindHullP(total1, S1, p, r1, r, threads1) :
                                 FindHull(total1, S1, p, r1, r);
                }
            }

            #pragma omp single nowait
            {
                #pragma omp task
                    rcount = (threads2 > 1) ?
                            FindHullP(total2, S2, r, r2, q, threads2) :
                            FindHull(total2, S2, r, r2, q);
            }
        }
    }

    /* Condense left and right hull into contiguous memory */
    memmove(S1.x + lcount + 1, S2.x, rcount * sizeof(double));
    memmove(S1.y + lcount + 1, S2.y, rcount * sizeof(double));
    P.x[lcount] = r.x;
    P.y[lcount] = r.y;

    return 1 + lcount + rcount;
}

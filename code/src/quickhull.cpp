#include <omp.h>
#include <string.h>

#include "quickhull.h"
#include "common.h"

/* Forward declaration */
size_t FindHull(size_t n, Points P, Point p, Point q, Point rk);


size_t Quickhull(size_t n, Points P)
{
    /* Find the points with left-most and right-most x-coordinate.
     * (In case of ties, with bottom and top y-coordinate.)
     * These are guaranteed to be on the convex hull, and will be our
     * first bisection. */
    size_t left, right;
    Point p, q;
    FindLeftRightV(n, P, &left, &right);

    swap(P, 0, left);
    swap(P, n - 1, right);

    Point r1, r2;
    size_t total1, total2;
    Points S1 = {P.x + 1, P.y + 1};
    TriPartitionV(n - 2, S1, p, q, p, &r1, &r2, &total1, &total2);

    Points S2 = {S1.x + total1, S1.y + total1};
    size_t lcount = FindHull(total1, S1, p, r1, q);
    size_t rcount = FindHull(total2, S2, q, r2, p);

    /* Condense left hull, right hull, q into contiguous memory */
    memmove(S1.x + lcount, S2.x, rcount * sizeof(double));
    memmove(S1.y + lcount, S2.y, rcount * sizeof(double));
    swap(P, 1 + lcount + rcount + 1, n - 1);

    return 2 + lcount + rcount;
}

size_t FindHull(size_t n, Points P, Point p, Point rk, Point q)
{
    if (n <= 1) return n;

    Point t1, t2;
    size_t total1, total2;
    TriPartitionV(n, P, p, rk, q, &t1, &t2, &total1, &total2);

    Points S1 = P;
    size_t lcount = FindHull(total1, S1, p, t1, rk);

    Points S2 = {P.x + total1, P.y + total1};
    size_t rcount = FindHull(total2, S2, rk, t2, q);

    /* Condense left and right hull into contiguous memory */
    memmove(S1.x + lcount + 1, S2.x, rcount * sizeof(double));
    memmove(S1.y + lcount + 1, S2.y, rcount * sizeof(double));
    P.x[lcount] = rk.x;
    P.y[lcount] = rk.y;

    return 1 + lcount + rcount;
}

size_t QuickhullP(size_t n, Points P)
{
    return n;
}

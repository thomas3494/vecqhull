#ifndef QUICKHULL_GUARD
#define QUICKHULL_GUARD

#include <cstddef>
#include <hwy/highway.h>

#include "typedefs.h"

using namespace hwy::HWY_NAMESPACE;

/* pbbs 2d sequence format */
void PrintPoints(size_t n, Points P);

/* pbbs 2d sequence format on stdin */
Points input(size_t *n /* out */);

/* Binary on stdin */
Points input_b(size_t *n /* out */);

/* Wallclock time in seconds */
double wtime(void);

inline size_t ceildiv(size_t a, size_t b)
{
    return (a + b - 1) / b;
}

inline size_t max(size_t a, size_t b)
{
    return (a > b) ? a : b;
}

inline size_t min(size_t a, size_t b)
{
    return (a < b) ? a : b;
}

inline void swap(Points P, size_t i, size_t j)
{
    double tempx = P.x[i];
    double tempy = P.y[i];
    P.x[i] = P.x[j];
    P.y[i] = P.y[j];
    P.x[j] = tempx;
    P.y[j] = tempy;
}

/* If this is negative, then p, q, u is a right-turn.
 * Returns precisely 0 for q = p or q = u.
 * For p = u, it mathematically returns 0, but may give
 * a round-off error. */
inline double orient(Point p, Point q, Point u)
{
    return (p.x - q.x) * (u.y - q.y) - (p.y - q.y) * (u.x - q.x);
}

inline Vec<ScalableTag<double>> orientV(Vec<ScalableTag<double>> px,
                                        Vec<ScalableTag<double>> py,
                                        Vec<ScalableTag<double>> qx,
                                        Vec<ScalableTag<double>> qy,
                                        Vec<ScalableTag<double>> ux,
                                        Vec<ScalableTag<double>> uy)
{
    return (px - qx) * (uy - qy) - (py - qy) * (ux - qx);
}

/**
 * True if u1 > u2
 **/
inline bool greater_orient(Point p, Point u1, Point u2, Point q)
{
    return ((p.y - q.y) * (u1.x - u2.x) - (p.x - q.x) * (u1.y - u2.y) > 0);
}

inline Mask<ScalableTag<double>> greater_orientV(Vec<ScalableTag<double>> px,
                                                 Vec<ScalableTag<double>> py,
                                                 Vec<ScalableTag<double>> u1x,
                                                 Vec<ScalableTag<double>> u1y,
                                                 Vec<ScalableTag<double>> u2x,
                                                 Vec<ScalableTag<double>> u2y,
                                                 Vec<ScalableTag<double>> qx,
                                                 Vec<ScalableTag<double>> qy)
{
    const ScalableTag<double> d;
    return ((py - qy) * (u1x - u2x) - (px - qx) * (u1y - u2y) > Zero(d));
}

/* Scalar */
void FindLeftRight(size_t n, Points P, size_t *left_out, size_t *right_out);

/* Vectorized */
void FindLeftRightV(size_t n, Points P, size_t *left_out, size_t *right_out);

/* Vectorized and parallelised with OpenMP */
void FindLeftRightVP(size_t n, Points P, size_t *left_out, size_t *right_out);

/* [0, c1_out[ contains S1, [c2_out, n[ contains S2 */
void TriPartitionV(size_t n, Points P, Point p, Point u, Point q,
                   Point *r1_out, Point *r2_out,
                   size_t *c1_out, size_t *c2_out);

/* Parallelised with OpenMP */
void TriPartitionP(size_t n, Points P, Point p, Point u, Point q,
                   Point *r1_out, Point *r2_out,
                   size_t *c1_out, size_t *c2_out, unsigned int nthreads);

#endif /* QUICKHULL_GUARD */

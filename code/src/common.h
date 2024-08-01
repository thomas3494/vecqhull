#ifndef QUICKHULL_GUARD
#define QUICKHULL_GUARD

#include <cstddef>

/* TODO make sure we can use this library from C code (and everything having
 * a C FFI). */

typedef struct {
    double x;
    double y;
} Point;

/* pbbs 2d sequence format */
void PrintPoints(size_t n, Point *P);

/* pbbs 2d sequence format on stdin */
Point *input(size_t *n /* out */);

/* Binary on stdin */
Point *input_b(size_t *n /* out */);

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

/* If this is negative, then p, q, u is a right-turn.
 * Returns precisely 0 for q = p or q = u.
 * For p = u, it mathematically returns 0, but may give
 * a round-off error. */
inline double orient(Point p, Point q, Point u)
{
    return (p.x - q.x) * (u.y - q.y) - (p.y - q.y) * (u.x - q.x);
}

/* Naive */
void FindLeftRight(size_t n, Point *P, Point *left_out, Point *right_out);

/* Vectorized */
void FindLeftRightV(size_t n, Point *P, Point *left_out, Point *right_out);

/* Wallclock time in seconds */
double wtime(void);

#endif /* QUICKHULL_GUARD */

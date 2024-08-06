#ifndef _QUICKHULL_GUARD_JLSDFKJSDJFDFJLKSDFJ
#define _QUICKHULL_GUARD_JLSDFKJSDJFDFJLKSDFJ

#include <stddef.h>

typedef struct {
    double *x;
    double *y;
} Points;

/**
 * Computes the convex hull of P in-place. Returns the number of elements in
 * this hull. Memory management, including a potential realloc is the
 * responsibility of the caller.
 *
 * The P version is parallelised with OpenMP.
 **/
size_t Quickhull(size_t n, Points P);
size_t QuickhullP(size_t n, Points P);

#endif /* _QUICKHULL_GUARD_JLSDFKJSDJFDFJLKSDFJ */

#ifndef VQHULL_GUARD_jalsdfjklasdjflasdjflasdghadf
#define VQHULL_GUARD_jalsdfjklasdjflasdjflasdghadf

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

typedef struct {
    double *x;
    double *y;
} Points;

/**
 * Takes a set of P of n points.
 * Replaces P in-place with its convex hull.
 * Returns the number of points in the convex hull.
 **/
size_t VecQuickhull(size_t n, Points P);

/**
 * Same as VecQuickhull, but parallelised with OpenMP 
 **/
size_t VecQuickhullP(size_t n, Points P);

#ifdef __cplusplus
}
#endif

#endif /* VQHULL_GUARD_jalsdfjklasdjflasdjflasdghadf */

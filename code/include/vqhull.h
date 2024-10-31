#ifndef VQHULL_GUARD_jalsdfjklasdjflasdjflasdghadf
#define VQHULL_GUARD_jalsdfjklasdjflasdjflasdghadf

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

/**
 * Takes a set of n points in the plane.
 * Replaces the set in-place with its convex hull.
 * Returns the number of points in the convex hull.
 **/
size_t VecQuickhull(size_t n, double *x_coor, double *y_coor);

/**
 * Same as VecQuickhull, but parallelised with OpenMP 
 * For the best performance, make sure x_coor = y_coor mod 64
 * (to avoid false sharing).
 **/
size_t VecQuickhullP(size_t n, double *x_coor, double *y_coor);

#ifdef __cplusplus
}
#endif

#endif /* VQHULL_GUARD_jalsdfjklasdjflasdjflasdghadf */

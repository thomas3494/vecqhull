/**
 * Taylor series around 0 give
 *
 * sin(x) = x + x^3 / 6 + O(x^5)
 * cos(x) = 1 + x^2 / 2 + O(x^4)
 *
 * If |x| < 2^(26), then fl(cos(x)) = 1
 *
 * So if we generate multiple theta < 2^26, we get colinear points.
 *
 * The same happens around 0.5 pi, pi, 1.5 pi.
 *
 * As orient is linear, we have orient(p, q, u(1 + delta)) =
 * (1 + delta)orient(p, q, u)
 *
 * |orient(p, q, u(1 + delta)) - orient(p, q, u)| / |orient(p, q, u)| =
 * delta, so this is well-conditioned. 
 *
 * However, when finding the maximum, we look at 
 *    orient(p, q, u) - orient(p, q, u')
 * Subtraction is ill-conditioned for values that are close to eachother,
 * so this is problematic.
 *
 * Alternatively, we could compare orient(p, q, u - u') to zero.
 * As u, u' are exact, it is not a problem that subtraction is ill-conditioned.
 **/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "quickhull.h"

int smaller(const void *x, const void *y)
{
    double xd = *(double *)x;
    double yd = *(double *)y;
    if (xd < yd) {
        return -1;
    } else if (xd == yd) {
        return 0;
    } else {
        return 1;
    }
}

int main(int argc, char **argv)
{
    if (argc != 1) {
        printf("Usage: %s takes binary points on stdin.\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t n;
    Points P = input_b(&n);

    qsort(P.x, n, sizeof(double), smaller);
    qsort(P.y, n, sizeof(double), smaller);

    size_t colinear_x = 0;
    size_t colinear_y = 0;
    for (size_t i = 2; i < n; i++) {
        colinear_x += (P.x[i] == P.x[i - 1] && P.x[i] == P.x[i - 2]);
        colinear_y += (P.y[i] == P.y[i - 1] && P.y[i] == P.y[i - 2]);
    }

    printf("Colinear (x): %zu\nColinear (y): %zu\n", colinear_x, colinear_y);

    free(P.x);
    free(P.y);

    return EXIT_SUCCESS;
}

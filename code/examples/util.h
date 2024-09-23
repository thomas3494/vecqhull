#include <stddef.h>
#include <vqhull.h>

typedef struct {
    double x;
    double y;
} Point;

/* pbbs 2d sequence format */
void PrintPoints(size_t n, Points P);

/* pbbs 2d sequence format on stdin */
Points input(size_t *n /* out */);

/* Binary on stdin */
Points input_b(size_t *n /* out */);

/* Wallclock time in seconds */
double wtime(void);

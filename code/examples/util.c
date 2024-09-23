#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <vqhull.h>
#include <string.h>

#include "util.h"

void PrintPoints(size_t n, Points P)
{
    printf("pbbs_sequencePoint2d\n");
    for (size_t i = 0; i < n; i++) {
        printf("%.17e %.17e\n", P.x[i], P.y[i]);
    }
}

Points input(size_t *n /* out param */)
{
    size_t line_length;
    char *file_header = NULL;

    if (getline(&file_header, &line_length, stdin) == -1) {
        perror("Error in getting file header");
    }

    if (strcmp(file_header, "pbbs_sequencePoint2d\n") != 0) {
        printf("Error: not in pbbs_sequencePoint2d format\n");
        abort();
    }

    free(file_header);

    size_t size = 1000;
    Points P;
    P.x = (double *)malloc(size * sizeof(double));
    P.y = (double *)malloc(size * sizeof(double));

    size_t i = 0;
    while (scanf("%lf %lf", &P.x[i], &P.y[i]) == 2) {
        i++;
        if (i >= size) {
            size *= 2;
            P.x = (double *)realloc(P.x, size * sizeof(Point));
            P.y = (double *)realloc(P.y, size * sizeof(Point));
        }
    }

    *n = i;
    P.x = (double *)realloc(P.x, i * sizeof(Point));
    P.y = (double *)realloc(P.y, i * sizeof(Point));

    return P;
}

Points input_b(size_t *n /* out param */)
{
    size_t size = 1000;
    Points P;
    P.x = (double *)malloc(size * sizeof(double));
    P.y = (double *)malloc(size * sizeof(double));

    size_t i = 0;
    double x, y;
    while ((fread(&x, sizeof(double), 1, stdin) == 1) &&
           (fread(&y, sizeof(double), 1, stdin) == 1))
    {
        P.x[i] = x;
        P.y[i] = y;
        i++;
        if (i >= size) {
            size *= 2;
            P.x = (double *)realloc(P.x, size * sizeof(Point));
            P.y = (double *)realloc(P.y, size * sizeof(Point));
        }
    }

    *n = i;
    P.x = (double *)realloc(P.x, i * sizeof(Point));
    P.y = (double *)realloc(P.y, i * sizeof(Point));

    return P;
}

double wtime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)(tv.tv_usec / 1e6 + tv.tv_sec);
}

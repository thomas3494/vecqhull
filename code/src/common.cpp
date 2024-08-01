#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cerrno>
#include <sys/time.h>
#include <hwy/highway.h>
#include "common.h"

using namespace hwy::HWY_NAMESPACE;

void PrintPoints(size_t n, Point *P)
{
    printf("pbbs_sequencePoint2d\n");
    for (size_t i = 0; i < n; i++) {
        printf("%.17e %.17e\n", P[i].x, P[i].y);
    }
}

Point *input(size_t *n /* out param */)
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
    Point *P = (Point *)malloc(size * sizeof(Point));

    size_t i = 0;
    while (scanf("%lf %lf", &P[i].x, &P[i].y) == 2) {
        i++;
        if (i >= size) {
            size *= 2;
            P = (Point *)realloc(P, size * sizeof(Point));
        }
    }

    *n = i;
    P = (Point *)realloc(P, i * sizeof(Point));

    return P;
}

Point *input_b(size_t *n /* out param */)
{
    size_t size = 1000;
    Point *P = (Point *)malloc(size * sizeof(Point));

    size_t i = 0;
    double x, y;
    while ((fread(&x, sizeof(double), 1, stdin) == 1) &&
           (fread(&y, sizeof(double), 1, stdin) == 1))
    {
        P[i].x = x;
        P[i].y = y;
        i++;
        if (i >= size) {
            size *= 2;
            P = (Point *)realloc(P, size * sizeof(Point));
        }
    }

    *n = i;
    P = (Point *)realloc(P, i * sizeof(Point));

    return P;
}

void FindLeftRight(size_t n, Point *P, Point *left_out, Point *right_out)
{
    Point p = P[0];
    Point q = P[0];

    for (size_t i = 1; i < n; i++) {
        if (P[i].x < p.x || (P[i].x == p.x && P[i].y < p.y)) {
            p = P[i];
        }
        if (P[i].x > q.x || (P[i].x == q.x && P[i].y > q.y)) {
            q = P[i];
        }
    }

    *left_out = p;
    *right_out = q;
}


void FindLeftRightV(size_t n, Point * HWY_RESTRICT P, 
                    Point *left_out, Point *right_out)
{
    const ScalableTag<double> d;

    Vec<ScalableTag<double>> leftx, lefty, rightx, righty, x_coor, y_coor;
    LoadInterleaved2(d, (double *)P, x_coor, y_coor);

    leftx  = x_coor;
    rightx = x_coor;
    lefty  = y_coor;
    righty = y_coor;

    size_t i = 0;
    for (; i + Lanes(d) <= n; i += Lanes(d)) {
        LoadInterleaved2(d, (double *)(P + i), x_coor, y_coor);
        auto mask1 = Or(x_coor < leftx, 
                        And((x_coor == leftx), (y_coor < lefty)));
        leftx = IfThenElse(mask1, x_coor, leftx);
        lefty = IfThenElse(mask1, y_coor, lefty);
        auto mask2 = Or(x_coor > rightx, 
                        And((x_coor == rightx), (y_coor > righty)));
        rightx = IfThenElse(mask2, x_coor, rightx);
        righty = IfThenElse(mask2, y_coor, righty);
    }
    if (i < n) {
        /* Still have < Lanes(d) elements left. It is not a problem to
         * check some elements twice, so we just take the last Lanes(d)
         * elements. */
        i = n - Lanes(d);
        LoadInterleaved2(d, (double *)(P + i), x_coor, y_coor);
        auto mask1 = Or(x_coor < leftx, 
                        And((x_coor == leftx), (y_coor < lefty)));
        leftx = IfThenElse(mask1, x_coor, leftx);
        lefty = IfThenElse(mask1, y_coor, lefty);
        auto mask2 = Or(x_coor > rightx, 
                        And((x_coor == rightx), (y_coor > righty)));
        rightx = IfThenElse(mask2, x_coor, rightx);
        righty = IfThenElse(mask2, y_coor, righty);
    }

    double leftx_arr[Lanes(d)];
    double lefty_arr[Lanes(d)];
    double rightx_arr[Lanes(d)];
    double righty_arr[Lanes(d)];

    Store(leftx,  d, leftx_arr);
    Store(rightx, d, rightx_arr);
    Store(lefty,  d, lefty_arr);
    Store(righty, d, righty_arr);

    size_t left_ind = 0;
    size_t right_ind = 0;
    for (size_t i = 1; i < Lanes(d); i++) {
        if ((leftx_arr[i] < leftx_arr[left_ind]) ||
                ((leftx_arr[i] == leftx_arr[left_ind]) &&
                    (lefty_arr[i] < lefty_arr[left_ind])))
        {
            left_ind = i;
        }
        if ((rightx_arr[i] > rightx_arr[right_ind]) ||
                ((rightx_arr[i] == rightx_arr[right_ind]) &&
                    (righty_arr[i] > righty_arr[right_ind])))
        {
            right_ind = i;
        }
    }

    left_out->x = leftx_arr[left_ind];
    left_out->y = lefty_arr[left_ind];
    right_out->x = rightx_arr[right_ind];
    right_out->y = righty_arr[right_ind];
}

double wtime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)(tv.tv_usec / 1e6 + tv.tv_sec);
}

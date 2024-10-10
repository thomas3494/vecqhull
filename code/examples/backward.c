#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <vqhull.h>
#include <assert.h>

#include "util.h"

bool right_turn(Point p, Point u, Point q)
{
    return (u.x - p.x) * (q.y - p.y) < (u.y - p.y) * (q.x - p.x);
}

bool greater_orient(Point p, Point u1, Point u2, Point q)
{
    return (q.y - p.y) * (u1.x - u2.x) < (q.x - p.x) * (u1.y - u2.y);
}

int main(int argc, char **argv)
{
    if (argc != 1) {
        printf("This program takes no arguments\n");
        return EXIT_FAILURE;
    }

    Point p = {0.50000000000002531, 0.5000000000000171};
    Point q = {24.00000000000005, 24.0000000000000517765};
    Point u = {17.300000000000001, 17.300000000000001};
    Point u2 = {u.x - 0x1p-49, u.y - 0x1p-49};


    assert(u.x != u2.x && u.y != u2.y);

    printf("test %.17e < %.17e\n", (u.x - p.x) * (q.y - p.y), 
                (u.y - p.y) * (q.x - p.x));
    printf("Test: %s\n", right_turn(p, u, q) ? "true" : "false");
    printf("u1 > u2?: %.17e < %.17e: %s\n", 
            (q.y - p.y) * (u.x - u2.x), (q.x - p.x) * (u.y - u2.y),
            greater_orient(p, u, u2, q) ? "true" : "false");

    return EXIT_SUCCESS;
}

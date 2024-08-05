#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cerrno>
#include <sys/time.h>
#include <cfloat>
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

Points input_b_soa(size_t *n /* out param */)
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

    size_t i = Lanes(d);
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
         * elements. 
         * TODO does not work if n < Lanes(d). Partial loads and stores? 
         * Then again, we are probably going to use Points instead of Point,
         * so may be a waste of time to handle this edge case. */
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

    StoreU(leftx,  d, leftx_arr);
    StoreU(rightx, d, rightx_arr);
    StoreU(lefty,  d, lefty_arr);
    StoreU(righty, d, righty_arr);

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

/* Other than eliding the unpack, it we also need only the x-coordinate
 * of most points, resulting in less bandwidth. But that is slightly more
 * involved because we store keep indices instead of points as we may need
 * the y-coordinate for the tie-break. Adepted from 
 * https://en.algorithmica.org/hpc/algorithms/argmin */
void FindLeftRightV_soa(size_t n, Points P, size_t *left_out, size_t *right_out)
{
    const ScalableTag<double> d;
    const ScalableTag<size_t> di;

    Vec<ScalableTag<double>> leftx, lefty, rightx, righty, x_coor, y_coor;
    /* On most amd64 archtitectures (Alderlake, Sapphire, Skylake, zen2-4) 
     * vtestpd has a throughput of 1, wheareas min, max, blend, vcmppd have a 
     * throughput of 0.5 or less. For this reason we unroll the loop. */
    Vec<ScalableTag<double>> x_coor2, y_coor2;
    Vec<ScalableTag<size_t>> l_i, r_i;
    l_i = Iota(di, 0);
    r_i = Iota(di, 0);

    /* Ideally, use LoadN to keep asan happy, but that is from very new
     * hwy version not present in most package managers. */
    // x_coor = LoadN(d, P.x, min(Lanes(d), n));
    // y_coor = LoadN(d, P.y, min(Lanes(d), n));
    auto last_mask = (Iota(d, 0) < Set(d, n));
    x_coor = MaskedLoad(last_mask, d, P.x);
    y_coor = MaskedLoad(last_mask, d, P.y);

    leftx  = x_coor;
    rightx = x_coor;
    lefty  = y_coor;
    righty = y_coor;

    Mask<ScalableTag<double>> mask1, mask2;
    size_t i = Lanes(d);
    for (; i + 2 * Lanes(d) <= n; i += 2 * Lanes(d)) {
        x_coor = LoadU(d, P.x + i);
        x_coor2 = LoadU(d, P.x + i + Lanes(d));
        mask1 = (Min(x_coor, x_coor2) <= leftx);
        mask2 = (Max(x_coor, x_coor2) >= rightx);
        /* Unlikely assuming no adverserial input. */
        if (HWY_UNLIKELY(!AllFalse(d, Or(mask1, mask2)))) {
            y_coor = LoadU(d, P.y + i);
            mask1 = Or((x_coor < leftx), 
                       And((x_coor == leftx), (y_coor < lefty)));
            mask2 = Or((x_coor > rightx), 
                       And((x_coor == rightx), (y_coor > righty)));
            leftx = IfThenElse(mask1,  x_coor, leftx);
            lefty = IfThenElse(mask1,  y_coor, lefty);
            rightx = IfThenElse(mask2, x_coor, rightx);
            righty = IfThenElse(mask2, y_coor, righty);
            l_i = IfThenElse(RebindMask(di, mask1), Iota(di, i), l_i);
            r_i = IfThenElse(RebindMask(di, mask2), Iota(di, i), r_i);

            y_coor2 = LoadU(d, P.y + i + Lanes(d));
            mask1 = Or((x_coor2 < leftx), 
                       And((x_coor2 == leftx), (y_coor2 < lefty)));
            mask2 = Or((x_coor2 > rightx), 
                       And((x_coor2 == rightx), (y_coor2 > righty)));
            leftx = IfThenElse(mask1,  x_coor2, leftx);
            lefty = IfThenElse(mask1,  y_coor2, lefty);
            rightx = IfThenElse(mask2, x_coor2, rightx);
            righty = IfThenElse(mask2, y_coor2, righty);
            l_i = IfThenElse(RebindMask(di, mask1), Iota(di, i + Lanes(d)), l_i);
            r_i = IfThenElse(RebindMask(di, mask2), Iota(di, i + Lanes(d)), r_i);
        }
    }
    for (; i + Lanes(d) <= n; i += Lanes(d)) {
        x_coor = LoadU(d, P.x + i);
        mask1 = (x_coor <= leftx);
        mask2 = (x_coor >= rightx);
        /* Unlikely assuming no adverserial input. */
        if (HWY_UNLIKELY(!AllFalse(d, Or(mask1, mask2)))) {
            y_coor = LoadU(d, P.y + i);
            mask1 = Or((x_coor < leftx), 
                       And((x_coor == leftx), (y_coor < lefty)));
            mask2 = Or((x_coor > rightx), 
                       And((x_coor == rightx), (y_coor > righty)));
            leftx = IfThenElse(mask1,  x_coor, leftx);
            lefty = IfThenElse(mask1,  y_coor, lefty);
            rightx = IfThenElse(mask2, x_coor, rightx);
            righty = IfThenElse(mask2, y_coor, righty);
            l_i = IfThenElse(RebindMask(di, mask1), Iota(di, i), l_i);
            r_i = IfThenElse(RebindMask(di, mask2), Iota(di, i), r_i);
        }
    }
    if (i < n) {
//        x_coor = LoadN(d, P.x + i, n - i);
        auto last_mask = (Iota(d, 0) < Set(d, n - i));
        x_coor = MaskedLoad(last_mask, d, P.x + i);
        mask1 = And((x_coor <= leftx), last_mask);
        mask2 = And((x_coor >= rightx), last_mask);
        /* Unlikely assuming no adverserial input. */
        if (HWY_UNLIKELY(!AllFalse(d, Or(mask1, mask2)))) {
            y_coor = LoadU(d, P.y + i);
            mask1 = And(Or((x_coor < leftx),
                           And((x_coor == leftx), (y_coor < lefty))),
                        last_mask);
            mask2 = And(Or((x_coor > rightx),
                           And((x_coor == rightx), (y_coor > righty))),
                        last_mask);
            leftx = IfThenElse(mask1,  x_coor, leftx);
            lefty = IfThenElse(mask1,  y_coor, lefty);
            rightx = IfThenElse(mask2, x_coor, rightx);
            righty = IfThenElse(mask2, y_coor, righty);
            l_i = IfThenElse(RebindMask(di, mask1), Iota(di, i), l_i);
            r_i = IfThenElse(RebindMask(di, mask2), Iota(di, i), r_i);
        }
    }

    double leftx_arr[Lanes(d)];
    double lefty_arr[Lanes(d)];
    double rightx_arr[Lanes(d)];
    double righty_arr[Lanes(d)];

    StoreU(leftx,  d, leftx_arr);
    StoreU(rightx, d, rightx_arr);
    StoreU(lefty,  d, lefty_arr);
    StoreU(righty, d, righty_arr);

    size_t left_ind = 0;
    size_t right_ind = 0;
    for (size_t i = 1; i < min(Lanes(d), n); i++) {
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

    *left_out  = ExtractLane(l_i, left_ind);
    *right_out = ExtractLane(r_i, right_ind);
}

void MinMaxV(size_t n, Point *P, Point p, Point u, Point q, 
             Point *max1_out, Point *max2_out)
{
    const ScalableTag<double> d;

    Vec<ScalableTag<double>> max1x, max1y, max2x, max2y, x_coor, y_coor, 
                             px, py, ux, uy, qx, qy, omax1, omax2;
    LoadInterleaved2(d, (double *)P, x_coor, y_coor);

    max1x = x_coor;
    max1y = x_coor;
    max2x = y_coor;
    max2y = y_coor;
    px = Set(d, p.x);
    py = Set(d, p.y);
    ux = Set(d, u.x);
    uy = Set(d, u.y);
    qx = Set(d, q.x);
    qy = Set(d, q.y);

    omax1 = orientV(px, py, x_coor, y_coor, ux, uy);
    omax2 = orientV(ux, uy, x_coor, y_coor, qx, qy);

    size_t i = Lanes(d);
    for (; i + Lanes(d) <= n; i += Lanes(d)) {
        LoadInterleaved2(d, (double *)(P + i), x_coor, y_coor);
        auto o1 = orientV(px, py, x_coor, y_coor, ux, uy);
        auto o2 = orientV(ux, uy, x_coor, y_coor, qx, qy);
        max1x = IfThenElse(o1 > omax1, x_coor, max1x);
        max1y = IfThenElse(o1 > omax1, y_coor, max1y);
        max2x = IfThenElse(o2 > omax2, x_coor, max2x);
        max2y = IfThenElse(o2 > omax2, y_coor, max2y);
        omax1 = Max(o1, omax1);
        omax2 = Max(o2, omax2);
    }
    if (i < n) {
        /* Still have < Lanes(d) elements left. It is not a problem to
         * check some elements twice, so we just take the last Lanes(d)
         * elements. */
        i = n - Lanes(d);
        auto o1 = orientV(px, py, x_coor, y_coor, ux, uy);
        auto o2 = orientV(ux, uy, x_coor, y_coor, qx, qy);
        max1x = IfThenElse(o1 > omax1, x_coor, max1x);
        max1y = IfThenElse(o1 > omax1, y_coor, max1y);
        max2x = IfThenElse(o2 > omax2, x_coor, max2x);
        max2y = IfThenElse(o2 > omax2, y_coor, max2y);
        omax1 = Max(o1, omax1);
        omax2 = Max(o2, omax2);
    }

    /* Horizontal reduction. */
    auto hmax1 = MaxOfLanes(d, omax1);
    auto hmax2 = MaxOfLanes(d, omax2);
    auto max1_ind = Iota(d, 0);
    auto max2_ind = Iota(d, 0);
    
    max1_ind = IfThenElseZero(omax1 == hmax1, max1_ind);
    max2_ind = IfThenElseZero(omax2 == hmax2, max2_ind);

    int i1 = GetLane(SumOfLanes(d, max1_ind));
    int i2 = GetLane(SumOfLanes(d, max2_ind));

    max1_out->x = ExtractLane(max1x, i1);
    max1_out->y = ExtractLane(max1y, i1);
    max2_out->x = ExtractLane(max2x, i2);
    max2_out->y = ExtractLane(max2y, i2);
}

void MinMaxV_soa(size_t n, Points P, Point p, Point u, Point q, 
                 Point *max1_out, Point *max2_out)
{
    const ScalableTag<double> d;

    Vec<ScalableTag<double>> max1x, max1y, max2x, max2y, x_coor, y_coor, 
                             px, py, ux, uy, qx, qy, omax1, omax2;

    auto last_mask = (Iota(d, 0) < Set(d, n));
    x_coor = MaskedLoad(last_mask, d, P.x);
    y_coor = MaskedLoad(last_mask, d, P.y);

    max1x = x_coor;
    max1y = x_coor;
    max2x = y_coor;
    max2y = y_coor;
    px = Set(d, p.x);
    py = Set(d, p.y);
    ux = Set(d, u.x);
    uy = Set(d, u.y);
    qx = Set(d, q.x);
    qy = Set(d, q.y);

    omax1 = orientV(px, py, x_coor, y_coor, ux, uy);
    omax2 = orientV(ux, uy, x_coor, y_coor, qx, qy);

    size_t i = Lanes(d);
    for (; i + Lanes(d) <= n; i += Lanes(d)) {
        x_coor = LoadU(d, P.x + i);
        y_coor = LoadU(d, P.y + i);
        auto o1 = orientV(px, py, x_coor, y_coor, ux, uy);
        auto o2 = orientV(ux, uy, x_coor, y_coor, qx, qy);
        max1x = IfThenElse(o1 > omax1, x_coor, max1x);
        max1y = IfThenElse(o1 > omax1, y_coor, max1y);
        max2x = IfThenElse(o2 > omax2, x_coor, max2x);
        max2y = IfThenElse(o2 > omax2, y_coor, max2y);
        omax1 = Max(o1, omax1);
        omax2 = Max(o2, omax2);
    }
    if (i < n) {
    }

    /* Horizontal reduction. */
    auto hmax1 = MaxOfLanes(d, omax1);
    auto hmax2 = MaxOfLanes(d, omax2);
    auto max1_ind = Iota(d, 0);
    auto max2_ind = Iota(d, 0);
    
    max1_ind = IfThenElseZero(omax1 == hmax1, max1_ind);
    max2_ind = IfThenElseZero(omax2 == hmax2, max2_ind);

    int i1 = GetLane(SumOfLanes(d, max1_ind));
    int i2 = GetLane(SumOfLanes(d, max2_ind));

    max1_out->x = ExtractLane(max1x, i1);
    max1_out->y = ExtractLane(max1y, i1);
    max2_out->x = ExtractLane(max2x, i2);
    max2_out->y = ExtractLane(max2y, i2);
}

/* Adepted from
 * https://arxiv.org/pdf/1704.08579
 * and 
 * https://github.com/google/highway/blob/master/hwy/contrib/sort/vqsort-inl.h
 */
void TriPartitionV(size_t n, Points P, Point p, Point u, Point q, 
                   Point *max1_out, Point *max2_out,
                   size_t *c1_out, size_t *c2_out)
{
    const ScalableTag<double> d;

    Vec<ScalableTag<double>> max1x, max1y, max2x, max2y, x_coor, y_coor, 
                             px, py, ux, uy, qx, qy, omax1, omax2,
                             vLx, vLy, vRx, vRy;
    Mask<ScalableTag<double>> maskL, maskR;

    if (n < 2 * Lanes(d)) {
        /* TODO implement scalar base case */
        return;
    }

    omax1 = Set(d, -DBL_MAX);
    omax2 = Set(d, -DBL_MAX);

    px = Set(d, p.x);
    py = Set(d, p.y);
    ux = Set(d, u.x);
    uy = Set(d, u.y);
    qx = Set(d, q.x);
    qy = Set(d, q.y);

    /**
     * Invariant
     *
     *
     *        writeL                                                 writeR
     *         \/                                                      \/
     *  |  S1  | undef | bufferL |   unpartitioned   | bufferR | undef | S2 |
     *                           \/                  \/                     \/
     *                          readL               readR                   num
     *
     **/
    size_t readL = Lanes(d);
    size_t readR = n - Lanes(d);
    size_t writeL = 0;
    size_t writeR = n;

    vLx = LoadU(d, P.x);
    vLy = LoadU(d, P.y);
    vRx = LoadU(d, P.x + readR);
    vRy = LoadU(d, P.y + readR);

    /* 5 gflops on 4 GHz avx2, 8 gflops on 2.6 GHz avx512. 
     *
     * Is this reasonable? It seems slow, but then there are a lot of 
     * instructions other than the actual orientation computation, so may
     * be limited by the decode / data movement operations. */
    while (readL + Lanes(d) <= readR) {
        size_t cap_left = readL - writeL;
        size_t cap_right = writeR - readR;
        if (cap_left <= cap_right) {
            x_coor = LoadU(d, P.x + readL);
            y_coor = LoadU(d, P.y + readL);
            readL += Lanes(d);
        } else {
            readR -= Lanes(d);
            x_coor = LoadU(d, P.x + readR);
            y_coor = LoadU(d, P.y + readR);
        }

        /* Finding r1, r2 */
        auto o1 = orientV(px, py, x_coor, y_coor, ux, uy);
        auto o2 = orientV(ux, uy, x_coor, y_coor, qx, qy);
        /* The if statement is optional, but seems to be a bit faster. */
        if (HWY_UNLIKELY(!AllFalse(d, Or((o1 > omax1), (o2 > omax2))))) {
            max1x = IfThenElse(o1 > omax1, x_coor, max1x);
            max1y = IfThenElse(o1 > omax1, y_coor, max1y);
            max2x = IfThenElse(o2 > omax2, x_coor, max2x);
            max2y = IfThenElse(o2 > omax2, y_coor, max2y);
            omax1 = Max(o1, omax1);
            omax2 = Max(o2, omax2);
        }

        /* Partition */
        maskL = (o1 > Zero(d));
        maskR = And((o2 > Zero(d)), Not(maskL));
        /* No blended store necessary because we write from left to right */
        size_t num_l = CompressStore(x_coor, maskL, d, P.x + writeL);
        CompressStore(y_coor, maskL, d, P.y + writeL);
        writeL += num_l;
        size_t num_r = CountTrue(d, maskR);
        writeR -= num_r;
        CompressBlendedStore(x_coor, maskR, d, P.x + writeR);
        CompressBlendedStore(y_coor, maskR, d, P.y + writeR);
    }

    /* TODO: [readL, readR[ */

    /* TODO: handle vLx, vLy, vRx, vRy */

    /* Horizontal reduction */
    auto hmax1 = MaxOfLanes(d, omax1);
    auto hmax2 = MaxOfLanes(d, omax2);
    auto max1_ind = Iota(d, 0);
    auto max2_ind = Iota(d, 0);
    
    max1_ind = IfThenElseZero(omax1 == hmax1, max1_ind);
    max2_ind = IfThenElseZero(omax2 == hmax2, max2_ind);

    int i1 = GetLane(SumOfLanes(d, max1_ind));
    int i2 = GetLane(SumOfLanes(d, max2_ind));

    max1_out->x = ExtractLane(max1x, i1);
    max1_out->y = ExtractLane(max1y, i1);
    max2_out->x = ExtractLane(max2x, i2);
    max2_out->y = ExtractLane(max2y, i2);

    *c1_out = writeL;
    *c2_out = n - writeR;

    /* Condense */
    memmove(P.x + writeL, P.x + writeR, n - writeR);
    memmove(P.y + writeL, P.y + writeR, n - writeR);
}

double wtime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)(tv.tv_usec / 1e6 + tv.tv_sec);
}

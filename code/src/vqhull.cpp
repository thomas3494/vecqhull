/**
 * vqhull - computes the convex hull
 * Copyright (C) 2024 Thomas Koopman and Jordy Aaldering 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 **/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <float.h>
#include <hwy/highway.h>

#include "vqhull.h"
#include "blockcyclic.h"

using namespace hwy::HWY_NAMESPACE;

/******************************************************************************
 * Type definitions
 *****************************************************************************/

typedef Vec<ScalableTag<double>> Vecd;
typedef Vec<ScalableTag<size_t>> Veci;

typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    double *x;
    double *y;
} Points;

/*****************************************************************************
 * Utility functions
 *****************************************************************************/

static inline size_t 
max(size_t a, size_t b)
{
    return (a > b) ? a : b;
}

static inline size_t 
min(size_t a, size_t b)
{
    return (a < b) ? a : b;
}

static inline void 
swap(Points P, size_t i, size_t j)
{
    double tempx = P.x[i];
    double tempy = P.y[i];
    P.x[i] = P.x[j];
    P.y[i] = P.y[j];
    P.x[j] = tempx;
    P.y[j] = tempy;
}

/*****************************************************************************
 * Geometric tests
 ****************************************************************************/

static inline bool
right_turn(Point p, Point u, Point q)
{
    return (u.x - p.x) * (q.y - p.y) < (u.y - p.y) * (q.x - p.x);
}

static inline Mask<ScalableTag<double>> 
right_turn(Vec<ScalableTag<double>> px,
           Vec<ScalableTag<double>> py,
           Vec<ScalableTag<double>> ux,
           Vec<ScalableTag<double>> uy,
           Vec<ScalableTag<double>> qx,
           Vec<ScalableTag<double>> qy)
{
    return (ux - px) * (qy - py) < (uy - py) * (qx - px);
}

/**
 * Returns true if orient(p, u1, q) > orient(p, u2, q) 
 * This is equivalent to
 *     (qx - px) * (u1y - u2y) > (qy - py) * (u1x - u2x)
 * We do the same technique as with right_turn.
 **/
static inline bool
greater_orient(Point p, Point u1, Point u2, Point q)
{
    return (q.y - p.y) * (u1.x - u2.x) < (q.x - p.x) * (u1.y - u2.y);
}

static inline Mask<ScalableTag<double>>
greater_orient(Vec<ScalableTag<double>> px,  Vec<ScalableTag<double>> py,
               Vec<ScalableTag<double>> u1x, Vec<ScalableTag<double>> u1y,
               Vec<ScalableTag<double>> u2x, Vec<ScalableTag<double>> u2y,
               Vec<ScalableTag<double>> qx,  Vec<ScalableTag<double>> qy)
{
    return (qx - px) * (u1y - u2y) > (qy - py) * (u1x - u2x);
}

/******************************************************************************
 * Finding initial points on the hull
 *****************************************************************************/

/* Adapted from https://en.algorithmica.org/hpc/algorithms/argmin */
static void 
FindLeftRightV(size_t n, Points P, size_t *left_out, size_t *right_out)
{
    const ScalableTag<double> d;
    const ScalableTag<size_t> di;

    Vecd leftx, lefty, rightx, righty, x_coor, y_coor;
    /* On most amd64 archtitectures (Alderlake, Sapphire, Skylake, zen2-4)
     * vtestpd has a throughput of 1, wheareas min, max, blend, vcmppd have a
     * throughput of 0.5 or less. For this reason we unroll the loop. */
    Vecd x_coor2, y_coor2;
    Veci l_i, r_i;
    Mask<ScalableTag<double>> mask1, mask2;
    l_i = Iota(di, 0);
    r_i = Iota(di, 0);

    size_t i = n % Lanes(d);
    x_coor = LoadN(d, P.x, i);
    y_coor = LoadN(d, P.y, i);

    leftx  = IfThenElse(FirstN(d, i), x_coor, Set(d,  DBL_MAX));
    lefty  = IfThenElse(FirstN(d, i), y_coor, Set(d,  DBL_MAX));
    rightx = IfThenElse(FirstN(d, i), x_coor, Set(d, -DBL_MAX));
    righty = IfThenElse(FirstN(d, i), y_coor, Set(d, -DBL_MAX));

    if ((n - i) % (2 * Lanes(d)) == Lanes(d)) {
        x_coor = LoadU(d, P.x + i);
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
        i += Lanes(d);
    }

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
            l_i = IfThenElse(RebindMask(di, mask1),
                             Iota(di, i + Lanes(d)), l_i);
            r_i = IfThenElse(RebindMask(di, mask2),
                             Iota(di, i + Lanes(d)), r_i);
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
    for (size_t j = 1; j < min(Lanes(d), n); j++) {
        if ((leftx_arr[j] < leftx_arr[left_ind]) ||
                ((leftx_arr[j] == leftx_arr[left_ind]) &&
                    (lefty_arr[j] < lefty_arr[left_ind])))
        {
            left_ind = j;
        }
        if ((rightx_arr[j] > rightx_arr[right_ind]) ||
                ((rightx_arr[j] == rightx_arr[right_ind]) &&
                    (righty_arr[j] > righty_arr[right_ind])))
        {
            right_ind = j;
        }
    }

    *left_out  = ExtractLane(l_i, left_ind);
    *right_out = ExtractLane(r_i, right_ind);
}

static void 
FindLeftRightVP(size_t n, Points P, size_t *left_out, size_t *right_out)
{
    unsigned int nthreads;
    size_t block;
    #pragma omp parallel master
    {
        nthreads = omp_get_num_threads();
        block = ceildiv(n, nthreads);
    }

    size_t lefts[nthreads];
    size_t rights[nthreads];

    #pragma omp parallel
    {
        unsigned int me = omp_get_thread_num();
        size_t start = me * block;
        size_t end = min((me + 1) * block, n);
        if (start < end) {
            Points Q = {P.x + start, P.y + start};
            FindLeftRightV(end - start, Q, lefts + me, rights + me);
            lefts[me] += start;
            rights[me] += start;
        }
    }

    size_t left = 0;
    size_t right = 0;
    for (unsigned int i = 1; i < min(n, nthreads); i++) {
        if ((P.x[lefts[i]] < P.x[lefts[left]]) ||
                ((P.x[lefts[i]] == P.x[lefts[left]]) &&
                    (P.y[lefts[i]] < P.y[lefts[left]])))
        {
            left = i;
        }
        if ((P.x[rights[i]] > P.x[rights[right]]) ||
                ((P.x[rights[i]] == P.x[rights[right]]) &&
                    (P.y[rights[i]] > P.y[rights[right]])))
        {
            right = i;
        }
    }

    *left_out = lefts[left];
    *right_out = rights[right];
}

/******************************************************************************
 * Partitioning
 *****************************************************************************/

static inline void
qhull_hmax(Vecd max1x, Vecd max1y,
           Vecd max2x, Vecd max2y,
           Vecd px, Vecd py,
           Vecd rx, Vecd ry,
           Vecd qx, Vecd qy,
           Point *max1_out, Point *max2_out)
{
    const ScalableTag<double> d;

    /**
     * For example, if Lanes(d) = 4
     * | A | B | C | D | ->
     * max(
     *    | A | B | C | D |
     *    | C | D | 0 | 0 |
     *    ) =
     * | max(A, C) | max(B, D) | ... | ... ->
     * max(
     *      | max(A, C) | max(B, D) | ... | ...
     *      | max(B, D) | ...       | ... | ...
     *    ) =
     * | max(max(A, C), max(B, D)) | ... | ... | ...
     **/
    for (int lanes = Lanes(d) / 2; lanes >= 1; lanes /= 2) {
        auto max1x_slide = SlideDownLanes(d, max1x, lanes);
        auto max1y_slide = SlideDownLanes(d, max1y, lanes);
        auto mask1 = greater_orient(px, py, max1x, max1y, 
                                    max1x_slide, max1y_slide, rx, ry);
        max1x = IfThenElse(mask1, max1x, max1x_slide);
        max1y = IfThenElse(mask1, max1y, max1y_slide);

        auto max2x_slide = SlideDownLanes(d, max2x, lanes);
        auto max2y_slide = SlideDownLanes(d, max2y, lanes);
        auto mask2 = greater_orient(rx, ry, max2x, max2y, 
                                    max2x_slide, max2y_slide, qx, qy);
        max2x = IfThenElse(mask2, max2x, max2x_slide);
        max2y = IfThenElse(mask2, max2y, max2y_slide);

    }
    max1_out->x = GetLane(max1x);
    max1_out->y = GetLane(max1y);
    max2_out->x = GetLane(max2x);
    max2_out->y = GetLane(max2y);
}

static inline void 
TriLoopBody(Vecd px, Vecd py,
            Vecd rx, Vecd ry,
            Vecd qx, Vecd qy,
            Vecd &max1x, Vecd &max1y,
            Vecd &max2x, Vecd &max2y,
            Points P, size_t &writeL, size_t &writeR,
            Vecd x_coor, Vecd y_coor)
{
    const ScalableTag<double> d;
    /* Finding r1, r2 */
    auto mask1 = greater_orient(px, py, x_coor, y_coor, max1x, max1y, rx, ry);
    auto mask2 = greater_orient(rx, ry, x_coor, y_coor, max2x, max2y, qx, qy);
    max1x = IfThenElse(mask1, x_coor, max1x);
    max1y = IfThenElse(mask1, y_coor, max1y);
    max2x = IfThenElse(mask2, x_coor, max2x);
    max2y = IfThenElse(mask2, y_coor, max2y);

    /* Partition. */
    auto maskL = right_turn(px, py, x_coor, y_coor, rx, ry);
    auto maskR = right_turn(rx, ry, x_coor, y_coor, qx, qy);
    /* Mathematically it is impossible to make a right turn pur and ruq,
     * but maskL and maskR can both evaluate to true due to rounding error.
     * To avoid out-of-bound writes in that case we do the following check. */
    maskR = And(maskR, Not(maskL));
    /* No blended store necessary because we write from left to right */
    size_t num_l = CompressStore(x_coor, maskL, d, P.x + writeL);
    CompressStore(y_coor, maskL, d, P.y + writeL);
    writeL += num_l;
    size_t num_r = CountTrue(d, maskR);
    writeR -= num_r;
    CompressBlendedStore(x_coor, maskR, d, P.x + writeR);
    CompressBlendedStore(y_coor, maskR, d, P.y + writeR);
}

static inline void 
TriLoopBodyPartial(Vecd px, Vecd py,
                   Vecd rx, Vecd ry,
                   Vecd qx, Vecd qy,
                   Vecd &max1x, Vecd &max1y,
                   Vecd &max2x, Vecd &max2y,
                   Points P, size_t n,
                   size_t &writeL, size_t &writeR,
                   Vecd x_coor, Vecd y_coor)
{
    const ScalableTag<double> d;
    /* Finding r1, r2 */
    auto mask1 = And(greater_orient(px, py, x_coor, y_coor, max1x, max1y,
                                    rx, ry),
                     FirstN(d, n));
    auto mask2 = And(greater_orient(rx, ry, x_coor, y_coor, max2x, max2y,
                                    qx, qy),
                     FirstN(d, n));
    max1x = IfThenElse(mask1, x_coor, max1x);
    max1y = IfThenElse(mask1, y_coor, max1y);
    max2x = IfThenElse(mask2, x_coor, max2x);
    max2y = IfThenElse(mask2, y_coor, max2y);

    /* Partition */
    auto maskL = And(right_turn(px, py, x_coor, y_coor, rx, ry), FirstN(d, n));
    auto maskR = And(right_turn(rx, ry, x_coor, y_coor, qx, qy), FirstN(d, n));
    maskR = And(maskR, Not(maskL));
    size_t num_l = CountTrue(d, maskL);
    CompressBlendedStore(x_coor, maskL, d, P.x + writeL);
    CompressBlendedStore(y_coor, maskL, d, P.y + writeL);
    writeL += num_l;
    size_t num_r = CountTrue(d, maskR);
    writeR -= num_r;
    CompressBlendedStore(x_coor, maskR, d, P.x + writeR);
    CompressBlendedStore(y_coor, maskR, d, P.y + writeR);
}

/* Adapted from
 * https://arxiv.org/pdf/1704.08579
 * and
 * https://github.com/google/highway/blob/master/hwy/contrib/sort/vqsort-inl.h
 */
static void 
TriPartitionV(size_t n, Points P, Point p, Point r, Point q,
              Point *r1_out, Point *r2_out,
              size_t *c1_out, size_t *c2_out)
{
    const ScalableTag<double> d;

    Vecd max1x, max1y, max2x, max2y, x_coor, y_coor,
         px, py, rx, ry, qx, qy,
         vLx, vLy, vRx, vRy;

    px = Set(d, p.x);
    py = Set(d, p.y);
    rx = Set(d, r.x);
    ry = Set(d, r.y);
    qx = Set(d, q.x);
    qy = Set(d, q.y);
    /* We do not have to make a case distinction n < Lanes(d) when
     * computing the horizontal maxes this way
     * because orient(p, r, r) = orient(r, r, q) = 0. */
    max1x = rx;
    max1y = ry;
    max2x = rx;
    max2y = ry;

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

    if (HWY_UNLIKELY((n < Lanes(d)))) {
        x_coor = LoadN(d, P.x, n);
        y_coor = LoadN(d, P.y, n);
        TriLoopBodyPartial(px, py, rx, ry, qx, qy,
                           max1x, max1y, max2x, max2y, P, n,
                           writeL, writeR, x_coor, y_coor);
        qhull_hmax(max1x, max1y, max2x, max2y, px, py, rx, ry, qx, qy,
                   r1_out, r2_out);
        *c1_out = writeL;
        *c2_out = writeR;
        return;
    } else if (HWY_UNLIKELY(n < 2 * Lanes(d))) {
        x_coor = LoadU(d, P.x);
        y_coor = LoadU(d, P.y);
        auto x_coor2 = LoadN(d, P.x + Lanes(d), n % Lanes(d));
        auto y_coor2 = LoadN(d, P.y + Lanes(d), n % Lanes(d));
        TriLoopBody(px, py, rx, ry, qx, qy,
                    max1x, max1y, max2x, max2y, P,
                    writeL, writeR, x_coor, y_coor);
        TriLoopBodyPartial(px, py, rx, ry, qx, qy,
                           max1x, max1y, max2x, max2y, P, n % Lanes(d),
                           writeL, writeR, x_coor2, y_coor2);

        qhull_hmax(max1x, max1y, max2x, max2y, px, py, rx, ry, qx, qy,
                   r1_out, r2_out);

        *c1_out = writeL;
        *c2_out = writeR;

        return;
    }

    vLx = LoadU(d, P.x);
    vLy = LoadU(d, P.y);
    vRx = LoadU(d, P.x + readR);
    vRy = LoadU(d, P.y + readR);

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

        TriLoopBody(px, py, rx, ry, qx, qy,
                    max1x, max1y, max2x, max2y,
                    P, writeL, writeR, x_coor, y_coor);
   }

    x_coor = LoadN(d, P.x + readL, readR - readL);
    y_coor = LoadN(d, P.y + readL, readR - readL);
    TriLoopBodyPartial(px, py, rx, ry, qx, qy,
                       max1x, max1y, max2x, max2y, P, readR - readL,
                       writeL, writeR, x_coor, y_coor);
    TriLoopBody(px, py, rx, ry, qx, qy,
                max1x, max1y, max2x, max2y,
                P, writeL, writeR, vLx, vLy);
    TriLoopBody(px, py, rx, ry, qx, qy,
                max1x, max1y, max2x, max2y,
                P, writeL, writeR, vRx, vRy);

    qhull_hmax(max1x, max1y, max2x, max2y, 
               px, py, rx, ry, qx, qy, r1_out, r2_out);

    *c1_out = writeL;
    *c2_out = writeR;
}

/**
 * Operates on block cyclic subarray of P described by block, nthreads, n.
 * See paper for a picture.
 * Writes num < Lanes(d) elements
 **/
static inline
void Blockcyc_Write(size_t num, BlockCycIndex i, 
                    Vecd x_coor, Vecd y_coor, Points P)
{
    const ScalableTag<double> d;

    if (i.j + num <= i.block) {
        StoreN(x_coor, d, P.x + i.i, num);
        StoreN(y_coor, d, P.y + i.i, num);
    } else {
        size_t countL = i.block - i.j;
        size_t countR = num - countL;
        StoreN(x_coor, d, P.x + i.i, countL);
        StoreN(y_coor, d, P.y + i.i, countL);
        i.k += i.p * i.block;
        x_coor = SlideDownLanes(d, x_coor, countL);
        y_coor = SlideDownLanes(d, y_coor, countL);
        StoreN(x_coor, d, P.x + i.k, countR);
        StoreN(y_coor, d, P.y + i.k, countR);
    }
}

static __attribute__((always_inline)) inline void
Blockcyc_TriLoopBody(Vecd px, Vecd py,
                     Vecd rx, Vecd ry,
                     Vecd qx, Vecd qy,
                     Vecd &max1x, Vecd &max1y,
                     Vecd &max2x, Vecd &max2y,
                     Points P,
                     BlockCycIndex &writeL, BlockCycIndex &writeR,
                     Vecd x_coor, Vecd y_coor)
{
    const ScalableTag<double> d;
    /* Finding r1, r2 */
    auto mask1 = greater_orient(px, py, x_coor, y_coor, max1x, max1y, rx, ry);
    auto mask2 = greater_orient(rx, ry, x_coor, y_coor, max2x, max2y, qx, qy);
    max1x = IfThenElse(mask1, x_coor, max1x);
    max1y = IfThenElse(mask1, y_coor, max1y);
    max2x = IfThenElse(mask2, x_coor, max2x);
    max2y = IfThenElse(mask2, y_coor, max2y);

    /* Partition */
    auto maskL = right_turn(px, py, x_coor, y_coor, rx, ry);
    auto maskR = right_turn(rx, ry, x_coor, y_coor, qx, qy);
    maskR = And(maskR, Not(maskL));
    size_t num_l = CountTrue(d, maskL);
    auto tempxL = Compress(x_coor, maskL);
    auto tempyL = Compress(y_coor, maskL);

    Blockcyc_Write(num_l, writeL, tempxL, tempyL, P);
    writeL += num_l;

    size_t num_r = CountTrue(d, maskR);
    writeR -= num_r;
    auto tempxR = Compress(x_coor, maskR);
    auto tempyR = Compress(y_coor, maskR);
    Blockcyc_Write(num_r, writeR, tempxR, tempyR, P);
}

/**
 * P is an array on n points. We do a block cyclic distribution determined
 * by nthreads and block. So we have indices start, ..., start + block - 1,
 * start + nthreads * block, ..., start + nthreads + block + block - 1, ...
 *
 * We assume we have at least 2 * Lanes(d) points, and start is a multiple
 * of block, n is a multiple of Lanes(d).
 */
static void 
TriPartititionBlockCyc(size_t n, Points P, Point p, Point r, Point q,
                       Point *r1_out, Point *r2_out,
                       size_t *c1_out, size_t *c2_out,
                       size_t *total1_out, size_t *total2_out,
                       size_t start, const size_t block,
                       unsigned int nthreads)
{
    const ScalableTag<double> d;

    Vecd max1x, max1y, max2x, max2y, x_coor, y_coor,
         px, py, rx, ry, qx, qy, vLx, vLy, vRx, vRy;

    px = Set(d, p.x);
    py = Set(d, p.y);
    rx = Set(d, r.x);
    ry = Set(d, r.y);
    qx = Set(d, q.x);
    qy = Set(d, q.y);
    max1x = rx;
    max1y = ry;
    max2x = rx;
    max2y = ry;

    BlockCycIndex writeL = BlockCycBegin(start / block, nthreads, block);
    BlockCycIndex writeR = BlockCycSup(start / block, block, nthreads, n);

    BlockCycIndex last_point = writeR;

    BlockCycIndex readL = writeL;
    BlockCycIndex readR = writeR;

    vLx = LoadU(d, P.x + readL.i);
    vLy = LoadU(d, P.y + readL.i);
    readL += Lanes(d);

    readR -= Lanes(d);
    vRx = LoadU(d, P.x + readR.i);
    vRy = LoadU(d, P.y + readR.i);

    assert(readR.i / block % nthreads == start / block);

    while (readL.i + Lanes(d) <= readR.i) {
        size_t cap_left = readL - writeL;
        size_t cap_right = writeR - readR;
        if (cap_left <= cap_right) {
            x_coor = LoadU(d, P.x + readL.i);
            y_coor = LoadU(d, P.y + readL.i);

            assert(readL.i / block % nthreads == start / block);
            readL += Lanes(d);
        } else {
            readR -= Lanes(d);
            assert(readR.i / block % nthreads == start / block);

            x_coor = LoadU(d, P.x + readR.i);
            y_coor = LoadU(d, P.y + readR.i);
        }

        /* Also advances the write pointer.
         * FIXME(?): less code that way, but you can't tell without inspecting
         * the function. Bad idea? Perhaps pass by pointer instead of pass
         * by reference. */
        assert(writeR.i / block % nthreads == start / block);
        Blockcyc_TriLoopBody(px, py, rx, ry, qx, qy,
                             max1x, max1y, max2x, max2y, P, writeL,
                             writeR, x_coor, y_coor);
        assert(writeR.i / block % nthreads == start / block);
    }

    /* [readL, readR[ is empty because they both start at something
     * divisible by Lanes(d) and are (in/de)creased by Lanes(d) at
     * a time. */
    assert(readL.i == readR.i);

    /* vL, vR */
    assert(writeR.i / block % nthreads == start / block);
    Blockcyc_TriLoopBody(px, py, rx, ry, qx, qy,
                         max1x, max1y, max2x, max2y, P, writeL,
                         writeR, vLx, vLy);

    assert(writeR.i / block % nthreads == start / block);
    Blockcyc_TriLoopBody(px, py, rx, ry, qx, qy,
                         max1x, max1y, max2x, max2y, P, writeL,
                         writeR, vRx, vRy);

    qhull_hmax(max1x, max1y, max2x, max2y, px, py, rx, ry, qx, qy,
               r1_out, r2_out);

    *c1_out = writeL.i;
    *c2_out = writeR.i;

    BlockCycIndex first_point = BlockCycBegin(start / block, nthreads, block);
    *total1_out = writeL - first_point;
    *total2_out = last_point - writeR;
}

/**
 * Invariant: [0, i[   subset S1
 *            [i, j[   undefined
 *            [j, k[   not classified
 *            [k, end[ subset S2
 * Assumes we have set undefined points to {-DBL_MAX, DBL_MAX}
 **/
static void dnf(Points P, size_t start, size_t end,
                Point p, Point r, Point q,
                size_t *i_out, size_t *k_out)
{
    (void)q;
    size_t i = start;
    size_t j = start;
    size_t k = end;

    while (j < k) {
        Point u = {P.x[j], P.y[j]};

        if (u.x == DBL_MAX) {
            j += 1;
        } else if (right_turn(p, u, r)) {
            swap(P, i, j);
            i += 1;
            j += 1;
        } else {
            assert(right_turn(r, u, q));

            k -= 1;
            swap(P, j, k);
        }
    }

    assert(i <= j);
    assert(j == k);
    assert(k <= end);

    *i_out = i;
    *k_out = k;
}

static void dnf_gap(Points P, size_t start, size_t end,
                    size_t gap_start, size_t gap_end,
                    Point p, Point r, Point q,
                    size_t *i_out, size_t *k_out)
{
    (void)q;
    // If this function is applied we are sure to have at least two threads,
    // each thread starts at an offset, so gap_start (c1_max) will always be
    // greater than 0 and gap_start.
    assert(start < gap_start);
    assert(gap_start < gap_end);

    size_t i = start;
    size_t j = start;
    size_t k = end;

    while (j < k) {
        Point u = {P.x[j], P.y[j]};

        if (u.x == DBL_MAX) {
            j += 1;
            if (j == gap_start) {
                j = gap_end;
            }
        } else if (right_turn(p, u, r)) {
            swap(P, i, j);

            i += 1;
            if (i == gap_start) {
                i = gap_end;
            }

            j += 1;
            if (j == gap_start) {
                j = gap_end;
            }
        } else {
            assert(right_turn(r, u, q));
            if (k == gap_end) {
                k = gap_start - 1;
            } else {
                k -= 1;
            }

            swap(P, j, k);
        }
    }

    assert(i <= j);
    assert(j == k);
    assert(k <= end);

    *i_out = i;
    *k_out = k;
}

/* Sets block-cyclic subset of P[start, end) to p. */
static void 
BlockCycSet(Points P, unsigned int t, size_t block, unsigned int nthreads,
            size_t start, size_t end, Point p)
{
    unsigned int start_t = (start / block) % nthreads;
    /* If necessary, round up to first index in I_t */
    size_t i = (start_t == t) ?
                    start :
                    t * block + ceildiv(start - t * block, nthreads * block) *
                                            nthreads * block;
    while (i < end) {
        P.x[i] = p.x;
        P.y[i] = p.y;
        i++;
        if (i % block == 0) {
            i += (nthreads - 1) * block;
        }
    }
}

static void 
TriPartitionP(size_t n, Points P, Point p, Point r, Point q,
              Point *r1_out, Point *r2_out,
              size_t *c1_out, size_t *c2_out,
              unsigned int nthreads)
{
    const ScalableTag<double> d;
    uintptr_t points_per_cacheline = 64 / sizeof(double);

    constexpr size_t block = 8192;
    if (nthreads <= 1 ||
            n + points_per_cacheline + Lanes(d) < block * nthreads)
    {
        TriPartitionV(n, P, p, r, q,
                      r1_out, r2_out, c1_out, c2_out);
        return;
    }

    /* To ensure blocks do not share cachelines */
    uintptr_t n_off = (points_per_cacheline - ((uintptr_t)P.x % 64) /
                                sizeof(double)) % points_per_cacheline;
    assert(n >= n_off);
    n -= n_off;
    P.x += n_off;
    P.y += n_off;

    /* To ensure readL and readR can always load contiguous memory */
    size_t n_end = n % Lanes(d);
    assert(n >= n_end);
    n -= n_end;
    assert(n % Lanes(d) == 0);

    assert((uintptr_t)P.x % 64 == 0);
    assert((uintptr_t)P.y % 64 == 0);

    /* Pad to avoid false sharing */
    size_t total1s[nthreads][8];
    size_t total2s[nthreads][8];
    size_t c1s[nthreads][8];
    size_t c2s[nthreads][8];
    Point r1s[nthreads][4];
    Point r2s[nthreads][4];

    #pragma omp parallel num_threads(nthreads)
    {
        unsigned int me = omp_get_thread_num();
        size_t c1, c2, total1, total2;
        Point r1, r2;

        size_t start = me * block;
        assert(n - me * block >= 2 * Lanes(d));
        TriPartititionBlockCyc(n, P, p, r, q,
                               &r1, &r2, &c1, &c2, &total1, &total2,
                               start, block, nthreads);

        assert(c1 >= start);
        assert(c1 <= n + block * nthreads);
        assert(c2 >= start);
        assert(c2 <= n * block * nthreads);

        c1s[me][0]     = c1;
        c2s[me][0]     = c2;
        total1s[me][0] = total1;
        total2s[me][0] = total2;
        r1s[me][0]     = r1;
        r2s[me][0]     = r2;
    }

    Point r1      = r1s[0][0];
    Point r2      = r2s[0][0];
    size_t total1 = total1s[0][0];
    size_t total2 = total2s[0][0];
    size_t c1_min = c1s[0][0], c1_max = c1s[0][0];
    size_t c2_min = c2s[0][0], c2_max = c2s[0][0];

    for (unsigned int t = 1; t < nthreads; t++) {
        total1 += total1s[t][0];
        total2 += total2s[t][0];

        /* Argmax over empty set is undefined */
        if ((total1s[t][0] != 0) && greater_orient(p, r1s[t][0], r1, r)) {
            r1 = r1s[t][0];
        }
        if ((total2s[t][0] != 0) && greater_orient(r, r2s[t][0], r2, q)) {
            r2 = r2s[t][0];
        }

        c1_min = min(c1_min, c1s[t][0]);
        c1_max = max(c1_max, c1s[t][0]);
        c2_min = min(c2_min, c2s[t][0]);
        c2_max = max(c2_max, c2s[t][0]);
    }

    assert(total1 + total2 <= n);

    /**
     * When working with index sets, say I = {0, ..., n - 1}, we usually
     * do not use the upper bound n - 1, but the smallest strict upper
     * bound n. This is not in I.
     *
     * We have extended this to the block-cyclic subsets I_t. So c2s[t][0]
     * is the smallest strict upper bound ub such that
     *   ub / block % nthreads == t.
     * This can be strictly larger than n.
     */
    c2_min = min(c2_min, n);
    c2_max = min(c2_max, n);
    c1_min = min(c1_min, n);
    c1_max = min(c1_max, n);

    assert(c1_min <= c1_max);
    assert(c2_min <= c2_max);
    assert(c1_min <= c2_max);
    assert(c2_max <= n);

    /* Change undefined points to recognize them later on. */
    #pragma omp parallel num_threads(nthreads)
    {
        unsigned int t = omp_get_thread_num();
        size_t c1 = c1s[t][0];
        size_t c2 = c2s[t][0];
        BlockCycSet(P, t, block, nthreads, c1, min(c1_max, c2),
                    {DBL_MAX, DBL_MAX});
        BlockCycSet(P, t, block, nthreads, max(c2_min, c1), c2,
                    {DBL_MAX, DBL_MAX});
    }

    {
        size_t i, k;
        if (c1_max >= c2_min) {
            dnf(P, c1_min, c2_max, p, r, q, &i, &k);
        } else {
            dnf_gap(P, c1_min, c2_max, c1_max, c2_min, p, r, q, &i, &k);
        }

        assert(total1 == i);
        assert(total2 == n - k);
    }

    /**
     * P now looks like:
     *
     * | S1 | undef | S2     | ... |
     * 0    total1  n-total2 n     n+n_end
     */
    if (n_end > 0) {
        Points LeftOver = {P.x + n, P.y + n};
        Point  r1_left_over, r2_left_over;
        size_t c1_left_over, c2_left_over;
        TriPartitionV(n_end, LeftOver, p, r, q,
                      &r1_left_over, &r2_left_over,
                      &c1_left_over, &c2_left_over);

        assert(c1_left_over <= c2_left_over);
        assert(c2_left_over <= n_end);

        if ((c1_left_over > 0) && greater_orient(p, r1_left_over, r1, r)) {
            r1 = r1_left_over;
        }
        if ((c2_left_over < n_end) && greater_orient(r, r2_left_over, r2, q)) {
            r2 = r2_left_over;
        }

        /**
         * P now looks like:
         *
         * | S1 | undef | S2     | S1 | undef        | S2           |
         * 0    total1  n-total2 n    n+c1_left_over n+c2_left_over n+n_end
         */
        if (total2 > 0) {
            size_t i = n;
            size_t k = n - total2;
            while (i < n + c2_left_over) {
                swap(P, i, k);
                i += 1;
                k += 1;
            }
        }

        /**
         * P now looks like:
         *
         * | S1 | -    | S1     | -          | S2  |
         * 0    total1 n-total2 n-total2+c1L       n+n_end
         */
        if (c1_left_over > 0) {
            size_t src = n - total2;
            size_t dest = total1;
            size_t len = c1_left_over;

            assert(src + len < n + n_end);
            assert(dest + len < n + n_end);

            memmove(P.x + dest, P.x + src, len * sizeof(double));
            memmove(P.y + dest, P.y + src, len * sizeof(double));
        }

        total1 += c1_left_over;
        total2 += (n_end - c2_left_over);
        n += n_end;
    }

    /* Insert the points cut off at the front */
    if (n_off > 0) {
        P.x -= n_off;
        P.y -= n_off;
        Point  r1_left_over, r2_left_over;
        size_t c1_left_over, c2_left_over;
        TriPartitionV(n_off, P, p, r, q,
                      &r1_left_over, &r2_left_over,
                      &c1_left_over, &c2_left_over);

        assert(c1_left_over <= c2_left_over);
        assert(c2_left_over <= n_off);

        if ((c1_left_over > 0) && greater_orient(p, r1_left_over, r1, r)) {
            r1 = r1_left_over;
        }
        if ((c2_left_over < n_off) && greater_orient(r, r2_left_over, r2, q)) {
            r2 = r2_left_over;
        }

        /**
         * P now looks like:
         *
         * | S1 | undef      | S2         | S1  | undef      | S2           |
         * 0    c1_left_over c2_left_over n_off n_off+total1 n+n_off-total2 n+n_off
         */
        if (total1 > 0) {
            size_t i = n_off;
            size_t k = n_off + total1;
            while (i > c1_left_over) {
                i -= 1;
                k -= 1;
                swap(P, i, k);
            }
        }

        /**
         * P now looks like:
         *
         * | S1 | -        | S2       | -          | S2           |
         * 0    total1+c1L total1+c2L total1+n_off n+n_off-total2 n+n_off
         */
        if (n_off > c2_left_over) {
            size_t src = total1 + c2_left_over;
            size_t dest = n - total2 + c2_left_over;
            size_t len = n_off - c2_left_over;

            assert(src + len < n + n_off);
            assert(dest + len < n + n_off);

            memmove(P.x + dest, P.x + src, len * sizeof(double));
            memmove(P.y + dest, P.y + src, len * sizeof(double));
        }

        total1 += c1_left_over;
        total2 += (n_off - c2_left_over);
        n += n_off;
    }

    assert(total1 <= n);
    assert(total2 <= n);

    assert(total1 <= n - total2);
    *c1_out = total1;
    *c2_out = n - total2;
    *r1_out = r1;
    *r2_out = r2;
}

/******************************************************************************
 * Quickhull Algorithm
 *****************************************************************************/

#ifdef MEASURE_BW
size_t read_bw     = 0;
size_t write_bw    = 0;
size_t condense_bw = 0;
#endif

static size_t FindHull(size_t n, Points P, Point p, Point r, Point q);
static size_t FindHullP(size_t n, Points P, Point p, Point r, Point q,
                        unsigned int nthreads);

size_t VecQuickhull(size_t n, double *x_coor, double *y_coor)
{
    Points P = {x_coor, y_coor};
#ifdef MEASURE_BW
    /* Most of the time, we read only x-coordinates */
    read_bw += n * sizeof(Point) / 2;
#endif
    size_t left, right;
    FindLeftRightV(n, P, &left, &right);
    Point p = {P.x[left], P.y[left]};
    Point q = {P.x[right], P.y[right]};

    if (right != 0) {
        swap(P, 0, left);
        swap(P, n - 1, right);
    } else {
        swap(P, n - 1, right);
        swap(P, 0, left);
    }

    Point r1, r2;
    size_t c1, c2;
    Points S1 = {P.x + 1, P.y + 1};
    TriPartitionV(n - 2, S1, p, q, p, &r1, &r2, &c1, &c2);
    size_t total1 = c1;
    size_t total2 = n - 2 - c2;

#ifdef MEASURE_BW
    read_bw  += (n - 2) * sizeof(Point);
    write_bw += (total1 + total2) * sizeof(Point);
#endif

    Points S2 = {S1.x + c2, S1.y + c2};
    assert(total1 <= c2);
    size_t lcount = FindHull(total1, S1, p, r1, q);
    size_t rcount = FindHull(total2, S2, q, r2, p);

    /* Put <p> <left hull> <q> <right hull> together */
    memmove(S1.x + lcount + 1, S2.x, rcount * sizeof(double));
    memmove(S1.y + lcount + 1, S2.y, rcount * sizeof(double));

#ifdef MEASURE_BW
    condense_bw  += 2 * rcount * sizeof(Point);
#endif

    P.x[lcount + 1] = q.x;
    P.y[lcount + 1] = q.y;

#ifdef MEASURE_BW
    printf("Total bw: %zu, read bw: %zu bytes, write bw: %zu, "
           "condese bw: %zu\n",
            read_bw + write_bw + condense_bw, read_bw, write_bw, condense_bw);
#endif

    return 2 + lcount + rcount;
}

//#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
static double wtime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)(tv.tv_usec / 1e6 + tv.tv_sec);
}
#endif

size_t VecQuickhullP(size_t n, double *x_coor, double *y_coor)
{
    Points P = {x_coor, y_coor};
#ifdef PROFILE
    double time1 = wtime();
#endif

    unsigned int nthreads;
    #pragma omp parallel master
    {
        nthreads = omp_get_num_threads();
    }

    omp_set_max_active_levels(nthreads);

    size_t left, right;
    FindLeftRightVP(n, P, &left, &right);
    Point p = {P.x[left], P.y[left]};
    Point q = {P.x[right], P.y[right]};

    if (right != 0) {
        swap(P, 0, left);
        swap(P, n - 1, right);
    } else {
        swap(P, n - 1, right);
        swap(P, 0, left);
    }

#ifdef PROFILE
    double time2 = wtime();
#endif

    Point r1, r2;
    size_t c1, c2;
    Points S1 = {P.x + 1, P.y + 1};
    TriPartitionP(n - 2, S1, p, q, p, &r1, &r2, &c1, &c2, nthreads);
    size_t total1 = c1;
    size_t total2 = n - 2 - c2;
    Points S2 = {S1.x + c2, S1.y + c2};

#ifdef PROFILE
    double time3 = wtime();
#endif

    unsigned int threads1 = roundf((double)total1 /
                                        (total1 + total2) * nthreads);
    unsigned int threads2 = nthreads - threads1;

    size_t lcount, rcount;
    if (threads1 == 0) {
        lcount = FindHull(total1, S1, p, r1, q);
        rcount = FindHullP(total2, S2, q, r2, p, threads2);
    } else if (threads2 == 0) {
        lcount = FindHullP(total1, S1, p, r1, q, threads1);
        rcount = FindHull(total2, S2, q, r2, p);
    } else {
        assert(total1 <= c2);
        #pragma omp parallel shared(threads1, threads2) num_threads(2)
        {
            #pragma omp single nowait
            {
                #pragma omp task
                lcount = FindHullP(total1, S1, p, r1, q, threads1);
            }

            #pragma omp single nowait
            {
                #pragma omp task
                rcount = FindHullP(total2, S2, q, r2, p, threads2);
            }
        }
    }

#ifdef PROFILE
    double time4 = wtime();
#endif

    /* Put <p> <left hull> <q> <right hull> together */
    memmove(S1.x + lcount + 1, S2.x, rcount * sizeof(double));
    memmove(S1.y + lcount + 1, S2.y, rcount * sizeof(double));
    P.x[lcount + 1] = q.x;
    P.y[lcount + 1] = q.y;

#ifdef PROFILE
    double time5 = wtime();
    fprintf(stderr, "Left right: %lf\nPartition: %lf\nRecursion: %lf\n"
                    "Condense: %lf\n", time2 - time1, time3 - time2, 
                    time4 - time3, time5 - time4);
#endif

    return 2 + lcount + rcount;
}

size_t FindHull(size_t n, Points P, Point p, Point r, Point q)
{
    if (n <= 1) return n;

    Point r1, r2;
    size_t c1, c2;
    TriPartitionV(n, P, p, r, q, &r1, &r2, &c1, &c2);
    size_t total1 = c1;
    size_t total2 = n - c2;

#ifdef MEASURE_BW
    read_bw  += n * sizeof(Point);
    write_bw += (total1 + total2) * sizeof(Point);
#endif

    Points S1 = P;
    size_t lcount = FindHull(total1, S1, p, r1, r);

    Points S2 = {P.x + c2, P.y + c2};
    size_t rcount = FindHull(total2, S2, r, r2, q);

    /* Condense left and right hull into contiguous memory */
    memmove(S1.x + lcount + 1, S2.x, rcount * sizeof(double));
    memmove(S1.y + lcount + 1, S2.y, rcount * sizeof(double));

#ifdef MEASURE_BW
    condense_bw  += 2 * rcount * sizeof(Point);
#endif

    P.x[lcount] = r.x;
    P.y[lcount] = r.y;

    return 1 + lcount + rcount;
}

size_t FindHullP(size_t n, Points P, Point p, Point r, Point q,
                 unsigned int nthreads)
{
    if (n <= 1) return n;

    Point r1, r2;
    size_t c1, c2;
    TriPartitionP(n, P, p, r, q, &r1, &r2, &c1, &c2, nthreads);
    size_t total1 = c1;
    size_t total2 = n - c2;

    Points S1 = P;
    Points S2 = {P.x + c2, P.y + c2};

    unsigned int threads1 = (total1 + total2 == 0) ?
                        0 :
                        roundf((double)total1 / (total1 + total2) * nthreads);
    unsigned int threads2 = nthreads - threads1;

    size_t lcount, rcount;
    if (threads1 == 0) {
        lcount = FindHull(total1, S1, p, r1, r);
        rcount = FindHullP(total2, S2, r, r2, q, threads2);
    } else if (threads2 == 0) {
        lcount = FindHullP(total1, S1, p, r1, r, threads1);
        rcount = FindHull(total2, S2, r, r2, q);
    } else {
        assert(total1 <= c2);
        #pragma omp parallel shared(threads1, threads2) num_threads(2)
        {
            #pragma omp single nowait
            {
                #pragma omp task
                {
                    lcount = (threads1 > 1) ?
                                 FindHullP(total1, S1, p, r1, r, threads1) :
                                 FindHull(total1, S1, p, r1, r);
                }
            }

            #pragma omp single nowait
            {
                #pragma omp task
                    rcount = (threads2 > 1) ?
                            FindHullP(total2, S2, r, r2, q, threads2) :
                            FindHull(total2, S2, r, r2, q);
            }
        }
    }

    /* Condense left and right hull into contiguous memory */
    memmove(S1.x + lcount + 1, S2.x, rcount * sizeof(double));
    memmove(S1.y + lcount + 1, S2.y, rcount * sizeof(double));
    P.x[lcount] = r.x;
    P.y[lcount] = r.y;

    return 1 + lcount + rcount;
}

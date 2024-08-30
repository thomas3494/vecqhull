#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cerrno>
#include <omp.h>
#include <sys/time.h>
#include <cfloat>
#include <hwy/highway.h>
#include <hwy/print-inl.h>
#include "common.h"

using namespace hwy::HWY_NAMESPACE;

typedef Vec<ScalableTag<double>> Vecd;
typedef Vec<ScalableTag<size_t>> Veci;

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

void FindLeftRight(size_t n, Points P, size_t *left_out, size_t *right_out)
{
    size_t left = 0;
    size_t right = 0;

    for (size_t i = 1; i < n; i++) {
        if (P.x[i] < P.x[left]
                || (P.x[i] == P.x[left] && P.y[i] < P.y[left]))
        {
            left = i;
        }
        if (P.x[i] > P.x[right] ||
                (P.x[i] == P.x[right] && P.y[i] > P.y[right]))
        {
            right = i;
        }
    }

    *left_out = left;
    *right_out = right;
}

/* Adepted from
 * https://en.algorithmica.org/hpc/algorithms/argmin */
void FindLeftRightV(size_t n, Points P, size_t *left_out, size_t *right_out)
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

void FindLeftRightVP(size_t n, Points P, size_t *left_out, size_t *right_out)
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

static void qhull_hmax(Vecd omax1, Vecd omax2,
                       Vecd max1x, Vecd max1y,
                       Vecd max2x, Vecd max2y,
                       Point *max1_out, Point *max2_out)
{
    const ScalableTag<double> d;

    double omax1_arr[Lanes(d)];
    double omax2_arr[Lanes(d)];
    size_t i1 = 0;
    size_t i2 = 0;
    StoreU(omax1, d, omax1_arr);
    StoreU(omax2, d, omax2_arr);
    for (size_t i = 1; i < Lanes(d); i++) {
        if (omax1_arr[i] > omax1_arr[i1]) {
            i1 = i;
        }
        if (omax2_arr[i] > omax2_arr[i2]) {
            i2 = i;
        }
    }

    max1_out->x = ExtractLane(max1x, i1);
    max1_out->y = ExtractLane(max1y, i1);
    max2_out->x = ExtractLane(max2x, i2);
    max2_out->y = ExtractLane(max2y, i2);
}

static inline void TriLoopBody(Vecd px, Vecd py,
                               Vecd rx, Vecd ry,
                               Vecd qx, Vecd qy,
                               Vecd &omax1, Vecd &omax2,
                               Vecd &max1x, Vecd &max1y,
                               Vecd &max2x, Vecd &max2y,
                               Points P, size_t &writeL, size_t &writeR,
                               Vecd x_coor, Vecd y_coor)
{
    const ScalableTag<double> d;
    /* Finding r1, r2 */
    auto o1 = orientV(px, py, x_coor, y_coor, rx, ry);
    auto o2 = orientV(rx, ry, x_coor, y_coor, qx, qy);
    /* The if statement is optional, but seems to be a bit (5%) faster. */
    if (HWY_UNLIKELY(!AllFalse(d, Or((o1 > omax1), (o2 > omax2))))) {
        max1x = IfThenElse(o1 > omax1, x_coor, max1x);
        max1y = IfThenElse(o1 > omax1, y_coor, max1y);
        max2x = IfThenElse(o2 > omax2, x_coor, max2x);
        max2y = IfThenElse(o2 > omax2, y_coor, max2y);
        omax1 = Max(o1, omax1);
        omax2 = Max(o2, omax2);
    }

    /* Partition */
    auto maskL = (o1 > Zero(d));
    auto maskR = And((o2 > Zero(d)), Not(maskL));
    /* No blended store necessary because we write from left to right */
    size_t num_l = CompressStore(x_coor, maskL, d, P.x + writeL);
    CompressStore(y_coor, maskL, d, P.y + writeL);
    writeL += num_l;
    size_t num_r = CountTrue(d, maskR);
    writeR -= num_r;
    CompressBlendedStore(x_coor, maskR, d, P.x + writeR);
    CompressBlendedStore(y_coor, maskR, d, P.y + writeR);
}

static inline void TriLoopBodyPartial(Vecd px, Vecd py,
                                      Vecd rx, Vecd ry,
                                      Vecd qx, Vecd qy,
                                      Vecd &omax1, Vecd &omax2,
                                      Vecd &max1x, Vecd &max1y,
                                      Vecd &max2x, Vecd &max2y,
                                      Points P, size_t n,
                                      size_t &writeL, size_t &writeR,
                                      Vecd x_coor, Vecd y_coor)
{
    const ScalableTag<double> d;
    /* Finding r1, r2 */
    auto o1 = orientV(px, py, x_coor, y_coor, rx, ry);
    auto o2 = orientV(rx, ry, x_coor, y_coor, qx, qy);
    o1 = IfThenElse(FirstN(d, n), o1, Set(d, -DBL_MAX));
    o2 = IfThenElse(FirstN(d, n), o2, Set(d, -DBL_MAX));
    /* The if statement is optional, but seems to be a bit (5%) faster. */
    if (HWY_UNLIKELY(!AllFalse(d, Or((o1 > omax1), (o2 > omax2))))) {
        max1x = IfThenElse(o1 > omax1, x_coor, max1x);
        max1y = IfThenElse(o1 > omax1, y_coor, max1y);
        max2x = IfThenElse(o2 > omax2, x_coor, max2x);
        max2y = IfThenElse(o2 > omax2, y_coor, max2y);
        omax1 = Max(o1, omax1);
        omax2 = Max(o2, omax2);
    }

    /* Partition */
    auto maskL = (o1 > Zero(d));
    size_t num_l = CountTrue(d, maskL);
    CompressBlendedStore(x_coor, maskL, d, P.x + writeL);
    CompressBlendedStore(y_coor, maskL, d, P.y + writeL);
    writeL += num_l;
    auto maskR = And((o2 > Zero(d)), Not(maskL));
    size_t num_r = CountTrue(d, maskR);
    writeR -= num_r;
    CompressBlendedStore(x_coor, maskR, d, P.x + writeR);
    CompressBlendedStore(y_coor, maskR, d, P.y + writeR);
}

/* Adepted from
 * https://arxiv.org/pdf/1704.08579
 * and
 * https://github.com/google/highway/blob/master/hwy/contrib/sort/vqsort-inl.h
 */
void TriPartitionV(size_t n, Points P, Point p, Point r, Point q,
                   Point *r1_out, Point *r2_out,
                   size_t *c1_out, size_t *c2_out)
{
    const ScalableTag<double> d;

    Vecd max1x, max1y, max2x, max2y, x_coor, y_coor,
         px, py, rx, ry, qx, qy, omax1, omax2,
         vLx, vLy, vRx, vRy;

    omax1 = Set(d, -DBL_MAX);
    omax2 = Set(d, -DBL_MAX);

    px = Set(d, p.x);
    py = Set(d, p.y);
    rx = Set(d, r.x);
    ry = Set(d, r.y);
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

    if (HWY_UNLIKELY((n < Lanes(d)))) {
        x_coor = LoadN(d, P.x, n);
        y_coor = LoadN(d, P.y, n);
        TriLoopBodyPartial(px, py, rx, ry, qx, qy,
                           omax1, omax2, max1x, max1y,
                           max2x, max2y, P, n,
                           writeL, writeR, x_coor, y_coor);
        qhull_hmax(omax1, omax2, max1x, max1y, max2x, max2y,
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
                    omax1, omax2, max1x, max1y,
                    max2x, max2y, P,
                    writeL, writeR, x_coor, y_coor);
        TriLoopBodyPartial(px, py, rx, ry, qx, qy,
                           omax1, omax2, max1x, max1y,
                           max2x, max2y, P, n % Lanes(d),
                           writeL, writeR, x_coor2, y_coor2);

        qhull_hmax(omax1, omax2, max1x, max1y, max2x, max2y,
                   r1_out, r2_out);

        *c1_out = writeL;
        *c2_out = writeR;

        return;
    }

    vLx = LoadU(d, P.x);
    vLy = LoadU(d, P.y);
    vRx = LoadU(d, P.x + readR);
    vRy = LoadU(d, P.y + readR);

    /**
     * 5 gflops on 4 GHz avx2, 8-9 gflops on 2.6 - 4.8 GHz avx512.
     *
     * Is this reasonable? It seems slow, but then there are a lot of
     * instructions other than the actual orientation computation, so may
     * be limited by data movement operations. There are also
     * register spills on avx2.
     *
     * Kuzmin being 10% faster than Circle, Disk is because the
     * (cap_left <= cap_right) branch is easier to predict. Doing it branchless
     * is not worth it on either zen2 or Skylake.
     **/
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

        TriLoopBody(px, py, rx, ry, qx, qy, omax1, omax2,
                    max1x, max1y, max2x, max2y,
                    P, writeL, writeR, x_coor, y_coor);
   }

    x_coor = LoadN(d, P.x + readL, readR - readL);
    y_coor = LoadN(d, P.y + readL, readR - readL);
    TriLoopBodyPartial(px, py, rx, ry, qx, qy,
                       omax1, omax2, max1x, max1y,
                       max2x, max2y, P, readR - readL,
                       writeL, writeR, x_coor, y_coor);
    TriLoopBody(px, py, rx, ry, qx, qy, omax1, omax2,
                max1x, max1y, max2x, max2y,
                P, writeL, writeR, vLx, vLy);
    TriLoopBody(px, py, rx, ry, qx, qy, omax1, omax2,
                max1x, max1y, max2x, max2y,
                P, writeL, writeR, vRx, vRy);

    qhull_hmax(omax1, omax2, max1x, max1y, max2x, max2y, r1_out, r2_out);

    *c1_out = writeL;
    *c2_out = writeR;
}

/**
 * Operates on block cyclic subarray of P described by block, nthreads, n.
 * See paper for a picture.
 * Writes block < Lanes(d) elements starting at i = k + j, where j < block.
 **/
static inline
void Blockcyc_Write(size_t num, size_t k, size_t j, Vecd x_coor, Vecd y_coor,
                    Points P, size_t block, unsigned int nthreads)
{
    const ScalableTag<double> d;

    if (j + num <= block) {
        StoreN(x_coor, d, P.x + k + j, num);
        StoreN(y_coor, d, P.y + k + j, num);
    } else {
        size_t countL = block - j;
        size_t countR = num - countL;
        StoreN(x_coor, d, P.x + k + j, countL);
        StoreN(y_coor, d, P.y + k + j, countL);
        k += nthreads * block;
        x_coor = SlideDownLanes(d, x_coor, countL);
        y_coor = SlideDownLanes(d, y_coor, countL);
        StoreN(x_coor, d, P.x + k, countR);
        StoreN(y_coor, d, P.y + k, countR);
    }
}

/* Assumes count < block */
static inline
void Blockcyc_Sub(size_t &k, size_t &j, size_t count,
                  size_t block, unsigned int nthreads)
{
    if (j >= count) {
        j -= count;
    } else {
        k -= nthreads * block;
        j += block - count;
    }
}

/* Assumes count < block */
static inline
void Blockcyc_Add(size_t &k, size_t &j, size_t count,
                  size_t block, unsigned int nthreads)
{
    if (j + count < block) {
        j += count;
    } else {
        k += nthreads * block;
        j = j + count - block;
    }
}

/**
 * Computes the distance between k1 + j1 and k2 + j2 as if
 * the block cyclic subarray was compacted.
 * Assumes k1 + j1 > k2 + j2.
 **/
static inline
size_t Blockcyc_Dist(size_t k1, size_t j1,
                     size_t k2, size_t j2,
                     unsigned int nthreads)
{
    return (k1 - k2) / nthreads + j1 - j2;
}

/**
 * Copies from block cyclic array src, indices k2 + j2 to k1 + j1,
 * to dense dest. Returns the number of elements copied.
 **/
static inline
size_t Blockcyc_Copy_From(double *dest, double *src,
                          size_t k1, size_t j1,
                          size_t k2, size_t j2,
                          size_t block,
                          unsigned int nthreads)
{
    if (k1 == k2) {
        memcpy(dest, src + k2 + j2, (j1 - j2) * sizeof(double));
        return j1 - j2;
    } else {
        /* Copy to the end of the k2 block */
        size_t copy_count = 0;
        memcpy(dest + copy_count, src + k2 + j2, (block - j2) * sizeof(double));
        copy_count += block - j2;
        k2 += nthreads * block;
        for (; k2 < k1; k2 += nthreads * block) {
            memcpy(dest + copy_count, src + k2, block * sizeof(double));
            copy_count += block;
        }
        memcpy(dest + copy_count, src + k2, j1 * sizeof(double));
        copy_count += j1;
        return copy_count;
    }
}

/**
 * Copies from dense array src to block cyclic array src, starting at k + j.
 **/
static inline
void Blockcyc_Copy_To(double *dest, double *src, size_t count,
                      size_t k, size_t j,
                      size_t block, unsigned int nthreads)
{
    size_t copy_count = 0;
    memcpy(dest + k + j, src + copy_count,
            min(block - j, count) * sizeof(double));
    copy_count += min(block - j, count);
    k += block * nthreads;
    while (copy_count + block <= count) {
        memcpy(dest + k, src + copy_count, block * sizeof(double));
        copy_count += block;
        k += block * nthreads;
    }
    /* copy_count + block > count iff count - copy_count < block */
    memcpy(dest + k, src + copy_count, (count - copy_count) * sizeof(double));
    copy_count += count - copy_count;
}

/**
 * Finds the 'supremum' of i in the subarray belonging to thread t.
 * That is, the smallest number greater or equal to i in
 * t's subarray. We assume this exists, so i >= t * block
 **/
static inline
void Blockcyc_Sup(unsigned int t, size_t i, size_t block, unsigned int nthreads,
                  size_t *k, size_t *j)
{
    if ((i / block) % nthreads == t) {
        /* i is in P(t), so sup is i */
        *j = i % block;
        *k = i - *j;
    } else if (i < t * block) {
        *k = t * block;
        *j = 0;
    } else {
        /* sup is smallest t * block + l * nthreads * block >= i. */
        size_t l = ceildiv(i - t * block, nthreads * block);
        *k = t * block + l * nthreads * block;
        *j = 0;
    }
}

/* writes do not seem disjoint */
static inline void Blockcyc_TriLoopBody(Vecd px, Vecd py,
                                        Vecd rx, Vecd ry,
                                        Vecd qx, Vecd qy,
                                        Vecd &omax1, Vecd &omax2,
                                        Vecd &max1x, Vecd &max1y,
                                        Vecd &max2x, Vecd &max2y,
                                        Points P,
                                        size_t &writeLk, size_t &writeLj,
                                        size_t &writeRk, size_t &writeRj,
                                        Vecd x_coor, Vecd y_coor,
                                        size_t block, unsigned int nthreads)
{
    const ScalableTag<double> d;
    /* Finding r1, r2 */
    auto o1 = orientV(px, py, x_coor, y_coor, rx, ry);
    auto o2 = orientV(rx, ry, x_coor, y_coor, qx, qy);
    /* The if statement is optional, but seems to be a bit (5%) faster. */
    if (HWY_UNLIKELY(!AllFalse(d, Or((o1 > omax1), (o2 > omax2))))) {
        max1x = IfThenElse(o1 > omax1, x_coor, max1x);
        max1y = IfThenElse(o1 > omax1, y_coor, max1y);
        max2x = IfThenElse(o2 > omax2, x_coor, max2x);
        max2y = IfThenElse(o2 > omax2, y_coor, max2y);
        omax1 = Max(o1, omax1);
        omax2 = Max(o2, omax2);
    }

    /* Partition */
    auto maskL = (o1 > Zero(d));
    size_t num_l = CountTrue(d, maskL);
    auto tempxL = Compress(x_coor, maskL);
    auto tempyL = Compress(y_coor, maskL);
//    printf("Write %zu\n", writeLk + writeLj);
    Blockcyc_Write(num_l, writeLk, writeLj, tempxL, tempyL, P, block, nthreads);
    Blockcyc_Add(writeLk, writeLj, num_l, block, nthreads);

    auto maskR = And((o2 > Zero(d)), Not(maskL));
    size_t num_r = CountTrue(d, maskR);
    Blockcyc_Sub(writeRk, writeRj, num_r, block, nthreads);
    auto tempxR = Compress(x_coor, maskR);
    auto tempyR = Compress(y_coor, maskR);
    Blockcyc_Write(num_r, writeRk, writeRj, tempxR, tempyR, P,
                   block, nthreads);
}

/**
 * TODO: this function is incorrect, it finds too little points. Smallest
 * reproducible example is 10000 points disk, second partition.
 *
 * I have checked that on the first partition, they produce the exact same
 * points.
 *
 * Also reproducible for n = 10007 and block = 4. When n is not divisible by
 * Lanes?
 *
 * P is an array on n points. We do a block cyclic distribution determined
 * by nthreads and block. So we have indices start, ..., start + block - 1,
 * start + nthreads * block, ..., start + nthreads + block + block - 1, ...
 *
 * We assume we have at least 2 * Lanes(d) points, and start is a multiple
 * of block, n is a multiple of Lanes(d).
 **/
void TriPartititionBlockCyc(size_t n, Points P, Point p, Point r, Point q,
                            Point *r1_out, Point *r2_out,
                            size_t *c1_out, size_t *c2_out,
                            size_t *total1_out, size_t *total2_out,
                            size_t start, const size_t block,
                            unsigned int nthreads)
{
    const ScalableTag<double> d;

    Vecd max1x, max1y, max2x, max2y, x_coor, y_coor,
         px, py, rx, ry, qx, qy, omax1, omax2,
         vLx, vLy, vRx, vRy;

    px = Set(d, p.x);
    py = Set(d, p.y);
    rx = Set(d, r.x);
    ry = Set(d, r.y);
    qx = Set(d, q.x);
    qy = Set(d, q.y);

    size_t writeLk = start;
    size_t writeLj = 0;

    size_t writeRk, writeRj;
    Blockcyc_Sup(start / block, n, block, nthreads, &writeRk, &writeRj);
//    printf("Last pointer is %zu\n", writeRk + writeRj);
//    Blockcyc_Sub(writeRk, writeRj, Lanes(d), block, nthreads);
//    printf("One to the left is %zu\n", writeRk + writeRj);
//    Blockcyc_Sub(writeRk, writeRj, Lanes(d), block, nthreads);
//    printf("One to the left of that is %zu\n", writeRk + writeRj);
//    Blockcyc_Add(writeRk, writeRj, Lanes(d), block, nthreads);
//    Blockcyc_Add(writeRk, writeRj, Lanes(d), block, nthreads);

    size_t last_pointk = writeRk;
    size_t last_pointj = writeRj;

    size_t readLk = writeLk;
    size_t readLj = writeLj;
    size_t readRk = writeRk;
    size_t readRj = writeRj;

    assert((readLk + readLj) / block % nthreads == start / block);
//    printf("Read %zu\n", readLk + readLj);
    vLx = LoadU(d, P.x + readLk + readLj);
    vLy = LoadU(d, P.y + readLk + readLj);
    Blockcyc_Add(readLk, readLj, Lanes(d), block, nthreads);

    Blockcyc_Sub(readRk, readRj, Lanes(d), block, nthreads);
    vRx = LoadU(d, P.x + readRk + readRj);
    vRy = LoadU(d, P.y + readRk + readRj);
    assert((readRk + readRj) / block % nthreads == start / block);
//    printf("Read %zu\n", readRk + readRj);

    while (readLk + readLj + Lanes(d) <= readRk + readRj) {
        size_t cap_left = Blockcyc_Dist(readLk, readLj, writeLk, writeLj,
                                        nthreads);
        size_t cap_right = Blockcyc_Dist(writeRk, writeRj, readRk, readRj,
                                         nthreads);
        if (cap_left <= cap_right) {
            x_coor = LoadU(d, P.x + readLk + readLj);
            y_coor = LoadU(d, P.y + readLk + readLj);
            assert((readLk + readLj) / block % nthreads == start / block);
//            printf("Read %zu\n", readLk + readLj);
            Blockcyc_Add(readLk, readLj, Lanes(d), block, nthreads);
        } else {
            Blockcyc_Sub(readRk, readRj, Lanes(d), block, nthreads);
            assert((readRk + readRj) / block % nthreads == start / block);
//            printf("Read %zu\n", readRk + readRj);
            x_coor = LoadU(d, P.x + readRk + readRj);
            y_coor = LoadU(d, P.y + readRk + readRj);
        }

        /* Also advances the write pointer.
         * FIXME(?): less code that way, but you can't tell without inspecting
         * the function. Bad idea? Perhaps pass by pointer instead of pass
         * by reference. */
        assert((writeRk + writeRj) / block % nthreads == start / block);
        Blockcyc_TriLoopBody(px, py, rx, ry, qx, qy,
                             omax1, omax2, max1x, max1y,
                             max2x, max2y, P, writeLk, writeLj,
                             writeRk, writeRj, x_coor, y_coor,
                             block, nthreads);
        assert((writeRk + writeRj) / block % nthreads == start / block);
    }

    /* [readL, readR[ is empty because they both start at something
     * divisible by Lanes(d) and are (in/de)creased by Lanes(d) at
     * a time. */
    assert(readLk + readLj == readRk + readRj);

    /* vL, vR */
    assert((writeRk + writeRj) / block % nthreads == start / block);
    Blockcyc_TriLoopBody(px, py, rx, ry, qx, qy,
                         omax1, omax2, max1x, max1y,
                         max2x, max2y, P, writeLk, writeLj,
                         writeRk, writeRj, vLx, vLy,
                         block, nthreads);

    assert((writeRk + writeRj) / block % nthreads == start / block);
    Blockcyc_TriLoopBody(px, py, rx, ry, qx, qy,
                         omax1, omax2, max1x, max1y,
                         max2x, max2y, P, writeLk, writeLj,
                         writeRk, writeRj, vRx, vRy,
                         block, nthreads);

    qhull_hmax(omax1, omax2, max1x, max1y, max2x, max2y, r1_out, r2_out);

    *c1_out = writeLk + writeLj;
    *c2_out = writeRk + writeRj;

    *total1_out = Blockcyc_Dist(writeLk, writeLj, start, 0, nthreads);
    *total2_out = Blockcyc_Dist(last_pointk, last_pointj, writeRk, writeRj,
                                nthreads);
}

static void dnf(Points P, size_t c1s[][8], size_t c2s[][8],
                unsigned int nthreads, const size_t block,
                size_t start, size_t end,
                size_t *i_out, size_t *j_out)
{
    // We might have a situation where (j=0) <= (k=0),
    // after which k -= 1 might occur, so we use int64_t instead of size_t.
    int64_t i = start;
    int64_t j = start;
    int64_t k = end - 1;

    while (j <= k) {
        unsigned int kt = (k / block) % nthreads;
        bool k_in_s1 = ((size_t)k < c1s[kt][0]);
        bool k_in_s2 = ((size_t)k >= c2s[kt][0]);

        while (j <= k && !(k_in_s1 || k_in_s2)) {
            k -= 1;

            if (k < 0) {
                break;
            }

            k_in_s1 = ((size_t)k < c1s[kt][0]);
            k_in_s2 = ((size_t)k >= c2s[kt][0]);
        }

        if (j > k) {
            break;
        }

        unsigned int jt = (j / block) % nthreads;
        bool j_in_s1 = ((size_t)j < c1s[jt][0]);
        bool j_in_s2 = ((size_t)j >= c2s[jt][0]);

        if (j_in_s1) {
            // Swap P[i] and P[j]
            double swap_x = P.x[j];
            double swap_y = P.y[j];
            P.x[j] = P.x[i];
            P.y[j] = P.y[i];
            P.x[i] = swap_x;
            P.y[i] = swap_y;

            i += 1;
            j += 1;
        } else if (j_in_s2) {
            j += 1;
        } else {
            if (k_in_s1) {
                // P[j] = P[i]
                P.x[j] = P.x[i];
                P.y[j] = P.y[i];
                // P[i] = P[k]
                P.x[i] = P.x[k];
                P.y[i] = P.y[k];

                i += 1;
                j += 1;
                k -= 1;
            } else {
                // P[j] = P[k]
                P.x[j] = P.x[k];
                P.y[j] = P.y[k];

                j += 1;
                k -= 1;
            }
        }
    }

    assert(i >= 0);
    assert(j >= i);
    assert(k == j - 1);

    printf("i = %zu, j = %zu\n", i, j);
    *i_out = (size_t)i;
    *j_out = (size_t)j;
}

/**
 * On the cluster it is better to use 8 threads instead of 16.
 *
 * The right-hand side is not correct
 **/
void TriPartitionP(size_t n, Points P, Point p, Point r, Point q,
                   Point *r1_out, Point *r2_out,
                   size_t *c1_out, size_t *c2_out,
                   unsigned int nthreads)
{
    constexpr size_t block = 4;
    if (nthreads <= 1 || ceildiv(n, block) < nthreads) {
        TriPartitionV(n, P, p, r, q,
                      r1_out, r2_out, c1_out, c2_out);
        return;
    }

    const ScalableTag<double> d;
    size_t n_end = n % Lanes(d);
    n -= n_end;

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

        TriPartititionBlockCyc(n, P, p, r, q,
                               &r1, &r2, &c1, &c2, &total1, &total2,
                               me * block, block, nthreads);

        printf("Thread %u, S1 = [%zu, %zu), S2 = [%zu, n), %zu, %zu elem\n",
               me, me * block, c1, c2, total1, total2);

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
    size_t c2_min = c1s[0][0], c2_max = c2s[0][0];

    for (unsigned int t = 1; t < nthreads; t++) {
        total1 += total1s[t][0];
        total2 += total2s[t][0];

        if (orient(p, r1s[t][0], r) > orient(p, r1, r)) {
            r1 = r1s[t][0];
        }
        if (orient(r, r2s[t][0], q) > orient(r, r2, q)) {
            r2 = r2s[t][0];
        }

        c1_min = min(c1_min, c1s[t][0]);
        c1_max = max(c1_max, c1s[t][0]);
        c2_min = min(c2_min, c2s[t][0]);
        c2_max = max(c2_max, c2s[t][0]);
    }

    assert(c1_min <= c1_max);
    assert(c2_min <= c2_max);

    size_t c1 = total1;
    size_t c2 = n - total2;
    printf("c1 = %zu; c2 = %zu\n", c1, c2);
    printf("c1_min = %zu; c1_max = %zu; c2_min = %zu; c2_max = %zu\n",
           c1_min, c1_max, c2_min, c2_max);

    if (c1_max > c2_min) {
        /**
         * P now looks like:
         *
         * | S1 | ?    | ?    | ?    | S2   | ... |
         * 0    c1_min c2_min c1_max c2_max n     n+n_end
         *
         * We apply dnf to [c1_min, c2_max)
         */
        size_t i, j;
        dnf(P, c1s, c2s, nthreads, block, c1_min, c2_max, &i, &j);

        /**
         * P now looks like:
         *
         * | S1 | S1   | S2 | undef | S2   | ... |
         * 0    c1_min i    j       c2_max n     n+n_end
         *
         * We move [i, j) to [c2_max - (j - i), c2_max)
         */
        if (j > i) {
            size_t len = j - i;
            memmove(P.x + i, P.x + c2_max - len, len * sizeof(double));
            memmove(P.y + i, P.y + c2_max - len, len * sizeof(double));
        }

        /**
         * P now looks like:
         *
         * | S1 | S1   | undef  | S2         | S2   | ... |
         * 0    c1_min i        c2_max-(j-i) c2_max n     n+n_end
         */
        // TODO: these assertions sometimes fail???
        //assert(c1 == i);
        //assert(c2 == c2_max - (j - i));
    } else /* c1_max <= c2_min */ {
        /**
         * P now looks like:
         *
         * | S1 | ?    | undef | ?    | S2   | ... |
         * 0    c1_min c1_max  c2_min c2_max n     n+n_end
         *
         * We apply dnf twice:
         * - once on [c1_min, c1_max),
         * - and once on [c2_min, c2_max)
         */
        size_t i1, j1;
        dnf(P, c1s, c2s, nthreads, block, c1_min, c1_max, &i1, &j1);

        /**
         * P now looks like:
         *
         * | S1 | S1   | S2 | undef | undef | ?    | S2   | ... |
         * 0    c1_min i1   j1      c1_max  c2_min c2_max n     n+n_end
         */
        size_t i2, j2;
        dnf(P, c1s, c2s, nthreads, block, c2_min, c2_max, &i2, &j2);

        /**
         * P now looks like:
         *
         * | S1 | S1   | S2 | undef | undef | S1   | S2 | undef | S2   | ... |
         * 0    c1_min i1   j1      c1_max  c2_min i2   j2      c2_max n     n+n_end
         *
         * - Buffer [c2_min, i2)
         * - Move [i2, j2) to [c2_max - (j2 - i2), c2_max)
         * - Move [i1, j1) to [c2_max - (j2 - i2) - (j1 - i1), c2_max - (j2 - i2))
         * - Move buffer to [i1, i1 + (i2 - c2_min))
         */
        double *buf_x, *buf_y;
        if (i2 > c2_min) {
            // Buffer [c2_min, i2)
            size_t len = i2 - c2_min;
            printf("Moving %zu elems to buffer\n", len);
            buf_x = (double *)malloc(len * sizeof(double));
            buf_y = (double *)malloc(len * sizeof(double));
            memmove(P.x + c2_min, buf_x, len * sizeof(double));
            memmove(P.y + c2_min, buf_y, len * sizeof(double));
        }

        if (j2 > i2) {
            // Move [i2, j2) to [c2_max - (j2 - i2), c2_max)
            size_t len = j2 - i2;
            printf("Moving %zu S2(r) elems to the right\n", len);
            memmove(P.x + i2, P.x + c2_max - len, len * sizeof(double));
            memmove(P.y + i2, P.y + c2_max - len, len * sizeof(double));
        }

        if (j1 > i1) {
            // Move [i1, j1) to [c2_max - (j2 - i2) - (j1 - i1), c2_max - (j2 - i2))
            size_t len = j1 - i1;
            printf("Moving %zu S2(l) elems to the right\n", len);
            memmove(P.x + i1, P.x + c2_max - (j2 - i2) - len, len * sizeof(double));
            memmove(P.y + i1, P.y + c2_max - (j2 - i2) - len, len * sizeof(double));
        }

        if (i2 > c2_min) {
            // Move buffer to [i1, i1 + (i2 - c2_min))
            size_t len = i2 - c2_min;
            printf("Moving %zu elems from buffer\n", len);
            memmove(buf_x, P.x + i1, len * sizeof(double));
            memmove(buf_y, P.y + i1, len * sizeof(double));
            free(buf_x);
            free(buf_y);
        }

        /**
         * P now looks like:
         *
         * | S1 | undef        | S2                   | ... |
         * 0    i1+(i2-c2_min) c2_max-(j2-i2)-(j1-i1) n     n+n_end
         */
        assert(c1 == i1 + (i2 - c2_min));
        assert(c2 == c2_max - (j2 - i2) - (j1 - i1));
    }

    /**
     * P now looks like:
     *
     * | S1 | undef | S2 | ... |
     * 0    c1      c2   n     n+n_end
     */

    // TODO: can you please do this bit @Thomas?

    /**
     * P now looks like:
     *
     * | S1 | undef | S2 | S1 | S2 | undef |
     * 0    c1      c2   n                 n+n_end
     */

    *c1_out = c1;
    *c2_out = c2;
}

double wtime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)(tv.tv_usec / 1e6 + tv.tv_sec);
}

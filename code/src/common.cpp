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

double wtime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)(tv.tv_usec / 1e6 + tv.tv_sec);
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

/* Adapted from
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

/* Adapted from
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

static bool in_s1(size_t idx, size_t c1s[][8], unsigned int nthreads, size_t block)
{
    assert(nthreads > 0);
    assert(block > 0);
    unsigned int thread_idx = (idx / block) % nthreads;
    assert(thread_idx < nthreads);
    size_t thread_s1_end = c1s[thread_idx][0];
    return idx < thread_s1_end;
}

static bool in_s2(size_t idx, size_t c2s[][8], unsigned int nthreads, size_t block)
{
    assert(nthreads > 0);
    assert(block > 0);
    unsigned int thread_idx = (idx / block) % nthreads;
    assert(thread_idx < nthreads);
    size_t thread_s2_start = c2s[thread_idx][0];
    return idx >= thread_s2_start;
}

static bool in_undef(bool in_s1, bool in_s2)
{
    return !(in_s1 || in_s2);
}

static void dnf(Points P, size_t c1s[][8], size_t c2s[][8],
                unsigned int nthreads, const size_t block,
                size_t start, size_t end,
                size_t *i_out, size_t *j_out)
{
    assert(start < end);

    size_t i = start;
    size_t j = start;
    size_t k = end - 1;

    while (j <= k) {
        bool k_in_s1 = in_s1(k, c1s, nthreads, block);
        bool k_in_s2 = in_s2(k, c2s, nthreads, block);
        assert(!(k_in_s1 && k_in_s2));
        bool k_undef = in_undef(k_in_s1, k_in_s2);

        while (k_undef && j <= k && k > 0) {
            k -= 1;

            k_in_s1 = in_s1(k, c1s, nthreads, block);
            k_in_s2 = in_s2(k, c2s, nthreads, block);
            k_undef = in_undef(k_in_s1, k_in_s2);
        }

        if (j > k) {
            break;
        }

        assert(k_in_s1 || k_in_s2);
        assert(!(k_in_s1 && k_in_s2));
        assert(i >= start);
        assert(j >= start);
        assert(j >= i);
        assert(k >= j);
        assert(k < end);

        bool j_in_s1 = in_s1(j, c1s, nthreads, block);
        bool j_in_s2 = in_s2(j, c2s, nthreads, block);
        assert(!(j_in_s1 && j_in_s2));

        if (j_in_s1) {
            // Swap P[i] and P[j]
            swap(P, i, j);
            i += 1;
            j += 1;
        } else if (j_in_s2) {
            j += 1;
        } else /* j is neither in s1 nor in s2 */ {
            assert(!j_in_s1 && !j_in_s2);

            if (k_in_s1) {
                // We might have j == k, so we need a swap
                double swap_x = P.x[k];
                double swap_y = P.y[k];
                // P[j] = P[i]
                P.x[j] = P.x[i];
                P.y[j] = P.y[i];
                // P[i] = P[k]
                P.x[i] = swap_x;
                P.y[i] = swap_y;

                i += 1;
                j += 1;
            } else /* k is in s2 */ {
                assert(k_in_s2);

                // P[j] = P[k]
                P.x[j] = P.x[k];
                P.y[j] = P.y[k];

                j += 1;
            }

            // Break early to ensure that k does not underflow
            if (k == 0) {
                break;
            }

            k -= 1;
        }
    }

    //printf("i = %zu, j = %zu\n", i, j);

    assert(i >= 0);
    assert(j >= i);
    assert((j == 0 && k == 0) ^ (k == j - 1));

    *i_out = i;
    *j_out = j;
}

static void MovePointsLocal(Points P, size_t src, size_t dest, size_t len)
{
    memmove(P.x + dest, P.x + src, len * sizeof(double));
    memmove(P.y + dest, P.y + src, len * sizeof(double));
}

static void CopyPointsToBuffer(Points P, size_t src, size_t len,
                               double **buf_x, double **buf_y)
{
    *buf_x = (double *)malloc(len * sizeof(double));
    *buf_y = (double *)malloc(len * sizeof(double));
    memcpy(*buf_x, P.x + src, len * sizeof(double));
    memcpy(*buf_y, P.y + src, len * sizeof(double));
}

static void CopyPointsFromBuffer(Points P, size_t dest, size_t len,
                                 double **buf_x, double **buf_y)
{
    memcpy(P.x + dest, *buf_x, len * sizeof(double));
    memcpy(P.y + dest, *buf_y, len * sizeof(double));
    free(*buf_x);
    free(*buf_y);
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
    constexpr size_t block = 8;
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

        //printf("Thread %2u, S1 = [%3zu, %zu), S2 = [%zu, n), %zu, %zu elem\n",
        //       me, me * block, c1, c2, total1, total2);

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

    //printf("n = %zu; n_end = %zu\n", n, n_end);
    //printf("total1 = %zu; total2 = %zu\n", total1, total2);
    //printf("c1_min = %zu; c1_max = %zu\n", c1_min, c1_max);
    //printf("c2_min = %zu; c2_max = %zu\n", c2_min, c2_max);
    assert(c1_min <= c1_max);
    assert(c2_min <= c2_max);
    assert(c1_min < c2_max);

    if (c1_max >= c2_min) {
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
        size_t len = j - i;

        /**
         * P now looks like:
         *
         * | S1 | S1   | S2 | undef | S2   | ... |
         * 0    c1_min i    j       c2_max n     n+n_end
         *
         * Namely:
         * - P[0, i) in S1
         * - P[i, j) in S2
         * - P[j, c2_max) undef
         *
         * We move [i, j) to [c2_max - (j - i), c2_max)
         */
        if (len > 0) {
            MovePointsLocal(P, i, c2_max - len, len);
        }

        /**
         * P now looks like:
         *
         * | S1 | S1   | undef  | S2         | S2   | ... |
         * 0    c1_min i        c2_max-(j-i) c2_max n     n+n_end
         */
        assert(total1 == i);
        assert(n - total2 == c2_max - len);
    } else /* c1_max < c2_min */ {
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
        size_t len1 = j1 - i1;

        /**
         * P now looks like:
         *
         * | S1 | S1   | S2 | undef | undef | ?    | S2   | ... |
         * 0    c1_min i1   j1      c1_max  c2_min c2_max n     n+n_end
         */
        size_t i2, j2;
        dnf(P, c1s, c2s, nthreads, block, c2_min, c2_max, &i2, &j2);
        size_t len2 = j2 - i2;

        /**
         * P now looks like:
         *
         * | S1 | S1   | S2 | undef | undef | S1   | S2 | undef | S2   | ... |
         * 0    c1_min i1   j1      c1_max  c2_min i2   j2      c2_max n     n+n_end
         *
         * - Move [c2_min, i2) to buffer
         * - Move [i2, j2) to [c2_max - (j2 - i2), c2_max)
         * - Move [i1, j1) to [c2_max - (j2 - i2) - (j1 - i1), c2_max - (j2 - i2))
         * - Move buffer to [i1, i1 + (i2 - c2_min))
         */
        double *buf_x, *buf_y;
        size_t buf_len = i2 - c2_min;
        if (buf_len > 0) {
            // Move [c2_min, i2) to buffer
            CopyPointsToBuffer(P, c2_min, buf_len, &buf_x, &buf_y);
        }

        if (len2 > 0) {
            // Move [i2, j2) to [c2_max - (j2 - i2), c2_max)
            MovePointsLocal(P, i2, c2_max - len2, len2);
        }

        if (len1 > 0) {
            // Move [i1, j1) to [c2_max - (j2 - i2) - (j1 - i1), c2_max - (j2 - i2))
            MovePointsLocal(P, i1, c2_max - len2 - len1, len1);
        }

        if (buf_len) {
            // Move buffer to [i1, i1 + (i2 - c2_min))
            CopyPointsFromBuffer(P, i1, buf_len, &buf_x, &buf_y);
        }

        /**
         * P now looks like:
         *
         * | S1 | undef        | S2                   | ... |
         * 0    i1+(i2-c2_min) c2_max-(j2-i2)-(j1-i1) n     n+n_end
         */
        assert(total1 == i1 + buf_len);
        assert(n - total2 == c2_max - len2 - len1);
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

        if ((c1_left_over > 0) && (orient(p, r1_left_over, r) > orient(p, r1, r))) {
            r1 = r1_left_over;
        }
        if ((c2_left_over < n_end) && (orient(r, r2_left_over, q) > orient(r, r2, q))) {
            r2 = r2_left_over;
        }

        //printf("c1_left_over = %zu; c2_left_over = %zu; n_end = %zu\n", c1_left_over, c2_left_over, n_end);
        assert(c1_left_over <= c2_left_over);
        assert(c2_left_over <= n_end);

        /**
         * P now looks like:
         *
         * | S1 | undef | S2     | S1 | undef        | S2           |
         * 0    total1  n-total2 n    n+c1_left_over n+c2_left_over n+n_end
         *
         * - Move [n,n+c1_left_over) to buffer
         * - Move [n-total2,n) to [n+c2_left_over-total2,n+c2_left_over)
         * - Move buffer to [total1,total1+c1_left_over)
         */
        double *buf_x, *buf_y;
        if (c1_left_over > 0) {
            // Move [n,n+c1_left_over) to buffer
            CopyPointsToBuffer(P, n, c1_left_over, &buf_x, &buf_y);
        }

        if (total2 > 0) {
            // Move [n-total2,n) to [n+c2_left_over-total2,n+c2_left_over)
            MovePointsLocal(P, n - total2, n + c2_left_over - total2, total2);
        }

        if (c1_left_over > 0) {
            // Move buffer to [total1,total1+c1_left_over)
            CopyPointsFromBuffer(P, total1, c1_left_over, &buf_x, &buf_y);
        }

        total1 += c1_left_over;
        total2 += n_end - c2_left_over;
    }

    //printf("S1: [0, %zu), S2: [%zu, %zu)\n", total1, n + n_end - total2, n + n_end);
    *c1_out = total1;
    *c2_out = n + n_end - total2;
    *r1_out = r1;
    *r2_out = r2;
}

/**
 * Routines to navigate an array P of size n, 
 * distributed block-cyclically with parameter b over nthreads threads.
 *
 * For thread t, we have
 *
 * P^t = {P[t * b + l * nthreads * b + j] | 
 *          0 <= t * b + l * nthreads * b + j < n, 
 *          0 <= j < block}
 *
 * FIXME: unfortunately, this is much slower than inlining the struct. 
 * Do not use for performance critical parts of the code.
 **/

#ifndef STRIDED_SUBARRAY_HSDFSDFDSF
#define STRIDED_SUBARRAY_HSDFSDFDSF

#include <stddef.h>
#include <string.h>

/**
 * BC_index i = k + j, where j < block 
 **/
typedef struct {
    size_t k;
    size_t j;
} BC_index;

/** 
 * Add/Sub assume count < block.
 * We do not need more than this, and handling that case comes at a performance
 * cost.
 **/
inline void 
BC_Add(BC_index *i, size_t count, size_t block, unsigned int nthreads)
{
    if (i->j + count < block) {
        i->j += count;
    } else {
        i->k += nthreads * block;
        i->j = i->j + count - block;
    }
}

inline void 
BC_Sub(BC_index *i, size_t count, size_t block, unsigned int nthreads)
{
    if (i->j >= count) {
        i->j -= count;
    } else {
        i->k -= nthreads * block;
        i->j = block + i->j - count;
    }
}

/**
 * Computes the number of elements between i1 and i2.
 * i1 > i2, and i1, i2 must belong to the same subarray.
 **/
inline size_t 
BC_Dist(BC_index i1, BC_index i2, unsigned int nthreads)
{
    return (i1.k - i2.k) /  nthreads + i1.j - i2.j;
}

/**
 * Copies count values from src, index i_s to dest, index i_d. 
 * src and dest do not need to belong to the same subarray!
 **/
inline void
BC_Copy(double *dest, BC_index i_d, double *src, BC_index i_s, size_t count,
        size_t block, unsigned int nthreads)
{
    while (count != 0) {
        /* Copy from the array having the least space in its block. */
        BC_index i = (i_d.j > i_s.j) ? i_d : i_s;
        size_t elems = (block - i.j < count) ? 
                             (block - i.j) :
                             count;  
        memcpy(dest + i_d.k + i_d.j, src + i_s.k + i_s.j,
               elems * sizeof(double));
        BC_Add(&i_d, elems, block, nthreads);
        BC_Add(&i_s, elems, block, nthreads);
        count -= elems;
    }
}

/**
 * Returns smallest i in P^t larger or equal to u
 **/
inline BC_index 
BC_Upper(unsigned int t, size_t u, size_t b, unsigned int nthreads)
{
    long ul = (long)u;
    /* We start with the largest index of the form
     *   t * b + l * b * nthreads < u,
     * or in other words
     *   t * b + l * b * nthreads <= u - 1.
     * We can obtain this by rounding down. */
    long kl = t * b + (ul - 1 - t * b) / (b * nthreads) * (b * nthreads);
    size_t j;
    /* Then we move up to u, or the start of the next block */
    if (u - kl < b) {
        j = u - kl;
    } else {
        kl += nthreads * b;
        j = 0;
    }
    size_t k = (long)kl;
    return {j, k};
}

#endif /* STRIDED_SUBARRAY_HSDFSDFDSF */

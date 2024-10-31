/******************************************************************************
 * Purpose of this file is to abstract iteration over block-cyclic subarrays.
 *****************************************************************************/

#include <stddef.h>

typedef struct {
    size_t block;   /* Block size */
    size_t t;       /* Processor index */
    unsigned int p; /* Number of processors */
    size_t i;       /* Index, i = k + j */
    size_t k;
    size_t j;       /* j < block */
} BlockCycIndex;

inline
void operator+=(BlockCycIndex& i, size_t count)
{
    assert(count < i.block);
    if (i.j + count < i.block) {
        i.j += count;
    } else {
        i.k += i.p * i.block;
        i.j = i.j + count - i.block;
    }
    i.i = i.k + i.j;
    assert(i.j < i.block);
}

inline
void operator-=(BlockCycIndex& i, size_t count)
{
    assert(count < i.block);
    if (i.j >= count) {
        i.j -= count;
    } else {
        i.k -= i.p * i.block;
        i.j += i.block - count;
    }
    i.i = i.k + i.j;
    assert(i.j < i.block);
}

/* Only defined for i1 > i2 in the same subarray. */
inline
size_t operator-(BlockCycIndex i1, BlockCycIndex i2)
{
    assert(i1.p == i2.p);
    assert(i1.t == i2.t);

    assert(i1.k >= i2.k);
    assert((i1.k - i2.k) / i1.p + i1.j >= i2.j);
    assert(i1.k + i1.j >= i2.k + i2.j);
    return (i1.k - i2.k) / i1.p + i1.j - i2.j;
}

static inline size_t 
ceildiv(size_t a, size_t b)
{
    return (a + b - 1) / b;
}

/**
 * Finds the 'supremum' of i in the subarray belonging to thread t.
 * That is, the smallest number greater or equal to i in t's subarray.
 **/
inline
BlockCycIndex BlockCycSup(unsigned int t, size_t block, unsigned int p,
                          size_t i)
{
    assert(i >= t * block);

    BlockCycIndex index;
    index.block = block;
    index.t     = t;
    index.p     = p;
    if ((i / block) % p == t) {
        /* i is in cyclic subarray, so sup is i */
        index.j = i % block;
        index.k = i - index.j;
    } else if (i < t * block) {
        index.k = t * block;
        index.j = 0;
    } else {
        /* sup is smallest t * block + l * p * block >= i. */
        size_t l = ceildiv(i - t * block, p * block);
        index.k = t * block + l * p * block;
        index.j = 0;
    }

    index.i = index.i + index.k;

    return index;
}

inline
BlockCycIndex BlockCycBegin(unsigned int t, unsigned int p, size_t block)
{
    BlockCycIndex index;

    index.block = block;
    index.t     = t;
    index.p     = p;
    index.i     = t * block;
    index.k     = t * block;
    index.j     = 0;

    return index;
}

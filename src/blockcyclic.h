/******************************************************************************
 * 'Iterator' for block-cyclic distribution, but very special-purpose for
 * maximum speed. See asserts for preconditions.
 *****************************************************************************/

#include <stddef.h>

template<int block>
struct BlockCycIndex {
    size_t t;       /* Processor index */
    unsigned int p; /* Number of processors */
    size_t k;       /* index is k + j */
    size_t j;       /* j < block */
};

template<int block>
inline size_t getIndex(BlockCycIndex<block> i)
{
    return i.k + i.j;
}

template<int block>
inline void operator+=(BlockCycIndex<block>& i, size_t count)
{
    assert(count < block);
    if (i.j + count < block) {
        i.j += count;
    } else {
        i.k += i.p * block;
        i.j = i.j + count - block;
    }
    assert(i.j < block);
}

template<int block>
inline void operator-=(BlockCycIndex<block>& i, size_t count)
{
    assert(count < block);
    if (i.j >= count) {
        i.j -= count;
    } else {
        i.k -= i.p * block;
        i.j += block - count;
    }
    assert(i.j < block);
}

template<int block>
inline size_t operator-(BlockCycIndex<block> i1, BlockCycIndex<block> i2)
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
template<int block>
inline BlockCycIndex<block> 
BlockCycSup(unsigned int t, unsigned int p, size_t n)
{
    assert(n >= t * block);
    assert(n % block == 0);

    BlockCycIndex<block> index;
    index.t     = t;
    index.p     = p;
    if ((n / block) % p == t) {
        /* i is in cyclic subarray, so sup is i */
        index.j = n % block;
        index.k = n - index.j;
    } else if (n < t * block) {
        index.k = t * block;
        index.j = 0;
    } else {
        /* sup is smallest t * block + l * p * block >= n. */
        size_t l = ceildiv(n - t * block, p * block);
        index.k = t * block + l * p * block;
        index.j = 0;
    }

    return index;
}

template<int block>
inline BlockCycIndex<block>
BlockCycBegin(unsigned int t, unsigned int p)
{
    BlockCycIndex<block> index;

    index.t     = t;
    index.p     = p;
    index.k     = t * block;
    index.j     = 0;

    return index;
}

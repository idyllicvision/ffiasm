#include <memory>
#include "msm.hpp"
#include "misc.hpp"

template <typename Curve, typename BaseField>
uint64_t MSM<Curve, BaseField>::estimateCost(uint8_t* _scalars,
                                             uint64_t _scalarSize,
                                             uint64_t _n,
                                             uint64_t _bitsPerChunk)
{
    const uint64_t nPoints = _n;
    if (nPoints == 0) return 0;

    ThreadPool &threadPool = ThreadPool::defaultPool();
    const uint64_t nThreads = threadPool.getThreadCount();

    scalars = _scalars;
    scalarSize = _scalarSize;

#ifdef MSM_BITS_PER_CHUNK
    bitsPerChunk = MSM_BITS_PER_CHUNK;
#else
    bitsPerChunk = _bitsPerChunk ? _bitsPerChunk : calcBitsPerChunk(nPoints, scalarSize);
#endif

    if (nPoints == 1) {
        return 1;
    }

    const uint64_t nChunks = calcChunkCount(scalarSize, bitsPerChunk);
    const uint64_t nBuckets = calcBucketCount(bitsPerChunk);

    std::unique_ptr<uint64_t[]> partial(new uint64_t[nThreads]);
    for (uint64_t t = 0; t < nThreads; t++) partial[t] = 0;

    threadPool.parallelFor(0, nPoints, [&] (int begin, int end, int numThread) {
        uint64_t localNonZero = 0;

        for (int i = begin; i < end; i++) {
            int carry = 0;

            for (uint64_t j = 0; j < nChunks; j++) {
                int bucketIndex = (int)getBucketIndex(i, j) + carry;

                if (bucketIndex >= (int)nBuckets) {
                    bucketIndex -= (int)(nBuckets * 2);
                    carry = 1;
                } else {
                    carry = 0;
                }

                if (bucketIndex != 0) {
                    localNonZero++;
                }
            }
        }

        partial[numThread] += localNonZero;
    });

    uint64_t nonZeroSlices = 0;
    for (uint64_t t = 0; t < nThreads; t++) {
        nonZeroSlices += partial[t];
    }

    const uint64_t bucketReduceAdds = nChunks * 2 * (nBuckets - 1);
    const uint64_t chunkMergeAdds = (nChunks > 0) ? (nChunks - 1) : 0;
    const uint64_t chunkMergeDbls = (nChunks > 0) ? ((nChunks - 1) * bitsPerChunk) : 0;

    return nonZeroSlices + bucketReduceAdds + chunkMergeAdds + chunkMergeDbls;
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::run(typename Curve::Point &r,
                                typename Curve::PointAffine *_bases,
                                uint8_t* _scalars,
                                uint64_t _scalarSize,
                                uint64_t _n,
                                uint64_t _nThreads)
{
    ThreadPool &threadPool = ThreadPool::defaultPool();

    const uint64_t nThreads = threadPool.getThreadCount();
    const uint64_t nPoints = _n;

    scalars = _scalars;
    scalarSize = _scalarSize;

#ifdef MSM_BITS_PER_CHUNK
    bitsPerChunk = MSM_BITS_PER_CHUNK;
#else
    bitsPerChunk = calcBitsPerChunk(nPoints, scalarSize);
#endif

    if (nPoints == 0) {
        g.copy(r, g.zero());
        return;
    }
    if (nPoints == 1) {
        g.mulByScalar(r, _bases[0], scalars, scalarSize);
        return;
    }

    const uint64_t nChunks = calcChunkCount(scalarSize, bitsPerChunk);
    const uint64_t nBuckets = calcBucketCount(bitsPerChunk);
    const uint64_t matrixSize = nThreads * nBuckets;
    const uint64_t nSlices = nChunks*nPoints;

    std::unique_ptr<typename Curve::Point[]> bucketMatrix(new typename Curve::Point[matrixSize]);
    std::unique_ptr<typename Curve::Point[]> chunks(new typename Curve::Point[nChunks]);
    std::unique_ptr<int32_t[]> slicedScalars(new int32_t[nSlices]);

    threadPool.parallelFor(0, nPoints, [&] (int begin, int end, int numThread) {

        for (int i = begin; i < end; i++) {
            int carry = 0;

            for (int j = 0; j < nChunks; j++) {
                int bucketIndex = getBucketIndex(i, j) + carry;

                if (bucketIndex >= nBuckets) {
                    bucketIndex -= nBuckets*2;
                    carry = 1;
                } else {
                    carry = 0;
                }

                slicedScalars[i*nChunks + j] = bucketIndex;
            }
        }
    });

    threadPool.parallelFor(0, nChunks, [&] (int begin, int end, int numThread) {

        for (int j = begin; j < end; j++) {

            typename Curve::Point *buckets = &bucketMatrix[numThread*nBuckets];

            for (int i = 0; i < nBuckets; i++) {
                g.copy(buckets[i], g.zero());
            }

            for (int i = 0; i < nPoints; i++) {
                const int bucketIndex = slicedScalars[i*nChunks + j];

                if (bucketIndex > 0) {
                    g.add(buckets[bucketIndex-1], buckets[bucketIndex-1], _bases[i]);

                } else if (bucketIndex < 0) {
                    g.sub(buckets[-bucketIndex-1], buckets[-bucketIndex-1], _bases[i]);
                }
            }

            typename Curve::Point t, tmp;

            g.copy(t, buckets[nBuckets - 1]);
            g.copy(tmp, t);

            for (int i = nBuckets - 2; i >= 0 ; i--) {
                g.add(tmp, tmp, buckets[i]);
                g.add(t, t, tmp);
            }

            chunks[j] = t;
        }
    });

    g.copy(r, chunks[nChunks - 1]);

    for (int j = nChunks - 2; j >= 0; j--) {
        for (int i = 0; i < bitsPerChunk; i++) {
            g.dbl(r, r);
        }
        g.add(r, r, chunks[j]);
    }
}

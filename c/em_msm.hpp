#pragma once

#include <cstdint>
#include <cstring>
#include <type_traits>
#include <memory>
#include <mutex>
#include <unordered_map>

#include "misc.hpp"

#include "msm.hpp"
#include "em.hpp"
#include "fq.hpp"

namespace em {

// inversion y
template <class PointAffine>
inline void invY(PointAffine &p)
{
    RawFq::Element tmp;
    RawFq::field.neg(tmp, p.y);
    p.y = tmp;
}

// msm::run w/ em
template <typename Curve, typename BaseField>
void run_msm_em(
    MSM<Curve, BaseField> &msm,
    Curve &g,
    typename Curve::Point &r,
    typename Curve::PointAffine *bases,
    uint8_t *scalars,
    uint64_t scalarSize,
    uint64_t n,
    uint64_t nThreads)
{
    // not G1 (BaseField != RawFq) -> common MSM w/o em
    if constexpr (!std::is_same<BaseField, RawFq>::value) {
        msm.run(r, bases, scalars, scalarSize, n, nThreads);
    } else {
        if (scalarSize != 32 || n == 0) {
            msm.run(r, bases, scalars, scalarSize, n, nThreads);
            return;
        }
        const uint64_t nPointsEM    = 2 * n;
        const uint64_t scalarSizeEM = 16; // k1,k2 <= 127 bits
        std::unique_ptr<typename Curve::PointAffine[]> bases2(
            new typename Curve::PointAffine[nPointsEM]
        );
        std::unique_ptr<uint8_t[]> scalars2(
            new uint8_t[nPointsEM * scalarSizeEM]
        );

        ThreadPool &tp = ThreadPool::defaultPool();
        tp.parallelFor(0, (int)n, [&] (int begin, int end, int /*tid*/) {
            for (int ii = begin; ii < end; ii++) {
                uint64_t i = (uint64_t)ii;
                const uint8_t *k = &scalars[i * 32];
                auto dec = em::decompose(k);

                // ---- k1 * Gi ----
                bases2[2*i] = bases[i];
                if (dec.neg1) {
                    invY(bases2[2*i]);
                }
                std::memcpy(&scalars2[(2*i) * scalarSizeEM], dec.k1, scalarSizeEM);

                // ---- k2 * phi(Gi) ----
                bases2[2*i + 1] = bases[i];
                em::phiP(bases2[2*i + 1]);
                if (dec.neg2) {
                    invY(bases2[2*i + 1]);
                }
                std::memcpy(&scalars2[(2*i + 1) * scalarSizeEM], dec.k2, scalarSizeEM);
            }
        });

        msm.run(
            r,
            bases2.get(),
            scalars2.get(),
            scalarSizeEM,
            nPointsEM,
            nThreads
        );
    }
}

} // namespace em

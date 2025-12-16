#pragma once

#include <cstdint>
#include <cstring>
#include <vector>
#include <type_traits>

#include "msm.hpp"
#include "glv_bn254_nogmp.hpp"
#include "fq.hpp"

namespace glv_bn254 {

// inversion y
template <class PointAffine>
inline void negate_g1_affine_inplace(PointAffine &p)
{
    RawFq::Element tmp;
    RawFq::field.neg(tmp, p.y);
    p.y = tmp;
}

// msm::run w/ glv
template <typename Curve, typename BaseField>
void run_msm_with_glv_if_g1_bn254(
    MSM<Curve, BaseField> &msm,
    Curve &g,
    typename Curve::Point &r,
    typename Curve::PointAffine *bases,
    uint8_t *scalars,
    uint64_t scalarSize,
    uint64_t n,
    uint64_t nThreads)
{
    // not G1 (BaseField != RawFq) -> common MSM w/o glv
    if constexpr (!std::is_same<BaseField, RawFq>::value) {
        msm.run(r, bases, scalars, scalarSize, n, nThreads);
    } else {
        if (scalarSize != 32 || n == 0) {
            msm.run(r, bases, scalars, scalarSize, n, nThreads);
            return;
        }
        const uint64_t nPointsGLV    = 2 * n;
        const uint64_t scalarSizeGLV = 16; // k1,k2 <= ~127 bits

        std::vector<typename Curve::PointAffine> bases2(nPointsGLV);
        std::vector<uint8_t> scalars2(nPointsGLV * scalarSizeGLV);

        for (uint64_t i = 0; i < n; i++) {
            const uint8_t *k_le = &scalars[i * 32];
            auto dec = glv_bn254::decompose_fr_le_32(k_le);

            // ---- k1 * Gi ----
            bases2[2*i] = bases[i];
            if (dec.neg1) {
                negate_g1_affine_inplace(bases2[2*i]);
            }
            std::memcpy(&scalars2[(2*i) * scalarSizeGLV], dec.k1, scalarSizeGLV);

            // ---- k2 * phi(Gi) ----
            bases2[2*i + 1] = bases[i];
            glv_bn254::apply_phi_inplace_g1(bases2[2*i + 1]);
            if (dec.neg2) {
                negate_g1_affine_inplace(bases2[2*i + 1]);
            }
            std::memcpy(&scalars2[(2*i + 1) * scalarSizeGLV], dec.k2, scalarSizeGLV);
        }

        msm.run(
            r,
            bases2.data(),
            scalars2.data(),
            scalarSizeGLV,
            nPointsGLV,
            nThreads
        );
    }
}

} // namespace glv_bn254

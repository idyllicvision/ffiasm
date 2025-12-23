#pragma once

#include <cstdint>
#include <cstring>
#include <gmp.h>
#include "fq.hpp"

namespace em {

static constexpr int N64 = 4;
static const mp_limb_t FR_N[N64] = {
    0x43e1f593f0000001ULL,
    0x2833e84879b97091ULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL
};

// -------- v1=(a1,b1), v2=(a2,b2) --------
static const mp_limb_t A1 = 0x89d3256894d213e3ULL; // 9931322734385697763
static const mp_limb_t B2 = 0x89d3256894d213e3ULL; // 9931322734385697763

static const mp_limb_t ABS_B1[2] = {
    0x8211bbeb7d4f1128ULL,
    0x6f4d8248eeb859fcULL
};

static const mp_limb_t A2[2] = {
    0x0be4e1541221250bULL,
    0x6f4d8248eeb859fdULL
};

static inline const RawFq::Element& beta_mont()
{
    static const RawFq::Element b = { {
        0x71930c11d782e155ULL,
        0xa6bb947cffbe3323ULL,
        0xaa303344d4741444ULL,
        0x2c3b3f0d26594943ULL
    } };
    return b;
}

static inline void bytes32_to_mpn4(mp_limb_t out[4], const uint8_t in[32])
{
    for (int i = 0; i < 4; i++) {
        uint64_t w = 0;
        for (int j = 0; j < 8; j++) w |= (uint64_t)in[i*8 + j] << (8*j);
        out[i] = (mp_limb_t)w;
    }
}

static inline void mpn_to_bytes(uint8_t* out, size_t outLen, const mp_limb_t* a, size_t an)
{
    std::memset(out, 0, outLen);
    for (size_t i = 0; i < an; i++) {
        uint64_t w = (uint64_t)a[i];
        for (int j = 0; j < 8; j++) {
            size_t idx = i*8 + (size_t)j;
            if (idx < outLen) out[idx] = (uint8_t)((w >> (8*j)) & 0xff);
        }
    }
}

static inline void round_div_pos(mp_limb_t* q, mp_limb_t* r, const mp_limb_t* num, size_t nn)
{
    mpn_tdiv_qr(q, r, 0, num, nn, FR_N, N64);
    mp_limb_t r2[N64];
    mp_limb_t carry = mpn_lshift(r2, r, N64, 1);
    int ge = carry ? 1 : (mpn_cmp(r2, FR_N, N64) >= 0);
    if (ge) {
        mpn_add_1(q, q, nn - N64 + 1, 1);
    }
}

struct k1k2 {
    bool neg1;
    bool neg2;
    uint8_t k1[16];
    uint8_t k2[16];
};

static inline k1k2 decompose(const uint8_t a[32])
{
    mp_limb_t k[4]; bytes32_to_mpn4(k, a);

    // q1 = round( k * b2 / N ) ; b2 = 64-bit
    mp_limb_t num1[5];
    num1[4] = mpn_mul_1(num1, k, 4, (mp_limb_t)B2);

    mp_limb_t q1[2] = {0,0};
    mp_limb_t r1[4] = {0,0,0,0};
    round_div_pos(q1, r1, num1, 5);

    // q2 = round( k * abs(b1) / N ) ; abs(b1)=2 limbs
    mp_limb_t num2[6];
    mpn_mul(num2, k, 4, ABS_B1, 2); // 6 limbs

    mp_limb_t q2[3] = {0,0,0};
    mp_limb_t r2[4] = {0,0,0,0};
    // quotient size = 6-4+1 = 3
    mpn_tdiv_qr(q2, r2, 0, num2, 6, FR_N, N64);
    // rounding:
    mp_limb_t r2x2[4];
    mp_limb_t c = mpn_lshift(r2x2, r2, 4, 1);
    int ge2 = c ? 1 : (mpn_cmp(r2x2, FR_N, 4) >= 0);
    if (ge2) mpn_add_1(q2, q2, 3, 1);

    // q2*a2 => 5 limbs (3*2)
    mp_limb_t q2a2[5] = {0,0,0,0,0};
    mpn_mul(q2a2, q2, 3, A2, 2);

    // q1*a1: (2 limbs * 1 limb) => 3 limbs
    mp_limb_t q1a1[3] = {0,0,0};
    q1a1[2] = mpn_mul_1(q1a1, q1, 2, (mp_limb_t)A1);

    // vx = q2a2 + q1a1
    mp_limb_t vx[5];
    std::memcpy(vx, q2a2, sizeof(vx));
    mpn_add(vx, vx, 5, q1a1, 3);

    // k1 = k - vx (signed)
    mp_limb_t k_ext[5] = {k[0],k[1],k[2],k[3],0};
    bool neg1 = (mpn_cmp(k_ext, vx, 5) < 0);

    mp_limb_t k1mag[5] = {0,0,0,0,0};
    if (!neg1) {
        mpn_sub_n(k1mag, k_ext, vx, 5);
    } else {
        mpn_sub_n(k1mag, vx, k_ext, 5);
    }

    // Compute k2 = q1*abs(b1) - q2*b2
    // q1*abs(b1): (2*2)->4 limbs
    mp_limb_t q1b1[4] = {0,0,0,0};
    mpn_mul(q1b1, q1, 2, ABS_B1, 2);

    // q2*b2: (3*1)->4 limbs via mul_1
    mp_limb_t q2b2[4] = {0,0,0,0};
    q2b2[3] = mpn_mul_1(q2b2, q2, 3, (mp_limb_t)B2);

    bool neg2 = (mpn_cmp(q1b1, q2b2, 4) < 0);
    mp_limb_t k2mag[4] = {0,0,0,0};
    if (!neg2) {
        mpn_sub_n(k2mag, q1b1, q2b2, 4);
    } else {
        mpn_sub_n(k2mag, q2b2, q1b1, 4);
    }

    k1k2 out;
    out.neg1 = neg1;
    out.neg2 = neg2;

    mpn_to_bytes(out.k1, 16, k1mag, 5);
    mpn_to_bytes(out.k2, 16, k2mag, 4);

    return out;
}

template <class PointAffine>
static inline void phiP(PointAffine &p)
{
    RawFq::Element tmp;
    RawFq::field.mul(tmp, p.x, beta_mont());
    p.x = tmp;
}

} // namespace em
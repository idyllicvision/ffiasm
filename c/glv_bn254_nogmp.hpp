#pragma once

#include <cstdint>
#include <cstring>
#include "fq.hpp"

namespace glv_bn254 {

// ---------------- BN254 Fr order N (u64 limbs) ----------------
static constexpr int N64 = 4;
static constexpr uint64_t FR_N[N64] = {
    0x43e1f593f0000001ULL,
    0x2833e84879b97091ULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL
};

// -------- v1=(a1,b1), v2=(a2,b2) constants --------
static constexpr uint64_t A1 = 0x89d3256894d213e3ULL; // 9931322734385697763
static constexpr uint64_t B2 = 0x89d3256894d213e3ULL; // 9931322734385697763

static constexpr uint64_t ABS_B1[2] = {
    0x8211bbeb7d4f1128ULL,
    0x6f4d8248eeb859fcULL
};

static constexpr uint64_t A2[2] = {
    0x0be4e1541221250bULL,
    0x6f4d8248eeb859fdULL
};

// beta in Montgomery for Fq (beta*R mod p) 4 limbs
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

// ---------------- small limb helpers (base B = 2^64) ----------------
static inline int cmp_u64(const uint64_t* a, const uint64_t* b, size_t n)
{
    for (size_t i = n; i-- > 0; ) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return  1;
    }
    return 0;
}

static inline uint64_t add_u64(uint64_t* r, const uint64_t* a, const uint64_t* b, size_t n)
{
    unsigned __int128 carry = 0;
    for (size_t i = 0; i < n; i++) {
        unsigned __int128 s = (unsigned __int128)a[i] + b[i] + carry;
        r[i] = (uint64_t)s;
        carry = s >> 64;
    }
    return (uint64_t)carry;
}

static inline uint64_t add_u64_inplace(uint64_t* r, const uint64_t* a, size_t n)
{
    unsigned __int128 carry = 0;
    for (size_t i = 0; i < n; i++) {
        unsigned __int128 s = (unsigned __int128)r[i] + a[i] + carry;
        r[i] = (uint64_t)s;
        carry = s >> 64;
    }
    return (uint64_t)carry;
}

static inline uint64_t add1_u64(uint64_t* r, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        uint64_t old = r[i];
        r[i] = old + 1;
        if (r[i] != 0) return 0;
    }
    return 1;
}

static inline uint64_t sub_u64(uint64_t* r, const uint64_t* a, const uint64_t* b, size_t n)
{
    uint64_t borrow = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t bi = b[i] + borrow;
        uint64_t new_borrow = (bi < b[i]) ? 1 : 0; // overflow in b[i]+borrow
        uint64_t ai = a[i];
        uint64_t out = ai - bi;
        borrow = (ai < bi) ? 1 : 0;
        borrow |= new_borrow;
        r[i] = out;
    }
    return borrow;
}

static inline uint64_t sub_u64_inplace(uint64_t* r, const uint64_t* b, size_t n)
{
    uint64_t borrow = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t bi = b[i] + borrow;
        uint64_t new_borrow = (bi < b[i]) ? 1 : 0;
        uint64_t ri = r[i];
        uint64_t out = ri - bi;
        borrow = (ri < bi) ? 1 : 0;
        borrow |= new_borrow;
        r[i] = out;
    }
    return borrow;
}

static inline uint64_t lshift_u64(uint64_t* out, const uint64_t* in, size_t n, unsigned s)
{
    if (s == 0) {
        std::memcpy(out, in, n * sizeof(uint64_t));
        return 0;
    }
    uint64_t carry = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t w = in[i];
        out[i] = (w << s) | carry;
        carry = (w >> (64 - s));
    }
    return carry;
}

static inline void rshift_u64(uint64_t* out, const uint64_t* in, size_t n, unsigned s)
{
    if (s == 0) {
        std::memcpy(out, in, n * sizeof(uint64_t));
        return;
    }
    uint64_t carry = 0;
    for (size_t i = n; i-- > 0; ) {
        uint64_t w = in[i];
        out[i] = (w >> s) | carry;
        carry = (w << (64 - s));
    }
}

static inline unsigned clz64(uint64_t x)
{
    return x ? (unsigned)__builtin_clzll(x) : 64u;
}

// num = a(na) * b(nb) -> out(na+nb), fixed-size small
static inline void mul_u64(uint64_t* out, const uint64_t* a, size_t na, const uint64_t* b, size_t nb)
{
    std::memset(out, 0, (na + nb) * sizeof(uint64_t));
    for (size_t i = 0; i < na; i++) {
        unsigned __int128 carry = 0;
        for (size_t j = 0; j < nb; j++) {
            unsigned __int128 cur = (unsigned __int128)a[i] * b[j]
                                  + out[i + j]
                                  + carry;
            out[i + j] = (uint64_t)cur;
            carry = cur >> 64;
        }
        out[i + nb] = (uint64_t)((unsigned __int128)out[i + nb] + carry);
    }
}

static inline void mul_u64_1(uint64_t* out, const uint64_t* a, size_t na, uint64_t b)
{
    unsigned __int128 carry = 0;
    for (size_t i = 0; i < na; i++) {
        unsigned __int128 cur = (unsigned __int128)a[i] * b + carry;
        out[i] = (uint64_t)cur;
        carry = cur >> 64;
    }
    out[na] = (uint64_t)carry;
}

// -------- Knuth division for divisor length n=4, numerator length n+m (m<=2) --------
// Returns q (m+1 limbs) and r (4 limbs)
static inline void divrem_4limb(uint64_t* q, size_t qn,
                               uint64_t r[4],
                               const uint64_t* num, size_t nn,
                               const uint64_t den[4])
{
    // nn must be 5 or 6, qn must be nn-4+1 (2 or 3)
    const size_t n = 4;
    const size_t m = nn - n; // 1 or 2

    // Normalize
    unsigned s = clz64(den[n - 1]);
    uint64_t v[4];
    lshift_u64(v, den, 4, s);

    // u needs nn+1 limbs for Knuth D
    uint64_t u[7] = {0,0,0,0,0,0,0};
    uint64_t carry = lshift_u64(u, num, nn, s);
    u[nn] = carry;

    // Main loop j = m .. 0
    for (size_t jj = 0; jj <= m; jj++) {
        size_t j = m - jj;

        // estimate qhat from top 2 limbs
        unsigned __int128 top = ((unsigned __int128)u[j + n] << 64) | u[j + n - 1];
        uint64_t qhat = (uint64_t)(top / v[n - 1]);
        uint64_t rhat = (uint64_t)(top % v[n - 1]);

        // clamp qhat to < B (it already is u64)
        // adjust while qhat*v[n-2] > B*rhat + u[j+n-2]
        while (true) {
            unsigned __int128 left  = (unsigned __int128)qhat * v[n - 2];
            unsigned __int128 right = ((unsigned __int128)rhat << 64) | u[j + n - 2];
            if (left <= right) break;
            qhat--;
            rhat = (uint64_t)(rhat + v[n - 1]);
            if (rhat < v[n - 1]) break; // overflow => stop
        }

        // subtract qhat*v from u[j..j+n]
        uint64_t borrow = 0;
        unsigned __int128 c2 = 0;
        for (size_t i = 0; i < n; i++) {
            unsigned __int128 p = (unsigned __int128)qhat * v[i] + c2;
            uint64_t p_lo = (uint64_t)p;
            c2 = (p >> 64);

            uint64_t ui = u[j + i];
            uint64_t sub = p_lo + borrow;
            uint64_t new_borrow = (sub < p_lo) ? 1 : 0;      // overflow in p_lo+borrow
            uint64_t out = ui - sub;
            borrow = (ui < sub) ? 1 : 0;
            borrow |= new_borrow;
            u[j + i] = out;
        }

        uint64_t ui_top = u[j + n];
        uint64_t sub_top = (uint64_t)c2 + borrow;
        uint64_t new_ui_top = ui_top - sub_top;
        uint64_t under = (ui_top < sub_top) ? 1 : 0;
        u[j + n] = new_ui_top;

        // if underflow: qhat-- and add back v
        if (under) {
            qhat--;

            unsigned __int128 carry3 = 0;
            for (size_t i = 0; i < n; i++) {
                unsigned __int128 ssum = (unsigned __int128)u[j + i] + v[i] + carry3;
                u[j + i] = (uint64_t)ssum;
                carry3 = ssum >> 64;
            }
            u[j + n] = (uint64_t)((unsigned __int128)u[j + n] + carry3);
        }

        q[j] = qhat;
    }

    // remainder = u[0..n-1] >> s
    uint64_t rem_norm[4] = { u[0], u[1], u[2], u[3] };
    rshift_u64(r, rem_norm, 4, s);

    // ensure q has exactly qn limbs; we wrote q[0..m] already
    (void)qn;
}

// round q = round(num / N) with remainder check: if 2r >= N => q++
static inline void round_div_N(uint64_t* q, size_t qn,
                               const uint64_t* num, size_t nn)
{
    uint64_t r[4] = {0,0,0,0};
    std::memset(q, 0, qn * sizeof(uint64_t));

    divrem_4limb(q, qn, r, num, nn, FR_N);

    // compare 2*r with N
    uint64_t r2[4];
    unsigned __int128 carry = 0;
    for (int i = 0; i < 4; i++) {
        unsigned __int128 s = (unsigned __int128)r[i] * 2u + carry;
        r2[i] = (uint64_t)s;
        carry = s >> 64;
    }
    bool ge = (carry != 0) || (cmp_u64(r2, FR_N, 4) >= 0);
    if (ge) {
        add1_u64(q, qn);
    }
}

static inline void bytes32_to_u64x4(uint64_t out[4], const uint8_t in[32])
{
    for (int i = 0; i < 4; i++) {
        uint64_t w = 0;
        for (int j = 0; j < 8; j++) w |= (uint64_t)in[i*8 + j] << (8*j);
        out[i] = w;
    }
}

static inline void u64_to_bytes(uint8_t* out, size_t outLen, const uint64_t* a, size_t an)
{
    std::memset(out, 0, outLen);
    for (size_t i = 0; i < an; i++) {
        uint64_t w = a[i];
        for (int j = 0; j < 8; j++) {
            size_t idx = i*8 + (size_t)j;
            if (idx < outLen) out[idx] = (uint8_t)((w >> (8*j)) & 0xff);
        }
    }
}

// ---------------- GLV decompose ----------------
struct GlvDecomp {
    bool neg1;
    bool neg2;
    uint8_t k1[16]; // abs(k1)
    uint8_t k2[16]; // abs(k2)
};

static inline GlvDecomp decompose_fr_le_32(const uint8_t a_le[32])
{
    uint64_t k[4]; bytes32_to_u64x4(k, a_le);

    // num1 = k * B2  => 5 limbs
    uint64_t num1[5];
    mul_u64_1(num1, k, 4, (uint64_t)B2);

    // q1 = round(num1 / N) => 2 limbs
    uint64_t q1[2] = {0,0};
    round_div_N(q1, 2, num1, 5);

    // num2 = k * abs(b1) => 6 limbs
    uint64_t num2[6];
    mul_u64(num2, k, 4, ABS_B1, 2);

    // q2 = round(num2 / N) => 3 limbs
    uint64_t q2[3] = {0,0,0};
    round_div_N(q2, 3, num2, 6);

    // q2*a2 => 5 limbs (3*2)
    uint64_t q2a2[5];
    mul_u64(q2a2, q2, 3, A2, 2);

    // q1*a1 => 3 limbs (2*1)
    uint64_t q1a1[3];
    mul_u64_1(q1a1, q1, 2, (uint64_t)A1);

    // vx = q2a2 + q1a1 (as 5 limbs)
    uint64_t vx[5] = { q2a2[0], q2a2[1], q2a2[2], q2a2[3], q2a2[4] };
    {
        uint64_t tmp3[5] = { q1a1[0], q1a1[1], q1a1[2], 0, 0 };
        add_u64_inplace(vx, tmp3, 5);
    }

    // k_ext = k (5 limbs)
    uint64_t k_ext[5] = { k[0], k[1], k[2], k[3], 0 };

    // k1 = k - vx (signed) => store abs + sign
    bool neg1 = (cmp_u64(k_ext, vx, 5) < 0);
    uint64_t k1mag[5];
    if (!neg1) sub_u64(k1mag, k_ext, vx, 5);
    else       sub_u64(k1mag, vx, k_ext, 5);

    // k2 = q1*abs(b1) - q2*b2 (signed)
    uint64_t q1b1[4];
    mul_u64(q1b1, q1, 2, ABS_B1, 2); // 4 limbs

    uint64_t q2b2[4];
    mul_u64_1(q2b2, q2, 3, (uint64_t)B2); // writes 4 limbs

    bool neg2 = (cmp_u64(q1b1, q2b2, 4) < 0);
    uint64_t k2mag[4];
    if (!neg2) sub_u64(k2mag, q1b1, q2b2, 4);
    else       sub_u64(k2mag, q2b2, q1b1, 4);

    GlvDecomp out{};
    out.neg1 = neg1;
    out.neg2 = neg2;

    // we only need 16 bytes (abs values)
    u64_to_bytes(out.k1, 16, k1mag, 5);
    u64_to_bytes(out.k2, 16, k2mag, 4);

    return out;
}

// Apply phi(P) = (beta*x, y) for affine point over RawFq
template <class PointAffine>
static inline void apply_phi_inplace_g1(PointAffine &p)
{
    RawFq::Element tmp;
    RawFq::field.mul(tmp, p.x, beta_mont());
    p.x = tmp;
}

} // namespace glv_bn254
#include <iostream>
#include <cstdint>
#include <cstring>

#include "gtest/gtest.h"
#include "alt_bn128.hpp"
#include "fft.hpp"
#include <chrono>
#include <iomanip>
#include <functional>
using namespace AltBn128;

namespace {

TEST(altBn128, f2_simpleMul) {
    F2Element e1;
    F2.fromString(e1, "(2,2)");

    F2Element e2;
    F2.fromString(e2, "(3,3)");

    F2Element e3;
    F2.mul(e3, e1, e2);

    F2Element e33;
    F2.fromString(e33, "(0,12)");

    ASSERT_TRUE(F2.eq(e3, e33));
}

TEST(altBn128, f6_mulDiv) {
    F6Element a;
    F6.fromString(a, "((1,2),(3,4),(5,6))");

    F6Element b;
    F6.fromString(b, "((12,11),(10,9),(8,7))");

    F6Element c,d;

    F6.mul(c,a,b);
    F6.div(d,c,b);

    ASSERT_TRUE(F6.eq(a,d));
}

TEST(altBn128, f6_inv) {
    F6Element a;
    F6.fromString(a, "((239487238491,2356249827341),"
                     "(082659782,182703523765),"
                     "(978236549263,64893242))");

    F6Element inv;
    F6.inv(inv, a);

    F6Element b;
    F6.mul(b, inv, a);

    ASSERT_TRUE(F6.isOne(b));
}

TEST(altBn128, f12_mulDiv) {
    F12Element a;
    F12.fromString(a, "(((1,2),(3,4),(5,6)),((7,8),(9,10),(11,12)))");

    F12Element b;
    F12.fromString(b, "(((12,11),(10,9),(8,7)),((6,5),(4,3),(2,1)))");

    F12Element c,d;

    F12.mul(c,a,b);
    F12.div(d,c,b);

    ASSERT_TRUE(F12.eq(a,d));
}

TEST(altBn128, f12_inv) {
    F12Element a;
    F12.fromString(a,
                   "(((239846234862342323958623,2359862352529835623),"
                     "(928836523,9856234),"
                     "(235635286,5628392833)),"
                    "((252936598265329856238956532167968,23596239865236954178968),"
                     "(95421692834,236548),"
                     "(924523,12954623)))" );

    F12Element inv;
    F12.inv(inv, a);

    F12Element b;
    F12.mul(b, inv, a);

    ASSERT_TRUE(F12.isOne(b));
}

TEST(altBn128, g1_PlusZero) {
    G1Point p1;
    G1.add(p1, G1.one(), G1.zero());
    ASSERT_TRUE(G1.eq(p1, G1.one()));
}

TEST(altBn128, g1_minus_g1) {
    G1Point p1;
    G1.sub(p1, G1.one(), G1.one());
    ASSERT_TRUE(G1.isZero(p1));
}

TEST(altBn128, g1_times_4) {
    G1Point p1;
    G1.add(p1, G1.one(), G1.one());
    G1.add(p1, p1, G1.one());
    G1.add(p1, p1, G1.one());

    G1Point p2;
    G1.dbl(p2, G1.one());
    G1.dbl(p2, p2);

    ASSERT_TRUE(G1.eq(p1,p2));
}

TEST(altBn128, g1_times_3) {
    G1Point p1;
    G1.add(p1, G1.one(), G1.one());
    G1.add(p1, p1, G1.one());

    G1Point p2;
    G1.dbl(p2, G1.one());
    G1.dbl(p2, p2);
    G1.sub(p2, p2, G1.one());

    ASSERT_TRUE(G1.eq(p1,p2));
}

TEST(altBn128, g1_times_3_exp) {
    G1Point p1;
    G1.add(p1, G1.one(), G1.one());
    G1.add(p1, p1, G1.one());

    mp_uint_t scalar;
    mp_uint_t x;
    mp_set(x, 3);
    mp_copy(scalar, x);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), (uint8_t*)scalar, MP_N);

    ASSERT_TRUE(G1.eq(p1,p2));
}

TEST(altBn128, g1_times_5) {
    G1Point p1;
    G1.dbl(p1, G1.one());
    G1.dbl(p1, p1);
    G1.add(p1, p1, p1);

    G1Point p2;
    G1Point p3;
    G1Point p4;
    G1Point p5;
    G1Point p6;
    G1.dbl(p2, G1.one());
    G1.dbl(p3, p2);
    G1.dbl(p4, G1.one());
    G1.dbl(p5, p4);
    G1.add(p6, p3, p5);

    ASSERT_TRUE(G1.eq(p1,p6));
}

TEST(altBn128, g1_times_65_exp) {
    G1Point p1;
    G1.dbl(p1, G1.one());
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.add(p1, p1, G1.one());

    mp_uint_t scalar;
    mp_uint_t x;
    mp_set(x, 65);
    mp_copy(scalar, x);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), (uint8_t*)scalar, MP_N);

    ASSERT_TRUE(G1.eq(p1,p2));
}

TEST(altBn128, g1_expToOrder) {
    mp_uint_t scalar;
    mp_uint_t x;
    ASSERT_TRUE(mp_set(x,
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10
    ));
    mp_copy(scalar, x);

    G1Point p1;
    G1.mulByScalar(p1, G1.one(), (uint8_t *)scalar, MP_N);

    ASSERT_TRUE(G1.isZero(p1));
}

TEST(altBn128, g2_expToOrder) {
    mp_uint_t scalar;
    mp_uint_t x;
    ASSERT_TRUE(mp_set(x,
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10
    ));
    mp_copy(scalar, x);

    Curve<F2Field<RawFq>>::Point p1;
    G2.mulByScalar(p1, G2.one(), (uint8_t *)scalar, MP_N);

    ASSERT_TRUE(G2.isZero(p1));
}

TEST(altBn128, multiExp_old) {
    int NMExp = 40000;

    typedef mp_uint_t Scalar;

    Scalar *scalars = new Scalar[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    uint64_t acc = 0;
    for (int i = 0; i < NMExp; i++) {
        if (i == 0) {
            G1.copy(bases[0], G1.one());
        } else {
            G1.add(bases[i], bases[i-1], G1.one());
        }

        mp_uint_t x;
        mp_set(x, (uint64_t)(i + 1));
        mp_copy(scalars[i], x);
        acc += (uint64_t)(i + 1) * (uint64_t)(i + 1);
    }

    G1Point p1;
    G1.multiMulByScalar(p1, bases, (uint8_t *)scalars, MP_N, NMExp);

    mp_uint_t sAcc;
    mp_uint_t x;
    mp_set(x, acc);
    mp_copy(sAcc, x);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), (uint8_t *)sAcc, MP_N);

    ASSERT_TRUE(G1.eq(p1, p2));
    //G1.printCounters();

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, multiExp) {
    int NMExp = 40000;

    typedef mp_uint_t Scalar;

    Scalar *scalars = new Scalar[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];
    G1Point *basesJac = new G1Point[NMExp];

    uint64_t acc = 0;

    G1.copy(basesJac[0], G1.one());
    for (int i = 1; i < NMExp; i++) {
        G1.add(basesJac[i], basesJac[i - 1], G1.one());
    }

    G1.batchToAffine(bases, basesJac, NMExp);

    for (int i = 0; i < NMExp; i++) {
        mp_uint_t x;
        mp_set(x, (uint64_t)(i + 1));
        mp_copy(scalars[i], x);
        acc += (uint64_t)(i + 1) * (uint64_t)(i + 1);
    }

    G1Point p1;
    G1.multiMulByScalar(p1, bases, (uint8_t *)scalars, MP_N, NMExp);

    mp_uint_t sAcc;
    mp_uint_t x;
    mp_set(x, acc);
    mp_copy(sAcc, x);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), (uint8_t *)sAcc, MP_N);

    ASSERT_TRUE(G1.eq(p1, p2));
    //G1.printCounters();

    delete[] basesJac;
    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, multiExpMSM_old) {
    int NMExp = 40000;

    typedef mp_uint_t Scalar;

    Scalar *scalars = new Scalar[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    uint64_t acc = 0;
    for (int i = 0; i < NMExp; i++) {
        if (i == 0) {
            G1.copy(bases[0], G1.one());
        } else {
            G1.add(bases[i], bases[i-1], G1.one());
        }

        mp_uint_t x;
        mp_set(x, (uint64_t)(i + 1));
        mp_copy(scalars[i], x);
        acc += (uint64_t)(i + 1) * (uint64_t)(i + 1);
    }

    G1Point p1;
    G1.multiMulByScalarMSM(p1, bases, (uint8_t *)scalars, MP_N, NMExp);

    mp_uint_t sAcc;
    mp_uint_t x;
    mp_set(x, acc);
    mp_copy(sAcc, x);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), (uint8_t *)sAcc, MP_N);

    ASSERT_TRUE(G1.eq(p1, p2));
    //G1.printCounters();

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, multiExpMSM) {
    int NMExp = 40000;

    typedef mp_uint_t Scalar;

    Scalar *scalars = new Scalar[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];
    G1Point *basesJac = new G1Point[NMExp];

    uint64_t acc = 0;

    G1.copy(basesJac[0], G1.one());
    for (int i = 1; i < NMExp; i++) {
        G1.add(basesJac[i], basesJac[i - 1], G1.one());
    }

    G1.batchToAffine(bases, basesJac, NMExp);

    for (int i = 0; i < NMExp; i++) {
        mp_uint_t x;
        mp_set(x, (uint64_t)(i + 1));
        mp_copy(scalars[i], x);
        acc += (uint64_t)(i + 1) * (uint64_t)(i + 1);
    }

    G1Point p1;
    G1.multiMulByScalarMSM(p1, bases, (uint8_t *)scalars, MP_N, NMExp);

    mp_uint_t sAcc;
    mp_uint_t x;
    mp_set(x, acc);
    mp_copy(sAcc, x);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), (uint8_t *)sAcc, MP_N);

    ASSERT_TRUE(G1.eq(p1, p2));

    delete[] basesJac;
    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, multiExp2) {
    int NMExp = 2;

    AltBn128::FrElement *scalars = new AltBn128::FrElement[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    F1.fromString(bases[0].x, "1626275109576878988287730541908027724405348106427831594181487487855202143055");
    F1.fromString(bases[0].y, "18706364085805828895917702468512381358405767972162700276238017959231481018884");
    F1.fromString(bases[1].x, "17245156998235704504461341147511350131061011207199931581281143511105381019978");
    F1.fromString(bases[1].y, "3858908536032228066651712470282632925312300188207189106507111128103204506804");

    Fr.fromString(scalars[0], "1");
    Fr.fromString(scalars[1], "20187316456970436521602619671088988952475789765726813868033071292105413408473");
    Fr.fromMontgomery(scalars[0], scalars[0]);
    Fr.fromMontgomery(scalars[1], scalars[1]);

    G1Point r;
    G1PointAffine ra;
    G1PointAffine ref;

    F1.fromString(ref.x, "9163953212624378696742080269971059027061360176019470242548968584908855004282");
    F1.fromString(ref.y, "20922060990592511838374895951081914567856345629513259026540392951012456141360");

    G1.multiMulByScalar(r, bases, (uint8_t *)scalars, MP_N, 2);
    G1.copy(ra, r);

    ASSERT_TRUE(G1.eq(ra, ref));

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, multiExp2MSM) {
    int NMExp = 2;

    AltBn128::FrElement *scalars = new AltBn128::FrElement[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    F1.fromString(bases[0].x, "1626275109576878988287730541908027724405348106427831594181487487855202143055");
    F1.fromString(bases[0].y, "18706364085805828895917702468512381358405767972162700276238017959231481018884");
    F1.fromString(bases[1].x, "17245156998235704504461341147511350131061011207199931581281143511105381019978");
    F1.fromString(bases[1].y, "3858908536032228066651712470282632925312300188207189106507111128103204506804");

    Fr.fromString(scalars[0], "1");
    Fr.fromString(scalars[1], "20187316456970436521602619671088988952475789765726813868033071292105413408473");
    Fr.fromMontgomery(scalars[0], scalars[0]);
    Fr.fromMontgomery(scalars[1], scalars[1]);

    G1Point r;
    G1PointAffine ra;
    G1PointAffine ref;

    F1.fromString(ref.x, "9163953212624378696742080269971059027061360176019470242548968584908855004282");
    F1.fromString(ref.y, "20922060990592511838374895951081914567856345629513259026540392951012456141360");

    G1.multiMulByScalarMSM(r, bases, (uint8_t *)scalars, MP_N, 2);
    G1.copy(ra, r);

    ASSERT_TRUE(G1.eq(ra, ref));

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, multiExp8MSM) {
    int NMExp = 8;

    AltBn128::FrElement *scalars = new AltBn128::FrElement[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    F1.fromString(bases[0].x, "1");
    F1.fromString(bases[0].y, "2");
    F1.fromString(bases[1].x, "1368015179489954701390400359078579693043519447331113978918064868415326638035");
    F1.fromString(bases[1].y, "9918110051302171585080402603319702774565515993150576347155970296011118125764");
    F1.fromString(bases[2].x, "3353031288059533942658390886683067124040920775575537747144343083137631628272");
    F1.fromString(bases[2].y, "19321533766552368860946552437480515441416830039777911637913418824951667761761");
    F1.fromString(bases[3].x, "10744596414106452074759370245733544594153395043370666422502510773307029471145");
    F1.fromString(bases[3].y, "848677436511517736191562425154572367705380862894644942948681172815252343932");
    F1.fromString(bases[4].x, "10415861484417082502655338383609494480414113902179649885744799961447382638712");
    F1.fromString(bases[4].y, "10196215078179488638353184030336251401353352596818396260819493263908881608606");
    F1.fromString(bases[5].x, "4444740815889402603535294170722302758225367627362056425101568584910268024244");
    F1.fromString(bases[5].y, "10537263096529483164618820017164668921386457028564663708352735080900270541420");
    F1.fromString(bases[6].x, "12852522211178622728088728121177131998585782282560100422041774753646305409836");
    F1.fromString(bases[6].y, "15918672909255108529698304535345707578139606904951176064731093256171019744261");
    F1.fromString(bases[7].x, "5841054468380737358126901601208759440531393618333695939860021859990434602332");
    F1.fromString(bases[7].y, "14496198492936799798866613472708085382638873795843716312315177042738122955744");

    Fr.fromString(scalars[0], "1");
    Fr.fromString(scalars[1], "4");
    Fr.fromString(scalars[2], "4");
    Fr.fromString(scalars[3], "4");
    Fr.fromString(scalars[4], "8");
    Fr.fromString(scalars[5], "10");
    Fr.fromString(scalars[6], "1");
    Fr.fromString(scalars[7], "5");
    Fr.fromMontgomery(scalars[0], scalars[0]);
    Fr.fromMontgomery(scalars[1], scalars[1]);
    Fr.fromMontgomery(scalars[2], scalars[2]);
    Fr.fromMontgomery(scalars[3], scalars[3]);
    Fr.fromMontgomery(scalars[4], scalars[4]);
    Fr.fromMontgomery(scalars[5], scalars[5]);
    Fr.fromMontgomery(scalars[6], scalars[6]);
    Fr.fromMontgomery(scalars[7], scalars[7]);

    G1Point r;
    G1PointAffine ra;
    G1PointAffine ref;

    F1.fromString(ref.x, "17747920359253913546551417160303297937542312574889904290131615776238588901697");
    F1.fromString(ref.y, "8815119438581789680513912776342567599606944899217792926373871775002956510503");

    G1.multiMulByScalarMSM(r, bases, (uint8_t *)scalars, MP_N, NMExp);
    G1.copy(ra, r);

    ASSERT_TRUE(G1.eq(ra, ref));

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, fft) {
    int NMExp = 1<<10;

    AltBn128::FrElement *a = new AltBn128::FrElement[NMExp];

    for (int i=0; i<NMExp; i++) {
        Fr.fromUI(a[i], i+1);
    }

    FFT<typename Engine::Fr> fft(NMExp);

    fft.fft(a, NMExp);
    fft.ifft(a, NMExp);

    AltBn128::FrElement aux;
    for (int i=0; i<NMExp; i++) {
        Fr.fromUI(aux, i+1);
        ASSERT_TRUE(Fr.eq(a[i], aux));
    }

    delete[] a;
}

using perf_clock = std::chrono::steady_clock;

static inline double bench_ms(std::function<void()> fn) {
    auto t0 = perf_clock::now();
    fn();
    auto t1 = perf_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

static inline double ns_per_op(double ms, size_t n) {
    return (ms * 1e6) / (double)n;
}

TEST(altBn128Perf, g1_curve_ops_timing_nogmp) {
    G1PointAffine a_aff, b_aff, c_aff;
    G1.copy(a_aff, G1.one());           //  P
    G1.add(b_aff, a_aff, G1.one());    // 2P
    G1.add(c_aff, b_aff, G1.one());    // 3P

    G1Point a_jac, b_jac, c_jac;
    G1.copy(a_jac, a_aff);
    G1.copy(b_jac, b_aff);
    G1.copy(c_jac, c_aff);

    const size_t N_MIXED = 300000;
    const size_t N_FULL  = 200000;
    const size_t N_DBL   = 300000;

    G1Point acc_mixed;
    G1.copy(acc_mixed, a_jac);
    double ms_mixed = bench_ms([&]() {
        for (size_t i = 0; i < N_MIXED; i++) {
            G1.add(acc_mixed, acc_mixed, b_aff);
        }
    });

    G1Point acc_full;
    G1.copy(acc_full, a_jac);
    double ms_full = bench_ms([&]() {
        for (size_t i = 0; i < N_FULL; i++) {
            G1.add(acc_full, acc_full, b_jac);
        }
    });

    G1Point acc_dbl;
    G1.copy(acc_dbl, c_jac);
    double ms_dbl = bench_ms([&]() {
        for (size_t i = 0; i < N_DBL; i++) {
            G1.dbl(acc_dbl, acc_dbl);
        }
    });

    ASSERT_FALSE(G1.isZero(acc_mixed));
    ASSERT_FALSE(G1.isZero(acc_full));
    ASSERT_FALSE(G1.isZero(acc_dbl));

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[perf][nogmp] G1.add mixed : " << ms_mixed << " ms (N=" << N_MIXED
              << ") => " << ns_per_op(ms_mixed, N_MIXED) << " ns/op\n";
    std::cout << "[perf][nogmp] G1.add full  : " << ms_full << " ms (N=" << N_FULL
              << ") => " << ns_per_op(ms_full, N_FULL) << " ns/op\n";
    std::cout << "[perf][nogmp] G1.dbl       : " << ms_dbl << " ms (N=" << N_DBL
              << ") => " << ns_per_op(ms_dbl, N_DBL) << " ns/op\n";
}

TEST(altBn128Perf, multiexp_phases_timing_nogmp) {
    const int NMExp = 40000;

    typedef mp_uint_t Scalar;

    Scalar *scalars = new Scalar[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    double ms_prepare = bench_ms([&]() {
        for (int i = 0; i < NMExp; i++) {
            if (i == 0) {
                G1.copy(bases[0], G1.one());
            } else {
                G1.add(bases[i], bases[i - 1], G1.one());
            }

            mp_uint_t x;
            mp_set(x, (uint64_t)(i + 1));
            mp_copy(scalars[i], x);
        }
    });

    G1Point r1;
    double ms_multiexp = bench_ms([&]() {
        G1.multiMulByScalar(r1, bases, (uint8_t *)scalars, MP_N, NMExp);
    });

    G1Point r2;
    double ms_msm = bench_ms([&]() {
        G1.multiMulByScalarMSM(r2, bases, (uint8_t *)scalars, MP_N, NMExp);
    });

    ASSERT_FALSE(G1.isZero(r1));
    ASSERT_FALSE(G1.isZero(r2));

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[perf][nogmp] prepare bases/scalars : " << ms_prepare << " ms (N=" << NMExp << ")\n";
    std::cout << "[perf][nogmp] G1.multiMulByScalar   : " << ms_multiexp << " ms (N=" << NMExp << ")\n";
    std::cout << "[perf][nogmp] G1.multiMulByScalarMSM: " << ms_msm << " ms (N=" << NMExp << ")\n";

    delete[] bases;
    delete[] scalars;
}



TEST(altBn128Perf, affine_build_path_timing_nogmp) {
    G1PointAffine p1_aff, p2_aff, out_aff;
    G1.copy(p1_aff, G1.one());
    G1.add(p2_aff, p1_aff, G1.one());   // 2P as affine

    G1Point out_jac;
    G1Point tmp_jac;

    const size_t N_ADD_AFFINE = 200000;
    const size_t N_TO_AFFINE  = 200000;
    const size_t N_BUILD_STEP = 100000;

    // 1) pure affine+affine -> Jacobian
    double ms_add_affine = bench_ms([&]() {
        for (size_t i = 0; i < N_ADD_AFFINE; i++) {
            G1.add(out_jac, p1_aff, p2_aff);
        }
    });

    // 2) pure Jacobian -> affine
    G1.add(tmp_jac, p1_aff, p2_aff);
    double ms_to_affine = bench_ms([&]() {
        for (size_t i = 0; i < N_TO_AFFINE; i++) {
            G1.copy(out_aff, tmp_jac);
        }
    });

    // 3) exact prepare-step shape:
    //    add(affine, affine) + copy(PointAffine <- Point)
    G1PointAffine cur_aff;
    G1.copy(cur_aff, p1_aff);

    double ms_build_step = bench_ms([&]() {
        for (size_t i = 0; i < N_BUILD_STEP; i++) {
            G1.add(tmp_jac, cur_aff, p1_aff);
            G1.copy(cur_aff, tmp_jac);
        }
    });

    ASSERT_FALSE(G1.isZero(out_aff));
    ASSERT_FALSE(G1.isZero(cur_aff));

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[perf][nogmp] G1.add affine+affine -> jac : "
              << ms_add_affine << " ms (N=" << N_ADD_AFFINE
              << ") => " << ns_per_op(ms_add_affine, N_ADD_AFFINE) << " ns/op\n";

    std::cout << "[perf][nogmp] G1.copy jac -> affine       : "
              << ms_to_affine << " ms (N=" << N_TO_AFFINE
              << ") => " << ns_per_op(ms_to_affine, N_TO_AFFINE) << " ns/op\n";

    std::cout << "[perf][nogmp] build step (add_aff + toAff): "
              << ms_build_step << " ms (N=" << N_BUILD_STEP
              << ") => " << ns_per_op(ms_build_step, N_BUILD_STEP) << " ns/op\n";
}

TEST(altBn128Perf, prepare_breakdown_nogmp) {
    const int NMExp = 40000;

    typedef mp_uint_t Scalar;

    Scalar *scalars = new Scalar[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    double ms_scalars = bench_ms([&]() {
        for (int i = 0; i < NMExp; i++) {
            mp_uint_t x;
            mp_set(x, (uint64_t)(i + 1));
            mp_copy(scalars[i], x);
        }
    });

    double ms_bases = bench_ms([&]() {
        for (int i = 0; i < NMExp; i++) {
            if (i == 0) {
                G1.copy(bases[0], G1.one());
            } else {
                G1.add(bases[i], bases[i - 1], G1.one());
            }
        }
    });

    ASSERT_FALSE(G1.isZero(bases[NMExp - 1]));

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[perf][nogmp] prepare scalars only : " << ms_scalars << " ms (N=" << NMExp << ")\n";
    std::cout << "[perf][nogmp] prepare bases only   : " << ms_bases   << " ms (N=" << NMExp << ")\n";

    delete[] bases;
    delete[] scalars;
}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/* Copyright (C) 2012-2017 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>
#include <FHE.h>
#include <timing.h>
#include <EncryptedArray.h>
#include <NTL/lzz_pXFactoring.h>

#include <cstdio>

#include "gtest/gtest.h"
#include "test_common.h"

#include "debugging.h"

/**************

  1. c1.multiplyBy(c0)
  2. c0 += random constant
  3. c2 *= random constant
  4. tmp = c1
  5. ea.shift(tmp, random amount in [-nSlots/2, nSlots/2])
  6. c2 += tmp
  7. ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
  8. c1.negate()
  9. c3.multiplyBy(c2)
  10. c0 -= c3

 **************/

namespace {

::testing::AssertionResult ciphertextMatches(const EncryptedArray &ea, const FHESecKey &sk, const PlaintextArray &p, const Ctxt &c)
{
    PlaintextArray pp(ea);
    ea.decrypt(c, sk, pp);
    if(equals(ea, pp, p)) {
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure()
            << "Ciphertext does not match plaintext:" << std::endl
            << "p = " << p << std::endl
            << "pp = " << pp << std::endl;
    }
}

struct Parameters {
    Parameters(long R, long p, long r, long d, long c, long k, long L, long s,long m, NTL::Vec<long> mvec, NTL::Vec<long> gens, NTL::Vec<long> ords, long seed, long nt) :
        R(R),
        p(p),
        r(r),
        d(d),
        c(c),
        k(k),
        L(L),
        s(s),
        m(m),
        mvec(mvec),
        gens(gens),
        ords(ords),
        seed(seed),
        nt(nt)
    {};

    long R; // number of rounds
    long p; // plaintext base
    long r; // lifting
    long d; // degree of the field extension
    // Note: d == 0 => factors[0] defines extension
    long c; // number of columns in the key-switching matrices
    long k; // security parameter
    long L; // # of bits in the modulus chain
    long s; // minimum number of slots
    long m; // use specified value as modulus
    NTL::Vec<long> mvec; // use product of the integers as modulus
    // Note: mvec takes priority over m
    NTL::Vec<long> gens; // use specified vector of generators
    NTL::Vec<long> ords; // use specified vector of orders
    // e.g., ords=[4 2 -4], negative means 'bad'
    long seed; // PRG seed
    long nt; // num threads

    // Let googletest know how to print the Parameters
    friend std::ostream& operator<<(std::ostream& os, const Parameters& params) {
        return os << "{"
            << "R=" << params.R << ","
            << "p=" << params.p << ","
            << "r=" << params.r << ","
            << "d=" << params.d << ","
            << "c=" << params.c << ","
            << "k=" << params.k << ","
            << "L=" << params.L << ","
            << "s=" << params.s << ","
            << "m=" << params.m << ","
            << "mvec=" << params.mvec << ","
            << "gens=" << params.gens << ","
            << "ords=" << params.ords << ","
            << "seed=" << params.seed << ","
            << "nt=" << params.nt
            << "}";
    };
};

class GTest_General : public ::testing::TestWithParam<Parameters> {
    protected:
        long R;
        long p;
        long r;
        long d;
        long c;
        long k;
        long w;
        long L;
        long m;

        std::vector<long> gens;
        std::vector<long> ords;

        FHEcontext context;

        FHESecKey secretKey;

        const FHEPubKey &publicKey;

        GTest_General () :
            R(GetParam().R),
            p(GetParam().p),
            r(GetParam().r),
            d(GetParam().d),
            c(GetParam().c),
            k(GetParam().k),
            w(64),
            L(GetParam().L),
            m(FindM(k, L, c, p, d, GetParam().s, GetParam().mvec.length() > 0 ? computeProd(GetParam().mvec) : GetParam().m, true)),
            gens(convert<std::vector<long>, NTL::Vec<long>>(GetParam().gens)),
            ords(convert<std::vector<long>, NTL::Vec<long>>(GetParam().ords)),
            context(m, p, r, gens, ords),
            secretKey((buildModChain(context, L, c), context)),
            publicKey(secretKey)
            {
            }

        virtual void SetUp() override
        {
            NTL::SetSeed(NTL::ZZ(GetParam().seed));
            NTL::SetNumThreads(GetParam().nt);
            secretKey.GenSecKey(w); // A Hamming-weight-w secret key
            addSome1DMatrices(secretKey); // compute key-switching matrices that we need
        };

        virtual void TearDown() override
        {
            cleanupGlobals();
        }
};

TEST_P(GTest_General, correctly_implements_mix_of_operations_over_four_ciphertexts)
{
    char buffer[32];
    if (!helib_test::noPrint) {
        std::cout << "\n\n******** " << (helib_test::dry ? "(dry run):" : ":");
        std::cout << " R=" << R
            << ", p=" << p
            << ", r=" << r
            << ", d=" << d
            << ", c=" << c
            << ", k=" << k
            << ", w=" << w
            << ", L=" << L
            << ", m=" << m
            << ", gens=" << gens
            << ", ords=" << ords
            << std::endl;
    }


    NTL::ZZX G;
    if (d == 0)
        G = context.alMod.getFactorsOverZZ()[0];
    else
        G = makeIrredPoly(p, d);

    if (!helib_test::noPrint) {
        context.zMStar.printout();
        std::cout << std::endl;

        std::cout << "security=" << context.securityLevel()<<std::endl;
        std::cout << "# small primes = " << context.smallPrimes.card() << "\n";
        std::cout << "# ctxt primes = " << context.ctxtPrimes.card() << "\n";
        std::cout << "# bits in ctxt primes = "
            << long(context.logOfProduct(context.ctxtPrimes)/log(2.0) + 0.5) << "\n";
        std::cout << "# special primes = " << context.specialPrimes.card() << "\n";
        std::cout << "# bits in special primes = "
            << long(context.logOfProduct(context.specialPrimes)/log(2.0) + 0.5) << "\n";
        std::cout << "G = " << G << "\n";
    }
    EncryptedArray ea(context, G);
    long nslots = ea.size();

    // Debugging additions
    dbgKey = &secretKey;
    dbgEa = &ea;


    PlaintextArray p0(ea);
    PlaintextArray p1(ea);
    PlaintextArray p2(ea);
    PlaintextArray p3(ea);

    random(ea, p0);
    random(ea, p1);
    random(ea, p2);
    random(ea, p3);

    Ctxt c0(publicKey), c1(publicKey), c2(publicKey), c3(publicKey);
    ea.encrypt(c0, publicKey, p0);
    // {ZZX ppp0; ea.encode(ppp0, p0); c0.DummyEncrypt(ppp0);} // dummy encryption
    ea.encrypt(c1, publicKey, p1); // real encryption
    ea.encrypt(c2, publicKey, p2); // real encryption
    ea.encrypt(c3, publicKey, p3); // real encryption

    resetAllTimers();

    FHE_NTIMER_START(Circuit);

    for (long i = 0; i < R; i++) {

        if (!helib_test::noPrint) std::cout << "*** round " << i << "..."<<std::endl;

        long shamt = NTL::RandomBnd(2*(nslots/2) + 1) - (nslots/2);
        // random number in [-nslots/2..nslots/2]
        long rotamt = NTL::RandomBnd(2*nslots - 1) - (nslots - 1);
        // random number in [-(nslots-1)..nslots-1]

        // two random constants
        PlaintextArray const1(ea);
        PlaintextArray const2(ea);
        random(ea, const1);
        random(ea, const2);

        NTL::ZZX const1_poly, const2_poly;
        ea.encode(const1_poly, const1);
        ea.encode(const2_poly, const2);

        mul(ea, p1, p0);     // c1.multiplyBy(c0)
        c1.multiplyBy(c0);
        if (!helib_test::noPrint) CheckCtxt(c1, "c1*=c0");
        EXPECT_TRUE(ciphertextMatches(ea,secretKey,p1,c1));

        add(ea, p0, const1); // c0 += random constant
        c0.addConstant(const1_poly);
        if (!helib_test::noPrint) CheckCtxt(c0, "c0+=k1");
        EXPECT_TRUE(ciphertextMatches(ea,secretKey,p0,c0));

        mul(ea, p2, const2); // c2 *= random constant
        c2.multByConstant(const2_poly);
        if (!helib_test::noPrint) CheckCtxt(c2, "c2*=k2");
        EXPECT_TRUE(ciphertextMatches(ea,secretKey,p2,c2));

        PlaintextArray tmp_p(p1); // tmp = c1
        Ctxt tmp(c1);
        sprintf(buffer, "tmp=c1>>=%d", (int)shamt);
        shift(ea, tmp_p, shamt); // ea.shift(tmp, random amount in [-nSlots/2,nSlots/2])
        ea.shift(tmp, shamt);
        if (!helib_test::noPrint) CheckCtxt(tmp, buffer);
        EXPECT_TRUE(ciphertextMatches(ea,secretKey,tmp_p,tmp));

        add(ea, p2, tmp_p);  // c2 += tmp
        c2 += tmp;
        if (!helib_test::noPrint) CheckCtxt(c2, "c2+=tmp");
        EXPECT_TRUE(ciphertextMatches(ea,secretKey,p2,c2));

        sprintf(buffer, "c2>>>=%d", (int)rotamt);
        rotate(ea, p2, rotamt); // ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
        ea.rotate(c2, rotamt);
        if (!helib_test::noPrint) CheckCtxt(c2, buffer);
        EXPECT_TRUE(ciphertextMatches(ea,secretKey,p2,c2));

        ::negate(ea, p1); // c1.negate()
        c1.negate();
        if (!helib_test::noPrint) CheckCtxt(c1, "c1=-c1");
        EXPECT_TRUE(ciphertextMatches(ea,secretKey,p1,c1));

        mul(ea, p3, p2); // c3.multiplyBy(c2)
        c3.multiplyBy(c2);
        if (!helib_test::noPrint) CheckCtxt(c3, "c3*=c2");
        EXPECT_TRUE(ciphertextMatches(ea,secretKey,p3,c3));

        sub(ea, p0, p3); // c0 -= c3
        c0 -= c3;
        if (!helib_test::noPrint) CheckCtxt(c0, "c0=-c3");
        EXPECT_TRUE(ciphertextMatches(ea,secretKey,p0,c0));
    }

    c0.cleanUp();
    c1.cleanUp();
    c2.cleanUp();
    c3.cleanUp();

    FHE_NTIMER_STOP(Circuit);

    if (!helib_test::noPrint) {
        std::cout << std::endl;
        printAllTimers();
        std::cout << std::endl;
    }
    resetAllTimers();
    FHE_NTIMER_START(Check);

    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p0, c0));
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p1, c1));
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p2, c2));
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p3, c3));

    FHE_NTIMER_STOP(Check);

    std::cout << std::endl;
    if (!helib_test::noPrint) {
        printAllTimers();
        std::cout << std::endl;
    }
};

INSTANTIATE_TEST_SUITE_P(various_parameters, GTest_General, ::testing::Values(
            //         R, p, r, d, c,  k,   L, s,  m,             mvec,             gens,             ords, seed, nt
            //FAST
            //Parameters(1, 2, 1, 1, 2, 80, 500, 0, 91, NTL::Vec<long>{}, NTL::Vec<long>{}, NTL::Vec<long>{},    0,  1)
            //DEEP
            Parameters(1, 2, 2, 1, 2, 10, 500, 0, 91, NTL::Vec<long>{}, NTL::Vec<long>{}, NTL::Vec<long>{},    0,  1),
            Parameters(1, 2, 1, 2, 2, 10, 500, 0, 91, NTL::Vec<long>{}, NTL::Vec<long>{}, NTL::Vec<long>{},    0,  1),
            Parameters(2, 7, 2, 1, 2, 10, 500, 0, 91, NTL::Vec<long>{}, NTL::Vec<long>{}, NTL::Vec<long>{},    0,  1)
            ));

// parameters to get an example where phi(m) is very
// close to m:
// m=18631 L=10 R=5

} // anonymous namespace

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

#include <NTL/BasicThreadPool.h>

#include "EvalMap.h"
#include "hypercube.h"
#include "powerful.h"
#include "debugging.h"

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

struct Parameters {
    Parameters(long p, long r, long c, long k, long L, long s, long seed, NTL::Vec<long> mvec, std::vector<long> gens, std::vector<long> ords, long nthreads, long useCache) :
        p(p),
        r(r),
        c(c),
        k(k),
        L(L),
        s(s),
        seed(seed),
        mvec(mvec),
        gens(gens),
        ords(ords),
        nthreads(nthreads),
        useCache(useCache)
    {};

    long p; // plaintext base
    long r; // lifting
    long c; // number of columns in the key-switching matrices
    long k; // security parameter
    long L; // # of bits in the modulus chain
    long s; // minimum number of slots
    long seed; // PRG seed
    NTL::Vec<long> mvec; // use specified factorization of m
    // e.g., mvec='[7 3 221]'
    std::vector<long> gens; // use specified vector of generators
    // e.g., gens='[3979 3095 3760]'
    std::vector<long> ords; // use specified vector of orders
    // e.g., ords='[6 2 -8]', negative means 'bad'
    long nthreads; // number of threads
    long useCache; // 0: zzX cache, 2: DCRT cache

    friend std::ostream& operator<<(std::ostream& os, const Parameters& params) {
        return os << "{" <<
            "p=" << params.p << "," <<
            "r=" << params.r << "," <<
            "c=" << params.c << "," <<
            "k=" << params.k << "," <<
            "L=" << params.L << "," <<
            "s=" << params.s << "," <<
            "seed=" << params.seed << "," <<
            "mvec=" << params.mvec << "," <<
            "gens=" << params.gens << "," <<
            "ords=" << params.ords << "," <<
            "nthreads=" << params.nthreads << "," <<
            "useCache=" << params.useCache <<
            "}";
    };
};

class GTest_EvalMap : public ::testing::TestWithParam<Parameters> {
    protected:
        long p;
        long r;
        long c;
        long k;
        long L;
        long s;
        long seed;
        NTL::Vec<long> mvec;
        long m;
        std::vector<long> gens;
        std::vector<long> ords;
        long nthreads;
        long useCache;
        FHEcontext context;
        long d;
        long phim;
        long nslots;
        FHESecKey secretKey;
        const FHEPubKey &publicKey;

        static NTL::Vec<long> getDefaultMvec()
        {
            NTL::Vec<long> defaultMvec;
            defaultMvec.SetLength(3);
            defaultMvec[0] = 7;
            defaultMvec[1] = 3;
            defaultMvec[2] = 221;
            return defaultMvec;
        };

        static std::vector<long> getDefaultGens()
        {
            return std::vector<long>{3979, 3095, 3760};
        };

        static std::vector<long> getDefaultOrds()
        {
            return std::vector<long>{6, 2, -8};
        };

        // Throws if mvec is invalid, returns mvec otherwise
        static NTL::Vec<long> checkMvec(NTL::Vec<long> mvec)
        {
            // mvec is supposed to include the prime-power factorization of m
            long nfactors = mvec.length();
            for (long i = 0; i < nfactors; i++)
                for (long j = i+1; j < nfactors; j++)
                    if(NTL::GCD(mvec[i], mvec[j]) != 1)
                        throw std::invalid_argument("mvec invalid: mvec["+std::to_string(i)+"]="+std::to_string(mvec[i]) +
                                ", mvec["+std::to_string(j)+"]="+std::to_string(mvec[j])+" are not coprime.");
            return mvec;
        };

        static long calculateM(const NTL::Vec<long> &mvec, long p)
        {
            // multiply all the prime powers to get m itself
            long m = computeProd(mvec);
            if(NTL::GCD(p, m) != 1)
                throw std::invalid_argument("mvec invalid: computeProd(mvec) and p are not coprime");
            return m;
        };

        static void prepareContext(FHEcontext &context, long L, long c)
        {
            buildModChain(context, L, c);

            if (!helib_test::noPrint) {
                context.zMStar.printout(); // print structure of Zm* /(p) to std::cout
                std::cout << std::endl;
            }
        };

        GTest_EvalMap() : 
            p(GetParam().p),
            r(GetParam().r),
            c(GetParam().c),
            k(GetParam().k),
            L(GetParam().L),
            s(GetParam().s),
            seed(GetParam().seed),
            mvec(checkMvec(
                        lsize(GetParam().mvec) < 1 ? getDefaultMvec() : GetParam().mvec)
                ),
            m(calculateM(mvec, p)),
            gens(lsize(GetParam().mvec) < 1 ? getDefaultGens() : GetParam().gens),
            ords(lsize(GetParam().mvec) < 1 ? getDefaultOrds() : GetParam().ords),
            nthreads(GetParam().nthreads),
            useCache(GetParam().useCache),
            context((setDryRun(false), m), p, r, gens, ords), // We need to get a 'real' context to test EvalMap.
            d((prepareContext(context, L, c), context.zMStar.getOrdP())),
            phim(context.zMStar.getPhiM()),
            nslots(phim/d),
            secretKey((setDryRun(helib_test::dry),context)),  // We can now switch to dry run if desired.
            publicKey(secretKey)
            {
            };

        virtual void SetUp()
        {
            NTL::SetNumThreads(nthreads);
            NTL::SetSeed(NTL::conv<NTL::ZZ>(seed));
            setTimersOn();

            secretKey.GenSecKey(); // A Hamming-weight-w secret key
            addSome1DMatrices(secretKey); // compute key-switching matrices that we need
            addFrbMatrices(secretKey); // compute key-switching matrices that we need
        };

        virtual void TearDown() override
        {
            cleanupGlobals();
        }
};


TEST_P(GTest_EvalMap, eval_map_behaves_correctly)
{
    // GG defines the plaintext space Z_p[X]/GG(X)
    NTL::ZZX GG;
    GG = context.alMod.getFactorsOverZZ()[0];
    EncryptedArray ea(context, GG);

    NTL::zz_p::init(context.alMod.getPPowR());
    NTL::zz_pX F;
    random(F, phim); // a random polynomial of degree phi(m)-1 modulo p

    // convert F to powerful representation: cube represents a multi-variate
    // polynomial with as many variables Xi as factors mi in mvec. cube has
    // degree phi(mi) in the variable Xi, and the coefficients are given
    // in lexicographic order.

    // compute tables for converting between powerful and zz_pX
    PowerfulTranslationIndexes ind(mvec); // indpendent of p
    PowerfulConversion pConv(ind);        // depends on p

    HyperCube<NTL::zz_p> cube(pConv.getShortSig());
    pConv.polyToPowerful(cube, F);

    // Sanity check: convert back and compare
    NTL::zz_pX F2;
    pConv.powerfulToPoly(F2, cube);
    EXPECT_EQ(F, F2) << " @@@ conversion error ):";

    // pack the coefficients from cube in the plaintext slots: the j'th
    // slot contains the polynomial pj(X) = \sum_{t=0}^{d-1} cube[jd+t] X^t
    std::vector<NTL::ZZX> val1;
    val1.resize(nslots);
    for (long i = 0; i < phim; i++) {
        val1[i/d] += NTL::conv<NTL::ZZX>(NTL::conv<NTL::ZZ>(cube[i])) << (i % d);
    }
    PlaintextArray pa1(ea);
    encode(ea, pa1, val1);

    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, pa1);

    resetAllTimers();
    FHE_NTIMER_START(ALL);

    // Compute homomorphically the transformation that takes the
    // coefficients packed in the slots and produces the polynomial
    // corresponding to cube

    if (!helib_test::noPrint) CheckCtxt(ctxt, "init");

    if (!helib_test::noPrint) std::cout << "build EvalMap\n";
    EvalMap map(ea, /*minimal=*/false, mvec, 
            /*invert=*/false, /*build_cache=*/false, /*normal_basis=*/false); 
    // compute the transformation to apply

    if (!helib_test::noPrint) std::cout << "apply EvalMap\n";
    if (useCache) map.upgrade();
    map.apply(ctxt); // apply the transformation to ctxt
    if (!helib_test::noPrint) CheckCtxt(ctxt, "EvalMap");
    if (!helib_test::noPrint) std::cout << "check results\n";

    NTL::ZZX FF1;
    secretKey.Decrypt(FF1, ctxt);
    NTL::zz_pX F1 = NTL::conv<NTL::zz_pX>(FF1);

    EXPECT_EQ(F1, F);

    publicKey.Encrypt(ctxt, FF1);
    if (!helib_test::noPrint) CheckCtxt(ctxt, "init");

    // Compute homomorphically the inverse transformation that takes the
    // polynomial corresponding to cube and produces the coefficients
    // packed in the slots

    if (!helib_test::noPrint) std::cout << "build EvalMap\n";
    EvalMap imap(ea, /*minimal=*/false, mvec, 
            /*invert=*/true, /*build_cache=*/false, /*normal_basis=*/false); 
    // compute the transformation to apply
    if (!helib_test::noPrint) std::cout << "apply EvalMap\n";
    if (useCache) imap.upgrade();
    imap.apply(ctxt); // apply the transformation to ctxt
    if (!helib_test::noPrint) {
        CheckCtxt(ctxt, "EvalMap");
        std::cout << "check results\n";
    }
    PlaintextArray pa2(ea);
    ea.decrypt(ctxt, secretKey, pa2);

    EXPECT_TRUE(equals(ea, pa1, pa2));

    FHE_NTIMER_STOP(ALL);

    if (!helib_test::noPrint) {
        std::cout << "\n*********\n";
        printAllTimers();
        std::cout << std::endl;
    }
};

INSTANTIATE_TEST_SUITE_P(some_parameters, GTest_EvalMap, ::testing::Values(
            //SLOW
            Parameters(2, 1, 2, 80, 300, 0, 0, convert<NTL::Vec<long>, std::vector<long>>(std::vector<long>{7, 3, 221}), std::vector<long>{3979, 3095, 3760}, std::vector<long>{6, 2, -8}, 1, 0)
            //FAST
            //Parameters(2, 1, 2, 80, 300, 0, 0, convert<NTL::Vec<long>, std::vector<long>>(std::vector<long>{3, 35}), std::vector<long>{71, 76}, std::vector<long>{2, 2}, 1, 0)
            ));

} // namespace

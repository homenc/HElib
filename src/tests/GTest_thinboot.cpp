/* Copyright (C) 2012-2018 IBM Corp.
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
#include "FHE.h"
#include "EncryptedArray.h"
#include "matmul.h"
#include "debugging.h"
#include "fhe_stats.h"


#include "gtest/gtest.h"
#include "test_common.h"

extern long fhe_force_chen_han;

namespace {
//NTL_CLIENT

/********************************************************************
 ********************************************************************/
    struct Parameters {
        long p; // plaintext base
        long r; // exponent
        // p^r is the plaintext-space modulus
        long c; // number of columns in the key-switching matrices
        long bits; // # of bits in the modulus chain
        long skHwt; // Hamming weight of recryption secret key
        long nthreads; // number of threads
        long seed; // random number seed
        int useCache; // 0: zzX cache, 1: DCRT cache
        int force_bsgs;
        int force_hoist;
        // int disable_intFactor // fhe_disable_intFactor (disabled by Victor)
        int chen_han;
        bool debug; // generate debugging output
        int scale; // scale parameter
        NTL::Vec<long> global_gens;
        NTL::Vec<long> global_ords;
        NTL::Vec<long> global_mvec;
        const int c_m; // = 100;
        const long iter;

        Parameters(long p, long r, long c, long bits, long skHwt, long nthreads, long seed, long useCache,
                   int c_m, int force_bsgs, int force_hoist, int chen_han, bool debug, int scale,
                   const std::vector<long> &global_gens, const std::vector<long> &global_ords,
                   const std::vector<long> &global_mvec, long iter) :
                p(p),
                r(r),
                c(c),
                bits(bits),
                skHwt(skHwt),
                nthreads(nthreads),
                seed(seed),
                useCache(useCache),
                force_bsgs(force_bsgs),
                force_hoist(force_hoist),
                chen_han(chen_han),
                debug(debug),
                scale(scale),
                global_gens(convert<NTL::Vec<long>>(global_gens)),
                global_ords(convert<NTL::Vec<long>>(global_ords)),
                global_mvec(convert<NTL::Vec<long>>(global_mvec)),
                c_m(c_m),
                iter(iter) {
            if (global_gens.empty() || global_ords.empty() || global_mvec.empty())
                throw helib::LogicError("gens, ords, and mvec must be non-empty");
        };

        friend std::ostream &operator<<(std::ostream &os, const Parameters &params) {
            return os << "{"
                      << "p" << params.p << ","
                      << "r" << params.r << ","
                      << "c" << params.c << ","
                      << "bits" << params.bits << ","
                      << "skHwt" << params.skHwt << ","
                      << "nthreads" << params.nthreads << ","
                      << "seed" << params.seed << ","
                      << "useCache" << params.useCache << ","
                      << "force_bsgs" << params.force_bsgs << ","
                      << "force_hoist" << params.force_hoist << ","
                      << "chen_han" << params.chen_han << ","
                      << "debug" << params.debug << ","
                      << "scale" << params.scale << ","
                      << "global_gens" << params.global_gens << ","
                      << "global_ords" << params.global_ords << ","
                      << "global_mvec" << params.global_mvec << ","
                      << "c_m" << params.c_m << ","
                      << "iter" << params.iter << ","
                      << "}";
        }

//
//
//  if (seed)
//    SetSeed(ZZ(seed));
//
//  SetNumThreads(nthreads);
//
//  TestIt(p,r,bits,c,t,useCache);

    };

    class GTest_thinboot : public ::testing::TestWithParam<Parameters> {
    private:


        void preContextSetup() {

            if (!helib_test::noPrint) {
                std::cout << "*** GTest_thinboot";
                if (isDryRun())
                    std::cout << " (dry run)";
                std::cout << ": p=" << p
                          << ", r=" << r
                          << ", bits=" << bits
                          << ", t=" << skHwt
                          << ", c=" << c
                          << ", m=" << m
                          << " mvec=" << mvec << ", gens=" << gens << ", ords=" << ords
                          << std::endl;
                std::cout << "Computing key-independent tables..." << std::flush;
            }
            setTimersOn();
            setDryRun(false); // Need to get a "real context" to test bootstrapping
            time = -NTL::GetTime();
        }

        void postContextSetup() {
            if (scale) {
                context.scale = scale;
            }

            context.zMStar.set_cM(c_m / 100.0);
        }

        static void setGlobals(int force_bsgs, int force_hoist, int chen_han) {
            fhe_test_force_bsgs = force_bsgs;
            fhe_test_force_hoist = force_hoist;
            fhe_force_chen_han = chen_han;
        };

        void cleanupBootstrappingGlobals() {
            fhe_test_force_bsgs = old_fhe_test_force_bsgs;
            fhe_test_force_hoist = old_fhe_test_force_hoist;
            fhe_force_chen_han = old_fhe_force_chen_han;
        }

        static void setSeedIfNeeded(const long seed) {
            if (seed)
                SetSeed(NTL::ZZ(seed));
        };

        static void checkPM(const long p, const long m) {
            helib::assertTrue(NTL::GCD(p, m) == 1, "GCD(p, m) == 1");
        }

    protected:
        const int old_fhe_test_force_bsgs;
        const int old_fhe_test_force_hoist;
        const int old_fhe_force_chen_han;
        const long p;
        const long r;
        const long c;
        const long bits;
        const long skHwt;
        const long nthreads;
        const long seed;
        const long useCache;
        const int force_bsgs;
        const int force_hoist;
        const int chen_han;
        const bool debug;
        const int scale;
        const NTL::Vec<long> mvec;
        const std::vector<long> gens;
        const std::vector<long> ords;
        const int c_m; // = 100;
        const long iter;

        const long m, phim;
        double time;
        FHEcontext context;

        GTest_thinboot() :
                old_fhe_test_force_bsgs(fhe_test_force_bsgs),
                old_fhe_test_force_hoist(fhe_test_force_hoist),
                old_fhe_force_chen_han(fhe_force_chen_han),
                p((setGlobals(GetParam().force_bsgs, GetParam().force_hoist, GetParam().chen_han),
                   GetParam().p)),
                r(GetParam().r),
                c(GetParam().c),
                bits(GetParam().bits),
                skHwt(GetParam().skHwt),
                nthreads((NTL::SetNumThreads(GetParam().nthreads), GetParam().nthreads)),
                seed((setSeedIfNeeded(GetParam().seed), GetParam().seed)),
                useCache(GetParam().useCache),
                force_bsgs(GetParam().force_bsgs),
                force_hoist(GetParam().force_hoist),
                chen_han(GetParam().chen_han),
                debug(GetParam().debug),
                scale(GetParam().scale),
                gens(convert<std::vector<long>>(GetParam().global_gens)),
                ords(convert<std::vector<long>>(GetParam().global_ords)),
                mvec(GetParam().global_mvec),
                c_m(GetParam().c_m),
                iter(GetParam().iter),
                m(computeProd(mvec)),
                phim((checkPM(p, m), phi_N(m))),
                time(0),
                context((preContextSetup(), m), p, r, gens, ords)
        { postContextSetup(); }

        void TearDown() override
        {
            if(!helib_test::noPrint) {
                printAllTimers();
            }
            if (fhe_stats)
                print_stats(std::cout);

            cleanupBootstrappingGlobals();
            cleanupGlobals();
        }

    };

    TEST_P(GTest_thinboot, correctly_performs_thinboot) {
        buildModChain(context, bits, c, /*willBeBootstrappable=*/true, /*t=*/skHwt);

        if (!helib_test::noPrint) {
            std::cout << "security=" << context.securityLevel() << std::endl;
            std::cout << "# small primes = " << context.smallPrimes.card() << std::endl;
            std::cout << "# ctxt primes = " << context.ctxtPrimes.card() << std::endl;
            std::cout << "# bits in ctxt primes = "
                      << long(context.logOfProduct(context.ctxtPrimes) / log(2.0) + 0.5) << std::endl;
            std::cout << "# special primes = " << context.specialPrimes.card() << std::endl;
            std::cout << "# bits in special primes = "
                      << long(context.logOfProduct(context.specialPrimes) / log(2.0) + 0.5) << std::endl;
            std::cout << "scale=" << context.scale << std::endl;
        }

        context.makeBootstrappable(mvec,/*t=*/skHwt, useCache,/*alsoThick=*/false);
        // save time...disable some fat boot precomputation

        time += NTL::GetTime();

        //if (skHwt>0) context.rcData.skHwt = skHwt;
        if (!helib_test::noPrint) {
            std::cout << " done in " << time << " seconds" << std::endl;
            std::cout << "  e=" << context.rcData.e
                 << ", e'=" << context.rcData.ePrime
                 << ", t=" << context.rcData.skHwt
                 << std::endl << "  ";
            context.zMStar.printout();
        }
        setDryRun(helib_test::dry); // Now we can set the dry-run flag if desired


        long p2r = context.alMod.getPPowR();

        if (!helib_test::noPrint) fhe_stats = true;

        for (long numkey = 0; numkey < iter; numkey++) { // test with 3 keys
            if (fhe_stats && numkey > 0 && numkey % 100 == 0)
                print_stats(std::cout);

            if (!helib_test::noPrint) std::cerr << "*********** iter=" << (numkey + 1) << std::endl;

            time = -NTL::GetTime();
            if (!helib_test::noPrint) std::cout << "Generating keys, " << std::flush;
            FHESecKey secretKey(context);
            secretKey.GenSecKey(skHwt);      // A Hamming-weight-64 secret key
            addSome1DMatrices(secretKey); // compute key-switching matrices that we need
            addFrbMatrices(secretKey);
            if (!helib_test::noPrint) std::cout << "computing key-dependent tables..." << std::flush;
            secretKey.genRecryptData();
            time += NTL::GetTime();
            if (!helib_test::noPrint) std::cout << " done in " << time << " seconds\n";

#ifdef DEBUG_PRINTOUT
            dbgEa = (EncryptedArray*) context.ea;
            dbgKey = &secretKey;
#endif

            FHEPubKey publicKey = secretKey;

            long d = context.zMStar.getOrdP();
            long phim = context.zMStar.getPhiM();
            long nslots = phim / d;

            // GG defines the plaintext space Z_p[X]/GG(X)
            NTL::ZZX GG;
            GG = context.alMod.getFactorsOverZZ()[0];
            EncryptedArray ea(context, GG);

            if (debug) {
                dbgKey = &secretKey;
                dbgEa = &ea;
            }

            NTL::zz_p::init(p2r);
            NTL::Vec<NTL::zz_p> val0(NTL::INIT_SIZE, nslots);
            for (auto &x: val0)
                random(x);

            std::vector<NTL::ZZX> val1;
            val1.resize(nslots);
            for (long i = 0; i < nslots; i++) {
                val1[i] = NTL::conv<NTL::ZZX>(NTL::conv<NTL::ZZ>(rep(val0[i])));
            }


            Ctxt c1(publicKey);
            ea.encrypt(c1, publicKey, val1);

            Ctxt c2(c1);

            if (!helib_test::noPrint) CheckCtxt(c2, "before squarings");

            long sqr_count = -1;
            Ctxt next_c2(c2);
            do {
                c2 = next_c2;
                next_c2.multiplyBy(next_c2);
                sqr_count++;
            } while (next_c2.bitCapacity() >= 100);

            if (!helib_test::noPrint) {
                std::cout << "sqr_count=" << sqr_count << std::endl;
                if (sqr_count > 0) {
                    std::cout << "BITS-PER-LEVEL: "
                         << ((c1.bitCapacity() - c2.bitCapacity()) / double(sqr_count)) << "\n";
                }
                CheckCtxt(c2, "before recryption");
            }

            resetAllTimers();

            {
                FHE_NTIMER_START(AAA_thinRecrypt);

                publicKey.thinReCrypt(c2);

            }


            if (!helib_test::noPrint) CheckCtxt(c2, "after recryption");

            for (auto &x: val0) {
                for (long i: range(sqr_count)) x = x * x;
            }

            for (long i = 0; i < nslots; i++) {
                val1[i] = NTL::conv<NTL::ZZX>(NTL::conv<NTL::ZZ>(rep(val0[i])));
            }

            std::vector<NTL::ZZX> val2;
            ea.decrypt(c2, secretKey, val2);

            EXPECT_EQ(val1, val2);
        }
    }

// LEGACY TEST DEFAULT PARAMETERS:
//long p=2;
//long r=1;
//long c=3;
//long bits=600;
//long t=64;
//long nthreads=1;
//long seed=0;
//long useCache=1;
//int c_m = 100;

    INSTANTIATE_TEST_SUITE_P(typical_parameters, GTest_thinboot, ::testing::Values(
            //SLOW
            Parameters(2, 1, 3, 600, 64, 1, 0, 1, 100, 0, 0, 0, 0, 0, {1026, 249}, {30, -2}, {31, 41}, 1),
            Parameters(17, 1, 3, 600, 64, 1, 0, 1, 100, 0, 0, 0, 0, 0, {556, 1037}, {6, 4}, {7, 5, 37}, 1)
            //FAST
            //Parameters(2, 1, 3, 600, 64, 1, 0, 1, 100, 0, 0, 0, 0, 0, {1026, 249}, {30, -2}, {31, 41}, 1)
    ));
} // namespace


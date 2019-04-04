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
// testPacking.cxx - testing uppack/repack functionality
#include "intraSlot.h"
#include "debugging.h"

#include "gtest/gtest.h"
#include "test_common.h"

namespace{

struct Parameters {

    Parameters(long p, long n, long r, long L, long m, long seed) :
        p(p),
        n(n),
        r(r),
        L(L),
        m(m),
        seed(seed)
    {};
    long p; // plaintext base
    long n; // number of packed ciphertexts
    long r; // lifting
    long L; // # of levels in the modulus chain
    long m; // use specified value as modulus
    long seed; // PRG seed
    
    // Let googletest know how to print the Parameters
    friend std::ostream& operator<<(std::ostream& os, const Parameters& params) {
        return os << "{"
            << "p=" << params.p << ","
            << "n=" << params.n << ","
            << "r=" << params.r << ","
            << "L=" << params.L << ","
            << "m=" << params.m << ","
            << "seed=" << params.seed
            << "}";
    };
};

class GTest_intraSlot : public ::testing::TestWithParam<Parameters> {

    static FHEcontext& setupContext(FHEcontext& context, long L)
    {
        if(helib_test::verbose) {
            context.zMStar.printout();
        }
        buildModChain(context, L, 3);
        return context;
    };

    protected:
        GTest_intraSlot () :
            p(GetParam().p),
            n(GetParam().n),
            r(GetParam().r),
            L(GetParam().L),
            m(GetParam().m),
            seed(GetParam().seed),
            context(m, p, r),
            secretKey(setupContext(context, L)),
            publicKey(secretKey)
    {};

        long p;
        long n;
        long r;
        long L;
        long m;
        long seed;
        FHEcontext context;
        FHESecKey secretKey;
        const FHEPubKey& publicKey;

        virtual void SetUp()
        {
            SetSeed(NTL::ZZ(seed));
            secretKey.GenSecKey(); // A +-1/0 secret key
            addSome1DMatrices(secretKey); // compute key-switching matrices that we need
            addFrbMatrices(secretKey);
        };

        virtual void TearDown() override
        {
            cleanupGlobals();
        }
};

TEST_P(GTest_intraSlot, packing_and_unpacking_works)
{
    NTL::ZZX G = context.alMod.getFactorsOverZZ()[0];
    EncryptedArray ea(context, G);

    long d = ea.getDegree(); // size of each slot

    std::vector<Ctxt> unpacked(d*n -1, Ctxt(publicKey));

    // generate (almost) d*n ciphertexts, with only integrs in the slots
    std::vector<PlaintextArray> p1(lsize(unpacked), PlaintextArray(ea));
    for (long i=0; i<lsize(unpacked); i++) {
        std::vector<long> slots;
        ea.random(slots);
        encode(ea, p1[i] ,slots);
        ea.encrypt(unpacked[i], publicKey, p1[i]);
    }

    // Pack (almost) d*n ciphetexts into only n of them
    std::vector<Ctxt> ct(n, Ctxt(publicKey));
    repack(CtPtrs_vectorCt(ct), CtPtrs_vectorCt(unpacked), ea);

    // Unpack them back
    std::vector<zzX> unpackSlotEncoding;
    buildUnpackSlotEncoding(unpackSlotEncoding, ea);
    unpack(CtPtrs_vectorCt(unpacked), CtPtrs_vectorCt(ct), ea, unpackSlotEncoding);

    PlaintextArray p2(ea);
    for (long i=0; i<lsize(unpacked); i++) {
        ea.decrypt(unpacked[i], secretKey, p2);
        ASSERT_TRUE(equals(ea, p1[i], p2)) <<
            "p2["<<i<<"]="<<p2;
    }
}
INSTANTIATE_TEST_SUITE_P(some_parameters, GTest_intraSlot, ::testing::Values(
            //FAST
            Parameters(2, 2, 1, 10, 91, 0)
            ));
} // namespace

//#include <helib/helib.h>

#include "gtest/gtest.h"
#include "test_common.h"

// only run if HEXL has been linked.
#ifdef USE_INTEL_HEXL
#include "../src/intelExt.h" // Private header

namespace {

TEST(TestHEXL, hexlInUse)
{

  //long N = 8; // This is phi_m
  //long modulus = 769;

  // For powers of 2, 2N == 2phi_m == m
  // FindPrimitiveRoot(ZZ_p(2 * N, modulus), );

  std::vector<long> args = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  auto expected_outputs = args;

  //intel::TinkyWinky(args.data(), args.data(), N, modulus);
  //intel::AltFFTRev1(args.data(), args.data(), N, modulus);

  EXPECT_TRUE(std::equal(args.begin(), args.end(), expected_outputs.begin()));
}

} // namespace

#else

namespace {

TEST(TestHEXL, noTestRequired) {}

} // namespace

#endif // USE_INTEL_HEXL

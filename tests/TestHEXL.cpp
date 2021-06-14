
// only run if HEXL has been linked.
#ifdef USE_INTEL_HEXL

#include "../src/intelExt.h" // Private header

#include <helib/helib.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

TEST(TestHEXL, hexlInUse)
{
   intel::initNTT(2,2,2);
}

} // namespace

#else

TEST(TestHEXL, noTestRequired)
{}

#endif // USE_INTEL_HEXL

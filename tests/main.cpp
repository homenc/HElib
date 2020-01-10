#include "gtest/gtest.h"
#include "test_common.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    // Now argc and argv will have had their gtest_* entries
    // processed and removed.  Parse for our own purposes.
    helib_test::parse_common_args(argc, argv);
    return RUN_ALL_TESTS();
}


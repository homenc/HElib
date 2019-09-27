#include "test_common.h"
#include <iostream>
#include <random>

namespace helib_test {
    char *path_of_executable = nullptr;
    bool noPrint = false;
    bool verbose = false;
    bool dry = false;
    unsigned int random_seed = 0U;
    
    void parse_common_args(int argc, char *argv[])
    {
        ArgMap amap;
        path_of_executable = argv[0];
        amap.arg("dry", dry, "dry=1 for a dry-run");
        amap.arg("noPrint", noPrint, "suppress printouts");
        amap.arg("verbose", verbose, "print more information");
        amap.arg("seed", random_seed, "specify random seed for test data");
        amap.parse(argc, argv);
        if(random_seed == 0U) // Not specified: use random seed
          random_seed = std::random_device{}();
        // TODO: change this printout so that the random_seed is simply another
        // parameter of parameterised tests.  This will only be possible once 
        // we parse the command-line args before gtest does.
        std::cout << "random seed: " << random_seed << std::endl;
    };
};


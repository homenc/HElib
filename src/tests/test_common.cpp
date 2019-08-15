#include "test_common.h"

namespace helib_test {
    char *path_of_executable = nullptr;
    bool noPrint = false;
    bool verbose = false;
    bool dry = false;

    
    void parse_common_args(int argc, char *argv[])
    {
        ArgMap amap;
        path_of_executable = argv[0];
        amap.arg("dry", dry, "dry=1 for a dry-run");
        amap.arg("noPrint", noPrint, "suppress printouts");
        amap.arg("verbose", verbose, "print more information");
        amap.parse(argc, argv);
    };
};


#ifndef TEST_COMMON_H
#define TEST_COMMON_H
#include "ArgMap.h" 

namespace helib_test {
    extern char *path_of_executable;
    extern bool noPrint;
    extern bool verbose;
    extern bool dry;

    static void parse_common_args(int argc, char *argv[])
    {
        ArgMap amap;
        path_of_executable = argv[0];
        amap.arg("dry", dry, "dry=1 for a dry-run");
        amap.arg("noPrint", noPrint, "suppress printouts");
        amap.arg("verbose", verbose, "print more information");
        amap.parse(argc, argv);
    };

};

#endif /* ifndef TEST_COMMON_H */


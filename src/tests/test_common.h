#ifndef TEST_COMMON_H
#define TEST_COMMON_H
#include "ArgMap.h" 

namespace helib_test {
    extern char *path_of_executable;
    extern bool noPrint;
    extern bool verbose;
    extern bool dry;

    void parse_common_args(int argc, char *argv[]);

};

#endif /* ifndef TEST_COMMON_H */


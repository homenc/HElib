#include <helib/version.h>

namespace helib {

// This is stored in the compiled library.
static constexpr char versionInLib[] = "v@PROJECT_VERSION@";

const char* version::libString() { return versionInLib; }

} // namespace helib

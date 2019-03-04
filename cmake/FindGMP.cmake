# Use cmake standard find_library package
include(FindPackageHandleStandardArgs)

if (GMP_DIR)
  # If user-specified folders: look there
  find_library(GMP_LIB NAMES gmp libgmp
               PATHS ${GMP_DIR}
               PATH_SUFFIXES lib
               NO_DEFAULT_PATH
               DOC "GMP library")

  find_path(GMP_HEADERS NAMES gmp.h
            PATHS ${GMP_DIR}
            PATH_SUFFIXES include
            NO_DEFAULT_PATH
            DOC "GMP headers")

else(GMP_DIR)
# Else: look in default paths
  find_library(GMP_LIB NAMES gmp libgmp
               PATH_SUFFIXES lib
               DOC "GMP library")

  find_path(GMP_HEADERS NAMES gmp.h
            PATH_SUFFIXES include
            DOC "GMP headers")
endif(GMP_DIR)

if (GMP_HEADERS AND GMP_LIB)
  # Find gmp version
  file(STRINGS "${GMP_HEADERS}/gmp.h" gmp_major REGEX "__GNU_MP_VERSION[ \t]+([0-9]+)")
  string(REGEX REPLACE "[^ \t]*[ \t]+__GNU_MP_VERSION[ \t]+([0-9]+)" "\\1" gmp_major "${gmp_major}")
  file(STRINGS "${GMP_HEADERS}/gmp.h" gmp_minor REGEX "__GNU_MP_VERSION_MINOR[ \t]+([0-9]+)")
  string(REGEX REPLACE "[^ \t]*[ \t]+__GNU_MP_VERSION_MINOR[ \t]+([0-9]+)" "\\1" gmp_minor "${gmp_minor}")
  file(STRINGS "${GMP_HEADERS}/gmp.h" gmp_patchlevel REGEX "__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+)")
  string(REGEX REPLACE "[^ \t]*[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+)" "\\1" gmp_patchlevel "${gmp_patchlevel}")
  if((gmp_major EQUAL "") OR
     (gmp_minor EQUAL "") OR
     (gmp_patchlevel EQUAL ""))
    # If the version encoding is wrong then it is set to "WRONG VERSION ENCODING" causing find_package_handle_standard_args to fail
    set(GMP_VERSION "GMP BAD VERSION ENCODING")
    set(GMP_VERSION_MAJOR "GMP BAD VERSION ENCODING")
    set(GMP_VERSION_MINOR "GMP BAD VERSION ENCODING")
    set(GMP_VERSION_PATCH "GMP BAD VERSION ENCODING")
  else()
    set(GMP_VERSION "${gmp_major}.${gmp_minor}.${gmp_patchlevel}")
    set(GMP_VERSION_MAJOR "${gmp_major}")
    set(GMP_VERSION_MINOR "${gmp_minor}")
    set(GMP_VERSION_PATCH "${gmp_patchlevel}")
  endif()

  unset(gmp_major)
  unset(gmp_minor)
  unset(gmp_patchlevel)
endif()

# Raising an error if gmp with required version has not been found
set(fail_msg "GMP required dynamic shared library has not been found. (Try cmake -DGMP_DIR=<GMP-root-path>).")
if (GMP_DIR)
  set(fail_msg "GMP required dynamic shared library has not been found in ${GMP_DIR}.")
endif(GMP_DIR)

# Check the library has been found, handle QUIET/REQUIRED arguments and set GMP_FOUND accordingly or raise the error
find_package_handle_standard_args(GMP
                                  REQUIRED_VARS GMP_LIB GMP_HEADERS
                                  VERSION_VAR GMP_VERSION
                                  FAIL_MESSAGE "${fail_msg}")

unset(fail_msg)

# If GMP has been found set the default variables
if (GMP_FOUND)
  set(GMP_LIBRARIES "${GMP_LIB}")
  set(GMP_INCLUDE_PATHS "${GMP_HEADERS}")
endif(GMP_FOUND)

# Mark GMP_LIBRARIES GMP_INCLUDE_PATHS as advanced variables
mark_as_advanced(GMP_LIBRARIES GMP_INCLUDE_PATHS)

# Use cmake standard find_library package
include(FindPackageHandleStandardArgs)

if (GMP_DIR)
  # If user-specified folders: look there
  find_library(GMP_LIB NAMES gmp libgmp
               PATHS ${GMP_DIR}
               PATH_SUFFIXES lib ${CMAKE_INSTALL_LIBDIR}
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
               PATH_SUFFIXES lib ${CMAKE_INSTALL_LIBDIR}
               DOC "GMP library")

  find_path(GMP_HEADERS NAMES gmp.h
            PATH_SUFFIXES include
            DOC "GMP headers")
endif(GMP_DIR)

if (GMP_HEADERS AND GMP_LIB)
  # Find gmp version
  try_run(gmp_ver_program_run_code gmp_ver_program_compile_code
    "${CMAKE_CURRENT_BINARY_DIR}/get_gmp_version"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/get_gmp_version.c"
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${GMP_HEADERS}"
    LINK_LIBRARIES "${GMP_LIB}"
    RUN_OUTPUT_VARIABLE gmp_version_string)

  if(NOT (gmp_ver_program_compile_code AND (gmp_ver_program_run_code EQUAL 0)))
    message(FATAL_ERROR "Failed to determine gmp version.")
  endif()

  STRING(REGEX REPLACE "([0-9]*)\.[0-9]*\.[0-9]*" "\\1" gmp_major "${gmp_version_string}")
  STRING(REGEX REPLACE "[0-9]*\.([0-9]*)\.[0-9]*" "\\1" gmp_minor "${gmp_version_string}")
  STRING(REGEX REPLACE "[0-9]*\.[0-9]*\.([0-9]*)" "\\1" gmp_patchlevel "${gmp_version_string}")

  if((gmp_major STREQUAL "") OR
     (gmp_minor STREQUAL "") OR
     (gmp_patchlevel STREQUAL ""))
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
  unset(gmp_version_string)
  unset(gmp_ver_program_compile_code)
  unset(gmp_ver_program_run_code)
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

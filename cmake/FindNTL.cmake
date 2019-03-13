# Use cmake standard find_library package
include(FindPackageHandleStandardArgs)

if (NTL_DIR)
  # If user-specified folders: look there
  find_library(NTL_LIB NAMES ntl libntl
               PATHS ${NTL_DIR}
               PATH_SUFFIXES lib
               NO_DEFAULT_PATH
               DOC "NTL library")

  find_path(NTL_HEADERS NAMES config.h
            PATHS ${NTL_DIR}
            PATH_SUFFIXES include/NTL
            NO_DEFAULT_PATH
            DOC "NTL headers")

else (NTL_DIR)
# Else: look in default paths
  find_library(NTL_LIB NAMES ntl libntl
               PATH_SUFFIXES lib
               DOC "NTL library")

  find_path(NTL_HEADERS NAMES config.h
            PATH_SUFFIXES include/NTL
            DOC "NTL headers")
endif(NTL_DIR)

if (NTL_HEADERS AND NTL_LIB)
  # Find ntl version
  try_run(ntl_ver_program_run_code ntl_ver_program_compile_code
    "${CMAKE_CURRENT_BINARY_DIR}/get_ntl_version"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/get_ntl_version.c"
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${NTL_HEADERS}"
    LINK_LIBRARIES "${NTL_LIB}"
    RUN_OUTPUT_VARIABLE ntl_version_string)

  if(NOT (${ntl_ver_program_compile_code} AND ${ntl_ver_program_run_code} EQUAL 0))
    message(FATAL_ERROR "Failed to determine ntl version.")
  endif()

  string(REGEX REPLACE "([0-9]+)\.([0-9]+)\.([0-9]+)" "\\1" ntl_major "${ntl_version_string}")
  string(REGEX REPLACE "([0-9]+)\.([0-9]+)\.([0-9]+)" "\\2" ntl_minor "${ntl_version_string}")
  string(REGEX REPLACE "([0-9]+)\.([0-9]+)\.([0-9]+)" "\\3" ntl_patchlevel "${ntl_version_string}")
  if((ntl_version_string EQUAL "") OR
     (ntl_major EQUAL "") OR
     (ntl_minor EQUAL "") OR
     (ntl_patchlevel EQUAL ""))
    # If the version encoding is wrong then it is set to "WRONG VERSION ENCODING" causing find_package_handle_standard_args to fail
    set(NTL_VERSION "WRONG VERSION ENCODING")
    set(NTL_VERSION_MAJOR "WRONG VERSION ENCODING")
    set(NTL_VERSION_MINOR "WRONG VERSION ENCODING")
    set(NTL_VERSION_PATCH "WRONG VERSION ENCODING")
  else()
    set(NTL_VERSION "${ntl_major}.${ntl_minor}.${ntl_patchlevel}")
    set(NTL_VERSION_MAJOR "${ntl_major}")
    set(NTL_VERSION_MINOR "${ntl_minor}")
    set(NTL_VERSION_PATCH "${ntl_patchlevel}")
  endif()

  unset(ntl_version_string)
  unset(ntl_major)
  unset(ntl_minor)
  unset(ntl_patchlevel)
endif(NTL_HEADERS AND NTL_LIB)

# Raising an error if ntl with required version has not been found
set(fail_msg "NTL required dynamic shared library has not been found. (Try cmake -DNTL_DIR=<NTL-root-path>).")
if (NTL_DIR)
  set(fail_msg "NTL required dynamic shared library has not been found in ${NTL_DIR}.")
endif(NTL_DIR)

# Check the library has been found, handle QUIET/REQUIRED arguments and set NTL_FOUND accordingly or raise the error
find_package_handle_standard_args(NTL
                                  REQUIRED_VARS NTL_LIB NTL_HEADERS
                                  VERSION_VAR NTL_VERSION
                                  FAIL_MESSAGE "${fail_msg}")

unset(fail_msg)

# If NTL has been found set the default variables
if (NTL_FOUND)
  set(NTL_LIBRARIES "${NTL_LIB}")
  set(NTL_INCLUDE_PATHS "${NTL_HEADERS}/..")
endif(NTL_FOUND)

# Mark NTL_LIBRARIES NTL_INCLUDE_PATHS as advanced variables
mark_as_advanced(NTL_LIBRARIES NTL_INCLUDE_PATHS)

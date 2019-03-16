# File cmake/rpath_fixer containing the function change_rpath
include(rpath_fixer)

# Adding meta-target dependent on gmp and ntl to delay the rpath fix
add_custom_target(dependencies_built)

# Adding meta-target upon which ntl and gmp will depend, in order
# to standardise installation directory (lib/lib64)
add_custom_target(setup_install_dir)

add_custom_command(TARGET "setup_install_dir"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory "${PACKAGE_DIR}/${CMAKE_INSTALL_LIBDIR}")
               
if (NOT CMAKE_INSTALL_LIBDIR STREQUAL "lib")
  add_custom_command(TARGET "setup_install_dir"
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_INSTALL_LIBDIR}" "lib"
                     WORKING_DIRECTORY "${PACKAGE_DIR}"
                     VERBATIM)
endif()

if (FETCH_GMP)
  ###########################################
  #                                         #
  #              Loading GMP                #
  #                                         #
  ###########################################
  add_subdirectory(gmp)
endif(FETCH_GMP)

###########################################
#                                         #
#              Loading NTL                #
#                                         #
###########################################
add_subdirectory(ntl)


if (FETCH_GMP)
  add_dependencies(gmp_external dependencies_built)
endif(FETCH_GMP)
add_dependencies(ntl_external dependencies_built)

# Fixing rpaths
# We are delaying the rpath fix after ntl is built since ntl in the configure step
# uses internal binaries with non relative rpath that links against gmp
# Getting gmp/ntl library locations
get_target_property(gmp_install_location gmp_external IMPORTED_LOCATION)
get_target_property(ntl_install_location ntl_external IMPORTED_LOCATION)
# Function definition
# change_rpath(lib_name_noext depend_target_name lib_path package_relative_rpath gmp_library_path)
if (FETCH_GMP)
  # Fixing gmp rpath
  change_rpath("gmp" dependencies_built "${gmp_install_location}" "${PACKAGE_RPATH}" "")
endif(FETCH_GMP)
# Fixing ntl rpath
change_rpath("ntl" dependencies_built "${ntl_install_location}" "${PACKAGE_RPATH}" "${gmp_install_location}")

# Unsetting temporary variables
unset(gmp_install_location)
unset(ntl_install_location)

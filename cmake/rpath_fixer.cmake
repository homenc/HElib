# Copyright (C) 2019-2020 IBM Corp.
#
# This program is Licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. See accompanying LICENSE file.

# Function to change external libraries rpath on mac and linux
function (change_rpath lib_name_noext
                       depend_target_name
                       lib_path
                       package_relative_rpath
                       gmp_library_path)
  # gmp_library_dir is the directory containing the libgmp.(so|dylib) file
  get_filename_component(gmp_library_dir "${gmp_library_path}" DIRECTORY)
  if (APPLE OR (CMAKE_CXX_PLATFORM_ID STREQUAL "Darwin"))
    if (NOT CMAKE_INSTALL_NAME_TOOL)
      message(FATAL_ERROR "CMAKE_INSTALL_NAME_TOOL is not set.")
    endif (NOT CMAKE_INSTALL_NAME_TOOL)

    # Find the readlink program to get the versioned name of the library
    find_program(READLINK_CMD NAMES readlink)
    if (NOT READLINK_CMD)
      message(FATAL_ERROR
              "Required readlink command not found.  Cannot fix rpaths.")
    endif (NOT READLINK_CMD)

    # set id as rpath
    set(change_id_command
        "LIBVERNAME=`${READLINK_CMD} ${lib_path}` && ${CMAKE_INSTALL_NAME_TOOL} -id @rpath/\${LIBVERNAME} ${lib_path}")
    add_custom_command(TARGET ${depend_target_name}
                       POST_BUILD
                       COMMAND eval
                       ARGS ${change_id_command}
                       VERBATIM)

    # set local relative rpath
    # hack to ignore install_name_tool weird errors when the
    # @loader_path/${package_relative_rpath} is already present in the binary
    # NOTE: This can be turned into an external script invoking it with
    # cmake -P to properly handle error
    set(add_rpath_command
        "${CMAKE_INSTALL_NAME_TOOL} -add_rpath @loader_path/${package_relative_rpath} ${lib_path}")
    add_custom_command(
        TARGET ${depend_target_name}
        POST_BUILD
        COMMAND eval
        ARGS "${add_rpath_command} 2> /dev/null || true"
        VERBATIM)

    if (NOT FETCH_GMP)
      # Add GMP location to NTL rpath if GMP is not fetched
      # NOTE: This can be turned into an external script invoking it with
      # cmake -P to properly handle error
      set(add_gmp_rpath_command
          "${CMAKE_INSTALL_NAME_TOOL} -add_rpath ${gmp_library_dir} ${lib_path}")
      add_custom_command(
          TARGET ${depend_target_name}
          POST_BUILD
          COMMAND eval
          ARGS "${add_gmp_rpath_command} 2> /dev/null || true"
          VERBATIM)
      unset(add_gmp_rpath_command)
    endif (NOT FETCH_GMP)

    # Change ID of GMP in the NTL library only if GMP is not the default one
    if (lib_name_noext MATCHES "ntl" AND FETCH_GMP)
      # NOTE: Here we assume the ID of the gmp library is equal to its absolute
      # path, including version name (i.e. not the symlink).
      # If GMP changes this convention, the following command may break.
      # Since gmp rpath fix has been delayed after ntl compilation due the
      # configuration step dependencies we have to change ntl rpath to gmp.
      set(change_gmp_rpath_cmd
          "LIBVERNAME=`${READLINK_CMD} ${gmp_library_path}` && LIB_FULL=\"${gmp_library_dir}/\${LIBVERNAME}\" && ${CMAKE_INSTALL_NAME_TOOL} -change \${LIB_FULL} @rpath/\${LIBVERNAME} ${lib_path}")
      add_custom_command(TARGET ${depend_target_name}
                         POST_BUILD
                         COMMAND eval
                         ARGS "${change_gmp_rpath_cmd}"
                         VERBATIM)
      unset(change_gmp_rpath_cmd)
    endif (lib_name_noext MATCHES "ntl" AND FETCH_GMP)
    unset(change_id_command)
    unset(add_rpath_command)

  elseif ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux"
          OR (CMAKE_CXX_PLATFORM_ID EQUAL "ELF"))
    # Find patchelf required to add a local rpath on linux
    find_program(LINUX_RPATH_TOOL NAMES patchelf)
    if (NOT LINUX_RPATH_TOOL)
      message(FATAL_ERROR
              "Cannot find patchelf, which is required for a package build.")
    endif (NOT LINUX_RPATH_TOOL)

    if (FETCH_GMP)
      set(rpath_patch_args "--set-rpath"
                           "$ORIGIN/${package_relative_rpath}"
                           "${lib_path}")
    else (FETCH_GMP)
      set(rpath_patch_args "--set-rpath"
                           "$ORIGIN/${package_relative_rpath}:${gmp_library_dir}"
                           "${lib_path}")
    endif (FETCH_GMP)
    # Adding $ORIGIN/${package_relative_rpath}
    add_custom_command(TARGET ${depend_target_name}
                       POST_BUILD
                       COMMAND ${LINUX_RPATH_TOOL}
                       ARGS ${rpath_patch_args}
                       VERBATIM)

    unset(rpath_patch_args)
  else ()
    message(FATAL_ERROR "OS not supported.")
  endif ()
endfunction (change_rpath)

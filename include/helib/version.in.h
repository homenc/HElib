/* Copyright (C) 2020 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

#ifndef HELIB_VERSION_H
#define HELIB_VERSION_H

namespace helib {

/**
 * @class version
 * @brief The class acts as a namespace with all members static.
 * Holds the version number for this code of HElib.
 **/
struct version
{

  // clang-format off
  /**
   * @brief The major number of this version of HElib.
   **/
  static constexpr long major = @PROJECT_VERSION_MAJOR@;
  /**
   * @brief The minor number of this version of HElib.
   **/
  static constexpr long minor = @PROJECT_VERSION_MINOR@;
  /**
   * @brief The patch number of this version of HElib.
   **/
  static constexpr long patch = @PROJECT_VERSION_PATCH@;
  // clang-format on

  /**
   * @brief The string representation of this version of HElib.
   **/
  static constexpr auto asString = "@PROJECT_VERSION@";

  /**
   * @brief Function that returns whether this version of HElib is equal to or
   * higher than a specified version.
   * @param major The major version number.
   * @param minor The minor version number.
   * @param patch The patch version number.
   * @return `true` if current HElib version is greater than or equal to the
   * specified version version.
   **/
  static inline constexpr bool greaterEquals(long major_,
                                             long minor_ = 0,
                                             long patch_ = 0)
  {
    if (major_ < 0 || minor_ < 0 || patch_ < 0)
      return false;
    long min_version = (((major_ << 8) ^ minor_) << 8) ^ patch_;
    long our_version = (((major << 8) ^ minor) << 8) ^ patch;
    return our_version >= min_version;
  }

  /**
   * @brief Get the string version from the HElib compiled library instead of
   * the one defined in the header.
   * @return A string representing the version of HElib stored in the compiled
   * library instead of from the one defined in the header.
   **/
  static const char* libString();

}; // struct version

} // namespace helib

#endif // HELIB_VERSION_H

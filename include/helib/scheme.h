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

#ifndef HELIB_SCHEME_H
#define HELIB_SCHEME_H

/**
 * @file scheme.h
 * @brief CKKS and BGV scheme tags contained as definitions of `CKKS` and `BGV`
 * structs.
 **/

namespace helib {

class PolyMod;

/**
 * @brief Type for CKKS scheme, to be used as template parameter.
 **/
struct CKKS
{
  /**
   * @brief Slot type used for CKKS plaintexts: `std::complex<double>`.
   **/
  using SlotType = std::complex<double>;

  /**
   * @brief Scheme label to be added to JSON serialization.
   */
  static constexpr std::string_view schemeName = "CKKS";
};

/**
 * @brief Type for BGV scheme, to be used as template parameter.
 **/
struct BGV
{
  /**
   * @brief Slot type used for BGV plaintexts: `helib::PolyMod` i.e. an integer
   * polynomial modulo p^r and G.
   **/
  using SlotType = PolyMod;

  /**
   * @brief Scheme label to be added to JSON serialization.
   */
  static constexpr std::string_view schemeName = "BGV";
};

} // namespace helib

#endif // HELIB_SHEME_H

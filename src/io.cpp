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

#include <helib/NumbTh.h>
#include <helib/exceptions.h>
#include "io.h"
#include <cxxabi.h>

namespace NTL {

void to_json(json& j, const NTL::xdouble& num)
{
  // e.g. 'stdev': {'mantissa': 10, 'exponent': 9 }
  j = {{"mantissa", num.mantissa()}, {"exponent", num.exponent()}};
}

void from_json(const json& j, NTL::xdouble& num)
{
  num.x = j.at("mantissa");
  num.e = j.at("exponent");
}

void to_json(json& j, const NTL::ZZ& num)
{
  std::stringstream str;
  str << num;
  j = {{"number", str.str()}};
}

void from_json(const json& j, NTL::ZZ& num)
{
  std::stringstream str;
  str << j.at("number").get<std::string>();
  str >> num;
}

void to_json(json& j, const NTL::Vec<long>& vec)
{
  j = helib::convert<std::vector<long>>(vec);
}

void from_json(const json& j, NTL::Vec<long>& vec)
{
  std::vector<long> repr = j;
  vec = helib::convert<NTL::Vec<long>>(repr);
}

void to_json(json& j, const NTL::ZZX& poly)
{
  if (poly == NTL::ZZX::zero()) {
    // Avoid string "[]" for zero ZZX.
    j = std::vector<int>(1, 0);
  } else {
    j = helib::convert<std::vector<long>>(poly.rep);
  }
}

void from_json(const json& j, NTL::ZZX& poly)
{
  if (j.is_number()) {
    poly = j.get<long>();
    return;
  } else {
    for (std::size_t i = 0; i < j.size(); ++i) {
      if (j[i].is_number_float()) {
        throw helib::IOError("Bad NTL::ZZX JSON serialization.  Expected an "
                             "integer number, got a floating-point.");
      }
      NTL::SetCoeff(poly, i, j[i].get<long>());
    }
  }
}

} // namespace NTL

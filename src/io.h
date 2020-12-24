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

#ifndef HELIB_IO_H
#define HELIB_IO_H
/**
 * @file IO.h
 * @brief - Internal header (not installed) containing convenience functions for
 * ASCII serialization (using JSON).
 *
 * Copyright IBM Corporation 2020 All rights reserved.
 */

#include <complex>

#include <json.hpp>
using json = ::nlohmann::json;

#include <NTL/vec_long.h>
#include <NTL/xdouble.h>

#include <helib/exceptions.h>
#include <helib/JsonWrapper.h>
#include <helib/version.h>

// For our convenience, this helps us.
namespace NTL {
void to_json(json& j, const NTL::xdouble& num);
void from_json(const json& j, NTL::xdouble& num);

void to_json(json& j, const NTL::ZZ& num);
void from_json(const json& j, NTL::ZZ& num);

void to_json(json& j, const NTL::Vec<long>& vec);
void from_json(const json& j, NTL::Vec<long>& vec);

void to_json(json& j, const NTL::ZZX& poly);
void from_json(const json& j, NTL::ZZX& poly);
} // namespace NTL

namespace std {
template <typename T>
inline void to_json(json& j, const std::complex<T>& num)
{
  // e.g. 'stdev': {'real': 10.1, 'imag': 9.3 }
  j = {num.real(), num.imag()};
}

// Permissive JSON to std::complex<T> function. I accepts single numbers,
// or a JSON object containing a real and/or imag part.
// Fails if the json value is not convertible to the type T.
template <typename T>
inline void from_json(const json& j, std::complex<T>& num)
{
  num.real(0);
  num.imag(0);

  if (j.is_number()) {
    // Set only the imaginary part
    num.real(j.get<T>());
    return;
  } else {
    if (j.size() > 2) {
      throw helib::IOError("Bad complex JSON serialization. Expected a maximum "
                           "of 2 elements, recieved " +
                           std::to_string(j.size()));
    }
    if (j.size() == 0) {
      return;
    }
    num.real(j[0].get<T>());
    if (j.size() == 2) {
      num.imag(j[1].get<T>());
    }
  }
}
} // namespace std

namespace helib {

class Context;

inline const std::string_view jsonSerializationVersion = "0.0.1";

inline JsonWrapper wrap(const json& j)
{
  return JsonWrapper(std::make_any<json>(j));
}

inline json unwrap(const JsonWrapper& jwrap)
{
  try {
    return std::any_cast<json>(jwrap.getJSONobj());
  } catch (const std::bad_any_cast& e) {
    throw LogicError(std::string("Cannot unwrap wrapper. Bad cast ") +
                     e.what());
  }
}

template <typename T>
inline std::vector<T> readVectorFromJSON(const json::array_t& j)
{
  std::vector<T> v;
  v.reserve(j.size());

  for (const auto& e : j) {
    v.emplace_back(T::readFromJSON(wrap(e)));
  }

  return v;
}

template <typename T, typename... TArgs>
inline std::vector<T> readVectorFromJSON(const json::array_t& j,
                                         const TArgs&... args)
{
  std::vector<T> v;
  v.reserve(j.size());

  for (const auto& e : j) {
    v.emplace_back(T::readFromJSON(wrap(e), args...));
  }

  return v;
}

template <typename T>
inline void readVectorFromJSON(const json::array_t& j,
                               std::vector<T>& v,
                               T& init)
{
  std::vector<json> jvec = j;

  v.resize(jvec.size(), init); // Make space in vector

  for (std::size_t i = 0; i < jvec.size(); i++) {
    v[i].readJSON(wrap(jvec[i]));
  }
}

template <typename T>
inline json writeVectorToJSON(const std::vector<T>& ts)
{
  std::vector<json> js;
  for (const auto& t : ts) {
    js.emplace_back(unwrap(t.writeToJSON()));
  }
  return js;
}

template <typename T>
static inline json toTypedJson(const json& tc)
{
  return {{"type", T::typeName},
          {"HElibVersion", version::asString},
          {"serializationVersion", jsonSerializationVersion},
          {"content", tc}};
}

template <typename T>
static inline json fromTypedJson(const json& j)
{
  std::string obj_ser_ver = j.at("serializationVersion").get<std::string>();
  if (obj_ser_ver != jsonSerializationVersion) {
    std::stringstream sstr;
    sstr << "Serialization version mismatch.  Expected: "
         << jsonSerializationVersion << " actual: " << obj_ser_ver;
    throw IOError(sstr.str());
  }

  std::string obj_helib_ver = j.at("HElibVersion").get<std::string>();
  if (obj_helib_ver != version::asString) {
    std::stringstream sstr;
    sstr << "HElib version mismatch.  Expected: " << version::asString
         << " actual: " << obj_helib_ver;
    throw IOError(sstr.str());
  }

  std::string obj_ty = j.at("type").get<std::string>();
  if (obj_ty != T::typeName) {
    std::stringstream fmt;
    fmt << "Type mismatch deserializing json object."
        << "  Expected: " << T::typeName << " actual: " << obj_ty;
    throw IOError(fmt.str());
  }
  return j.at("content");
}

template <typename T, typename TCALL>
inline T executeRedirectJsonError(const TCALL& f)
{
  try {
    return f();
  } catch (const nlohmann::detail::exception& e) {
    throw IOError(std::string("Error with JSON IO. ") + e.what());
  }
}

} // namespace helib

#endif // HELIB_IO_H

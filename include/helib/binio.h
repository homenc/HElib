/* Copyright (C) 2012-2020 IBM Corp.
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
#ifndef HELIB_BINIO_H
#define HELIB_BINIO_H
#include <iostream>
#include <vector>
#include <type_traits>
#include <NTL/xdouble.h>
#include <NTL/vec_long.h>

#define BINIO_32BIT 4
#define BINIO_48BIT 6
#define BINIO_64BIT 8

#define BINIO_EYE_SIZE 4

// clang-format off
#define BINIO_EYE_CONTEXTBASE_BEGIN "|BS["
#define BINIO_EYE_CONTEXTBASE_END   "]BS|"
#define BINIO_EYE_CONTEXT_BEGIN     "|CN["
#define BINIO_EYE_CONTEXT_END       "]CN|"
#define BINIO_EYE_CTXT_BEGIN        "|CX["
#define BINIO_EYE_CTXT_END          "]CX|"
#define BINIO_EYE_PK_BEGIN          "|PK["
#define BINIO_EYE_PK_END            "]PK|"
#define BINIO_EYE_SK_BEGIN          "|SK["
#define BINIO_EYE_SK_END            "]SK|"
#define BINIO_EYE_SKM_BEGIN         "|KM["
#define BINIO_EYE_SKM_END           "]KM|"
// clang-format on

namespace helib {

/* This struct (or similar) is a nice to have not used at the moment. */
// struct BinaryHeader {
//  uint8_t structId[4];
//  uint8_t version[4] = {0, 0, 0, 1};
//  uint64_t id;
//  uint64_t payloadSize;
//};

/* Some utility functions for binary IO */

int readEyeCatcher(std::istream& str, const char* expect);
void writeEyeCatcher(std::ostream& str, const char* eye);

void write_ntl_vec_long(std::ostream& str,
                        const NTL::vec_long& vl,
                        long intSize = BINIO_64BIT);
void read_ntl_vec_long(std::istream& str, NTL::vec_long& vl);

long read_raw_int(std::istream& str);
int read_raw_int32(std::istream& str);
void write_raw_int(std::ostream& str, long num);
void write_raw_int32(std::ostream& str, int num);

void write_raw_double(std::ostream& str, const double d);
double read_raw_double(std::istream& str);

void write_raw_xdouble(std::ostream& str, const NTL::xdouble xd);
NTL::xdouble read_raw_xdouble(std::istream& str);

void write_raw_ZZ(std::ostream& str, const NTL::ZZ& zz);
void read_raw_ZZ(std::istream& str, NTL::ZZ& zz);

template <typename T>
void write_raw_vector(std::ostream& str, const std::vector<T>& v)
{
  write_raw_int(str, v.size());

  for (const T& n : v) {
    n.write(str);
  }
}

// vector<long> has a different implementation, since long.write does not work
template <>
void write_raw_vector<long>(std::ostream& str, const std::vector<long>& v);

// vector<double> has a different implementation, since double.write does not
// work
template <>
void write_raw_vector<double>(std::ostream& str, const std::vector<double>& v);

template <typename T>
void read_raw_vector(std::istream& str, std::vector<T>& v, T& init)
{
  long sz = read_raw_int(str);
  v.resize(sz, init); // Make space in vector

  for (auto& n : v) {
    n.read(str);
  }
}

template <typename T>
void read_raw_vector(std::istream& str, std::vector<T>& v)
{
  read_raw_vector<T>(str, v, T());
}

// vector<long> has a different implementation, since long.read does not work
template <>
void read_raw_vector<long>(std::istream& str, std::vector<long>& v);

// vector<double> has a different implementation, since double.read does not
// work
template <>
void read_raw_vector<double>(std::istream& str, std::vector<double>& v);

// KeySwitch::read(...) (in keySwitching.cpp) requires the context.
class Context;
template <typename T>
void read_raw_vector(std::istream& str,
                     std::vector<T>& v,
                     const Context& context)
{
  long sz = read_raw_int(str);
  v.resize(sz); // Make space in vector

  for (auto& n : v) {
    n.read(str, context);
  }
}

} // namespace helib
#endif // ifndef HELIB_BINIO_H

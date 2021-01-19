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
#include <array>
#include <vector>
#include <type_traits>
#include <cstdint>
#include <sstream>
#include <helib/assertions.h>
#include <helib/version.h>

#include <NTL/xdouble.h>
#include <NTL/vec_long.h>

namespace helib {

struct Binio
{
  static constexpr int BIT32 = 4;
  static constexpr int BIT64 = 8;

  static constexpr std::array<char, 4> VERSION_0_0_1_0 = {0, 0, 1, 0};
};

struct EyeCatcher
{
  static constexpr int SIZE = 4;
  // clang-format off
  static constexpr std::array<char, SIZE> HEADER_BEGIN  = {'|','H','E','['};
  static constexpr std::array<char, SIZE> HEADER_END    = {']','H','E','|'};
  static constexpr std::array<char, SIZE> CONTEXT_BEGIN = {'|','C','N','['};
  static constexpr std::array<char, SIZE> CONTEXT_END   = {']','C','N','|'};
  static constexpr std::array<char, SIZE> CTXT_BEGIN    = {'|','C','X','['};
  static constexpr std::array<char, SIZE> CTXT_END      = {']','C','X','|'};
  static constexpr std::array<char, SIZE> PK_BEGIN      = {'|','P','K','['};
  static constexpr std::array<char, SIZE> PK_END        = {']','P','K','|'};
  static constexpr std::array<char, SIZE> SK_BEGIN      = {'|','S','K','['};
  static constexpr std::array<char, SIZE> SK_END        = {']','S','K','|'};
  static constexpr std::array<char, SIZE> SKM_BEGIN     = {'|','K','M','['};
  static constexpr std::array<char, SIZE> SKM_END       = {']','K','M','|'};
  // clang-format on
};

template <typename T>
inline constexpr char nameToStructId()
{
  static_assert(true, "Type without a struct id.");
  return 0; // Should not reach.
}

class Context;
class PubKey;
class SecKey;
class Ctxt;

template <>
inline constexpr char nameToStructId<Context>()
{
  return 5;
}
template <>
inline constexpr char nameToStructId<PubKey>()
{
  return 10;
}
template <>
inline constexpr char nameToStructId<SecKey>()
{
  return 15;
}
template <>
inline constexpr char nameToStructId<Ctxt>()
{
  return 20;
}

// Already broken into bytes, thus should be the same written and read in bog
// or little endian.
template <typename T>
struct SerializeHeader
{
  // Header eye catcher
  const std::array<char, EyeCatcher::SIZE> beginCatcher =
      EyeCatcher::HEADER_BEGIN;
  // 32 bit number: 8 bits for major, minor, patch, fix
  const std::array<char, 4> version = Binio::VERSION_0_0_1_0;
  // The helib version that output this header.
  const std::array<char, 4> helibVersion = {version::major,
                                            version::minor,
                                            version::patch,
                                            0};
  // ObjectType
  char structId = nameToStructId<T>();
  // Reserved for future use
  const char reserved[7] = {0, 0, 0, 0, 0, 0, 0};
  // End
  const std::array<char, EyeCatcher::SIZE> endCatcher = EyeCatcher::HEADER_END;

  void writeTo(std::ostream& os)
  {
    os.write(reinterpret_cast<const char*>(this), sizeof(*this));
  }

  static SerializeHeader readFrom(std::istream& is)
  {
    SerializeHeader<T> header;

    // Broken in bytes, should be the same written and read in bog or little
    // endian.
    is.read(reinterpret_cast<char*>(&header), sizeof(header));

    // Checks
    if (header.beginCatcher != EyeCatcher::HEADER_BEGIN ||
        header.endCatcher != EyeCatcher::HEADER_END) {

      std::ostringstream oss;
      oss << "Eye catchers for header mismatch '";
      oss.write(header.beginCatcher.data(), EyeCatcher::SIZE);
      oss << ", ";
      oss.write(header.endCatcher.data(), EyeCatcher::SIZE);
      oss << "' (begin, end).";
      throw IOError(oss.str());
    }

    return header;
  }

  std::string versionString() const
  {
    // version has only 4 numbers.
    return std::to_string(version[0]) + "." + std::to_string(version[1]) + "." +
           std::to_string(version[2]) + "." + std::to_string(version[3]);
  }
};

/* Some utility functions for binary IO */

bool readEyeCatcher(std::istream& str,
                    const std::array<char, EyeCatcher::SIZE>& expect);
void writeEyeCatcher(std::ostream& str,
                     const std::array<char, EyeCatcher::SIZE>& eye);

void write_ntl_vec_long(std::ostream& str,
                        const NTL::vec_long& vl,
                        long intSize = Binio::BIT64);
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
    n.writeTo(str);
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

// FIXME: Change other method adding _inplace or this to read_return to
// distinguish them
template <typename T, typename CTy>
std::vector<T> read_raw_vector(std::istream& str, const CTy& ctx)
{
  std::vector<T> v;
  long sz = read_raw_int(str);
  v.reserve(sz); // Make space in vector

  for (long i = 0; i < sz; i++) {
    v.emplace_back(T::readFrom(str, ctx));
  }

  return v;
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

} // namespace helib
#endif // ifndef HELIB_BINIO_H

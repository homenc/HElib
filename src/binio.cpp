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
#include "binio.h"
#include <helib/assertions.h>
#include <sys/types.h> // byte order macros in a platform-independent way.

namespace helib {

/* Some utility functions for binary IO */
bool readEyeCatcher(std::istream& str,
                    const std::array<char, EyeCatcher::SIZE>& expect)
{
  std::array<char, EyeCatcher::SIZE> eye;
  str.read(eye.data(), EyeCatcher::SIZE);
  return eye == expect;
}

void writeEyeCatcher(std::ostream& str,
                     const std::array<char, EyeCatcher::SIZE>& eye)
{
  str.write(eye.data(), eye.size());
}

// compile only 64-bit (-m64) therefore long must be at least 64-bit
long read_raw_int(std::istream& str)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  long result = 0;
  str.read(reinterpret_cast<char*>(&result), Binio::BIT64);
  return result;
#else
  long result = 0;
  char byte;

  for (long i = 0; i < Binio::BIT64; i++) {
    str.read(&byte, 1); // read a byte
    result |= (static_cast<long>(byte) & 0xff)
              << i * 8; // must be in little endian
  }

  return result;
#endif
}

int read_raw_int32(std::istream& str)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  int result = 0;
  str.read(reinterpret_cast<char*>(&result), Binio::BIT32);
  return result;
#else
  int result = 0;
  char byte;

  for (long i = 0; i < Binio::BIT32; i++) {
    str.read(&byte, 1); // read a byte
    result |= (static_cast<long>(byte) & 0xff)
              << i * 8; // must be in little endian
  }

  return result;
#endif
}

// compile only 64-bit (-m64) therefore long must be at least 64-bit
void write_raw_int(std::ostream& str, long num)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  str.write(reinterpret_cast<const char*>(&num), Binio::BIT64);
#else
  char byte;

  for (long i = 0; i < Binio::BIT64; i++) {
    byte = num >> 8 * i; // serializing in little endian
    str.write(&byte, 1); // write byte out
  }
#endif
}

void write_raw_int32(std::ostream& str, int num)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  str.write(reinterpret_cast<const char*>(&num), Binio::BIT32);
#else
  char byte;

  for (long i = 0; i < Binio::BIT32; i++) {
    byte = num >> 8 * i; // serializing in little endian
    str.write(&byte, 1); // write byte out
  }
#endif
}

void write_ntl_vec_long(std::ostream& str,
                        const NTL::vec_long& vl,
                        long intSize)
{
  assertTrue<InvalidArgument>(intSize == Binio::BIT64 ||
                                  intSize == Binio::BIT32,
                              "intSize must be 32 or 64 bit for binary IO");
  write_raw_int32(str, vl.length());
  write_raw_int32(str, intSize);

  if (intSize == Binio::BIT64) {
    for (long i = 0; i < vl.length(); i++) {
      write_raw_int(str, vl[i]);
    }
  } else {
    for (long i = 0; i < vl.length(); i++) {
      write_raw_int32(str, vl[i]);
    }
  }
}

void read_ntl_vec_long(std::istream& str, NTL::vec_long& vl)
{
  int sizeOfVL = read_raw_int32(str);
  int intSize = read_raw_int32(str);
  assertTrue<InvalidArgument>(intSize == Binio::BIT64 ||
                                  intSize == Binio::BIT32,
                              "intSize must be 32 or 64 bit for binary IO");

  // Remember to check and increase Vec before trying to fill it.
  if (vl.length() < sizeOfVL) {
    vl.SetLength(sizeOfVL);
  }

  if (intSize == Binio::BIT64) {
    for (long i = 0; i < sizeOfVL; i++) {
      vl[i] = read_raw_int(str);
    }
  } else {
    for (long i = 0; i < sizeOfVL; i++) {
      vl[i] = read_raw_int32(str);
    }
  }
}

void write_raw_double(std::ostream& str, const double d)
{
  // FIXME: this is not portable:
  //  * we might have sizeof(long) < sizeof(double)
  //  * we also don't know if the bit layout is really compatible
  const long* pd = reinterpret_cast<const long*>(&d);
  write_raw_int(str, *pd);
}

double read_raw_double(std::istream& str)
{
  // FIXME: see FIXME for write_raw_double
  long d = read_raw_int(str);
  double* pd = reinterpret_cast<double*>(&d);
  return *pd;
}

void write_raw_xdouble(std::ostream& str, const NTL::xdouble xd)
{
  double m = xd.mantissa();
  long e = xd.exponent();
  write_raw_double(str, m);
  write_raw_int(str, e);
}

NTL::xdouble read_raw_xdouble(std::istream& str)
{
  double m = read_raw_double(str);
  long e = read_raw_int(str);
  return NTL::xdouble(m, e);
}

void write_raw_ZZ(std::ostream& str, const NTL::ZZ& zz)
{
  long noBytes = NumBytes(zz);
  assertTrue<InvalidArgument>(noBytes > 0,
                              "Number of bytes to write must be non-negative");
  unsigned char* zzBytes = new unsigned char[noBytes];
  BytesFromZZ(zzBytes, zz, noBytes); // From ZZ.h
  write_raw_int(str, noBytes);
  // TODO - ZZ appears to be endian agnostic
  str.write(reinterpret_cast<char*>(zzBytes), noBytes);
  delete[] zzBytes;
}

void read_raw_ZZ(std::istream& str, NTL::ZZ& zz)
{
  long noBytes = read_raw_int(str);
  assertTrue<InvalidArgument>(noBytes > 0,
                              "Number of bytes to write must be non-negative");
  unsigned char* zzBytes = new unsigned char[noBytes];
  // TODO - ZZ appears to be endian agnostic
  str.read(reinterpret_cast<char*>(zzBytes), noBytes);
  zz = NTL::ZZFromBytes(zzBytes, noBytes);
  delete[] zzBytes;
}

// FIXME: there is some repetitive code here.
// We should think about a better overloading strategy for
// read/write_raw_vector to avoid this.

template <>
void read_raw_vector<long>(std::istream& str, std::vector<long>& v)
{

  long sz = read_raw_int(str);
  v.resize(sz); // Make space in vector

  for (long i = 0; i < sz; i++)
    v[i] = read_raw_int(str);
}

template <>
void write_raw_vector<long>(std::ostream& str, const std::vector<long>& v)
{
  long sz = v.size();
  write_raw_int(str, sz);

  for (long n : v)
    write_raw_int(str, n);
}

template <>
void read_raw_vector<double>(std::istream& str, std::vector<double>& v)
{

  long sz = read_raw_int(str);
  v.resize(sz); // Make space in vector

  for (long i = 0; i < sz; i++)
    v[i] = read_raw_double(str);
}

template <>
void write_raw_vector<double>(std::ostream& str, const std::vector<double>& v)
{
  long sz = v.size();
  write_raw_int(str, sz);

  for (double n : v)
    write_raw_double(str, n);
}

} // namespace helib

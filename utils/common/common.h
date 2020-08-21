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

#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <fstream>

#include <helib/helib.h>

std::string stripExtension(const std::string& s)
{
  std::size_t dotPos = s.find_last_of(".");
  return (dotPos == std::string::npos) ? s : s.substr(0, dotPos);
}

std::string readline(std::istream& is)
{
  std::string s;
  getline(is, s);
  return s;
}

void readKeyBinary(std::istream& keyFile, helib::PubKey& pk)
{
  helib::readPubKeyBinary(keyFile, pk);
}

void readKeyBinary(std::istream& keyFile, helib::SecKey& sk)
{
  helib::readSecKeyBinary(keyFile, sk);
}

template <typename T1, typename T2>
using uniq_pair = std::pair<std::unique_ptr<T1>, std::unique_ptr<T2>>;

template <typename KEY>
uniq_pair<helib::Context, KEY> loadContextAndKey(const std::string& keyFilePath)
{
  std::ifstream keyFile(keyFilePath, std::ios::binary);
  if (!keyFile.is_open())
    throw std::runtime_error("Cannot open Public Key file '" + keyFilePath +
                             "'.");

  unsigned long m, p, r;
  std::vector<long> gens, ords;

  helib::readContextBaseBinary(keyFile, m, p, r, gens, ords);
  std::unique_ptr<helib::Context> contextp =
      std::make_unique<helib::Context>(m, p, r, gens, ords);
  helib::readContextBinary(keyFile, *contextp);

  std::unique_ptr<KEY> keyp = std::make_unique<KEY>(*contextp);
  readKeyBinary(keyFile, *keyp);

  return {std::move(contextp), std::move(keyp)};
}

inline long estimateCtxtSize(const helib::Context& context, long offset)
{
  // Return in bytes.

  // We assume that the size of each element in the DCRT is BINIO_64BIT

  // sizeof(BINIO_EYE_CTXT_BEGIN) = 4;
  // BINIO_32BIT = 4
  // sizeof(long) = BINIO_64BIT = 8
  // xdouble = s * sizeof(long) = 2 * BINIO_64BIT = 16

  // We assume that primeSet after encryption is context.ctxtPrimes
  // We assume we have exactly 2 parts after encryption
  // We assume that the DCRT prime set is the same as the ctxt one

  long size = 0;

  size += 4;

  // Begin Ctxt metadata
  // 64 = header_size = ptxtSpace (long) + intFactor (long) + ptxtMag (xdouble)
  //                    + ratFactor (xdouble) + noiseBound (xdouble)
  size += 64;

  // primeSet.write(str);
  // size of set (long) + each prime (long)
  size += 8 + context.ctxtPrimes.card() * 8;

  // Begin Ctxt content size
  // write_raw_vector(str, parts);
  // Size of the parts vector (long)
  size += 8;

  long part_size = 0;
  // Begin CtxtPart size

  // skHandle.write(str);
  // powerOfS (long) + powerOfX (long) + secretKeyID (long)
  part_size += 24;

  // Begin DCRT size computation

  // this->DoubleCRT::write(str);
  // map.getIndexSet().write(str);
  // size of set (long) + each prime (long)
  part_size += 8 + context.ctxtPrimes.card() * 8;

  // DCRT data write as write_ntl_vec_long(str, map[i]);
  // For each prime in the ctxt modulus chain
  //    size of DCRT column (long) + size of each element (long) +
  //    size of all the slots (column in DCRT) (PhiM long elements)
  long dcrt_size =
      (8 + 8 * context.zMStar.getPhiM()) * context.ctxtPrimes.card();

  part_size += dcrt_size;

  // End DCRT size
  // End CtxtPart size

  size += 2 * part_size; // 2 * because we assumed 2 parts
  // End Ctxt content size

  // End eye-catcher
  size += 4;

  return size + offset;
}

#endif // COMMON_H

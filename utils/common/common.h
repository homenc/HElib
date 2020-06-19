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
     throw std::runtime_error(
       "Cannot open Public Key file '" + keyFilePath + "'.");

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

#endif // COMMON_H

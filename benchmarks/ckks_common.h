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

#ifndef HELIB_CKKS_COMMON_H
#define HELIB_CKKS_COMMON_H

#include "bench_common.h" // HE_BENCH_CAPTURE

#include <helib/helib.h>

#include <memory>

struct Params
{
  const long m, r, L;
  Params(long _m, long _r, long _L) : m(_m), r(_r), L(_L) {}
  Params(const Params& other) : Params(other.m, other.r, other.L) {}
  bool operator!=(Params& other) const { return !(*this == other); }
  bool operator==(Params& other) const
  {
    return m == other.m && r == other.r && L == other.L;
  }
};

struct ContextAndKeys
{
  const Params params;

  helib::Context context;
  helib::SecKey secretKey;
  const helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  ContextAndKeys(Params& _params) :
      params(_params),
      context(helib::ContextBuilder<helib::CKKS>()
                  .m(params.m)
                  .precision(params.r)
                  .bits(params.L)
                  .scale(10)
                  .build()),
      secretKey(context),
      publicKey((secretKey.GenSecKey(),
                 helib::addSome1DMatrices(secretKey),
                 secretKey)),
      ea(context.getEA())
  {
    context.printout();
  }
};

struct Meta
{
  std::unique_ptr<ContextAndKeys> data;
  Meta& operator()(Params& params)
  {
    // Only change if nullptr or different.
    if (data == nullptr || data->params != params)
      data = std::make_unique<ContextAndKeys>(params);
    return *this;
  }
};

#endif // HELIB_CKKS_COMMON_H

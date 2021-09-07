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

#ifndef HELIB_BGV_COMMON_H
#define HELIB_BGV_COMMON_H

#include "bench_common.h" // HE_BENCH_CAPTURE

#include <helib/helib.h>

#include <memory>

struct Params
{
  const long m, p, r, qbits;
  const std::vector<long> gens;
  const std::vector<long> ords;
  const std::vector<long> mvec;
  Params(long _m,
         long _p,
         long _r,
         long _qbits,
         const std::vector<long>& _gens = {},
         const std::vector<long>& _ords = {},
         const std::vector<long>& _mvec = {}) :
      m(_m), p(_p), r(_r), qbits(_qbits), gens(_gens), ords(_ords), mvec(_mvec)
  {}
  Params(const Params& other) :
      Params(other.m,
             other.p,
             other.r,
             other.qbits,
             other.gens,
             other.ords,
             other.mvec)
  {}
  bool operator!=(Params& other) const { return !(*this == other); }
  bool operator==(Params& other) const
  {
    return m == other.m && p == other.p && r == other.r &&
           qbits == other.qbits && gens == other.gens && ords == other.ords &&
           mvec == other.mvec;
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
      context(helib::ContextBuilder<helib::BGV>()
                  .m(params.m)
                  .p(params.p)
                  .r(params.r)
                  .bits(params.qbits)
                  .gens(params.gens)
                  .ords(params.ords)
                  .bootstrappable(!params.mvec.empty())
                  .mvec(params.mvec)
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

#endif // HELIB_BGV_COMMON_H

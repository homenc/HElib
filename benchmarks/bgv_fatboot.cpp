/* Copyright (C) 2019-2020 IBM Corp.
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

#include <benchmark/benchmark.h>
#include <iostream>

#include <helib/helib.h>
#include <helib/debugging.h>

namespace {

void squareWithFatBoot(helib::PubKey& pk, helib::Ctxt& c)
{
  if (c.bitCapacity() <= 50) {
    pk.reCrypt(c);
  }
  c.square();
}

static void BM_fatboot(benchmark::State& state,
                       long m,
                       long p,
                       long r,
                       long c,
                       long bits,
                       long t,
                       int c_m,
                       std::vector<long> mvec,
                       std::vector<long> gens,
                       std::vector<long> ords)
{
  // clang-format off
  std::cout << "m=" << m
            << ", p=" << p
            << ", r=" << r
            << ", bits=" << bits
            << ", c=" << c
            << ", skHwt=" << t
            << ", c_m=" << c_m
            << ", mvec=" << helib::vecToStr(mvec)
            << ", gens=" << helib::vecToStr(gens)
            << ", ords=" << helib::vecToStr(ords)
            << std::endl;
  // clang-format on
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .gens(gens)
                               .ords(ords)
                               .bits(bits)
                               .c(c)
                               .bootstrappable(true)
                               .skHwt(t)
                               .mvec(mvec)
                               .thickboot()
                               .build();

  // Print the context
  context.printout();
  std::cout << std::endl;
  std::cout << "Security: " << context.securityLevel() << std::endl;

  std::cout << "Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  std::cout << "Generating key-switching matrices..." << std::endl;
  addSome1DMatrices(secret_key);
  addFrbMatrices(secret_key);

  // Generate bootstrapping data
  secret_key.genRecryptData();

  // NOTE: For some reason the reCrypt method is not marked const so
  //       I had to remove the const from the public key
  helib::PubKey& public_key = secret_key;
  const helib::EncryptedArray& ea = context.getEA();

  long nslots = ea.size();
  std::cout << "Number of slots: " << nslots << std::endl;

  std::vector<long> ptxt(nslots);
  for (int i = 0; i < nslots; ++i) {
    ptxt[i] = std::rand() % 2; // Random 0s and 1s
  }

  helib::Ctxt ctxt(public_key);
  ea.encrypt(ctxt, public_key, ptxt);
  for (auto _ : state)
    squareWithFatBoot(public_key, ctxt);
  std::cout << "Multiplications performed = " << state.iterations()
            << std::endl;
}

BENCHMARK_CAPTURE(BM_fatboot,
                  tiny_params,
                  /*m =*/31 * 41,
                  /*p =*/2,
                  /*r =*/1,
                  /*c =*/2,
                  /*bits =*/580,
                  /*t =*/64,
                  /*c_m =*/100,
                  /*mvec =*/std::vector<long>{31, 41},
                  /*gens =*/std::vector<long>{1026, 249},
                  /*ords =*/std::vector<long>{30, -2})
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

BENCHMARK_CAPTURE(BM_fatboot,
                  small_params,
                  /*m =*/31775,
                  /*p =*/2,
                  /*r =*/1,
                  /*c =*/2,
                  /*bits =*/580,
                  /*t =*/64,
                  /*c_m =*/100,
                  /*mvec =*/std::vector<long>{41, 775},
                  /*gens =*/std::vector<long>{6976, 24806},
                  /*ords =*/std::vector<long>{40, 30})
    ->Unit(benchmark::kMillisecond)
    ->MinTime(200);

BENCHMARK_CAPTURE(BM_fatboot,
                  big_params,
                  /*m =*/35113,
                  /*p =*/2,
                  /*r =*/1,
                  /*c =*/2,
                  /*bits =*/580,
                  /*t =*/64,
                  /*c_m =*/100,
                  /*mvec =*/std::vector<long>{37, 949},
                  /*gens =*/std::vector<long>{16134, 8548},
                  /*ords =*/std::vector<long>{36, 24})
    ->Unit(benchmark::kMillisecond)
    ->MinTime(200);

} // namespace

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

#include <benchmark/benchmark.h>
#include <iostream>

#include <helib/helib.h>
#include <helib/debugging.h>
#include "bgv_common.h"

namespace {

static void benchContextBinaryIO(benchmark::State& state, Meta& meta)
{
  std::stringstream ss;

  for (auto _ : state) {
    meta.data->context.writeTo(ss);
    helib::Context newContext = helib::Context::readFrom(ss);
    ::benchmark::DoNotOptimize(newContext);
  }
}

static void benchContextJSONIO(benchmark::State& state, Meta& meta)
{
  std::stringstream ss;

  for (auto _ : state) {
    meta.data->context.writeToJSON(ss);
    helib::Context newContext = helib::Context::readFromJSON(ss);
    ::benchmark::DoNotOptimize(newContext);
  }
}

static void benchPublicKeyBinaryIO(benchmark::State& state, Meta& meta)
{
  std::stringstream ss;

  for (auto _ : state) {
    meta.data->publicKey.writeTo(ss);
    helib::PubKey newPublicKey =
        helib::PubKey::readFrom(ss, meta.data->context);
    ::benchmark::DoNotOptimize(newPublicKey);
  }
}

static void benchPublicKeyJSONIO(benchmark::State& state, Meta& meta)
{
  std::stringstream ss;

  for (auto _ : state) {
    meta.data->publicKey.writeToJSON(ss);
    helib::PubKey newPublicKey =
        helib::PubKey::readFromJSON(ss, meta.data->context);
    ::benchmark::DoNotOptimize(newPublicKey);
  }
}

static void benchSecretKeyBinaryIO(benchmark::State& state, Meta& meta)
{
  std::stringstream ss;

  for (auto _ : state) {
    meta.data->secretKey.writeTo(ss);
    helib::SecKey newSecretKey =
        helib::SecKey::readFrom(ss, meta.data->context);
    ::benchmark::DoNotOptimize(newSecretKey);
  }
}

static void benchSecretKeyJSONIO(benchmark::State& state, Meta& meta)
{
  std::stringstream ss;

  for (auto _ : state) {
    meta.data->secretKey.writeToJSON(ss);
    helib::SecKey newSecretKey =
        helib::SecKey::readFromJSON(ss, meta.data->context);
    ::benchmark::DoNotOptimize(newSecretKey);
  }
}

static void benchCiphertextBinaryIO(benchmark::State& state, Meta& meta)
{
  std::stringstream ss;
  helib::Ctxt ctxt(meta.data->publicKey);

  for (auto _ : state) {
    ctxt.writeTo(ss);
    helib::Ctxt newCtxt = helib::Ctxt::readFrom(ss, meta.data->publicKey);
    ::benchmark::DoNotOptimize(newCtxt);
  }
}

static void benchCiphertextJSONIO(benchmark::State& state, Meta& meta)
{
  std::stringstream ss;
  helib::Ctxt ctxt(meta.data->publicKey);

  for (auto _ : state) {
    ctxt.writeToJSON(ss);
    helib::Ctxt newCtxt = helib::Ctxt::readFromJSON(ss, meta.data->publicKey);
    ::benchmark::DoNotOptimize(newCtxt);
  }
}

Meta fn;
Params no_boot_params(/*m =*/45,
                      /*p =*/19,
                      /*r =*/1,
                      /*bits =*/30);
// Binary IO benchmarks
BENCHMARK_CAPTURE(benchContextBinaryIO, no_boot_params, fn(no_boot_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(benchPublicKeyBinaryIO, no_boot_params, fn(no_boot_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(benchSecretKeyBinaryIO, no_boot_params, fn(no_boot_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(benchCiphertextBinaryIO, no_boot_params, fn(no_boot_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

// JSON IO benchmarks
BENCHMARK_CAPTURE(benchContextJSONIO, no_boot_params, fn(no_boot_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(benchPublicKeyJSONIO, no_boot_params, fn(no_boot_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(benchSecretKeyJSONIO, no_boot_params, fn(no_boot_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(benchCiphertextJSONIO, no_boot_params, fn(no_boot_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

Params tiny_params(/*m =*/31 * 41,
                   /*p =*/2,
                   /*r =*/1,
                   /*bits =*/580,
                   /*gens =*/std::vector<long>{1026, 249},
                   /*ords =*/std::vector<long>{30, -2},
                   /*mvec =*/std::vector<long>{31, 41});
// Binary IO benchmarks
BENCHMARK_CAPTURE(benchContextBinaryIO, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(10);
BENCHMARK_CAPTURE(benchPublicKeyBinaryIO, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(10);
BENCHMARK_CAPTURE(benchSecretKeyBinaryIO, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(10);
BENCHMARK_CAPTURE(benchCiphertextBinaryIO, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

// JSON IO benchmarks
BENCHMARK_CAPTURE(benchContextJSONIO, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(10);
BENCHMARK_CAPTURE(benchPublicKeyJSONIO, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(10);
BENCHMARK_CAPTURE(benchSecretKeyJSONIO, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(10);
BENCHMARK_CAPTURE(benchCiphertextJSONIO, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

Params small_params(/*m =*/31775,
                    /*p =*/2,
                    /*r =*/1,
                    /*bits =*/580,
                    /*gens =*/std::vector<long>{6976, 24806},
                    /*ords =*/std::vector<long>{40, 30},
                    /*mvec =*/std::vector<long>{41, 775});
// Binary IO benchmarks
BENCHMARK_CAPTURE(benchContextBinaryIO, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchPublicKeyBinaryIO, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchSecretKeyBinaryIO, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchCiphertextBinaryIO, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

// JSON IO benchmarks
BENCHMARK_CAPTURE(benchContextJSONIO, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchPublicKeyJSONIO, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchSecretKeyJSONIO, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchCiphertextJSONIO, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

Params big_params(/*m =*/35113,
                  /*p =*/2,
                  /*r =*/1,
                  /*bits =*/580,
                  /*gens =*/std::vector<long>{16134, 8548},
                  /*ords =*/std::vector<long>{36, 24},
                  /*mvec =*/std::vector<long>{37, 949});
// Binary IO benchmarks
BENCHMARK_CAPTURE(benchContextBinaryIO, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchPublicKeyBinaryIO, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchSecretKeyBinaryIO, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchCiphertextBinaryIO, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

// JSON IO benchmarks
BENCHMARK_CAPTURE(benchContextJSONIO, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchPublicKeyJSONIO, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchSecretKeyJSONIO, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
BENCHMARK_CAPTURE(benchCiphertextJSONIO, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

} // namespace

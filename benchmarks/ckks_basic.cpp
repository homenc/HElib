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
#include <memory>

#include <helib/helib.h>
#include <helib/debugging.h>

#include "ckks_common.h"

namespace {

static void adding_two_ciphertexts(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::CKKS> ptxt1(meta.data->context);
  helib::Ptxt<helib::CKKS> ptxt2(meta.data->context);
  ptxt1.random();
  ptxt2.random();

  helib::Ctxt ctxt1(meta.data->publicKey);
  helib::Ctxt ctxt2(meta.data->publicKey);
  meta.data->publicKey.Encrypt(ctxt1, ptxt1);
  meta.data->publicKey.Encrypt(ctxt2, ptxt2);
  // Benchmark adding ciphertexts
  for (auto _ : state)
    ctxt1 += ctxt2;
  std::cout << "Additions performed = " << state.iterations() << std::endl;
}

static void multiplying_two_ciphertexts(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::CKKS> ptxt1(meta.data->context);
  helib::Ptxt<helib::CKKS> ptxt2(meta.data->context);
  ptxt1.random();
  ptxt2.random();

  helib::Ctxt ctxt1(meta.data->publicKey);
  helib::Ctxt ctxt2(meta.data->publicKey);
  meta.data->publicKey.Encrypt(ctxt1, ptxt1);
  meta.data->publicKey.Encrypt(ctxt2, ptxt2);
  // Benchmark multiplying ciphertexts
  for (auto _ : state) {
    ctxt1.multiplyBy(ctxt2);
    benchmark::DoNotOptimize(ctxt1);
  }
  std::cout << "Multiplications performed = " << state.iterations()
            << std::endl;
}

static void encrypting_ciphertexts(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::CKKS> ptxt(meta.data->context);
  ptxt.random();

  helib::Ctxt ctxt(meta.data->publicKey);

  // Benchmark encrypting ciphertexts
  for (auto _ : state) {
    meta.data->publicKey.Encrypt(ctxt, ptxt);
    benchmark::DoNotOptimize(ctxt);
  }
  std::cout << "Encryptions performed = " << state.iterations() << std::endl;
}

static void decrypting_ciphertexts(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::CKKS> ptxt(meta.data->context);
  ptxt.random();

  helib::Ctxt ctxt(meta.data->publicKey);
  meta.data->publicKey.Encrypt(ctxt, ptxt);
  helib::Ptxt<helib::CKKS> decrypted_result(meta.data->context);

  // Benchmark decrypting ciphertexts
  for (auto _ : state) {
    meta.data->secretKey.Decrypt(decrypted_result, ctxt);
    benchmark::DoNotOptimize(decrypted_result);
  }
  std::cout << "Decryptions performed = " << state.iterations() << std::endl;
}

static void multiply_and_add_two_ciphertexts(benchmark::State& state,
                                             Meta& meta)
{
  helib::Ptxt<helib::CKKS> ptxt1(meta.data->context);
  helib::Ptxt<helib::CKKS> ptxt2(meta.data->context);
  ptxt1.random();
  ptxt2.random();

  helib::Ctxt ctxt1(meta.data->publicKey);
  helib::Ctxt ctxt2(meta.data->publicKey);

  meta.data->publicKey.Encrypt(ctxt1, ptxt1);
  meta.data->publicKey.Encrypt(ctxt2, ptxt2);
  // Benchmark multiplying and adding ciphertexts
  for (auto _ : state) {
    ctxt1.multiplyBy(ctxt2);
    benchmark::DoNotOptimize(ctxt1 += ctxt2);
  }
  std::cout << "Operations performed = " << state.iterations() << std::endl;
}

Meta fn;
Params tiny_params(1024, 1, 5800);
BENCHMARK_CAPTURE(adding_two_ciphertexts, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200)
    ->Complexity();
BENCHMARK_CAPTURE(multiplying_two_ciphertexts, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(encrypting_ciphertexts, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(decrypting_ciphertexts, tiny_params, fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(multiply_and_add_two_ciphertexts,
                  tiny_params,
                  fn(tiny_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

Params small_params(16384, 1, 5800);
BENCHMARK_CAPTURE(adding_two_ciphertexts, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->MinTime(200);
BENCHMARK_CAPTURE(multiplying_two_ciphertexts, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->MinTime(200);
BENCHMARK_CAPTURE(encrypting_ciphertexts, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(decrypting_ciphertexts, small_params, fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(multiply_and_add_two_ciphertexts,
                  small_params,
                  fn(small_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

Params big_params(65536, 1, 5800);
BENCHMARK_CAPTURE(adding_two_ciphertexts, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->MinTime(200);
BENCHMARK_CAPTURE(multiplying_two_ciphertexts, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->MinTime(200);
BENCHMARK_CAPTURE(encrypting_ciphertexts, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(decrypting_ciphertexts, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);
BENCHMARK_CAPTURE(multiply_and_add_two_ciphertexts, big_params, fn(big_params))
    ->Unit(benchmark::kMillisecond)
    ->Iterations(200);

} // namespace

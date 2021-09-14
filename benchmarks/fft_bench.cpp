/* Copyright (C) 2021 Intel Corporation
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "bgv_common.h"

#include <NTL/BasicThreadPool.h>
#include <helib/helib.h>
#include "../src/PrimeGenerator.h" // Private header

#include <benchmark/benchmark.h>
#include <iostream>

namespace {

static void helib_fft_forward(benchmark::State& state, Meta& meta)
{
  NTL::SetNumThreads(1);

  long N = meta.data->ea.size();
  long m = meta.data->context.getM();
  auto zms = meta.data->context.getZMStar();

  helib::PrimeGenerator prime_generator(49, m);

  long q = prime_generator.next();
  helib::Cmodulus cmod(zms, q, 0);

  NTL::zz_pX poly(N, 1);

  NTL::vec_long transformed(NTL::INIT_SIZE, N);

  for (auto _ : state) {
    cmod.FFT(transformed, poly);
  }
}

static void helib_fft_inverse(benchmark::State& state, Meta& meta)
{
  NTL::SetNumThreads(1);

  long N = meta.data->ea.size();
  long m = meta.data->context.getM();
  auto zms = meta.data->context.getZMStar();

  helib::PrimeGenerator prime_generator(49, m);

  long q = prime_generator.next();
  helib::Cmodulus cmod(zms, q, 0);

  NTL::zz_pX inverse;
  inverse.SetLength(N);

  NTL::vec_long transformed(NTL::INIT_SIZE, N);
  for (long i = 0; i < N; ++i)
    transformed[i] = i;

  for (auto _ : state)
    cmod.iFFT(inverse, transformed);
}

Meta fn;
Params hexl_F4_params(/*m=*/16384, /*p=*/65537, /*r=*/1, /*qbits=*/5800);
HE_BENCH_CAPTURE(helib_fft_forward, hexl_F4_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(helib_fft_inverse, hexl_F4_params, fn); //->Iterations(200);

Params hexl_F3_params(/*m=*/16, /*p=*/257, /*r=*/1, /*qbits=*/5800);
HE_BENCH_CAPTURE(helib_fft_forward, hexl_F3_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(helib_fft_inverse, hexl_F3_params, fn); //->Iterations(200);

Params hexl_F3d2_params(/*m=*/512, /*p=*/257, /*r=*/1, /*qbits=*/5800);
HE_BENCH_CAPTURE(helib_fft_forward, hexl_F3d2_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(helib_fft_inverse, hexl_F3d2_params, fn); //->Iterations(200);

} // namespace

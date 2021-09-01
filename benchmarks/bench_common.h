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

#ifndef HELIB_BENCH_COMMON_H
#define HELIB_BENCH_COMMON_H

#define HE_BENCH_CAPTURE(adding_two_ciphertexts, tiny_params, fn)              \
  BENCHMARK_CAPTURE(adding_two_ciphertexts, tiny_params, fn(tiny_params))

#endif // HELIB_BENCH_COMMON_H

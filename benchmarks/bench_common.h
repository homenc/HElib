// Copyright (C) 2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#ifndef HELIB_BENCH_COMMON_H
#define HELIB_BENCH_COMMON_H

#define HE_BENCH_CAPTURE(adding_two_ciphertexts, tiny_params, fn)              \
  BENCHMARK_CAPTURE(adding_two_ciphertexts, tiny_params, fn(tiny_params))      


#endif // HELIB_BENCH_COMMON_H

#ifndef HELIB_BENCH_COMMON_H
#define HELIB_BENCH_COMMON_H

#define HE_BENCH_CAPTURE(adding_two_ciphertexts, tiny_params, fn)              \
  BENCHMARK_CAPTURE(adding_two_ciphertexts, tiny_params, fn(tiny_params))      \
      ->Unit(benchmark::kMillisecond)

#endif // HELIB_BENCH_COMMON_H

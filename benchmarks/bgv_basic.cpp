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

#include "bgv_common.h"

#include <NTL/BasicThreadPool.h>
#include <helib/helib.h>

#include <benchmark/benchmark.h>
#include <iostream>

namespace {

static void adding_two_ciphertexts(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::BGV> ptxt1(meta.data->context);
  helib::Ptxt<helib::BGV> ptxt2(meta.data->context);

  ptxt1.random();
  ptxt2.random();

  helib::Ctxt ctxt1(meta.data->publicKey);
  helib::Ctxt ctxt2(meta.data->publicKey);

  meta.data->publicKey.Encrypt(ctxt1, ptxt1);
  meta.data->publicKey.Encrypt(ctxt2, ptxt2);
  // Benchmark adding ciphertexts
  for (auto _ : state)
    ctxt1 += ctxt2;
}

static void subtracting_two_ciphertexts(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::BGV> ptxt1(meta.data->context);
  helib::Ptxt<helib::BGV> ptxt2(meta.data->context);

  ptxt1.random();
  ptxt2.random();

  helib::Ctxt ctxt1(meta.data->publicKey);
  helib::Ctxt ctxt2(meta.data->publicKey);

  meta.data->publicKey.Encrypt(ctxt1, ptxt1);
  meta.data->publicKey.Encrypt(ctxt2, ptxt2);
  // Benchmark adding ciphertexts
  for (auto _ : state)
    ctxt1 -= ctxt2;
}

static void negating_a_ciphertext(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::BGV> ptxt(meta.data->context);

  ptxt.random();

  helib::Ctxt ctxt(meta.data->publicKey);

  meta.data->publicKey.Encrypt(ctxt, ptxt);
  // Benchmark adding ciphertexts
  for (auto _ : state)
    ctxt.negate();
}

static void square_a_ciphertext(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::BGV> ptxt(meta.data->context);

  ptxt.random();

  helib::Ctxt ctxt(meta.data->publicKey);

  meta.data->publicKey.Encrypt(ctxt, ptxt);
  // Benchmark adding ciphertexts
  for (auto _ : state)
    ctxt.square();
}

static void multiplying_two_ciphertexts_no_relin(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::BGV> ptxt1(meta.data->context);
  helib::Ptxt<helib::BGV> ptxt2(meta.data->context);

  ptxt1.random();
  ptxt2.random();
  
  helib::Ctxt ctxt1(meta.data->publicKey);
  helib::Ctxt ctxt2(meta.data->publicKey);

  meta.data->publicKey.Encrypt(ctxt1, ptxt1);
  meta.data->publicKey.Encrypt(ctxt2, ptxt2);
  // Benchmark adding ciphertexts
  for (auto _ : state)
    ctxt1.multLowLvl(ctxt2);
}

static void multiplying_two_ciphertexts(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::BGV> ptxt1(meta.data->context);
  helib::Ptxt<helib::BGV> ptxt2(meta.data->context);

  ptxt1.random();
  ptxt2.random();

  helib::Ctxt ctxt1(meta.data->publicKey);
  helib::Ctxt ctxt2(meta.data->publicKey);

  meta.data->publicKey.Encrypt(ctxt1, ptxt1);
  meta.data->publicKey.Encrypt(ctxt2, ptxt2);
  // Benchmark adding ciphertexts
  for (auto _ : state)
    ctxt1.multiplyBy(ctxt2);
}

static void rotate_a_ciphertext_by1(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::BGV> ptxt(meta.data->context);

  ptxt.random();

  helib::Ctxt ctxt(meta.data->publicKey);

  meta.data->publicKey.Encrypt(ctxt, ptxt);
  // Benchmark adding ciphertexts
  for (auto _ : state)
    meta.data->ea.rotate(ctxt, 1);
}

static void encrypting_ciphertexts(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::BGV> ptxt(meta.data->context);
  ptxt.random();

  helib::Ctxt ctxt(meta.data->publicKey);

  // Benchmark encrypting ciphertexts
  for (auto _ : state)
    meta.data->publicKey.Encrypt(ctxt, ptxt);
}

static void decrypting_ciphertexts(benchmark::State& state, Meta& meta)
{
  helib::Ptxt<helib::BGV> ptxt(meta.data->context);
  ptxt.random();

  helib::Ctxt ctxt(meta.data->publicKey);

  meta.data->publicKey.Encrypt(ctxt, ptxt);
  helib::Ptxt<helib::BGV> decrypted_result(meta.data->context);

  // Benchmark decrypting ciphertexts
  for (auto _ : state)
    meta.data->secretKey.Decrypt(decrypted_result, ctxt);
}

Meta fn;
Params tiny_params(/*m=*/257, /*p=*/2, /*r=*/1, /*qbits=*/5800);
HE_BENCH_CAPTURE(adding_two_ciphertexts, tiny_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(subtracting_two_ciphertexts,
                 tiny_params,
                 fn);                                     //->Iterations(200);
HE_BENCH_CAPTURE(negating_a_ciphertext, tiny_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(square_a_ciphertext, tiny_params, fn);   //->Iterations(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts,
                 tiny_params,
                 fn);                                       //->Iterations(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts_no_relin,
                 tiny_params,
                 fn);                                       //->Iterations(200);
HE_BENCH_CAPTURE(rotate_a_ciphertext_by1, tiny_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(encrypting_ciphertexts, tiny_params, fn);  //->Iterations(200);
HE_BENCH_CAPTURE(decrypting_ciphertexts, tiny_params, fn);  //->Iterations(200);

Params small_params(/*m=*/8009, /*p=*/2, /*r=*/1, /*qbits=*/5800);
HE_BENCH_CAPTURE(adding_two_ciphertexts, small_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(subtracting_two_ciphertexts,
                 small_params,
                 fn);                                      //->Iterations(200);
HE_BENCH_CAPTURE(negating_a_ciphertext, small_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(square_a_ciphertext, small_params, fn);   //->MinTime(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts,
                 small_params,
                 fn);                                        //->MinTime(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts_no_relin,
                 small_params,
                 fn);                                        //->MinTime(200);
HE_BENCH_CAPTURE(rotate_a_ciphertext_by1, small_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(encrypting_ciphertexts, small_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(decrypting_ciphertexts, small_params, fn); //->Iterations(200);

Params big_params(/*m=*/32003, /*p=*/2, /*r=*/1, /*qbits=*/5800);
HE_BENCH_CAPTURE(adding_two_ciphertexts, big_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(subtracting_two_ciphertexts,
                 big_params,
                 fn);                                    //->Iterations(200);
HE_BENCH_CAPTURE(negating_a_ciphertext, big_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(square_a_ciphertext, big_params, fn);   //->MinTime(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts, big_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts_no_relin, big_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(rotate_a_ciphertext_by1, big_params, fn);     //->MinTime(200);
HE_BENCH_CAPTURE(encrypting_ciphertexts, big_params, fn);      //->MinTime(200);
HE_BENCH_CAPTURE(decrypting_ciphertexts, big_params, fn);      //->MinTime(200);

Params hexl_F4_params(/*m=*/32768, /*p=*/65537, /*r=*/1, /*qbits=*/5800);
HE_BENCH_CAPTURE(adding_two_ciphertexts, hexl_F4_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(subtracting_two_ciphertexts,
                 hexl_F4_params,
                 fn);                                    //->Iterations(200);
HE_BENCH_CAPTURE(negating_a_ciphertext, hexl_F4_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(square_a_ciphertext, hexl_F4_params, fn);   //->MinTime(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts, hexl_F4_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts_no_relin, hexl_F4_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(rotate_a_ciphertext_by1, hexl_F4_params, fn);     //->MinTime(200);
HE_BENCH_CAPTURE(encrypting_ciphertexts, hexl_F4_params, fn);      //->MinTime(200);
HE_BENCH_CAPTURE(decrypting_ciphertexts, hexl_F4_params, fn);      //->MinTime(200);


Params hexl_F3_params(/*m=*/16, /*p=*/257, /*r=*/1, /*qbits=*/5800);
HE_BENCH_CAPTURE(adding_two_ciphertexts, hexl_F3_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(subtracting_two_ciphertexts,
                 hexl_F3_params,
                 fn);                                    //->Iterations(200);
HE_BENCH_CAPTURE(negating_a_ciphertext, hexl_F3_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(square_a_ciphertext, hexl_F3_params, fn);   //->MinTime(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts, hexl_F3_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts_no_relin, hexl_F3_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(rotate_a_ciphertext_by1, hexl_F3_params, fn);     //->MinTime(200);
HE_BENCH_CAPTURE(encrypting_ciphertexts, hexl_F3_params, fn);      //->MinTime(200);
HE_BENCH_CAPTURE(decrypting_ciphertexts, hexl_F3_params, fn);      //->MinTime(200);

Params hexl_F3d2_params(/*m=*/512, /*p=*/257, /*r=*/1, /*qbits=*/5800);
HE_BENCH_CAPTURE(adding_two_ciphertexts, hexl_F3d2_params, fn); //->Iterations(200);
HE_BENCH_CAPTURE(subtracting_two_ciphertexts,
                 hexl_F3d2_params,
                 fn);                                    //->Iterations(200);
HE_BENCH_CAPTURE(negating_a_ciphertext, hexl_F3d2_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(square_a_ciphertext, hexl_F3d2_params, fn);   //->MinTime(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts, hexl_F3d2_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(multiplying_two_ciphertexts_no_relin, hexl_F3d2_params, fn); //->MinTime(200);
HE_BENCH_CAPTURE(rotate_a_ciphertext_by1, hexl_F3d2_params, fn);     //->MinTime(200);
HE_BENCH_CAPTURE(encrypting_ciphertexts, hexl_F3d2_params, fn);      //->MinTime(200);
HE_BENCH_CAPTURE(decrypting_ciphertexts, hexl_F3d2_params, fn);      //->MinTime(200);

} // namespace

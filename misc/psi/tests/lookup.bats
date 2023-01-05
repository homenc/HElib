#!/usr/bin/env bats

# Copyright (C) 2020 IBM Corp.
# This program is Licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. See accompanying LICENSE file.

# Copyright (C) 2022 Intel Corporation
# SPDX-License-Identifier: Apache-2.0

utils_dir="../../../../utils"
load "../../../utils/tests/std"

nslots=3
modulus=13
#nslots=13860
#modulus=257
#nslots=79872
#modulus=1278209
datadir="data_and_params"
keygen_datadir="keygen_data_and_params"
sk_keygen_datadir="sk_keygen_data_and_params"
data_prefix="../$datadir"
keygen_data_prefix="../$keygen_datadir"
sk_keygen_data_prefix="../$sk_keygen_datadir"

db_encoded="../db${nslots}.enc"
query_encoded="../query${nslots}.enc"
lookup="../../build/bin/lookup"

# Run once per test
function setup {
  mkdir -p $tmp_folder
  cd $tmp_folder
  print-info-location
}

# Run once per test
function teardown {
  cd -
  remove-test-directory "$tmp_folder"
}

@test "lookup 64 threads" {
  skip
  echo "lookup 64 threads" > README
  $lookup ${data_prefix}/${prefix_bgv}.pk $data_prefix/db.ctxt $data_prefix/query.ctxt result.ctxt -n=64

  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
#  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
#  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
#  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
#  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
#  diff "result.ptxt_or" "expected.mask_or"
#  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "lookup 32 threads" {
  skip
  echo "lookup 32 threads" > README
  $lookup ${data_prefix}/${prefix_bgv}.pk $data_prefix/db.ctxt $data_prefix/query.ctxt result.ctxt -n=32

  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "lookup 16 threads" {
  skip
  echo "lookup 16 threads" > README
  $lookup ${data_prefix}/${prefix_bgv}.pk $data_prefix/db.ctxt $data_prefix/query.ctxt result.ctxt -n=16

  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "lookup 8 threads" {
  skip
  echo "lookup 8 threads" > README
  $lookup ${data_prefix}/${prefix_bgv}.pk $data_prefix/db.ctxt $data_prefix/query.ctxt result.ctxt -n=8

  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "lookup 4 threads" {
  skip
  echo "lookup 4 threads" > README
  $lookup ${data_prefix}/${prefix_bgv}.pk $data_prefix/db.ctxt $data_prefix/query.ctxt result.ctxt -n=4

  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "lookup 2 threads" {
  skip
  echo "lookup 2 threads" > README
  $lookup ${data_prefix}/${prefix_bgv}.pk $data_prefix/db.ctxt $data_prefix/query.ctxt result.ctxt -n=2

  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "lookup 1 thread" {
  skip
  echo "lookup 1 thread" > README
  $lookup ${data_prefix}/${prefix_bgv}.pk $data_prefix/db.ctxt $data_prefix/query.ctxt result.ctxt -n=1

  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen lookup 64 threads" {
  skip
  echo "keygen lookup 64 threads" > README
  $lookup ${keygen_data_prefix}/${prefix_bgv}Eval.pk $keygen_data_prefix/db.ctxt $keygen_data_prefix/query.ctxt result.ctxt -n=64

  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
#  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
#  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
#  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
#  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
#  diff "result.ptxt_or" "expected.mask_or"
#  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen lookup 32 threads" {
  skip
  echo "keygen lookup 32 threads" > README
  $lookup ${keygen_data_prefix}/${prefix_bgv}Eval.pk $keygen_data_prefix/db.ctxt $keygen_data_prefix/query.ctxt result.ctxt -n=32

  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen lookup 16 threads" {
  skip
  echo "keygen lookup 16 threads" > README
  $lookup ${keygen_data_prefix}/${prefix_bgv}Eval.pk $keygen_data_prefix/db.ctxt $keygen_data_prefix/query.ctxt result.ctxt -n=16

  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen lookup 8 threads" {
  skip
  echo "keygen lookup 8 threads" > README
  $lookup ${keygen_data_prefix}/${prefix_bgv}Eval.pk $keygen_data_prefix/db.ctxt $keygen_data_prefix/query.ctxt result.ctxt -n=8

  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen lookup 4 threads" {
  skip
  echo "keygen lookup 4 threads" > README
  $lookup ${keygen_data_prefix}/${prefix_bgv}Eval.pk $keygen_data_prefix/db.ctxt $keygen_data_prefix/query.ctxt result.ctxt -n=4

  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen lookup 2 threads" {
  skip
  echo "keygen lookup 2 threads" > README
  $lookup ${keygen_data_prefix}/${prefix_bgv}Eval.pk $keygen_data_prefix/db.ctxt $keygen_data_prefix/query.ctxt result.ctxt -n=2

  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen lookup 1 thread" {
  skip
  echo "keygen lookup 1 thread" > README
  $lookup ${keygen_data_prefix}/${prefix_bgv}Eval.pk $keygen_data_prefix/db.ctxt $keygen_data_prefix/query.ctxt result.ctxt -n=1

  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or"
  $decrypt ${keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand"

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen, sk only lookup 64 threads" {
  skip
  echo "keygen, sk only lookup 64 threads" > README
  $lookup ${sk_keygen_data_prefix}/${prefix_bgv}Eval.pk $sk_keygen_data_prefix/db.ctxt $sk_keygen_data_prefix/query.ctxt result.ctxt -n=64

  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and" -s
#  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or" -s
#  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand" -s

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
#  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
#  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
#  diff "result.ptxt_or" "expected.mask_or"
# diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen, sk only lookup 32 threads" {
  skip
  echo "keygen, sk only lookup 32 threads" > README
  $lookup ${sk_keygen_data_prefix}/${prefix_bgv}Eval.pk $sk_keygen_data_prefix/db.ctxt $sk_keygen_data_prefix/query.ctxt result.ctxt -n=32

  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand" -s

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen, sk only lookup 16 threads" {
  skip
  echo "keygen, sk only lookup 16 threads" > README
  $lookup ${sk_keygen_data_prefix}/${prefix_bgv}Eval.pk $sk_keygen_data_prefix/db.ctxt $sk_keygen_data_prefix/query.ctxt result.ctxt -n=16

  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand" -s

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen, sk only lookup 8 threads" {
  skip
  echo "keygen, sk only lookup 8 threads" > README
  $lookup ${sk_keygen_data_prefix}/${prefix_bgv}Eval.pk $sk_keygen_data_prefix/db.ctxt $sk_keygen_data_prefix/query.ctxt result.ctxt -n=8

  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand" -s

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen, sk only lookup 4 threads" {
  skip
  echo "keygen, sk only lookup 4 threads" > README
  $lookup ${sk_keygen_data_prefix}/${prefix_bgv}Eval.pk $sk_keygen_data_prefix/db.ctxt $sk_keygen_data_prefix/query.ctxt result.ctxt -n=4

  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand" -s

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen, sk only lookup 2 threads" {
  skip
  echo "keygen, sk only lookup 2 threads" > README
  $lookup ${sk_keygen_data_prefix}/${prefix_bgv}Eval.pk $sk_keygen_data_prefix/db.ctxt $sk_keygen_data_prefix/query.ctxt result.ctxt -n=2

  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand" -s

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

@test "keygen, sk only lookup 1 threads" {
  skip
  echo "keygen, sk only lookup 1 threads" > README
  $lookup ${sk_keygen_data_prefix}/${prefix_bgv}Eval.pk $sk_keygen_data_prefix/db.ctxt $sk_keygen_data_prefix/query.ctxt result.ctxt -n=1

  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt -o "result.ptxt" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_and -o "result.ptxt_and" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_or -o "result.ptxt_or" -s
  $decrypt ${sk_keygen_data_prefix}/${prefix_bgv}.sk result.ctxt_expand -o "result.ptxt_expand" -s

  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test SAME > "expected.mask"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test AND > "expected.mask_and"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test OR > "expected.mask_or"
  ../gen-expected-mask.py ${query_encoded} ${db_encoded} --mod-p $modulus --test EXPAND > "expected.mask_expand"

  diff "result.ptxt" "expected.mask"
  diff "result.ptxt_and" "expected.mask_and"
  diff "result.ptxt_or" "expected.mask_or"
  diff "result.ptxt_expand" "expected.mask_expand"
}

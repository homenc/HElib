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

utils_dir="../../../../utils"
load "../../../utils/tests/std"

modulus=257
nslots=13860
#nslots=79872
#modulus=1278209
datadir="data_and_params"
data_prefix="../$datadir"
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

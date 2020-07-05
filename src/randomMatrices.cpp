/* Copyright (C) 2012-2020 IBM Corp.
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

/**
 * @file randomMatrices.cpp
 * @brief implementation of random matrices of various forms build functions,
 * used for testing
 */
#include <helib/randomMatrices.h>

namespace helib {

MatMul1D* buildRandomMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
  case PA_GF2_tag: {
    return new RandomMatrix<PA_GF2>(ea, dim);
  }
  case PA_zz_p_tag: {
    return new RandomMatrix<PA_zz_p>(ea, dim);
  }
  default:
    return 0;
  }
}

MatMul1D* buildRandomMultiMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
  case PA_GF2_tag: {
    return new RandomMultiMatrix<PA_GF2>(ea, dim);
  }
  case PA_zz_p_tag: {
    return new RandomMultiMatrix<PA_zz_p>(ea, dim);
  }
  default:
    return 0;
  }
}

//********************************

BlockMatMul1D* buildRandomBlockMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
  case PA_GF2_tag: {
    return new RandomBlockMatrix<PA_GF2>(ea, dim);
  }
  case PA_zz_p_tag: {
    return new RandomBlockMatrix<PA_zz_p>(ea, dim);
  }
  default:
    return 0;
  }
}

//********************************

BlockMatMul1D* buildRandomMultiBlockMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
  case PA_GF2_tag: {
    return new RandomMultiBlockMatrix<PA_GF2>(ea, dim);
  }
  case PA_zz_p_tag: {
    return new RandomMultiBlockMatrix<PA_zz_p>(ea, dim);
  }
  default:
    return 0;
  }
}

MatMulFull* buildRandomFullMatrix(const EncryptedArray& ea)
{
  switch (ea.getTag()) {
  case PA_GF2_tag: {
    return new RandomFullMatrix<PA_GF2>(ea);
  }
  case PA_zz_p_tag: {
    return new RandomFullMatrix<PA_zz_p>(ea);
  }
  default:
    return nullptr;
  }
}

BlockMatMulFull* buildRandomFullBlockMatrix(const EncryptedArray& ea)
{
  switch (ea.getTag()) {
  case PA_GF2_tag: {
    return new RandomFullBlockMatrix<PA_GF2>(ea);
  }
  case PA_zz_p_tag: {
    return new RandomFullBlockMatrix<PA_zz_p>(ea);
  }
  default:
    return nullptr;
  }
}

} // namespace helib

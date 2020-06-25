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
#ifndef HELIB_CTPTRS_H
#define HELIB_CTPTRS_H
/**
 * @file CtPtrs.h
 * @brief Unified interface for vector of pointers to ciphertexts
 **/
#include <initializer_list>
#include <helib/Ctxt.h>
#include <helib/PtrVector.h>
#include <helib/PtrMatrix.h>

namespace helib {

typedef PtrVector<Ctxt> CtPtrs;
// CtPtrs_VecCt(NTL::Vec<Ctxt>)
typedef PtrVector_VecT<Ctxt> CtPtrs_VecCt;
// CtPtrs_vectorCt(std::vector<Ctxt>)
typedef PtrVector_vectorT<Ctxt> CtPtrs_vectorCt;
// CtPtrs_VecPt(NTL::Vec<Ctxt*>)
typedef PtrVector_VecPt<Ctxt> CtPtrs_VecPt;
// CtPtrs_vectorPt(std::vector<Ctxt*>)
typedef PtrVector_vectorPt<Ctxt> CtPtrs_vectorPt;

// A slice of CtPtrs
typedef PtrVector_slice<Ctxt> CtPtrs_slice;

typedef PtrMatrix<Ctxt> CtPtrMat;
typedef PtrMatrix_Vec<Ctxt> CtPtrMat_VecCt;
typedef PtrMatrix_vector<Ctxt> CtPtrMat_vectorCt;
typedef PtrMatrix_ptVec<Ctxt> CtPtrMat_ptVecCt;
typedef PtrMatrix_ptvector<Ctxt> CtPtrMat_ptvectorCt;

// Use packed bootstrapping, so we can bootstrap all in just a few calls
void packedRecrypt(const CtPtrs& cPtrs,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea);

// recrypt all ctxt below level 'belowLvl'
void packedRecrypt(const CtPtrs& array, // vector of Ctxts
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea,
                   long belowLvl);

void packedRecrypt(const CtPtrMat& m, // matrix of Ctxts
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea,
                   long belowLvl = LONG_MAX);

// Find the lowest level among many ciphertexts
// FIXME: using bitCapacity isn't really the right thing.
// this could break some code
inline long findMinBitCapacity(const CtPtrs& v)
{
  long lvl = LONG_MAX;
  for (long i = 0; i < v.size(); i++)
    if (v.isSet(i) && !v[i]->isEmpty())
      lvl = std::min(lvl, v[i]->bitCapacity());
  return lvl;
}

inline long findMinBitCapacity(const CtPtrMat& m)
{
  long lvl = LONG_MAX;
  for (long i = 0; i < m.size(); i++)
    lvl = std::min(lvl, findMinBitCapacity(m[i]));
  return lvl;
}

inline long findMinBitCapacity(std::initializer_list<const CtPtrs*> list)
{
  long lvl = LONG_MAX;
  for (auto elem : list)
    lvl = std::min(lvl, findMinBitCapacity(*elem));
  return lvl;
}

void innerProduct(Ctxt& result, const CtPtrs& v1, const CtPtrs& v2);
inline Ctxt innerProduct(const CtPtrs& v1, const CtPtrs& v2)
{
  Ctxt ret(ZeroCtxtLike, *v1[0]);
  innerProduct(ret, v1, v2);
  return ret;
}

} // namespace helib

#endif // ifndef HELIB_CTPTRS_H

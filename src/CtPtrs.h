#ifndef _CTPTRS_H
#define _CTPTRS_H
/** CtxtPtrs.h: unified interface for vector of pointers to ciphertexts */
#include "Ctxt.h"
#include "PtrVector.h"
#include "PtrMatrix.h"

typedef PtrVector<Ctxt> CtPtrs;
typedef PtrVector_VecT<Ctxt> CtPtrs_VecCt;      // CtPtrs_VecCt(NTL::Vec<Ctxt>)
typedef PtrVector_vectorT<Ctxt> CtPtrs_vectorCt;//CtPtrs_vectorCt(std::vector<Ctxt>)
typedef PtrVector_VecPt<Ctxt> CtPtrs_VecPt;     // CtPtrs_VecPt(NTL::Vec<Ctxt*>)
typedef PtrVector_vectorPt<Ctxt> CtPtrs_vectorPt;//CtPtrs_vectorPt(std::vector<Ctxt*>)

typedef PtrVector_slice<Ctxt>  CtPtrs_slice;    // A slice of CtPtrs

typedef PtrMatrix<Ctxt> CtPtrMat;
typedef PtrMatrix_Vec<Ctxt> CtPtrMat_VecCt;
typedef PtrMatrix_vector<Ctxt> CtPtrMat_vectorCt;
typedef PtrMatrix_ptVec<Ctxt> CtPtrMat_ptVecCt;
typedef PtrMatrix_ptvector<Ctxt> CtPtrMat_ptvectorCt;

// Use packed bootstrapping, so we can bootstrap all in just a few calls
void packedRecrypt(const CtPtrs& cPtrs,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea);

// recrypt all ctxt at level < belowLvl
void packedRecrypt(const CtPtrs& array,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea, long belowLvl);
void packedRecrypt(const CtPtrMat& m,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea, long belowLvl=LONG_MAX);

// Find the lowest level among many ciphertexts
inline long findMinLevel(const CtPtrs& v)
{
  long lvl = LONG_MAX;
  for (long i=0; i<v.size(); i++)
    if (v.isSet(i) && !v[i]->isEmpty())
      lvl = std::min(lvl, v[i]->findBaseLevel());
  return lvl;
}
inline long findMinLevel(const CtPtrMat& m)
{
  long lvl = LONG_MAX;
  for (long i=0; i<m.size(); i++)
    lvl = std::min(lvl, findMinLevel(m[i]));
  return lvl;
}

#include <initializer_list>
inline long findMinLevel(std::initializer_list<const CtPtrs*> list)
{
  long lvl = LONG_MAX;
  for (auto elem : list)
    lvl = std::min(lvl, findMinLevel(*elem));
  return lvl;
}

#endif // _CTPTRS_H

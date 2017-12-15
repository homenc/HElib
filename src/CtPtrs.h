#ifndef _CTPTRS_H
#define _CTPTRS_H
/** CtxtPtrs.h: unified interface for vector of pointers to ciphertexts */
#include "Ctxt.h"
#include "PtrVector.h"

typedef PtrVector<Ctxt> CtPtrs;
typedef PtrVector_VecT<Ctxt> CtPtrs_VecCt;    // CtPtrs_VecCt(NTL::Vec<Ctxt>)
typedef PtrVector_VecPt<Ctxt> CtPtrs_VecPt;   // CtPtrs_VecPt(NTL::Vec<Ctxt*>)
typedef PtrVector_vectorT<Ctxt> CtPtrs_vectorCt;//CtPtrs_vectorCt(std::vector<Ctxt>)
typedef PtrVector_vectorPt<Ctxt> CtPtrs_vectorPt;//CtPtrs_vectorPt(std::vector<Ctxt*>)

typedef PtrVector_slice<Ctxt>  CtPtrs_slice;    // A slice of CtPtrs
#endif // _CTPTRS_H


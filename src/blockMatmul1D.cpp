/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
/* blockMatmul1D.cpp - Data-movement operations on arrays of slots
 */
#include <algorithm>
#include <NTL/BasicThreadPool.h>
#include "matmul.h"

// A useful utility function from matmul.cpp, returns true if cache[i] is zero
extern bool
getDataFromCache(CachedConstants& cache, long i,
		 CachedConstants::CacheTag tag, const FHEcontext& context,
		 NTL::ZZX*& zzxPtr, DoubleCRT*& dcrtPtr);

// A useful utility function from blockMatmul.cpp, converts cache formats
extern void CachedBlockMatrixConvert(CachedDCRTPtxtBlockMatrix& v, 
            const CachedPtxtBlockMatrix& w, const FHEcontext& context);


// Extract one "column" from a matrix that was built with buildLinPolyCoeffs
// Extract one "column" from a matrix that was built with buildLinPolyCoeffs
template<class RX, class EAtype>
static bool shiftedColumInDiag(NTL::ZZX& zpoly, long f,
			const std::vector< std::vector<RX> >& diag,
			const EAtype& ea)
{
  // extract f "column" and store it in cvev
  bool zero = true;
  long nslots = ea.size();
  vector<RX> cvec(nslots);
  for (long j = 0; j < nslots; j++) {
    cvec[j] = diag[j][f];
    if (!NTL::IsZero(cvec[j])) zero = false;
  }
  if (zero) { // all are zeros
    clear(zpoly);
    return true;
  }
  ea.encode(zpoly, cvec);

  if (f>0) {
    long p = ea.getContext().zMStar.getP();
    long m = ea.getContext().zMStar.getM();
    long d = ea.getDegree();
    long exp = PowerMod(mcMod(p, m), d-f, m); // apply inverse automorphism
    const auto& F = ea.getTab().getPhimXMod();
    RX rpoly1, rpoly2;

    conv(rpoly1, zpoly);
    plaintextAutomorph(rpoly2,rpoly1, exp, m, F);
    conv(zpoly, rpoly2);
  }
  return false;
}

// A class that implements the basic (sparse) 1D matrix-vector functions
template<class type>
class blockMat_mul1D_impl{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, Ctxt& ctxt,
    const PlaintextBlockMatrixBaseInterface& mat, long dim) 
  {
    const PAlgebra& zMStar = ea.getContext().zMStar;

    assert(&ea == &mat.getEA().getDerived(type()));
    assert(&ea.getContext() == &ctxt.getContext());
    assert(dim >= 0 && dim <=  LONG(zMStar.numOfGens()));

    long p = zMStar.getP(); 
    long m = zMStar.getM();
    const RXModulus& F = ea.getTab().getPhimXMod();

    // special case for the extra dimension
    bool special = (dim == LONG(zMStar.numOfGens()));
    long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator
    bool bad = !special && !zMStar.SameOrd(dim);

    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup the NTL modulus

    // Get the derived type
    const PlaintextBlockMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

    ctxt.cleanUp(); // not sure, but this may be a good idea

    long nslots = ea.size();
    long d = ea.getDegree();

    Vec< shared_ptr<Ctxt> > acc;
    acc.SetLength(d);
    for (long k = 0; k < d; k++)
      acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

    mat_R entry;
    entry.SetDims(d, d);

    vector<RX> entry1;
    entry1.resize(d);
    
    vector< vector<RX> > diag;
    diag.resize(D);
    for (long j = 0; j < D; j++) diag[j].resize(d);

    // Process the diagonals one at a time
    for (long i = 0; i < D; i++) { // process diagonal i
      bool zDiag = true; // is this a zero diagonal?
      long nzLast = -1;  // index of last non-zero entry

      // Process the entries in this diagonal one at a time
      for (long j = 0; j < D; j++) { // process entry j
        bool zEntry = mat1.get(entry, mcMod(j-i, D), j); // entry [i,j-i mod D]
        assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));
          // get(...) returns true if the entry is empty, false otherwise

        if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too

        if (!zEntry) {    // non-empty entry
          zDiag = false;  // mark diagonal as non-empty

          // clear entries between last nonzero entry and this one
          for (long jj = nzLast+1; jj < j; jj++) {
            for (long k = 0; k < d; k++)
              clear(diag[jj][k]);
          }
          nzLast = j;

          // recode entry as a vector of polynomials
          for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

          // compute the lin poly coeffs
          ea.buildLinPolyCoeffs(diag[j], entry1);
        }
      }
      if (zDiag) continue; // zero diagonal, continue

      // clear trailing zero entries    
      for (long jj = nzLast+1; jj < D; jj++) {
        for (long k = 0; k < d; k++)
          clear(diag[jj][k]);
      }

      // now diag[j] contains the lin poly coeffs

      vector<Ctxt> shCtxt;
      vector<RX> shMask;
      if (i == 0) {
        shCtxt.resize(1, ctxt);
        shMask.resize(1, conv<RX>(1));
      }
      else if (!bad) {
        shCtxt.resize(1, ctxt);
        shMask.resize(1, conv<RX>(1));
        ea.rotate1D(shCtxt[0], dim, i);
        shCtxt[0].cleanUp();
      }
      else {
        // we fold the masking constants into the linearized polynomial
        // constants to save a level. We lift some code out of rotate1D
        // to do this.

        shCtxt.resize(2, ctxt);
        shMask.resize(2);

        long val = PowerMod(zMStar.ZmStarGen(dim), i, m);
        long ival = PowerMod(zMStar.ZmStarGen(dim), i-D, m);
        const RX& mask = ea.getTab().getMaskTable()[dim][D-i];

        shCtxt[0].smartAutomorph(val); 
        shCtxt[0].cleanUp();

        shCtxt[1].smartAutomorph(ival); 
        shCtxt[1].cleanUp();

        plaintextAutomorph(shMask[0], 1 - mask, val, zMStar.getM(), F);
        plaintextAutomorph(shMask[1], mask, ival, zMStar.getM(), F);
      }
      
      RX cpoly1, cpoly2, cpoly3;
      ZZX cpoly;

      // apply the linearlized polynomial
      for (long k = 0; k < d; k++) {

        // compute the constant
        bool zConst = true;
        vector<RX> cvec;
        cvec.resize(nslots);
        for (long j = 0; j < nslots; j++) {
          cvec[j] = diag[ special ? 0 : zMStar.coordinate(dim, j) ][k];
          if (!IsZero(cvec[j])) zConst = false;
        }

        if (zConst) continue;

        ea.encode(cpoly, cvec);
        conv(cpoly1, cpoly);

        // apply inverse automorphism to constant
        plaintextAutomorph(cpoly2,cpoly1, PowerMod(mcMod(p, m), mcMod(-k,d), m), zMStar.getM(), F);

        for (long j = 0; j < LONG(shCtxt.size()); j++) {
          MulMod(cpoly3, cpoly2, shMask[j], F);
          conv(cpoly, cpoly3);
          Ctxt shCtxt1 = shCtxt[j];;
          shCtxt1.multByConstant(cpoly);
          *acc[k] += shCtxt1;
        }
      }
    }

    Ctxt res(ZeroCtxtLike, ctxt);

    for (long k = 0; k < d; k++) {
      acc[k]->frobeniusAutomorph(k);
      res += *acc[k];
    }

    ctxt = res;
  }
};


void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextBlockMatrixBaseInterface& mat, long dim) 
{
  FHE_TIMER_START;
  ea.dispatch<blockMat_mul1D_impl>(Fwd(ctxt), mat, dim);
}




// A class that implements the basic caching functionality for (sparse)
// 1D matrix-vector products
template<class type> class compBlockMat1D_impl{
public:
  PA_INJECT(type)

  // This "apply" function only computes the cached block matrix cmat
  static void apply(const EncryptedArrayDerived<type>& ea, 
    CachedPtxtBlockMatrix& cmat, const PlaintextBlockMatrixBaseInterface& mat, 
    long dim) 
  {
    const PAlgebra& zMStar = ea.getContext().zMStar;

    assert(&ea == &mat.getEA().getDerived(type()));
    assert(dim >= 0 && dim <=  LONG(zMStar.numOfGens()));

    long p = zMStar.getP(); 
    long m = zMStar.getM();
    const RXModulus& F = ea.getTab().getPhimXMod();

    // special case for the extra dimension
    bool special = (dim == LONG(zMStar.numOfGens()));
    long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator
    bool bad = !special && !zMStar.SameOrd(dim);

    long nslots = ea.size();
    long d = ea.getDegree();

    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup the NTL modulus

    // Get the derived type
    const PlaintextBlockMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

    mat_R entry;
    entry.SetDims(d, d);

    vector<RX> entry1;
    entry1.resize(d);
    
    if (bad)
      cmat.SetDims(D, 2*d);
    else
      cmat.SetDims(D, d);

    vector< vector<RX> > diag;
    diag.resize(D);
    for (long j = 0; j < D; j++) diag[j].resize(d);

    // Process the diagonals one at a time
    for (long i = 0; i < D; i++) { // process diagonal i
      bool zDiag = true; // is this a zero diagonal?
      long nzLast = -1;  // index of last non-zero entry

      // Process the entries in this diagonal one at a time
      for (long j = 0; j < D; j++) { // process entry j
        bool zEntry = mat1.get(entry, mcMod(j-i, D), j); // entry [i,j-i mod D]
        assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));
        // get(...) returns true if the entry is empty, false otherwise

        if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too

        if (!zEntry) {    // non-empty entry
          zDiag = false;  // mark diagonal as non-empty

          // clear entries between last nonzero entry and this one
          for (long jj = nzLast+1; jj < j; jj++) {
            for (long k = 0; k < d; k++)
              clear(diag[jj][k]);
          }
          nzLast = j;

          // recode entry as a vector of polynomials
          for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

          // compute the lin poly coeffs
          ea.buildLinPolyCoeffs(diag[j], entry1);
        }
      }
      if (zDiag) continue; // zero diagonal, continue

      // clear trailing zero entries    
      for (long jj = nzLast+1; jj < D; jj++) {
        for (long k = 0; k < d; k++)
          clear(diag[jj][k]);
      }

      // now diag[j] contains the lin poly coeffs

      vector<RX> shMask;
      if (i == 0) {
        shMask.resize(1, conv<RX>(1));
      }
      else if (!bad) {
        shMask.resize(1, conv<RX>(1));
      }
      else {
        // we fold the masking constants into the linearized polynomial
        // constants to save a level. We lift some code out of rotate1D
        // to do this.

        shMask.resize(2);

        long val = PowerMod(zMStar.ZmStarGen(dim), i, m);
        long ival = PowerMod(zMStar.ZmStarGen(dim), i-D, m);
        const RX& mask = ea.getTab().getMaskTable()[dim][D-i];

        plaintextAutomorph(shMask[0], 1 - mask, val, zMStar.getM(), F);
        plaintextAutomorph(shMask[1], mask, ival, zMStar.getM(), F);
      }
      
      RX cpoly1, cpoly2, cpoly3;
      ZZX cpoly;

      // apply the linearlized polynomial
      for (long k = 0; k < d; k++) {

        // compute the constant
        bool zConst = true;
        vector<RX> cvec;
        cvec.resize(nslots);
        for (long j = 0; j < nslots; j++) {
          cvec[j] = diag[ special ? 0 : zMStar.coordinate(dim, j) ][k];
          if (!IsZero(cvec[j])) zConst = false;
        }

        if (zConst) continue;

        ea.encode(cpoly, cvec);
        conv(cpoly1, cpoly);

        // apply inverse automorphism to constant
        plaintextAutomorph(cpoly2,cpoly1, PowerMod(mcMod(p, m), mcMod(-k,d), m), zMStar.getM(), F);

        for (long j = 0; j < LONG(shMask.size()); j++) {
          MulMod(cpoly3, cpoly2, shMask[j], F);
          conv(cpoly, cpoly3);
          cmat[i][k + d*j] = ZZXptr(new ZZX(cpoly));
        }
      }
    }
  }
};


void compMat1D(const EncryptedArray& ea, CachedPtxtBlockMatrix& cmat, 
	       const PlaintextBlockMatrixBaseInterface& mat, long dim) 
{
  FHE_TIMER_START;
  ea.dispatch<compBlockMat1D_impl>(Fwd(cmat), mat, dim);
}


void compMat1D(const EncryptedArray& ea, CachedDCRTPtxtBlockMatrix& cmat, 
	       const PlaintextBlockMatrixBaseInterface& mat, long dim)
{
  FHE_TIMER_START;
  CachedPtxtBlockMatrix zzxMat;
  compMat1D(ea, zzxMat, mat, dim);
  CachedBlockMatrixConvert(cmat, zzxMat, ea.getContext());
}



#ifdef FHE_BOOT_THREADS

template<class CachedMatrix>
static void blockMat_mul1D_tmpl(const EncryptedArray& ea, Ctxt& ctxt, 
				const CachedMatrix& cmat, long dim)
{
  FHE_TIMER_START;
  const FHEcontext& context = ctxt.getContext();
  const PAlgebra& zMStar = context.zMStar;
  assert(dim >= 0 && dim <=  LONG(zMStar.numOfGens()));

  long m = zMStar.getM();

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator
  bool bad = !special && !zMStar.SameOrd(dim);
  long d = ea.getDegree();

  Vec< shared_ptr<Ctxt> > acc;
  acc.SetLength(d);
  ctxt.cleanUp(); // not sure, but this may be a good idea
  for (long k = 0; k < d; k++)
    acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  FHE_NTIMER_START(blockMat1);
  vector< vector<Ctxt> > shCtxt;

  shCtxt.resize(D);

  // Process the diagonals one at a time
  NTL_EXEC_RANGE(D, first, last)
      for (long i = first; i < last; i++) { // process diagonal i
        if (i == 0) {
          shCtxt[i].resize(1, ctxt);
        }
        else if (!bad) {
          shCtxt[i].resize(1, ctxt);
          ea.rotate1D(shCtxt[i][0], dim, i);
          shCtxt[i][0].cleanUp();
        }
        else {
          // we fold the masking constants into the linearized polynomial
          // constants to save a level. We lift some code out of rotate1D
          // to do this.
    
          shCtxt[i].resize(2, ctxt);
    
          long val = PowerMod(zMStar.ZmStarGen(dim), i, m);
          long ival = PowerMod(zMStar.ZmStarGen(dim), i-D, m);
    
          shCtxt[i][0].smartAutomorph(val); 
          shCtxt[i][0].cleanUp();
    
          shCtxt[i][1].smartAutomorph(ival); 
          shCtxt[i][1].cleanUp();
        }
      }
  NTL_EXEC_RANGE_END

  FHE_NTIMER_STOP(blockMat1);


  FHE_NTIMER_START(blockMat3);

  NTL_EXEC_RANGE(d, first, last)
      for (long k = first; k < last; k++) {
        for (long i = 0; i < D; i++) {
          for (long j = 0; j < LONG(shCtxt[i].size()); j++) {
            if (!cmat[i][k + d*j]) continue; // zero constant
    
            Ctxt shCtxt1 = shCtxt[i][j];
            shCtxt1.multByConstant(*cmat[i][k + d*j]);
            *acc[k] += shCtxt1;
          }
        }
        acc[k]->frobeniusAutomorph(k);
      }
  NTL_EXEC_RANGE_END

  FHE_NTIMER_STOP(blockMat3);

  FHE_NTIMER_START(blockMat4);
  Ctxt res(ZeroCtxtLike, ctxt);
  for (long k = 0; k < d; k++) {
    res += *acc[k];
  }
  FHE_NTIMER_STOP(blockMat4);

  ctxt = res;
}


#else


template<class CachedMatrix>
static void blockMat_mul1D_tmpl(const EncryptedArray& ea, Ctxt& ctxt, 
				const CachedMatrix& cmat, long dim)
{
  FHE_TIMER_START;
  const FHEcontext& context = ctxt.getContext();
  const PAlgebra& zMStar = context.zMStar;
  assert(dim >= 0 && dim <=  LONG(zMStar.numOfGens()));

  long m = zMStar.getM();

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator
  bool bad = !special && !zMStar.SameOrd(dim);
  long d = ea.getDegree();

  Vec< shared_ptr<Ctxt> > acc;
  acc.SetLength(d);
  ctxt.cleanUp(); // not sure, but this may be a good idea
  for (long k = 0; k < d; k++)
    acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  FHE_NTIMER_START(blockMat1);
  // Process the diagonals one at a time
  for (long i = 0; i < D; i++) { // process diagonal i
    vector<Ctxt> shCtxt;
    if (i == 0) {
      shCtxt.resize(1, ctxt);
    }
    else if (!bad) {
      shCtxt.resize(1, ctxt);
      ea.rotate1D(shCtxt[0], dim, i);
      shCtxt[0].cleanUp();
    }
    else {
      // we fold the masking constants into the linearized polynomial
      // constants to save a level. We lift some code out of rotate1D
      // to do this.

      shCtxt.resize(2, ctxt);

      long val = PowerMod(zMStar.ZmStarGen(dim), i, m);
      long ival = PowerMod(zMStar.ZmStarGen(dim), i-D, m);

      shCtxt[0].smartAutomorph(val); 
      shCtxt[0].cleanUp();

      shCtxt[1].smartAutomorph(ival); 
      shCtxt[1].cleanUp();
    }

    // apply the linearlized polynomial
    for (long k = 0; k < d; k++) {
      for (long j = 0; j < LONG(shCtxt.size()); j++) {
        if (!cmat[i][k + d*j]) continue; // zero constant

        Ctxt shCtxt1 = shCtxt[j];
        shCtxt1.multByConstant(*cmat[i][k + d*j]);
        *acc[k] += shCtxt1;
      }
    }
  }
  FHE_NTIMER_STOP(blockMat1);

  FHE_NTIMER_START(blockMat2);
  Ctxt res(ZeroCtxtLike, ctxt);
  for (long k = 0; k < d; k++) {
    acc[k]->frobeniusAutomorph(k);
    res += *acc[k];
  }
  FHE_NTIMER_STOP(blockMat2);

  ctxt = res;
}
#endif

void mat_mul1D( const EncryptedArray& ea, Ctxt& ctxt, 
		const CachedPtxtBlockMatrix& cmat, long dim)
{ blockMat_mul1D_tmpl(ea, ctxt, cmat, dim); }


void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
	       const CachedDCRTPtxtBlockMatrix& cmat, long dim)
{ blockMat_mul1D_tmpl(ea, ctxt, cmat, dim); }




/*********************************************************************/
/* mat_multi1D: Similar to mat_mul1D but different submatrices have
 * different transformations. 
 */

// This implementation uses one set of procedures to handle all of the
// caching options (none/ZZX/DCRT), at some point we should migrate
// all the stuff above to the same format.


template<class type>
class mat_multi1D_block_impl {
public:
  PA_INJECT(type)

  // Process a single block diagonal with index idx along dimenssion dim,
  // making calls to the get(i,j,k) method and storing the result in the
  // diag matrix (# of vectors = extenssion-degree). If all the calls to
  // get(i,j,k) return empty entries then return true and do not touch
  // the diag vector (in that case we may need to keep in it the last
  // non-zero diagonal)
  static bool
  processDiagonal(std::vector< std::vector<RX> >& diag, long dim, long idx,
		  PlaintextMultiBlockMatrixInterface<type>& mats)
  {
    const EncryptedArray& ea = mats.getEA();
    const PAlgebra& zMStar = ea.getContext().zMStar;
    long extDeg = ea.getDegree();

    mat_R entry;
    entry.SetDims(extDeg, extDeg);

    std::vector<RX> entry1;
    entry1.resize(extDeg);
    
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry
    long D = (dim == ea.dimension())? 1 : ea.sizeOfDimension(dim);

    // Get the slots in this diagonal one at a time
    long blockIdx, rowIdx, colIdx;
    for (long j = 0; j < ea.size(); j++) { // process entry j
      if (dim == ea.dimension()) { // "special" last dimenssion of size 1
	rowIdx = colIdx = 0; blockIdx=j;
      } else {
	std::pair<long,long> idxes = zMStar.breakIndexByDim(j, dim);
	blockIdx = idxes.first;  // index into mats,
	colIdx = idxes.second;   // index along diemnssion dim
	rowIdx = mcMod(colIdx-idx,D);
      }
      // process entry j
      bool zEntry = mats.get(entry,rowIdx,colIdx,blockIdx);
      // entry [i,j-i mod D] in the block corresponding to blockIdx

      assert(zEntry ||
             (entry.NumRows() == extDeg && entry.NumCols() == extDeg));
      // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too

      if (!zEntry) {    // non-empty entry
	zDiag = false;  // mark diagonal as non-empty

	// clear entries between last nonzero entry and this one
	for (long jj = nzLast+1; jj < j; jj++) {
	  for (long k = 0; k < extDeg; k++)
	    clear(diag[jj][k]);
	}
	nzLast = j;

	// recode entry as a vector of polynomials
	for (long k = 0; k < extDeg; k++) conv(entry1[k], entry[k]);

	// compute the lin poly coeffs
	ea.getDerived(type()).buildLinPolyCoeffs(diag[j], entry1);
      }
    }    
    if (zDiag) return true; // zero diagonal, nothing to do

    // clear trailing zero entries
    for (long jj = nzLast+1; jj < ea.size(); jj++)
      for (long k = 0; k < extDeg; k++)
	clear(diag[jj][k]);

    return false; // a nonzero diagonal
  }

  /************************************************************
   * The transformation encoded in mats can be described as
   *
   * \sum_{e=0}^{D-1}\sum_{f=0}^{d-1} \lambda_{e,f}*\sigma^f(\rho^e(x))
   *
   * and currently diag contains the d coefficients \lambda_{i,f}
   * for f=0...d-1. But our optimization computes the d sums
   *
   *     y_f = \sum_{e=0}^{D-1} \sigma^{-f}(\lambda_{e,f})*\rho^e(x)
   *
   * and then gets the result as \sum_{f=0}^{d-1} \sigma_f(y_f).
   *
   * Hence we need to transform the constants \lambda_{e,f} in diag,
   * into \lambda'_{e,f} = \sigma^{-f}(\lambda_{e,f}). To do that, we
   * encode the \lambda_{i,f}'s in d polynomials q_0,...,q_{d-1} and
   * apply the automorphisms \sigma^{-f}(q_f) to get an encoding of
   * the \lambda'_{i,f}'s.
   **/

  static void apply(const EncryptedArrayDerived<type> &ea, Ctxt &ctxt, long dim,
                    const PlaintextBlockMatrixBaseInterface& mats,
		    CachedConstants::CacheTag tag)
  {
    assert(&ea == &mats.getEA().getDerived(type()));
    assert(&ea.getContext() == &ctxt.getContext());
    assert(dim >= 0 && dim <= ea.dimension());

    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup NTL modulus

    long nslots = ea.size();
    long d = ea.getDegree();

    // special case for the extra dimension
    bool special = (dim == ea.dimension());
    long D = special ? 1 : ea.sizeOfDimension(dim);
    bool bad = !special && !ea.nativeDimension(dim);

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long p = zMStar.getP(); 
    long m = zMStar.getM();

    // Get the derived type (DIRT: need to get rid of const modifier)
    PlaintextMultiBlockMatrixInterface<type>& mats1 =
      (PlaintextMultiBlockMatrixInterface<type>&)
      dynamic_cast< const PlaintextMultiBlockMatrixInterface<type>& >( mats );

    // It is assumed that a cache of the right size is ready to use
    CachedConstants& cache = mats1.getCache();
    long cacheEntrySize = ea.nativeDimension(dim)? d : (2*d);
    bool cacheAvailable = (cache.size() == D*cacheEntrySize);
    if (!cacheAvailable) {
      cache.clear(); // ensure that there is no stale cache
      if (tag == CachedConstants::tagZZX || tag == CachedConstants::tagDCRT)
	cache.resize(D*cacheEntrySize);
    }

    ctxt.cleanUp(); // not sure, but this may be a good idea
    Ctxt tmp(ZeroCtxtLike, ctxt);
    std::vector<Ctxt> acc(d, tmp);

    // diag is a scratch space for calculations
    std::vector< std::vector<RX> > diag;
    diag.resize(nslots);
    for (long j = 0; j < nslots; j++) diag[j].resize(d);

    // Process the diagonals one at a time

    tmp = ctxt;
    long lastShift = 0;
    for (long e = 0; e < D; e++) { // process diagonal e
      bool zDiag = true;
      long procDiag = false; // did we call processDiagonal

      // For each diagonal e, we update the d accumulators y_0,..,y_{d-1}
      // with y_f += \sigma^{-f}(\lambda_{e,f}) * \rho^e(x)

      std::vector<bool> zero(d, false);
      std::vector<NTL::ZZX> zpoly(d, NTL::ZZX::zero());
      std::vector<NTL::ZZX*> zzxPtr(d,NULL);
      std::vector<DoubleCRT*> dcrtPtr(d,NULL);

      for (long f=0; f<d; f++) {
	long i = e*d + f; // index into cache

	// See if data is available in cache
        if (cacheAvailable && !cache.isEmpty(i)) {
	  zero[f] = getDataFromCache(cache,i,tag,ctxt.getContext(),
				     zzxPtr[f],dcrtPtr[f]);
	  if (!zero[f]) zDiag = false;
	}
	else { // no cache, need to compute this constant
	  if (!procDiag) { // diag is not set yet
	    zDiag = processDiagonal(diag, dim, e, mats1);
	    procDiag = true; // mark diag as set
            if (zDiag && tag != CachedConstants::tagEmpty)  // update cache
		for (long ii=e*d; ii<(e+1)*d; ii++) cache.setZero(ii);
	  }

	  // extract f'th "column" from diag and encode it in zpoly
	  if (!zDiag) {
	    zero[f] = shiftedColumInDiag(zpoly[f], f, diag, ea);
	    if (!zero[f]) {
	      if (tag == CachedConstants::tagDCRT) // allocate a new DoubleCRT
		dcrtPtr[f] = new DoubleCRT(zpoly[f], ctxt.getContext());
	      else if (tag != CachedConstants::tagEmpty) // allocate a new ZZX
		zzxPtr[f] = new NTL::ZZX(zpoly[f]);
	      else                                   // just use temporary ZZX
		zzxPtr[f] = &(zpoly[f]);
	    }

	    // update the cache if needed
	    if (tag != CachedConstants::tagEmpty) {
	      if (zero[f])                 cache.setZero(i);
	      else if (dcrtPtr[f] != NULL) cache.setAt(i,dcrtPtr[f]);
	      else if (tag == CachedConstants::tagZZX && zzxPtr[f] != NULL)
		cache.setAt(i,zzxPtr[f]);
	    }
	  }
	} // end of "else" case (no cache)
      }
      // done preparing all the zero, zzxPtr, dcrtPtr variables

      if (zDiag) break; // nothing to do for this diagonal

      // Rotate the ciphertext to position corresponding to diagonal e.
      // The code below uses only rotate-by-one operations if this is
      // a good dimenssion and the matrix is dense, hence it may require
      // fewer key-switching matrices

      if (e > 0) { // rotate the ciphertext
	if (ea.nativeDimension(dim)) // rotate the previous version
	  ea.rotate1D(tmp, dim, e-lastShift);
	else {                       // rotate the original ciphertext
	  tmp = ctxt;
	  ea.rotate1D(tmp, dim, e);
	}
	lastShift = e;
      } // if (e>0)

      // Depending on zero, zzxPtr, dcrtPtr, update the accumulated sums
      for (long f=0; f<d; f++) if (!zero[f]) {
          Ctxt tmp1(tmp);
	  if (dcrtPtr[f] != NULL) tmp1.multByConstant(*(dcrtPtr[f]));
	  else                    tmp1.multByConstant(*(zzxPtr[f]));
	  acc[f] += tmp1;
	}
	// The implementation above incurs an extra mult-by-constant due
	// to the masks in rotate1D when applied in a "bad dimension".
	// These masks can be folded into the constants here, but then
	// we would need to store two constants for each (e,f) rather
	// than one, namely const*mask and const*(1-mask).
	// We should implement that optimization at some point.

    } // end of this diagonal

    Ctxt res(ZeroCtxtLike, ctxt);

    // Finally, compute the result as \sum_{f=0}^{d-1} \sigma_f(y_f)
    ctxt = acc[0];
    for (long f = 1; f < d; f++) {
      acc[f].frobeniusAutomorph(f);
      ctxt += acc[f];
    }
  } // end of apply(...)
}; // end of class mat_multi1D_block_impl


//! @brief Multiply ctx by plaintext matrix over the base field/ring
void mat_multi1D_block(Ctxt& ctxt, const EncryptedArray& ea, long dim,
		       const PlaintextBlockMatrixBaseInterface& mats,
		       CachedConstants::CacheTag tag) 
{
  FHE_TIMER_START;
  ea.dispatch<mat_multi1D_block_impl>(Fwd(ctxt), dim, mats, tag);
}




// A useful helper function to get information from cache
bool
getDataFromCache(CachedConstants& cache, long i,
		 CachedConstants::CacheTag tag, const FHEcontext& context,
		 NTL::ZZX*& zzxPtr, DoubleCRT*& dcrtPtr)
{
  if (cache.isZero(i)) return true; // zero constant

  if (cache.isDCRT(i)) dcrtPtr = cache.getDCRT(i);
  else if (cache.isZZX(i)) {
    zzxPtr = cache.getZZX(i);
    if (tag == CachedConstants::tagDCRT) { // upgrade cache to DoubleCRT
      // DIRT: this "upgrade" logic may not be thread safe
      dcrtPtr = new DoubleCRT(*zzxPtr, context);
      cache.setAt(i,dcrtPtr);
      zzxPtr = NULL;
    }
  }
  else throw std::logic_error("cached constant is NULL");
  return false;
}

// helper routines

void CachedMatrixConvert(CachedDCRTPtxtMatrix& v, 
			 const CachedPtxtMatrix& w, const FHEcontext& context)
{
  long n = w.length();
  v.SetLength(n);
  for (long i = 0; i < n; i++)
    if (w[i]) v[i] = DCRTptr(new DoubleCRT(*w[i], context));
    // DoubleCRT defined relative to all primes, even the "special" ones
}


// Applying matmul to plaintext, useful for debugging
template<class type>
class mat_mul_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const PlaintextMatrixBaseInterface& mat)
  {
    PA_BOILER

    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    vector<RX> res;
    res.resize(n);
    for (long j = 0; j < n; j++) {
      RX acc, val, tmp; 
      acc = 0;
      for (long i = 0; i < n; i++) {
        if (!mat1.get(val, i, j)) {
          NTL::mul(tmp, data[i], val);
          NTL::add(acc, acc, tmp);
        }
      }
      rem(acc, acc, G);
      res[j] = acc;
    }

    data = res;
  }
}; 
void mat_mul(const EncryptedArray& ea, NewPlaintextArray& pa, 
	     const PlaintextMatrixBaseInterface& mat)
{
  ea.dispatch<mat_mul_pa_impl>(Fwd(pa), mat); 
}

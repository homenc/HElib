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
/* matmul.cpp - Implementation of homomorphic linear transforms
 */
#include <algorithm>
#include <NTL/ZZ.h>
NTL_CLIENT
#include "EncryptedArray.h"
//#include "timing.h"
//#include "cloned_ptr.h"
//#include "NumbTh.h"

// plaintextAutomorph: an auxilliary routine...maybe palce in NumbTh?
// result is calclated in the output x "in place", so a should
// not alias x

template <class RX, class RXModulus>
void plaintextAutomorph(RX& x, const RX& a, long k, const PAlgebra& zMStar, 
                       const RXModulus& F)
{
  long m  = zMStar.getM();

  assert(zMStar.inZmStar(k));

  x.SetLength(m);
  for (long j = 0; j < m; j++) x[j] = 0;

  long d = deg(a);
  
  mulmod_precon_t precon = PrepMulModPrecon(k, m, 1/(double)m);
  for (long j = 0; j <= d; j++) 
    x[MulModPrecon(j, k, m, precon)] = a[j];

  x.normalize();

  rem(x, x, F);
}

// A recursive matrix-by-vector multiply routine, used by the dense matrix
// code. This routine is optimized to use only the rotate1D routine rather
// than the more expensive linear-array rotations.
template<class type>
void EncryptedArrayDerived<type>::
  rec_mul(long dim, 
          Ctxt& res, 
          const Ctxt& pdata, 
          const vector<long>& idx,
          const PlaintextMatrixInterface<type>& mat,
          const vector<long>& dimx) const
{
  long ndims = dimension();
  long nslots = size();

  if (dim >= ndims) {
    vector<RX> pmat;
    pmat.resize(nslots);
    for (long j = 0; j < nslots; j++) {
      long i = idx[j];
      RX val;
      if (!mat.get(val, i, j))
        pmat[j] = val;
      else
        clear(pmat[j]);
    }

    ZZX epmat;
    encode(epmat, pmat);

    Ctxt tmp = pdata;
    tmp.multByConstant(epmat);
    res += tmp;
  }
  else {
    long sdim = sizeOfDimension(dimx[dim]);

    for (long offset = 0; offset < sdim; offset++) {
      Ctxt pdata1 = pdata;
      vector<long> idx1;
      rotate1D(pdata1, dimx[dim], offset);
      this->EncryptedArrayBase::rotate1D(idx1, idx, dimx[dim], offset);
      rec_mul(dim+1, res, pdata1, idx1, mat, dimx);
    }
  }
}


// helper class to sort dimensions, so that
//    - bad dimensions come before good dimensions (primary sort key)
//    - small dimensions come before large dimesnions (secondary sort key)
// this is a good order to process the dimensions in the recursive mat_mul_dense
// routine: it ensures that the work done at the work done at the
// leaves of the recursion is minimized, and that the work done
// at the non-leaves is dominated by the work done at the leaves.

template<class type>
struct MatMulDimComp {
  const EncryptedArrayDerived<type> *ea;
  MatMulDimComp(const EncryptedArrayDerived<type> *_ea) : ea(_ea) {}

  bool operator()(long i, long j) { 
    return (!ea->nativeDimension(i) && ea->nativeDimension(j)) ||
           (  (ea->nativeDimension(i) == ea->nativeDimension(j)) &&
              (ea->sizeOfDimension(i) < ea->sizeOfDimension(j))  );
  }
};


template<class type>
void EncryptedArrayDerived<type>::mat_mul_dense(Ctxt& ctxt, const PlaintextMatrixBaseInterface& mat) const
{
  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());

  RBak bak; bak.save(); tab.restoreContext();

  const PlaintextMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace()); // fresh encryption of zero
  vector<long> idx;
  idx.resize(size());
  for (long i = 0; i < size(); i++)
     idx[i] = i;

  vector<long> dimx;
  dimx.resize(dimension());
  for (long i = 0; i < dimension(); i++)
    dimx[i] = i;

  sort(dimx.begin(), dimx.end(), MatMulDimComp<type>(this));
  // sort the dimenesions so that bad ones come before good,
  // and then small ones come before large

  rec_mul(0, res, ctxt, idx, mat1, dimx);

  ctxt = res;
}


// this mat_mul is optimized for diagonally sparse matrices

template<class type>
void EncryptedArrayDerived<type>::mat_mul(Ctxt& ctxt, const PlaintextMatrixBaseInterface& mat) const
{
  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());

  RBak bak; bak.save(); tab.restoreContext();

  const PlaintextMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace());
  // a new ciphertext, encrypting zero
  

  long nslots = size();
  long d = getDegree();

  RX entry;
  vector<RX> diag;
  diag.resize(nslots);

  for (long i = 0; i < nslots; i++) {
    // process diagonal i


    bool zDiag = true;
    long nzLast = -1;

    for (long j = 0; j < nslots; j++) {
      bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j);
      assert(zEntry || deg(entry) < d);

      if (!zEntry && IsZero(entry)) zEntry = true;

      if (!zEntry) {

        zDiag = false;

        // clear entries between last nonzero entry and this one

        for (long jj = nzLast+1; jj < j; jj++) clear(diag[jj]);
        nzLast = j;

        diag[j] = entry;
      }
    }

    
    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries
    for (long jj = nzLast+1; jj < nslots; jj++) clear(diag[jj]);
    
    Ctxt shCtxt = ctxt;
    rotate(shCtxt, i); 

    ZZX cpoly;
    encode(cpoly, diag);

    shCtxt.multByConstant(cpoly);
    res += shCtxt;
  }

  ctxt = res;
}

template<class type> void EncryptedArrayDerived<type>::
mat_mul1D(Ctxt& ctxt, const PlaintextMatrixBaseInterface& mat, long dim) const
{
  FHE_TIMER_START;

  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());
  
  const PAlgebra& zMStar = context.zMStar;

  assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim);

  RBak bak; bak.save(); tab.restoreContext();

  const PlaintextMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  Ctxt res(ZeroCtxtLike, ctxt);

  long nslots = size();
  long d = getDegree();

  RX entry;
  vector<RX> diag, diag1;
  diag.resize(D);
  diag1.resize(nslots);
  ZZX cpoly;
  

  for (long i = 0; i < D; i++) {
    // process diagonal i

    bool zDiag = true;
    long nzLast = -1;

    for (long j = 0; j < D; j++) {
      bool zEntry = mat1.get(entry, mcMod(j-i, D), j);
      assert(zEntry || deg(entry) < d);

      if (!zEntry && IsZero(entry)) zEntry = true;

      if (!zEntry) {

        zDiag = false;

        // clear entries between last nonzero entry and this one

        for (long jj = nzLast+1; jj < j; jj++) clear(diag[jj]);
        nzLast = j;

        diag[j] = entry;
      }
    }

    
    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries
    for (long jj = nzLast+1; jj < D; jj++) clear(diag[jj]);
    
    Ctxt shCtxt = ctxt;
    if (i != 0) rotate1D(shCtxt, dim, i); 

    
    for (long j = 0; j < nslots; j++)
      diag1[j] = diag[ special ? 0 : zMStar.coordinate(dim, j) ];

    encode(cpoly, diag1);

    shCtxt.multByConstant(cpoly);
    res += shCtxt;
  }

  ctxt = res;
}


// This code has a complexity of N+d (instead of N*d) where N is the number of
// nonzero diagonal blocks. However, it requires space for d extra ciphertexts
template<class type>
void EncryptedArrayDerived<type>::mat_mul(Ctxt& ctxt, const PlaintextBlockMatrixBaseInterface& mat) const
{
  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());

  const PAlgebra& zMStar = context.zMStar;
  long p = zMStar.getP(); 
  long m = zMStar.getM();
  const RXModulus& F = tab.getPhimXMod();

  RBak bak; bak.save(); tab.restoreContext();

  const PlaintextBlockMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  long nslots = size();
  long d = getDegree();

  Vec< shared_ptr<Ctxt> > acc;
  acc.SetLength(d);
  for (long k = 0; k < d; k++)
    acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  mat_R entry;
  entry.SetDims(d, d);

  vector<RX> entry1;
  entry1.resize(d);
  
  vector< vector<RX> > diag;
  diag.resize(nslots);
  for (long j = 0; j < nslots; j++) diag[j].resize(d);

  for (long i = 0; i < nslots; i++) {
    // process diagonal i


    bool zDiag = true;
    long nzLast = -1;

    for (long j = 0; j < nslots; j++) {
      bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j);
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
        buildLinPolyCoeffs(diag[j], entry1);
      }
    }

    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries    
    for (long jj = nzLast+1; jj < nslots; jj++) {
      for (long k = 0; k < d; k++)
        clear(diag[jj][k]);
    }

    // now diag[j] contains the lin poly coeffs

    Ctxt shCtxt = ctxt;
    rotate(shCtxt, i); 
    shCtxt.cleanUp();

    RX cpoly1, cpoly2;
    ZZX cpoly;

    // apply the linearlized polynomial
    for (long k = 0; k < d; k++) {

      // compute the constant
      bool zConst = true;
      vector<RX> cvec;
      cvec.resize(nslots);
      for (long j = 0; j < nslots; j++) {
        cvec[j] = diag[j][k];
        if (!IsZero(cvec[j])) zConst = false;
      }

      if (zConst) continue;

      encode(cpoly, cvec);
      conv(cpoly1, cpoly);

      // apply inverse automorphism to constant
      plaintextAutomorph(cpoly2, cpoly1, PowerMod(p, mcMod(-k, d), m), zMStar, F);
      conv(cpoly, cpoly2);
      Ctxt shCtxt1 = shCtxt;
      shCtxt1.multByConstant(cpoly);
      *acc[k] += shCtxt1;
    }
  }

  Ctxt res(ZeroCtxtLike, ctxt);

  for (long k = 0; k < d; k++) {
    acc[k]->frobeniusAutomorph(k);
    res += *acc[k];
  }

  ctxt = res;
}


template<class type> void EncryptedArrayDerived<type>::
mat_mul1D(Ctxt& ctxt, const PlaintextBlockMatrixBaseInterface& mat, long dim) const
{
  FHE_TIMER_START;

  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());

  const PAlgebra& zMStar = context.zMStar;
  long p = zMStar.getP(); 
  long m = zMStar.getM();
  const RXModulus& F = tab.getPhimXMod();

  assert(dim >= 0 && dim <=  LONG(zMStar.numOfGens()));

  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim);

  RBak bak; bak.save(); tab.restoreContext();

  const PlaintextBlockMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  long nslots = size();
  long d = getDegree();

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

  for (long i = 0; i < D; i++) {
    // process diagonal i


    bool zDiag = true;
    long nzLast = -1;

    for (long j = 0; j < D; j++) {
      bool zEntry = mat1.get(entry, mcMod(j-i, D), j);
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
        buildLinPolyCoeffs(diag[j], entry1);
      }
    }

    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries    
    for (long jj = nzLast+1; jj < D; jj++) {
      for (long k = 0; k < d; k++)
        clear(diag[jj][k]);
    }

    // now diag[j] contains the lin poly coeffs

    Ctxt shCtxt = ctxt;
    if (i != 0) rotate1D(shCtxt, dim, i); 
    shCtxt.cleanUp();

    RX cpoly1, cpoly2;
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

      encode(cpoly, cvec);
      conv(cpoly1, cpoly);

      // apply inverse automorphism to constant
      plaintextAutomorph(cpoly2, cpoly1, PowerMod(p, mcMod(-k, d), m), zMStar, F);
      conv(cpoly, cpoly2);
      Ctxt shCtxt1 = shCtxt;
      shCtxt1.multByConstant(cpoly);
      *acc[k] += shCtxt1;
    }
  }

  Ctxt res(ZeroCtxtLike, ctxt);

  for (long k = 0; k < d; k++) {
    acc[k]->frobeniusAutomorph(k);
    res += *acc[k];
  }

  ctxt = res;
}


// Linearized polynomials.
// L describes a linear map M by describing its action on the standard
// power basis: M(x^j mod G) = (L[j] mod G), for j = 0..d-1.  
// The result is a coefficient vector C for the linearized polynomial
// representing M: a polynoamial h in Z/(p^r)[X] of degree < d is sent to
//
//    M(h(X) \bmod G)= \sum_{i=0}^{d-1}(C[j] \cdot h(X^{p^j}))\bmod G).
template<class type> void
EncryptedArrayDerived<type>::buildLinPolyCoeffs(vector<ZZX>& C, 
						const vector<ZZX>& L) const
{
  RBak bak; bak.save(); restoreContext();
  vector<RX> CC, LL;
  convert(LL, L);
  buildLinPolyCoeffs(CC, LL);
  convert(C, CC);
}

template<class type> void
EncryptedArrayDerived<type>::buildLinPolyCoeffs(vector<RX>& C, 
						const vector<RX>& L) const
{
  FHE_TIMER_START;

  RBak bak; bak.save(); restoreContext();
  REBak ebak; ebak.save(); restoreContextForG();

  if (!linPolyMatrix) {
    FHE_NTIMER_START(buildLinPolyCoeffs_buildMatrix);

    long p = tab.getZMStar().getP();
    long r = tab.getR();

    Mat<RE> M1;
    buildLinPolyMatrix(M1, p);
    Mat<RE> M2;
    ppInvert(M2, M1, p, r);

    linPolyMatrix = shared_ptr< Mat<RE> >(new Mat<RE>(M2) );
  }

  Vec<RE> CC, LL;
  convert(LL, L);
  mul(CC, LL, *linPolyMatrix);
  convert(C, CC);
}


// Apply the same linear transformation to all the slots.
// C is the output of ea.buildLinPolyCoeffs
void applyLinPoly1(const EncryptedArray& ea, Ctxt& ctxt, const vector<ZZX>& C)
{
  assert(&ea.getContext() == &ctxt.getContext());
  long d = ea.getDegree();
  assert(d == lsize(C));

  long nslots = ea.size();

  vector<ZZX> encodedC(d);
  for (long j = 0; j < d; j++) {
    vector<ZZX> v(nslots);
    for (long i = 0; i < nslots; i++) v[i] = C[j];
    ea.encode(encodedC[j], v);
  }

  applyLinPolyLL(ea, ctxt, encodedC);
}


// Apply different transformations to different slots. Cvec is a vector of
// length ea.size(), with each entry the output of ea.buildLinPolyCoeffs
void applyLinPolyMany(const EncryptedArray& ea, Ctxt& ctxt, 
                      const vector< vector<ZZX> >& Cvec)
{
  assert(&ea.getContext() == &ctxt.getContext());
  long d = ea.getDegree();
  long nslots = ea.size();

  assert(nslots == lsize(Cvec));
  for (long i = 0; i < nslots; i++)
    assert(d == lsize(Cvec[i]));

  vector<ZZX> encodedC(d);
  for (long j = 0; j < d; j++) {
    vector<ZZX> v(nslots);
    for (long i = 0; i < nslots; i++) v[i] = Cvec[i][j];
    ea.encode(encodedC[j], v);
  }

  applyLinPolyLL(ea, ctxt, encodedC);
}

// A low-level variant: encodedCoeffs has all the linPoly coeffs encoded
// in slots; different transformations can be encoded in different slots
void applyLinPolyLL(const EncryptedArray& ea, 
                    Ctxt& ctxt, const vector<ZZX>& encodedC)
{
  long d = ea.getDegree();
  assert(d == lsize(encodedC));

  ctxt.cleanUp();  // not sure, but this may be a good idea

  Ctxt tmp(ctxt);

  ctxt.multByConstant(encodedC[0]);
  for (long j = 1; j < d; j++) {
    Ctxt tmp1(tmp);
    tmp1.frobeniusAutomorph(j);
    tmp1.multByConstant(encodedC[j]);
    ctxt += tmp1;
  }
}

// Explicit instantiation

template class EncryptedArrayDerived<PA_GF2>;
template class EncryptedArrayDerived<PA_zz_p>;

#if 0
/************************* OLD UNUSED CODE *************************/
template<class type>
void EncryptedArrayDerived<type>::mat_mul(Ctxt& ctxt, const PlaintextBlockMatrixBaseInterface& mat) const
{
  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());

  RBak bak; bak.save(); tab.restoreContext();

  const PlaintextBlockMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace());
  // a new ciphertext, encrypting zero
  

  long nslots = size();
  long d = getDegree();

  mat_R entry;
  entry.SetDims(d, d);

  vector<RX> entry1;
  entry1.resize(d);
  
  vector< vector<RX> > diag;
  diag.resize(nslots);
  for (long j = 0; j < nslots; j++) diag[j].resize(d);

  for (long i = 0; i < nslots; i++) {
    // process diagonal i


    bool zDiag = true;
    long nzLast = -1;

    for (long j = 0; j < nslots; j++) {
      bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j);
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
        buildLinPolyCoeffs(diag[j], entry1);
      }
    }

    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries    
    for (long jj = nzLast+1; jj < nslots; jj++) {
      for (long k = 0; k < d; k++)
        clear(diag[jj][k]);
    }

    // now diag[j] contains the lin poly coeffs

    Ctxt shCtxt = ctxt;
    rotate(shCtxt, i); 

    // apply the linearlized polynomial
    for (long k = 0; k < d; k++) {

      // compute the constant
      bool zConst = true;
      vector<RX> cvec;
      cvec.resize(nslots);
      for (long j = 0; j < nslots; j++) {
        cvec[j] = diag[j][k];
        if (!IsZero(cvec[j])) zConst = false;
      }

      if (zConst) continue;

      ZZX cpoly;
      encode(cpoly, cvec);
      // FIXME: record the encoded polynomial for future use

      Ctxt shCtxt1 = shCtxt;
      shCtxt1.frobeniusAutomorph(k);
      shCtxt1.multByConstant(cpoly);
      res += shCtxt1;
    }
  }
  ctxt = res;
}
#endif

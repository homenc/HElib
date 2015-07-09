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
/* matrix.cpp - Data-movement operations on arrays of slots
 */

#include "matrix.h"
#include <algorithm>

// Copied from EncryptedArray.cpp -- it shoudn't be there.

// plaintextAutomorph: an auxilliary routine...maybe palce in NumbTh?
// Compute b(X) = a(X^k) mod Phi_m(X). Result is calclated in the output b
// "in place", so a should not alias b.
template <class RX, class RXModulus> static
void plaintextAutomorph(RX& b, const RX& a, long k, const PAlgebra& zMStar, 
  const RXModulus& PhimX)
{
  long m  = zMStar.getM();

  assert(zMStar.inZmStar(k));

  b.SetLength(m);
  for (long j = 0; j < m; j++) b[j] = 0;

  long d = deg(a);

  // compute b(X) = a(X^k) mod (X^m-1)
  mulmod_precon_t precon = PrepMulModPrecon(k, m);
  for (long j = 0; j <= d; j++) 
    b[MulModPrecon(j, k, m, precon)] = a[j]; // b[j*k mod m] = a[j]
  b.normalize();

  rem(b, b, PhimX); // reduce modulo the m'th cyclotomic
}



// This code has a complexity of N+d (instead of N*d) where N is the number of
// nonzero diagonal blocks. However, it requires space for d extra ciphertexts
template<class type>
class mat_mul_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, 
    Ctxt& ctxt, const PlaintextMatrixBaseInterface& mat) 
  {
    assert(&ea == &mat.getEA().getDerived(type()));
    assert(&ea.getContext() == &ctxt.getContext());

    RBak bak; bak.save(); ea.getTab().restoreContext();

    // Get the derived type
    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    ctxt.cleanUp(); // not sure, but this may be a good idea

    Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace()); // fresh encryption of zero

    long nslots = ea.size();
    long d = ea.getDegree();

    RX entry;
    vector<RX> diag;
    diag.resize(nslots);

    // Process the diagonals one at a time
    for (long i = 0; i < nslots; i++) {  // process diagonal i
      bool zDiag = true; // is this a zero diagonal?
      long nzLast = -1;  // index of last non-zero entry on this diagonal

      // Compute constants for each entry on this diagonal
      for (long j = 0; j < nslots; j++) { // process entry j
        bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j); // callback
        assert(zEntry || deg(entry) < d);

        if (!zEntry && IsZero(entry)) zEntry = true; // check for zero

        if (!zEntry) { // non-zero diagonal entry

          zDiag = false; // diagonal is non-zero

          // clear entries between last nonzero entry and this one
          for (long jj = nzLast+1; jj < j; jj++) clear(diag[jj]);
          nzLast = j;

          diag[j] = entry;
        }
      }
      
      if (zDiag) continue; // zero diagonal, continue

      // clear trailing zero entries
      for (long jj = nzLast+1; jj < nslots; jj++) clear(diag[jj]);

      // Now we have the constants for all the diagonal entries, encode the
      // diagonal as a single polynomial with these constants in the slots
      ZZX cpoly;
      ea.encode(cpoly, diag);

      // rotate by i, multiply by the polynomial, then add to the result
      Ctxt shCtxt = ctxt;
      ea.rotate(shCtxt, i); // rotate by i
      shCtxt.multByConstant(cpoly);
      res += shCtxt;
    }
    ctxt = res;
  }


  static void apply(const EncryptedArrayDerived<type>& ea, Ctxt& ctxt, 
    const PlaintextBlockMatrixBaseInterface& mat) 
  {
    assert(&ea == &mat.getEA().getDerived(type()));
    assert(&ea.getContext() == &ctxt.getContext());

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long p = zMStar.getP(); 
    long m = zMStar.getM();
    const RXModulus& F = ea.getTab().getPhimXMod();

    RBak bak; bak.save(); ea.getTab().restoreContext();

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
    diag.resize(nslots);
    for (long j = 0; j < nslots; j++) diag[j].resize(d);

    for (long i = 0; i < nslots; i++) { // process diagonal i
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
          ea.buildLinPolyCoeffs(diag[j], entry1);
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
      ea.rotate(shCtxt, i); 
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

        ea.encode(cpoly, cvec);
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



};


void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextBlockMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  ea.dispatch<mat_mul_impl>(Fwd(ctxt), mat);
}

void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  ea.dispatch<mat_mul_impl>(Fwd(ctxt), mat);
}






template<class type>
class mat_mul_dense_impl {
public:
  PA_INJECT(type)

  // A recursive matrix-by-vector multiply, used by the dense matrix code.
  // This routine is optimized to use only the rotate1D routine rather
  // than the more expensive linear-array rotations.
  static void rec_mul(const EncryptedArrayDerived<type>& ea,
    long dim, Ctxt& res, const Ctxt& pdata, const vector<long>& idx,
    const PlaintextMatrixInterface<type>& mat,
    const vector<long>& dimx) 
  {
    long ndims = ea.dimension();
    long nslots = ea.size();

    if (dim >= ndims) { // Last dimension (recursion edge condition)
      vector<RX> pmat;  // the plaintext diagonal
      pmat.resize(nslots);
      bool zDiag = true; // is this a zero diagonal
      for (long j = 0; j < nslots; j++) {
        long i = idx[j];
        RX val;
        if (mat.get(val, i, j)) // returns true if the entry is zero
          clear(pmat[j]);
        else {           // not a zero entry
          pmat[j] = val;
          zDiag = false; // not a zero diagonal
        }
      }
      if (zDiag) return; // zero diagonal, nothing to do

      // Now we have the constants for all the diagonal entries, encode the
      // diagonal as a single polynomial with these constants in the slots
      ZZX epmat;
      ea.encode(epmat, pmat);

      // multiply by the polynomial, then add to the result
      Ctxt tmp = pdata;
      tmp.multByConstant(epmat);
      res += tmp;
    }
    else { // not the last dimension, make a recursive call
      long sdim = ea.sizeOfDimension(dimx[dim]);

      // compute "in spirit" sum_i (pdata >> i) * i'th-diagonal, but
      // adjust the indexes so that we only need to rotate the cipehrtext
      // along the different dimensions separately
      for (long offset = 0; offset < sdim; offset++) {
        Ctxt pdata1 = pdata;
        vector<long> idx1;
        ea.rotate1D(pdata1, dimx[dim], offset);
        ea.EncryptedArrayBase::rotate1D(idx1, idx, dimx[dim], offset);
        // indexes adjusted, make the recursive call
        rec_mul(ea, dim+1, res, pdata1, idx1, mat, dimx);
      }
    }
  }


  // helper class to sort dimensions, so that
  //    - bad dimensions come before good dimensions (primary sort key)
  //    - small dimensions come before large dimesnions (secondary sort key)
  // this is a good order to process the dimensions in the recursive mat_mul_dense
  // routine: it ensures that the work done at the leaves of the recursion is
  // minimized, and that the work done at the non-leaves is dominated by the
  // work done at the leaves.

  struct MatMulDimComp {
    const EncryptedArrayDerived<type> *ea;
    MatMulDimComp(const EncryptedArrayDerived<type> *_ea) : ea(_ea) {}

    bool operator()(long i, long j) { 
      return (!ea->nativeDimension(i) && ea->nativeDimension(j)) ||
             (  (ea->nativeDimension(i) == ea->nativeDimension(j)) &&
                (ea->sizeOfDimension(i) < ea->sizeOfDimension(j))  );
    }
  };

  // Multiply a ciphertext vector by a plaintext dense matrix
  static void apply(const EncryptedArrayDerived<type>& ea, Ctxt& ctxt, 
    const PlaintextMatrixBaseInterface& mat) 
  {
    assert(&ea == &mat.getEA().getDerived(type()));
    assert(&ea.getContext() == &ctxt.getContext());

    RBak bak; bak.save(); ea.getTab().restoreContext();

    // Get the derived type
    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    ctxt.cleanUp(); // not sure, but this may be a good idea

    Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace()); // fresh encryption of zero
    vector<long> idx;
    idx.resize(ea.size());
    for (long i = 0; i < ea.size(); i++)
       idx[i] = i;

    vector<long> dimx;
    dimx.resize(ea.dimension());
    for (long i = 0; i < ea.dimension(); i++)
      dimx[i] = i;

    sort(dimx.begin(), dimx.end(), MatMulDimComp(&ea));
    // sort the dimenesions so that bad ones come before good,
    // and then small ones come before large

    // call the recursive procedure to do the actual work
    rec_mul(ea, 0, res, ctxt, idx, mat1, dimx);

    ctxt = res; // copy the result back to ctxt
  }


};



void mat_mul_dense(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  ea.dispatch<mat_mul_dense_impl>(Fwd(ctxt), mat);
}



template<class type>
class compMat_impl{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, 
    CachedPtxtMatrix& cmat, const PlaintextMatrixBaseInterface& mat) 
  {
    assert(&ea == &mat.getEA().getDerived(type()));

    RBak bak; bak.save(); ea.getTab().restoreContext();

    // Get the derived type
    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    long nslots = ea.size();
    long d = ea.getDegree();
    RX entry;
    vector<RX> diag;
    diag.resize(nslots);
    cmat.SetLength(nslots);

    // Process the diagonals one at a time
    for (long i = 0; i < nslots; i++) {  // process diagonal i
      bool zDiag = true; // is this a zero diagonal?
      long nzLast = -1;  // index of last non-zero entry on this diagonal

      // Compute constants for each entry on this diagonal
      for (long j = 0; j < nslots; j++) { // process entry j
        bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j); // callback
        assert(zEntry || deg(entry) < d);

        if (!zEntry && IsZero(entry)) zEntry = true; // check for zero

        if (!zEntry) { // non-zero diagonal entry

          zDiag = false; // diagonal is non-zero

          // clear entries between last nonzero entry and this one
          for (long jj = nzLast+1; jj < j; jj++) clear(diag[jj]);
          nzLast = j;

          diag[j] = entry;
        }
      }
      
      if (zDiag) continue; // zero diagonal, continue

      // clear trailing zero entries
      for (long jj = nzLast+1; jj < nslots; jj++) clear(diag[jj]);

      // Now we have the constants for all the diagonal entries, encode the
      // diagonal as a single polynomial with these constants in the slots
      ZZX cpoly;
      ea.encode(cpoly, diag);
      cmat[i] = ZZXptr(new ZZX(cpoly));
    }
  }


  static void apply(const EncryptedArrayDerived<type>& ea, 
    CachedPtxtBlockMatrix& cmat, const PlaintextBlockMatrixBaseInterface& mat) 
  {
    assert(&ea == &mat.getEA().getDerived(type()));
    const PAlgebra& zMStar = ea.getContext().zMStar;
    long p = zMStar.getP(); 
    long m = zMStar.getM();
    const RXModulus& F = ea.getTab().getPhimXMod();

    RBak bak; bak.save(); ea.getTab().restoreContext();

    // Get the derived type
    const PlaintextBlockMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

    long nslots = ea.size();
    long d = ea.getDegree();

    mat_R entry;
    entry.SetDims(d, d);
    vector<RX> entry1;
    entry1.resize(d);
    
    vector< vector<RX> > diag;
    diag.resize(nslots);
    for (long j = 0; j < nslots; j++) diag[j].resize(d);
    cmat.SetDims(nslots, d);

    for (long i = 0; i < nslots; i++) { // process diagonal i
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
          ea.buildLinPolyCoeffs(diag[j], entry1);
        }
      }
      if (zDiag) continue; // zero diagonal, continue

      // clear trailing zero entries    
      for (long jj = nzLast+1; jj < nslots; jj++) {
        for (long k = 0; k < d; k++)
          clear(diag[jj][k]);
      }
      // now diag[j] contains the lin poly coeffs

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

        ea.encode(cpoly, cvec);
        conv(cpoly1, cpoly);

        // apply inverse automorphism to constant
        plaintextAutomorph(cpoly2, cpoly1, PowerMod(p, mcMod(-k, d), m), zMStar, F);
        conv(cpoly, cpoly2);
        cmat[i][k] = ZZXptr(new ZZX(cpoly));
      }
    }
  }




};


// helper routines


static void CachedMatrixConvert(CachedDCRTPtxtMatrix& v, 
  const CachedPtxtMatrix& w, const FHEcontext& context)
{
  long n = w.length();
  v.SetLength(n);
  for (long i = 0; i < n; i++)
    if (w[i]) v[i] = DCRTptr(new DoubleCRT(*w[i], context));
    // DoubleCRT defined relative to all primes, even the "special" ones
}

static void CachedBlockMatrixConvert(CachedDCRTPtxtBlockMatrix& v, 
  const CachedPtxtBlockMatrix& w, const FHEcontext& context)
{
  long m = w.NumRows();
  long n = w.NumCols();

  v.SetDims(m, n);
  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++)
    if (w[i][j]) v[i][j] = DCRTptr(new DoubleCRT(*w[i][j], context));
    // DoubleCRT defined relative to all primes, even the "special" ones
    // NOTE: Fixes potential bug in original version
}




void compMat(const EncryptedArray& ea, CachedPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  ea.dispatch<compMat_impl>(Fwd(cmat), mat);
}




void compMat(const EncryptedArray& ea, CachedDCRTPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  CachedPtxtMatrix zzxMat;
  compMat(ea, zzxMat, mat);
  CachedMatrixConvert(cmat, zzxMat, ea.getContext());
}



void compMat(const EncryptedArray& ea, CachedPtxtBlockMatrix& cmat, 
  const PlaintextBlockMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  ea.dispatch<compMat_impl>(Fwd(cmat), mat);
}




void compMat(const EncryptedArray& ea, CachedDCRTPtxtBlockMatrix& cmat, 
  const PlaintextBlockMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  CachedPtxtBlockMatrix zzxMat;
  compMat(ea, zzxMat, mat);
  CachedBlockMatrixConvert(cmat, zzxMat, ea.getContext());
}





template<class type, class RX>
static bool processDiagonal(vector<RX>& diag, long dim, long i,
  vector<RX>& tmpDiag, const PlaintextMatrixInterface<type>& mat,
  const PAlgebra& zMStar, long d, bool special)
{
  long D = tmpDiag.size();
  long nslots = diag.size();
  bool zDiag = true; // is this a zero diagonal?
  long nzLast = -1;  // index of last non-zero entry
  RX entry;

  // Process the entries in this diagonal one at a time
  for (long j = 0; j < D; j++) { // process entry j
    bool zEntry = mat.get(entry, mcMod(j-i, D), j); // entry [i,j-i mod D]
    assert(zEntry || deg(entry) < d);
    // get(...) returns true if the entry is empty, false otherwise

    if (!zEntry && IsZero(entry)) zEntry = true; // zero is an empty entry too

    if (!zEntry) {   // not a zero entry
      zDiag = false; // mark diagonal as non-empty

      // clear entries between last nonzero entry and this one
      for (long jj = nzLast+1; jj < j; jj++) clear(tmpDiag[jj]);
      nzLast = j;

      tmpDiag[j] = entry;
    }
  }    
  if (zDiag) return true; // zero diagonal, nothing to do

  // clear trailing zero entries
  for (long jj = nzLast+1; jj < D; jj++) clear(tmpDiag[jj]);

  if (special) diag.assign(nslots, tmpDiag[0]); // order-1 dimension
  else for (long j = 0; j < nslots; j++)
    diag[j] = tmpDiag[ zMStar.coordinate(dim,j) ];
    // rearrange the indexes based on the current dimension

  return false; // a nonzero diagonal
}




template<class type>
class mat_mul1D_impl{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, Ctxt& ctxt,
    const PlaintextMatrixBaseInterface& mat, long dim) 
  {
    const PAlgebra& zMStar = ea.getContext().zMStar;

    assert(&ea == &mat.getEA().getDerived(type()));
    assert(&ea.getContext() == &ctxt.getContext());
    assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

    // special case for the extra dimension
    bool special = (dim == LONG(zMStar.numOfGens()));
    long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup the NTL modulus

    // Get the derived type
    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    ctxt.cleanUp();  // not sure, but this may be a good idea
    Ctxt res(ZeroCtxtLike, ctxt); // fresh encryption of zero

    long nslots = ea.size();
    long d = ea.getDegree();
    vector<RX> diag, diag1;
    diag.resize(D);
    diag1.resize(nslots);
    ZZX cpoly;

    // Process the diagonals one at a time
    for (long i = 0; i < D; i++) { // process diagonal i
      if (processDiagonal(diag1, dim, i, diag, mat1, zMStar, d, special))
        continue; // zero diagonal

      // encode as a polynomial, then multiply and add
      ea.encode(cpoly, diag1);
      Ctxt shCtxt = ctxt;
      if (i != 0) ea.rotate1D(shCtxt, dim, i);   
      shCtxt.multByConstant(cpoly);
      res += shCtxt;
    }

    ctxt = res;
  }




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

        plaintextAutomorph(shMask[0], 1 - mask, val, zMStar, F);
        plaintextAutomorph(shMask[1], mask, ival, zMStar, F);
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
        plaintextAutomorph(cpoly2,cpoly1, PowerMod(p, mcMod(-k,d), m), zMStar, F);

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
  const PlaintextMatrixBaseInterface& mat, long dim) 
{
  FHE_TIMER_START;
  ea.dispatch<mat_mul1D_impl>(Fwd(ctxt), mat, dim);
}



void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextBlockMatrixBaseInterface& mat, long dim) 
{
  FHE_TIMER_START;
  ea.dispatch<mat_mul1D_impl>(Fwd(ctxt), mat, dim);
}






template<class type>
class compMat1D_impl{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, 
    CachedPtxtMatrix& cmat, const PlaintextMatrixBaseInterface& mat, long dim) 
  {
    const PAlgebra& zMStar = ea.getContext().zMStar;

    assert(&ea == &mat.getEA().getDerived(type()));
    assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

    // special case for the extra dimension
    bool special = (dim == LONG(zMStar.numOfGens()));
    long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup the NTL modulus

    // Get the derived type
    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    long nslots = ea.size();
    long d = ea.getDegree();
    vector<RX> diag, diag1;
    diag.resize(D);
    diag1.resize(nslots);
    cmat.SetLength(D);
    ZZX cpoly;

    // Process the diagonals one at a time
    for (long i = 0; i < D; i++) { // process diagonal i
      if (processDiagonal(diag1, dim, i, diag, mat1, zMStar, d, special))
        continue; // zero diagonal

      // encode as a polynomial, then multiply and add
      ea.encode(cpoly, diag1);
      cmat[i] = ZZXptr(new ZZX(cpoly));
    }
  }



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

        plaintextAutomorph(shMask[0], 1 - mask, val, zMStar, F);
        plaintextAutomorph(shMask[1], mask, ival, zMStar, F);
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
        plaintextAutomorph(cpoly2,cpoly1, PowerMod(p, mcMod(-k,d), m), zMStar, F);

        for (long j = 0; j < LONG(shMask.size()); j++) {
          MulMod(cpoly3, cpoly2, shMask[j], F);
          conv(cpoly, cpoly3);
          cmat[i][k + d*j] = ZZXptr(new ZZX(cpoly));
        }
      }
    }
  }



};


void compMat1D(const EncryptedArray& ea, CachedPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat, long dim) 
{
  FHE_TIMER_START;
  ea.dispatch<compMat1D_impl>(Fwd(cmat), mat, dim);
}

void compMat1D(const EncryptedArray& ea, CachedDCRTPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat, long dim)
{
  FHE_TIMER_START;
  CachedPtxtMatrix zzxMat;
  compMat1D(ea, zzxMat, mat, dim);
  CachedMatrixConvert(cmat, zzxMat, ea.getContext());
}


void compMat1D(const EncryptedArray& ea, CachedPtxtBlockMatrix& cmat, 
  const PlaintextBlockMatrixBaseInterface& mat, long dim) 
{
  FHE_TIMER_START;
  ea.dispatch<compMat1D_impl>(Fwd(cmat), mat, dim);
}

void compMat1D(const EncryptedArray& ea, CachedDCRTPtxtBlockMatrix& cmat, 
  const PlaintextBlockMatrixBaseInterface& mat, long dim)
{
  FHE_TIMER_START;
  CachedPtxtBlockMatrix zzxMat;
  compMat1D(ea, zzxMat, mat, dim);
  CachedBlockMatrixConvert(cmat, zzxMat, ea.getContext());
}







template<class CachedMatrix>
void mat_mul_tmpl(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedMatrix& cmat)
{
  FHE_TIMER_START;
  ctxt.cleanUp(); // not sure, but this may be a good idea
  Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace()); // fresh encryption of zero

  // Process the diagonals one at a time
  long nslots = ea.size();
  for (long i = 0; i < nslots; i++) {  // process diagonal i
    if (!cmat[i]) continue; // a zero diagonal

    // rotate by i, multiply, and add to the result
    Ctxt shCtxt = ctxt;
    ea.rotate(shCtxt, i); // rotate by i
    shCtxt.multByConstant(*cmat[i]);
    res += shCtxt;
  }
  ctxt = res;
}
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtMatrix& cmat)
{
  mat_mul_tmpl(ea, ctxt, cmat);
}
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtMatrix& cmat)
{
  mat_mul_tmpl(ea, ctxt, cmat);
}




#ifdef FHE_BOOT_THREADS

template<class CachedMatrix>
static void mat_mul1D_tmpl(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedMatrix& cmat, long dim)
{
  FHE_TIMER_START;
  const FHEcontext& context = ctxt.getContext();
  const PAlgebra& zMStar = context.zMStar;
  assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

  ctxt.cleanUp();  // not sure, but this may be a good idea
  Ctxt res(ZeroCtxtLike, ctxt); // fresh encryption of zero

  Vec< shared_ptr<Ctxt> > tvec;
  tvec.SetLength(D);
  for (long i = 0; i < D; i++)
    tvec[i] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  bootTask->exec1(D,
    [&](long first, long last) {
      for (long i = first; i < last; i++) { // process diagonal i
        if (!cmat[i]) continue;      // zero diagonal
    
        // rotate and multiply
        (*tvec[i]) = ctxt;
        if (i != 0) ea.rotate1D(*tvec[i], dim, i);   
        tvec[i]->multByConstant(*cmat[i]);
      }
    }
  );

  for (long i = 0; i < D; i++)
    res += *tvec[i];

  ctxt = res;
}

#else

template<class CachedMatrix>
static void mat_mul1D_tmpl(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedMatrix& cmat, long dim)
{
  FHE_TIMER_START;
  const FHEcontext& context = ctxt.getContext();
  const PAlgebra& zMStar = context.zMStar;
  assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

  ctxt.cleanUp();  // not sure, but this may be a good idea
  Ctxt res(ZeroCtxtLike, ctxt); // fresh encryption of zero

  // Process the diagonals one at a time
  for (long i = 0; i < D; i++) { // process diagonal i
    if (!cmat[i]) continue;      // zero diagonal

    // rotate, multiply and add
    Ctxt shCtxt = ctxt;
    if (i != 0) ea.rotate1D(shCtxt, dim, i);   
    shCtxt.multByConstant(*cmat[i]);
    res += shCtxt;
  }
  ctxt = res;
}


#endif



void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtMatrix& cmat, long dim)
{ mat_mul1D_tmpl(ea, ctxt, cmat, dim); }

void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtMatrix& cmat, long dim)
{ mat_mul1D_tmpl(ea, ctxt, cmat, dim); }



template<class CachedMatrix>
void blockMat_mul_tmpl(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedMatrix& cmat)
{
  FHE_TIMER_START;
  ctxt.cleanUp(); // not sure, but this may be a good idea

  long nslots = ea.size();
  long d = ea.getDegree();

  Vec< shared_ptr<Ctxt> > acc;
  acc.SetLength(d);
  for (long k = 0; k < d; k++)
    acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  for (long i = 0; i < nslots; i++) { // process diagonal i
    // apply the linearlized polynomial
    Ctxt shCtxt = ctxt;
    ea.rotate(shCtxt, i); 
    shCtxt.cleanUp();

    for (long k = 0; k < d; k++) {
      if (!cmat[i][k]) continue; // a zero constant
      Ctxt shCtxt1 = shCtxt;
      shCtxt1.multByConstant(*cmat[i][k]);
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
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtBlockMatrix& cmat)
{
  blockMat_mul_tmpl(ea, ctxt, cmat);
}
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtBlockMatrix& cmat)
{
  blockMat_mul_tmpl(ea, ctxt, cmat);
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
  bootTask->exec1(D,
    [&](long first, long last) {
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
    }
  );
  FHE_NTIMER_STOP(blockMat1);


  FHE_NTIMER_START(blockMat3);
  bootTask->exec1(d,
    [&](long first, long last) {
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
    }
  );
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








//=======================================================================================




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



  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const PlaintextBlockMatrixBaseInterface& mat)
  {
    PA_BOILER

    const PlaintextBlockMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

    vector<RX> res;
    res.resize(n);
    for (long j = 0; j < n; j++) {
      vec_R acc, tmp, tmp1;
      mat_R val;

      acc.SetLength(d);
      for (long i = 0; i < n; i++) {
         if (!mat1.get(val, i, j)) {
            VectorCopy(tmp1, data[i], d);
            mul(tmp, tmp1, val);
            add(acc, acc, tmp);
         }
      }
      conv(res[j], acc);
    }

    data = res;
  }

}; 


void mat_mul(const EncryptedArray& ea, NewPlaintextArray& pa, 
  const PlaintextMatrixBaseInterface& mat)
{
  ea.dispatch<mat_mul_pa_impl>(Fwd(pa), mat); 
}



void mat_mul(const EncryptedArray& ea, NewPlaintextArray& pa, 
  const PlaintextBlockMatrixBaseInterface& mat)
{
  ea.dispatch<mat_mul_pa_impl>(Fwd(pa), mat); 
}




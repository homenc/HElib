/* Copyright (C) 2012-2014 IBM Corp.
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
#include "EvalMap.h"
#include "powerful.h"
#include "matrix.h"

// The callback interface for the matrix-multiplication routines.

//! \cond FALSE (make doxygen ignore these classes)
template<class type>
class AltStep2Matrix : public PlaintextMatrixInterface<type> 
{
public:
  PA_INJECT(type)

private:
  const EncryptedArray& ea;
  shared_ptr<CubeSignature> sig;
  long dim;

  Mat<RX> A;

public:
  // constructor
  AltStep2Matrix(const EncryptedArray& _ea, 
              shared_ptr<CubeSignature> _sig,
              const Vec<long>& reps,
              long _dim,
              long cofactor,
              bool invert = false);

  virtual const EncryptedArray& getEA() const { return ea; }

  virtual bool get(RX& out, long i, long j) const;
              
};

template<class type>
bool AltStep2Matrix<type>::get(RX& out, long i, long j) const
{
  out = A[i][j];
  return false;
}

template<class type>
AltStep2Matrix<type>::AltStep2Matrix(const EncryptedArray& _ea, 
                               shared_ptr<CubeSignature> _sig,
                               const Vec<long>& reps,
                               long _dim,
                               long cofactor,
                               bool invert)
: ea(_ea), sig(_sig), dim(_dim)
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  const RX& G = ea.getDerived(type()).getG();

  long sz = sig->getDim(dim);
  assert(sz == reps.length());


  Vec<RX> points;
  points.SetLength(sz);
  for (long j = 0; j < sz; j++) 
    points[j] = RX(reps[j]*cofactor, 1) % G;

  A.SetDims(sz, sz);
  for (long j = 0; j < sz; j++)
    A[0][j] = 1;

  for (long i = 1; i < sz; i++)
    for (long j = 0; j < sz; j++)
      A[i][j] = (A[i-1][j] * points[j]) % G;

  if (invert) {
    REBak ebak; ebak.save(); ea.getDerived(type()).restoreContextForG();

    mat_RE A1, A2;
    conv(A1, A);

    long p = ea.getAlMod().getZMStar().getP();
    long r = ea.getAlMod().getR();

    ppInvert(A2, A1, p, r);
    conv(A, A2);
 }
}


PlaintextMatrixBaseInterface*
buildAltStep2Matrix(const EncryptedArray& ea, 
                 shared_ptr<CubeSignature> sig,
                 const Vec<long>& reps,
                 long dim,
                 long cofactor,
                 bool invert = false)
{
  switch (ea.getAlMod().getTag()) {
  case PA_GF2_tag: 
    return new AltStep2Matrix<PA_GF2>(ea, sig, reps, dim, cofactor, invert);

  case PA_zz_p_tag: 
    return new AltStep2Matrix<PA_zz_p>(ea, sig, reps, dim, cofactor, invert);

  default: return 0;
  }
}

template<class type>
class AltStep1Matrix : public PlaintextBlockMatrixInterface<type> 
{
public:
  PA_INJECT(type)

private:
  const EncryptedArray& ea;

  shared_ptr<CubeSignature> sig;
  long dim;

  Mat< mat_R > A;

public:
  // constructor
  AltStep1Matrix(const EncryptedArray& _ea, 
              shared_ptr<CubeSignature> _sig,
              const Vec<long>& reps,
              long _dim,
              long cofactor,
              bool invert,
              bool normal_basis);

  virtual const EncryptedArray& getEA() const { return ea; }

  virtual bool get(mat_R& out, long i, long j) const;
              
};

template<class type>
bool AltStep1Matrix<type>::get(mat_R& out, long i, long j) const
{
  out = A[i][j];
  return false;
}

template<class type>
AltStep1Matrix<type>::AltStep1Matrix(const EncryptedArray& _ea, 
                               shared_ptr<CubeSignature> _sig,
                               const Vec<long>& reps,
                               long _dim,
                               long cofactor,
                               bool invert,
                               bool normal_basis)
: ea(_ea), sig(_sig), dim(_dim)
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  const RX& G = ea.getDerived(type()).getG();

  assert(dim == sig->getNumDims() - 1);
  assert(sig->getSize() == ea.size());

  long sz = sig->getDim(dim);
  assert(sz == reps.length());

  long d = deg(G);

  // so sz == phi(m_last)/d, where d = deg(G) = order of p mod m

  Vec<RX> points;
  points.SetLength(sz);
  for (long j = 0; j < sz; j++) 
    points[j] = RX(reps[j]*cofactor, 1) % G;

  Mat<RX> AA;

  AA.SetDims(sz*d, sz);
  for (long j = 0; j < sz; j++)
    AA[0][j] = 1;

  for (long i = 1; i < sz*d; i++)
    for (long j = 0; j < sz; j++)
      AA[i][j] = (AA[i-1][j] * points[j]) % G;

  A.SetDims(sz, sz);
  for (long i = 0; i < sz; i++)
    for (long j = 0; j < sz; j++) {
      A[i][j].SetDims(d, d);
      for (long k = 0; k < d; k++)
        VectorCopy(A[i][j][k], AA[i*d + k][j], d);
    }


  if (invert) {

    mat_R A1, A2;

    A1.SetDims(sz*d, sz*d);
    for (long i = 0; i < sz*d; i++)
      for (long j = 0; j < sz*d; j++)
        A1[i][j] = A[i/d][j/d][i%d][j%d];


    long p = ea.getAlMod().getZMStar().getP();
    long r = ea.getAlMod().getR();

    ppInvert(A2, A1, p, r);

    for (long i = 0; i < sz*d; i++)
      for (long j = 0; j < sz*d; j++)
        A[i/d][j/d][i%d][j%d] = A2[i][j];
    
    
    if (normal_basis) {
      const Mat<R>& CB = ea.getDerived(type()).getNormalBasisMatrix();

      // multiply each entry of A on the right by CB
      for (long i = 0; i < sz; i++)
        for (long j = 0; j < sz; j++)
          A[i][j] =  A[i][j] * CB;
    }
  }
}


PlaintextBlockMatrixBaseInterface*
buildAltStep1Matrix(const EncryptedArray& ea, 
                 shared_ptr<CubeSignature> sig,
                 const Vec<long>& reps,
                 long dim,
                 long cofactor,
                 bool invert,
                 bool normal_basis)
{
  switch (ea.getAlMod().getTag()) {
  case PA_GF2_tag: 
    return new AltStep1Matrix<PA_GF2>(ea, sig, reps, dim, cofactor, invert, normal_basis);

  case PA_zz_p_tag: 
    return new AltStep1Matrix<PA_zz_p>(ea, sig, reps, dim, cofactor, invert, normal_basis);

  default: return 0;
  }
}
//! \endcond


void init_representatives(Vec<long>& representatives, long dim, 
                          const Vec<long>& mvec, const PAlgebra& zMStar)
{
  assert(dim >= 0 && dim < mvec.length());

  // special case
  if (dim >= LONG(zMStar.numOfGens())) {
    representatives.SetLength(1);
    representatives[0] = 1;
    return;
  }
  
  long m = mvec[dim];
  long D = zMStar.OrderOf(dim);
  long g = InvMod(zMStar.ZmStarGen(dim) % m, m);

  representatives.SetLength(D);
  for (long i = 0; i < D; i++)
    representatives[i] = PowerMod(g, i, m);
}


EvalMap::EvalMap(const EncryptedArray& _ea, const Vec<long>& mvec, bool _invert,
                 bool normal_basis)

  : ea(_ea), invert(_invert)
{
  const FHEcontext& context = ea.getContext();
  const PAlgebra& zMStar = context.zMStar;
  
  long p = zMStar.getP();
  long d = zMStar.getOrdP();

  // FIXME: we should check that ea was initilized with 
  // G == factors[0], but this is a slight pain to check
  // currently

  // NOTE: this code is derived from a more general setting, and
  // could certainly be greatly simplified

  nfactors = mvec.length();

  assert(nfactors > 0);

  for (long i = 0; i < nfactors; i++)
    for (long j = i+1; j < nfactors; j++)
      assert(GCD(mvec[i], mvec[j]) == 1);

  long m = computeProd(mvec);
  assert(m == long(zMStar.getM()));

  Vec<long> phivec(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)  phivec[i] = phi_N(mvec[i]);
  long phim = computeProd(phivec);

  Vec<long> dprodvec(INIT_SIZE, nfactors+1);
  dprodvec[nfactors] = 1;
  
  for (long i = nfactors-1; i >= 0; i--)
    dprodvec[i] = dprodvec[i+1] *
      multOrd(PowerMod(p % mvec[i], dprodvec[i+1], mvec[i]), mvec[i]);

  Vec<long> dvec(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    dvec[i] = dprodvec[i] / dprodvec[i+1];

  long nslots = phim/d;
  assert(d == dprodvec[0]);
  assert(nslots == long(zMStar.getNSlots()));

  long inertPrefix = 0;
  for (long i = 0; i < nfactors && dvec[i] == 1; i++) {
    inertPrefix++;
  }

  if (inertPrefix != nfactors-1)
    Error("EvalMap: case not handled: bad inertPrefix");

  Vec< Vec<long> > local_reps(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    init_representatives(local_reps[i], i, mvec, zMStar);

  Vec<long> crtvec(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++) 
    crtvec[i] = (m/mvec[i]) * InvMod((m/mvec[i]) % mvec[i], mvec[i]);

  Vec<long> redphivec(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    redphivec[i] = phivec[i]/dvec[i];

  CubeSignature redphisig(redphivec);

  Vec< shared_ptr<CubeSignature> > sig_sequence;
  sig_sequence.SetLength(nfactors+1);
  sig_sequence[nfactors] = shared_ptr<CubeSignature>(new CubeSignature(phivec));

  Vec<long> reduced_phivec = phivec;

  for (long dim = nfactors-1; dim >= 0; dim--) {
    reduced_phivec[dim] /= dvec[dim];
    sig_sequence[dim] = 
      shared_ptr<CubeSignature>(new CubeSignature(reduced_phivec));
  }

  {long dim = nfactors - 1;
  shared_ptr<PlaintextBlockMatrixBaseInterface> blockMat(
        buildAltStep1Matrix(ea, sig_sequence[dim],
			    local_reps[dim], dim, m/mvec[dim], invert, normal_basis));
#ifdef EVALMAP_CACHED // cache the matrix of constants
  compMat1D(ea, mat1, *blockMat, dim);
#else
  mat1 = blockMat;
#endif
  }
  matvec.SetLength(nfactors-1);
  for (long dim=nfactors-2; dim>=0; --dim) {
    shared_ptr<PlaintextMatrixBaseInterface> mat_dim(
          buildAltStep2Matrix(ea, sig_sequence[dim],
                              local_reps[dim], dim, m/mvec[dim], invert));
#ifdef EVALMAP_CACHED // cache the matrix of constants
    compMat1D(ea, matvec[dim], *mat_dim, dim);
#else
    matvec[dim] = mat_dim;
#endif
  }
}

void EvalMap::apply(Ctxt& ctxt) const
{
  if (!invert) {
    // forward direction

#ifdef EVALMAP_CACHED // cached matrix of constants
    mat_mul1D(ea, ctxt, mat1, nfactors-1);
#else
    mat_mul1D(ea, ctxt, *mat1, nfactors-1);
#endif

    for (long i = matvec.length()-1; i >= 0; i--) {
#ifdef EVALMAP_CACHED // cached matrix of constants
      mat_mul1D(ea, ctxt, matvec[i], i);
#else
      mat_mul1D(ea, ctxt, *matvec[i], i);
#endif
    }
  }
  else {
    for (long i = 0; i < matvec.length(); i++) {
#ifdef EVALMAP_CACHED // cached matrix of constants
      mat_mul1D(ea, ctxt, matvec[i], i);
#else
      mat_mul1D(ea, ctxt, *matvec[i], i);
#endif
    }

#ifdef EVALMAP_CACHED // cached matrix of constants
    mat_mul1D(ea, ctxt, mat1, nfactors-1);
#else
    mat_mul1D(ea, ctxt, *mat1, nfactors-1);
#endif
  }
}

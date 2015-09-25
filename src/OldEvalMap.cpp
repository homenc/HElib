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
#include "OldEvalMap.h"
#include "powerful.h"

//! \cond FALSE (make doxygen ignore these classes)
template<class type>
class Step1Matrix : public PlaintextBlockMatrixInterface<type> 
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
  Step1Matrix(const EncryptedArray& _ea, 
              shared_ptr<CubeSignature> _sig,
              const Vec<long>& reps,
              long _dim,
              long cofactor,
              bool invert = false);

  virtual const EncryptedArray& getEA() const { return ea; }

  virtual bool get(mat_R& out, long i, long j) const;
              
};

template<class type>
bool Step1Matrix<type>::get(mat_R& out, long i, long j) const
{
  long i1 = sig->getCoord(i, dim);
  long j1 = sig->getCoord(j, dim);

  if (sig->addCoord(i, dim, -i1) != sig->addCoord(j, dim, -j1)) 
    return true;

  out = A[i1][j1];
  return false;
}

template<class type>
Step1Matrix<type>::Step1Matrix(const EncryptedArray& _ea, 
                               shared_ptr<CubeSignature> _sig,
                               const Vec<long>& reps,
                               long _dim,
                               long cofactor,
                               bool invert)
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
    REBak ebak; ebak.save(); ea.getDerived(type()).restoreContextForG();

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
 }
}


PlaintextBlockMatrixBaseInterface*
buildStep1Matrix(const EncryptedArray& ea, 
                 shared_ptr<CubeSignature> sig,
                 const Vec<long>& reps,
                 long dim,
                 long cofactor,
                 bool invert = false)
{
  switch (ea.getAlMod().getTag()) {
  case PA_GF2_tag: 
    return new Step1Matrix<PA_GF2>(ea, sig, reps, dim, cofactor, invert);

  case PA_zz_p_tag: 
    return new Step1Matrix<PA_zz_p>(ea, sig, reps, dim, cofactor, invert);

  default: return 0;
  }
}

/***** END Step1 stuff *****/



template<class type>
class Step2Matrix : public PlaintextMatrixInterface<type> 
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
  Step2Matrix(const EncryptedArray& _ea, 
              shared_ptr<CubeSignature> _sig,
              const Vec<long>& reps,
              long _dim,
              long cofactor,
              bool invert = false);

  virtual const EncryptedArray& getEA() const { return ea; }

  virtual bool get(RX& out, long i, long j) const;
              
};

template<class type>
bool Step2Matrix<type>::get(RX& out, long i, long j) const
{
  long i1 = sig->getCoord(i, dim);
  long j1 = sig->getCoord(j, dim);

  if (sig->addCoord(i, dim, -i1) != sig->addCoord(j, dim, -j1)) 
    return true;

  out = A[i1][j1];
  return false;
}

template<class type>
Step2Matrix<type>::Step2Matrix(const EncryptedArray& _ea, 
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
buildStep2Matrix(const EncryptedArray& ea, 
                 shared_ptr<CubeSignature> sig,
                 const Vec<long>& reps,
                 long dim,
                 long cofactor,
                 bool invert = false)
{
  switch (ea.getAlMod().getTag()) {
  case PA_GF2_tag: 
    return new Step2Matrix<PA_GF2>(ea, sig, reps, dim, cofactor, invert);

  case PA_zz_p_tag: 
    return new Step2Matrix<PA_zz_p>(ea, sig, reps, dim, cofactor, invert);

  default: return 0;
  }
}



// two-step tower stuff...

class TowerBase { 
public:
  virtual ~TowerBase() { }
};

template<class type>
class Tower : public TowerBase {
public:
  PA_INJECT(type)

  long cofactor, d1, d2, p, r;

  long d;
  RE zeta; // = [X^cofactor mod G]

  RX H;  // = the min poly of zeta over R

  Mat<R> M2, M2i;
  // M2 is the matrix that takes us from the two-step tower
  // to the one-step tower, and M2i is its inverse.

  mutable shared_ptr< Mat<RE> > linPolyMatrix;


  Tower(long _cofactor, long _d1, long _d2, long _p, long _r)
    : cofactor(_cofactor), d1(_d1), d2(_d2), p(_p), r(_r) 
  {
    d = RE::degree();
    assert(d == d1*d2);

    const RXModulus& G = RE::modulus();

    zeta = conv<RE>(RX(cofactor, 1));  // zeta = [X^cofactor mod G]

    // compute H = min poly of zeta over R

    Mat<R> M1;
    M1.SetDims(d1, d);

    for (long i = 0; i < d1; i++) {
      VectorCopy(M1[i], rep(power(zeta, i)), d);
    }

    Vec<R> V1;
    VectorCopy(V1, rep(power(zeta, d1)), d);

    Mat<R> M1sq;

    Mat<R> R1;
    R1.SetDims(d, d1);

    for (;;) {
      for (long i = 0; i < d; i++)
        for (long j = 0; j < d1; j++)
          random(R1[i][j]);

      M1sq = M1*R1;
    
      Mat<long> M1sqInt = conv< Mat<long> >(M1sq);
      {
         RBak bak; bak.save();
         GenericModulus<R>::init(p);
         Mat<R> M1sq_modp = conv< Mat<R> >(M1sqInt);
         if (determinant(M1sq_modp) != 0) break;
      }
   
    }

    Vec<R> V1sq = V1*R1;

    Mat<R> M1sqi;
    ppInvert(M1sqi, M1sq, p, r);

    Vec<R> W1 = V1sq * M1sqi;

    assert(W1*M1 == V1);

    H = RX(d1, 1) - conv<RX>(W1);
    // H is the min poly of zeta

    assert(eval(H, zeta) == 0);

    // compute matrices M2 and M2i

    M2.SetDims(d, d);

    for (long i = 0; i < d2; i++) {
      // construct rows [i..i+d1)
      for (long j = 0; j < d1; j++) {
         VectorCopy(M2[i*d1+j], (rep(power(zeta, j)) << i) % G, d);
      }
    }

    ppInvert(M2i, M2, p, r);
  }


  // converts an object represented in the two-step
  // tower representation to the one-step representation

  RE convert2to1(const Vec<RX>& v) const
  {
    assert(v.length() <= d2);

    Vec<R> w;
    w.SetLength(d);

    for (long i = 0; i < v.length(); i++) {
      assert(deg(v[i]) < d1);
      for (long j = 0; j <= deg(v[i]); j++) 
        w[i*d1 + j] = v[i][j];
    }

    Vec<R> z = w*M2;
    return conv<RE>( conv<RX>( z ) );
  }

  // performs the reverse conversion

  Vec<RX> convert1to2(const RE& beta) const
  {
    Vec<R> z = VectorCopy(rep(beta), d);
    Vec<R> w = z*M2i;

    Vec<RX> res;
    res.SetLength(d2);

    for (long i = 0; i < d2; i++) {
      Vec<R> tmp;
      tmp.SetLength(d1);
      for (long j = 0; j < d1; j++)
        tmp[j] = w[i*d1+j];

      res[i] = conv<RX>(tmp);
    }

    return res;
  }

  void buildLinPolyMatrix(Mat<RE>& M) const
  {
     ZZ q = power_ZZ(p, d1);

     M.SetDims(d2, d2);

     for (long j = 0; j < d2; j++) 
        conv(M[0][j], RX(j, 1));

     for (long i = 1; i < d2; i++)
        for (long j = 0; j < d2; j++)
           M[i][j] = power(M[i-1][j], q);
  }



  void buildLinPolyCoeffs(Vec<RE>& C_out, const Vec<RE>& L) const
  {
     FHE_TIMER_START;

     if (!linPolyMatrix) {
       FHE_NTIMER_START(buildLinPolyCoeffs_buildMatrix);

       Mat<RE> M;
       buildLinPolyMatrix(M);
       Mat<RE> Minv;
       ppInvert(Minv, M, p, r);
       linPolyMatrix = shared_ptr< Mat<RE> >(new Mat<RE>(Minv) );
     }

     Vec<RE> C;
     mul(C, L, *linPolyMatrix);

     C_out = C;
  }

  void applyLinPoly(RE& beta, const Vec<RE>& C, const RE& alpha) const
  {
     assert(d2 == C.length());
     ZZ q = power_ZZ(p, d1);

     RE gamma, res;

     gamma = conv<RE>(RX(1, 1));
     res = C[0]*alpha;
     for (long i = 1; i < d2; i++) {
        gamma = power(gamma, q);
        res += C[i]*conv<RE>(CompMod(rep(alpha), rep(gamma), RE::modulus()));
     }

     beta = res;
  }


  void print(const EncryptedArray& ea, ostream& s, const NewPlaintextArray& v, long nrows) const
  {
    vector<ZZX> v1;
    decode(ea, v1, v); 
    Vec<RX> v2;
    convert(v2, v1);
    for (long i = 0; i < v2.length(); i++) {
      if (i % nrows == 0) s << "\n";
      s << convert1to2(conv<RE>(v2[i])) << "\n";
    }
  }

};

template class Tower<PA_GF2>;
template class Tower<PA_zz_p>;

TowerBase* buildTowerBase(const EncryptedArray& ea, 
                      long cofactor, long d1, long d2)
{
  long p = ea.getAlMod().getZMStar().getP();
  long r = ea.getAlMod().getR();

  switch (ea.getAlMod().getTag()) {
    case PA_GF2_tag: {
      GF2EBak ebak; ebak.save(); 
      ea.restoreContextForG();
      return new Tower<PA_GF2>(cofactor, d1, d2, p, r);
    }
    case PA_zz_p_tag: {
      zz_pBak bak; bak.save(); zz_pEBak ebak; ebak.save();
      ea.restoreContext(); ea.restoreContextForG();
      return new Tower<PA_zz_p>(cofactor, d1, d2, p, r);
    }
    default: return 0;
  }
}


template<class type>
class Step1aMatrix : public PlaintextBlockMatrixInterface<type> 
{
public:
  PA_INJECT(type)

private:
  const EncryptedArray& ea;
  long cofactor, d1, d2, phim1;
  shared_ptr<TowerBase> towerBase;
  bool invert;



  Mat< mat_R > A;

public:
  // constructor
  Step1aMatrix(const EncryptedArray& _ea, 
              const Vec<long>& reps,
              long _cofactor, long _d1, long _d2, long _phim1,
              shared_ptr<TowerBase> _towerBase,
              bool _invert = false);

  virtual const EncryptedArray& getEA() const { return ea; }

  virtual bool get(mat_R& out, long i, long j) const;
              
};

template<class type>
bool Step1aMatrix<type>::get(mat_R& out, long i, long j) const
{
  long sz = phim1/d1;

  long i_lo = i*d2;
  long i_hi = i_lo + d2 - 1;
  long j_lo = j*d2;
  long j_hi = j_lo + d2 - 1;

  if (i_hi/sz < j_lo/sz || j_hi/sz < i_lo/sz) return true;

  long d = d1*d2;

  Mat<R> tmp;

  tmp.SetDims(d, d);
  clear(tmp);

  for (long i1 = i_lo; i1 <= i_hi; i1++)
    for (long j1 = j_lo; j1 <= j_hi; j1++) 
      if (i1/sz == j1/sz) {
        long i2 = i1 % sz;
        long j2 = j1 % sz;
        for (long i3 = 0; i3 < d1; i3++)
          for (long j3 = 0; j3 < d1; j3++)
            tmp[(i1-i_lo)*d1 + i3][(j1-j_lo)*d1 + j3] = A[i2][j2][i3][j3];
      }

  Tower<type> *tower = dynamic_cast<Tower<type> *>(towerBase.get());

  if (invert)
    mul(out, tower->M2i, tmp);
  else
    mul(out, tmp, tower->M2);

  return false;
}

template<class type>
Step1aMatrix<type>::Step1aMatrix(const EncryptedArray& _ea, 
                               const Vec<long>& reps,
                               long _cofactor, long _d1, long _d2, long _phim1,
                               shared_ptr<TowerBase> _towerBase,
                               bool _invert)
: ea(_ea), cofactor(_cofactor), d1(_d1), d2(_d2), phim1(_phim1), 
  towerBase(_towerBase), invert(_invert)
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  REBak ebak; ebak.save(); ea.getDerived(type()).restoreContextForG();

  long p = ea.getAlMod().getZMStar().getP();
  long r = ea.getAlMod().getR();

  Tower<type> *tower = dynamic_cast<Tower<type> *>(towerBase.get());

  const RX& H = tower->H;

  assert(phim1 % d1 == 0);

  long sz = phim1/d1;

  Vec<RX> points;
  points.SetLength(sz);
  for (long j = 0; j < sz; j++)
    points[j] = RX(reps[j], 1) % H;

  Mat<RX> AA;
  AA.SetDims(sz*d1, sz);
  for (long j = 0; j < sz; j++)
    AA[0][j] = 1;

  for (long i = 1; i < sz*d1; i++)
    for (long j = 0; j < sz; j++)
      AA[i][j] = (AA[i-1][j] * points[j]) % H;

  A.SetDims(sz, sz);
  for (long i = 0; i < sz; i++)
    for (long j = 0; j < sz; j++) {
      A[i][j].SetDims(d1, d1);
      for (long k = 0; k < d1; k++)
        VectorCopy(A[i][j][k], AA[i*d1 + k][j], d1);
    }

  if (invert) {
    Mat<R> A1, A2;
    A1.SetDims(sz*d1, sz*d1);
    for (long i = 0; i < sz*d1; i++)
      for (long j = 0; j < sz*d1; j++)
        A1[i][j] = A[i/d1][j/d1][i%d1][j%d1];


    ppInvert(A2, A1, p, r);

    for (long i = 0; i < sz*d1; i++)
      for (long j = 0; j < sz*d1; j++)
        A[i/d1][j/d1][i%d1][j%d1] = A2[i][j];
  }
}


PlaintextBlockMatrixBaseInterface*
buildStep1aMatrix(const EncryptedArray& ea, 
                 const Vec<long>& reps,
                 long cofactor, long d1, long d2, long phim1,
                 shared_ptr<TowerBase> towerBase,
                 bool invert = false)
{
  switch (ea.getAlMod().getTag()) {
  case PA_GF2_tag: 
    return new Step1aMatrix<PA_GF2>(ea, reps, cofactor, d1, d2, phim1, towerBase, invert);

  case PA_zz_p_tag: 
    return new Step1aMatrix<PA_zz_p>(ea, reps, cofactor, d1, d2, phim1, towerBase, invert);

  default: return 0;
  }
}

/***** END Step1a stuff *****/


class Step2aShuffleBase {
public:

Vec<long> new_order;

virtual void apply(NewPlaintextArray& v) const = 0;
virtual void apply(Ctxt& v) const = 0;

};

template<class type>
class Step2aShuffle : public Step2aShuffleBase
{
public:
  PA_INJECT(type)

  virtual void apply(NewPlaintextArray& v) const
  {
    if (invert) 
      applyBack(v);
    else
      applyFwd(v);
  }

  virtual void apply(Ctxt& v) const
  {
    if (invert) 
      applyBack(v);
    else
      applyFwd(v);
  }



private:
  const EncryptedArray& ea;
  shared_ptr<CubeSignature> sig;
  Vec<long> reps;
  long dim;
  long cofactor;
  shared_ptr<TowerBase> towerBase;
  bool invert;
  

  long p, r, d, d1, d2, phim1, phim2, nrows;

  long hfactor;
  Vec<long> cshift;
  Mat<long> intraSlotPerm;
  Mat<long> eval_reordering;
  Mat<RE> eval_mat;
  Mat<RX> inv_mat;

  bool get(Vec<RE>& entry, long i, long j) const;
  bool iget(Vec<RE>& entry, long i, long j) const;

  typedef bool (Step2aShuffle<type>::*get_type)(Vec<RE>&, long, long) const; 

  void mat_mul(NewPlaintextArray& ctxt, get_type) const;
  void mat_mul(Ctxt& ctxt, get_type) const;

  void applyBack(NewPlaintextArray& v) const;
  void applyBack(Ctxt& v) const;

  void applyFwd(NewPlaintextArray& v) const;
  void applyFwd(Ctxt& v) const;

public:
  // constructor
  Step2aShuffle(const EncryptedArray& _ea, 
              shared_ptr<CubeSignature> _sig,
              const Vec<long>& _reps,
              long _dim,
              long _cofactor,
              shared_ptr<TowerBase> _towerBase,
              bool _invert = false);

};


template<class type>
Step2aShuffle<type>::Step2aShuffle(const EncryptedArray& _ea, 
                                   shared_ptr<CubeSignature> _sig,
                                   const Vec<long>& _reps,
                                   long _dim,
                                   long _cofactor,
                                   shared_ptr<TowerBase> _towerBase,
                                   bool _invert)

: ea(_ea), sig(_sig), reps(_reps), dim(_dim), cofactor(_cofactor),
  towerBase(_towerBase), invert(_invert)
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  REBak ebak; ebak.save(); ea.getDerived(type()).restoreContextForG();

  Tower<type> *tower = dynamic_cast<Tower<type> *>(towerBase.get());

  p = tower->p;
  r = tower->r;
  d = tower->d;
  d1 = tower->d1;
  d2 = tower->d2;

  phim1 = sig->getDim(dim+1); // actually, phim1/d1
  phim2 = sig->getDim(dim) * d2;

  // cout << "phim1=" << phim1 << ", phim2=" << phim2 << ", d2=" << d2 << "\n";

  nrows = phim1*phim2/d2;

  hfactor = GCD(d2, phim1);

  cshift.SetLength(d2);

  Mat< Pair<long, long> > mapping;
  mapping.SetDims(nrows, d2);

  for (long i = 0; i < phim1*phim2; i++)
    mapping[i/d2][i%d2] = Pair<long,long>(i%phim1, i/phim1);

  // cout << "mapping:\n";
  // cout << mapping << "\n";

  Mat<long> hshift;
  
  hshift.SetDims(nrows, d2/hfactor); 

  for (long j = 0 ; j < d2/hfactor; j++) {
    hshift[0][j] = 0;

    for (long i = 1; i < nrows; i++) 
      if (mapping[i][j*hfactor].a != 0)
        hshift[i][j] = hshift[i-1][j];
      else
        hshift[i][j] = (hshift[i-1][j] + 1) % hfactor;
  }

  // cout << "hshift:\n";
  // cout << hshift << "\n";

  // apply the hshift's to mapping

  for (long i = 0; i < nrows; i++) { 
    for (long j = 0; j < d2/hfactor; j++) {
      // rotate subarray mapping[i][j*hfactor..j*hfactor+hfactor-1]
      // by hshift[i][j]

      long amt = hshift[i][j];

      Vec< Pair<long, long> > tmp1, tmp2;
      tmp1.SetLength(hfactor);
      tmp2.SetLength(hfactor);
 
      for (long k = 0; k < hfactor; k++) tmp1[k] = mapping[i][j*hfactor+k];
      for (long k = 0; k < hfactor; k++) tmp2[(k+amt)%hfactor] = tmp1[k];
      for (long k = 0; k < hfactor; k++) mapping[i][j*hfactor+k] = tmp2[k];
    }
  }

  
  // cout << "mapping:\n";
  // cout << mapping << "\n";

  for (long j = 0; j < d2; j++) {
    long amt = 0;

    while (mapping[0][j].a != 0) {
      amt++;

      // rotate column j of mapping mapping by 1
      Vec< Pair<long, long> > tmp1, tmp2;
      tmp1.SetLength(nrows);
      tmp2.SetLength(nrows);

      for (long i = 0; i < nrows; i++) tmp1[i] = mapping[i][j];
      for (long i = 0; i < nrows; i++) tmp2[(i+1)%nrows] = tmp1[i];
      for (long i = 0; i < nrows; i++) mapping[i][j] = tmp2[i];
    } 

    cshift[j] = amt;
  }

  // cout << "mapping:\n";
  // cout << mapping << "\n";

  new_order.SetLength(phim1);
  for (long i = 0; i < phim1; i++)
    new_order[i] = mapping[i][0].a;

  // cout << new_order << "\n";


  intraSlotPerm.SetDims(nrows, d2);

  for (long i = 0; i < nrows; i++)
    for (long j = 0; j < d2; j++)
      intraSlotPerm[i][j] = (j/hfactor)*hfactor + mcMod(j - hshift[i][j/hfactor], hfactor);

  eval_reordering.SetDims(phim1, phim2);

  for (long i = 0; i < nrows; i++)
    for (long j = 0; j < d2; j++) {
      eval_reordering[i % phim1][(i / phim1)*d2 + j] = mapping[i][j].b;
      assert(mapping[i][j].a == new_order[i % phim1]);
    }

  // cout << "eval_reordering: \n";
  // cout << eval_reordering << "\n";

  eval_mat.SetDims(phim2, phim2/d2);

  Vec<RE> points;
  points.SetLength(phim2/d2);
  for (long j = 0; j < phim2/d2; j++)
    points[j] = conv<RE>(RX(reps[j]*cofactor, 1));

  for (long j = 0; j < phim2/d2; j++)
    eval_mat[0][j] = 1;

  for (long i = 1; i < phim2; i++)
    for (long j = 0; j < phim2/d2; j++)
      eval_mat[i][j] = eval_mat[i-1][j] * points[j];

  if (invert) {
    Mat<RX> inv_mat1;
    inv_mat1.SetDims(phim2, phim2);
    for (long i = 0; i < phim2; i++) {
      for (long j = 0; j < phim2/d2; j++) {
        Vec<RX> tmp1 = tower->convert1to2(eval_mat[i][j]);
        for (long k = 0; k < d2; k++)
          inv_mat1[i][j*d2+k] = tmp1[k];
      }
    }

    eval_mat.kill(); // we no longer need it

    { // temporarily switch RE::modulus to the minpoly of the subring
      REBak ebak1; ebak1.save();
      RE::init(tower->H);
      Mat<RE> inv_mat2, inv_mat3;
      conv(inv_mat2, inv_mat1);
      ppInvert(inv_mat3, inv_mat2, p, r);
      conv(inv_mat, inv_mat3);
    }
  }

}

template<class type>
bool Step2aShuffle<type>::get(Vec<RE>& entry, long i, long j) const 
{
  long i1 = sig->getCoord(i, dim);
  long j1 = sig->getCoord(j, dim);

  if (sig->addCoord(i, dim, -i1) != sig->addCoord(j, dim, -j1)) 
    return true;

  long i2 = sig->getCoord(i, dim+1);

  for (long i3 = 0; i3 < d2; i3++) {
    entry[i3] = eval_mat[eval_reordering[i2][i1*d2+i3]][j1];
  }

  return false;
}


template<class type>
bool Step2aShuffle<type>::iget(Vec<RE>& entry, long i, long j) const 
{
  long i1 = sig->getCoord(i, dim);
  long j1 = sig->getCoord(j, dim);

  if (sig->addCoord(i, dim, -i1) != sig->addCoord(j, dim, -j1)) 
    return true;

  long j2 = sig->getCoord(j, dim+1);

  Mat<RX> tmp;
  tmp.SetDims(d2, d2);
  for (long i3 = 0; i3 < d2; i3++)
    for (long j3 = 0; j3 < d2; j3++)
      tmp[i3][j3] = inv_mat[i1*d2+i3][eval_reordering[j2][j1*d2+j3]];

  Tower<type> *tower = dynamic_cast<Tower<type> *>(towerBase.get());

  for (long i3 = 0; i3 < d2; i3++) {
    entry[i3] = tower->convert2to1(tmp[i3]);
  }

  return false;
}


template<class type>
void Step2aShuffle<type>::mat_mul(NewPlaintextArray& ctxt, get_type get_fn) const
{
  Tower<type> *tower = dynamic_cast<Tower<type> *>(towerBase.get());
  long nslots = ea.size();

  // ctxt.cleanUp(); 

  NewPlaintextArray res(ea);

  Vec<RE> entry;
  entry.SetLength(d2);

  Vec<RE> C;
  C.SetLength(d2);

  
  Vec< Vec<RX> > diag;
  diag.SetLength(nslots);
  for (long j = 0; j < nslots; j++) diag[j].SetLength(d2);

  for (long i = 0; i < nslots; i++) {
    // process diagonal i


    bool zDiag = true;
    long nzLast = -1;

    for (long j = 0; j < nslots; j++) {
      bool zEntry = (this->*get_fn)(entry, mcMod(j-i, nslots), j);

      if (!zEntry) {    // non-empty entry

        zDiag = false;  // mark diagonal as non-empty

        // clear entries between last nonzero entry and this one

        for (long jj = nzLast+1; jj < j; jj++) {
          for (long k = 0; k < d2; k++)
            clear(diag[jj][k]);
        }

        nzLast = j;

        // compute the lin poly coeffs
        tower->buildLinPolyCoeffs(C, entry);
        conv(diag[j], C);
      }
    }

    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries    
    for (long jj = nzLast+1; jj < nslots; jj++) {
      for (long k = 0; k < d2; k++)
        clear(diag[jj][k]);
    }

    // now diag[j] contains the lin poly coeffs

    NewPlaintextArray shCtxt = ctxt;
    rotate(ea, shCtxt, i); 

    // apply the linearlized polynomial
    for (long k = 0; k < d2; k++) {

      // compute the constant
      bool zConst = true;
      vector<ZZX> cvec;
      cvec.resize(nslots);
      for (long j = 0; j < nslots; j++) {
        convert(cvec[j], diag[j][k]);
        if (!IsZero(cvec[j])) zConst = false;
      }

      if (zConst) continue;

      NewPlaintextArray cpoly(ea);
      encode(ea, cpoly, cvec);
      // FIXME: record the encoded polynomial for future use

      NewPlaintextArray shCtxt1 = shCtxt;
      frobeniusAutomorph(ea, shCtxt1, k*d1);
      mul(ea, shCtxt1, cpoly);
      add(ea, res, shCtxt1);
    }
  }
  ctxt = res;
}

// FIXME: this mat_mul can probably be replaced by one that 
// is more like the mat_mul1D routine in EncryptedArray
template<class type>
void Step2aShuffle<type>::mat_mul(Ctxt& ctxt, get_type get_fn) const
{
  Tower<type> *tower = dynamic_cast<Tower<type> *>(towerBase.get());
  long nslots = ea.size();

  ctxt.cleanUp(); 

  Ctxt res(ZeroCtxtLike, ctxt);

  Vec<RE> entry;
  entry.SetLength(d2);

  Vec<RE> C;
  C.SetLength(d2);

  
  Vec< Vec<RX> > diag;
  diag.SetLength(nslots);
  for (long j = 0; j < nslots; j++) diag[j].SetLength(d2);

  for (long i = 0; i < nslots; i++) {
    // process diagonal i

    bool zDiag = true;
    long nzLast = -1;

    for (long j = 0; j < nslots; j++) {
      bool zEntry = (this->*get_fn)(entry, mcMod(j-i, nslots), j);

      if (!zEntry) {    // non-empty entry

        zDiag = false;  // mark diagonal as non-empty

        // clear entries between last nonzero entry and this one

        for (long jj = nzLast+1; jj < j; jj++) {
          for (long k = 0; k < d2; k++)
            clear(diag[jj][k]);
        }

        nzLast = j;

        // compute the lin poly coeffs
        tower->buildLinPolyCoeffs(C, entry);
        conv(diag[j], C);
      }
    }

    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries    
    for (long jj = nzLast+1; jj < nslots; jj++) {
      for (long k = 0; k < d2; k++)
        clear(diag[jj][k]);
    }

    // now diag[j] contains the lin poly coeffs

    Ctxt shCtxt = ctxt;
    ea.rotate(shCtxt, i); 

    // apply the linearlized polynomial
    for (long k = 0; k < d2; k++) {

      // compute the constant
      bool zConst = true;
      vector<ZZX> cvec;
      cvec.resize(nslots);
      for (long j = 0; j < nslots; j++) {
        convert(cvec[j], diag[j][k]);
        if (!IsZero(cvec[j])) zConst = false;
      }

      if (zConst) continue;

      ZZX cpoly;
      ea.encode(cpoly, cvec);
      // FIXME: record the encoded polynomial for future use

      Ctxt shCtxt1 = shCtxt;
      shCtxt1.frobeniusAutomorph(k*d1);
      shCtxt1.multByConstant(cpoly);
      res += shCtxt1;
    }
  }
  ctxt = res;
}






template<class type>
void Step2aShuffle<type>::applyFwd(NewPlaintextArray& v) const
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  REBak ebak; ebak.save(); ea.getDerived(type()).restoreContextForG();

  Tower<type> *tower = dynamic_cast<Tower<type> *>(towerBase.get());
  long nslots = ea.size();

  // cout << "starting shuffle...\n";

  // tower->print(ea, cout, v, nrows);

  // build linPolyCoeffs

  Mat< Vec<ZZX> > C;

  C.SetDims(d2, d2);
  for (long i = 0; i < d2; i++)
    for (long j = 0; j < d2; j++)
      C[i][j].SetLength(nrows);

  // C[i][j][k] is the j-th lin-poly coefficient
  // of the map that projects subslot intraSlotPerm[k][i]
  // onto subslot i

  for (long k = 0; k < nrows; k++) {
    for (long i = 0; i < d2; i++) {
      long idx_in = intraSlotPerm[k][i];
      long idx_out = i;

      Vec< Vec<RX> > map2;
      map2.SetLength(d2);
      map2[idx_in].SetLength(idx_out+1);
      map2[idx_in][idx_out] = 1;
      // map2 projects idx_in ontot idx_out

      Vec<RE> map1;
      map1.SetLength(d2);
      for (long j = 0; j < d2; j++)
        map1[j] = tower->convert2to1(map2[j]);

      Vec<RE> C1;
      tower->buildLinPolyCoeffs(C1, map1);

      for (long j = 0; j < d2; j++)
        C[i][j][k] = conv<ZZX>(rep(C1[j]));
    }
  }

  // mask each sub-slot

  Vec< shared_ptr<NewPlaintextArray> > frobvec; 
  frobvec.SetLength(d2);
  for (long j = 0; j < d2; j++) {
    shared_ptr<NewPlaintextArray> ptr(new NewPlaintextArray(v));
    frobeniusAutomorph(ea, *ptr, j*d1);
    frobvec[j] = ptr;
  }

  Vec< shared_ptr<NewPlaintextArray> > colvec;
  colvec.SetLength(d2);
  for (long i = 0; i < d2; i++) {
    shared_ptr<NewPlaintextArray> acc(new NewPlaintextArray(ea));

    for (long j = 0; j < d2; j++) {
      NewPlaintextArray const1(ea);

      vector<ZZX> vec1;
      vec1.resize(nslots);
      for (long k = 0; k < nslots; k++)
        vec1[k] = C[i][j][k % nrows];
      encode(ea, const1, vec1);

      NewPlaintextArray ctxt1(*frobvec[j]);

      mul(ea, ctxt1, const1);
      add(ea, *acc, ctxt1);
    }

    colvec[i] = acc;
  }

  // for (long i = 0; i < d2; i++) {
    // cout << "column " << i << "\n";
    // tower->print(ea, cout, *colvec[i], nrows);
  // }

  // rotate each subslot 

  for (long i = 0; i < d2; i++) {
    if (cshift[i] == 0) continue;

    if (nrows == nslots) {
      // simple rotation

      rotate(ea, *colvec[i], cshift[i]);

    }
    else {
      // synthetic rotation 

      vector<long> mask;
      mask.resize(nslots);

      for (long j = 0; j < nslots; j++) 
        mask[j] = ((j % nrows) < (nrows - cshift[i]));

      NewPlaintextArray emask(ea);
      encode(ea, emask, mask);

      NewPlaintextArray tmp1(*colvec[i]), tmp2(*colvec[i]);

      mul(ea, tmp1, emask);
      sub(ea, tmp2, tmp1);

      rotate(ea, tmp1, cshift[i]);
      rotate(ea, tmp2, -(nrows-cshift[i]));
      
      add(ea, tmp1, tmp2);
      *colvec[i] = tmp1;
    }
  }

  // for (long i = 0; i < d2; i++) {
    // cout << "column " << i << "\n";
    // tower->print(ea, cout, *colvec[i], nrows);
  // }

  // combine columns

  NewPlaintextArray v1(ea);
  for (long i = 0; i < d2; i++) 
    add(ea, v1, *colvec[i]);


  // apply the matrix

  mat_mul(v1, &Step2aShuffle<type>::get);

  v = v1;
}


template<class type>
void Step2aShuffle<type>::applyBack(NewPlaintextArray& v) const
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  REBak ebak; ebak.save(); ea.getDerived(type()).restoreContextForG();

  Tower<type> *tower = dynamic_cast<Tower<type> *>(towerBase.get());
  long nslots = ea.size();


  mat_mul(v, &Step2aShuffle<type>::iget);

  Mat< Vec<ZZX> > C;

  C.SetDims(d2, d2);
  for (long i = 0; i < d2; i++)
    for (long j = 0; j < d2; j++)
      C[i][j].SetLength(nrows);

  // C[i][j][k] is the j-th lin-poly coefficient
  // of the map that projects subslot i
  // onto subslot i if hfactor == 1, or
  // onto subslot 0 if hfactor != 1

  for (long k = 0; k < nrows; k++) {
    for (long i = 0; i < d2; i++) {
      long idx_in = i;
      long idx_out = (hfactor == 1) ? i : 0;

      Vec< Vec<RX> > map2;
      map2.SetLength(d2);
      map2[idx_in].SetLength(idx_out+1);
      map2[idx_in][idx_out] = 1;
      // map2 projects idx_in ontot idx_out

      Vec<RE> map1;
      map1.SetLength(d2);
      for (long j = 0; j < d2; j++)
        map1[j] = tower->convert2to1(map2[j]);

      Vec<RE> C1;
      tower->buildLinPolyCoeffs(C1, map1);

      for (long j = 0; j < d2; j++)
        C[i][j][k] = conv<ZZX>(rep(C1[j]));
    }
  }

  // mask each sub-slot

  Vec< shared_ptr<NewPlaintextArray> > frobvec; 
  frobvec.SetLength(d2);
  for (long j = 0; j < d2; j++) {
    shared_ptr<NewPlaintextArray> ptr(new NewPlaintextArray(v));
    frobeniusAutomorph(ea, *ptr, j*d1);
    frobvec[j] = ptr;
  }

  Vec< shared_ptr<NewPlaintextArray> > colvec;
  colvec.SetLength(d2);
  for (long i = 0; i < d2; i++) {
    shared_ptr<NewPlaintextArray> acc(new NewPlaintextArray(ea));

    for (long j = 0; j < d2; j++) {
      NewPlaintextArray const1(ea);

      vector<ZZX> vec1;
      vec1.resize(nslots);
      for (long k = 0; k < nslots; k++)
        vec1[k] = C[i][j][k % nrows];
      encode(ea, const1, vec1);

      NewPlaintextArray ctxt1(*frobvec[j]);

      mul(ea, ctxt1, const1);
      add(ea, *acc, ctxt1);
    }

    colvec[i] = acc;
  }
  

  // rotate each subslot 

  for (long i = 0; i < d2; i++) {
    long shamt = mcMod(-cshift[i], nrows);

    if (shamt == 0) continue;

    if (nrows == nslots) {
      // simple rotation

      rotate(ea, *colvec[i], shamt);

    }
    else {
      // synthetic rotation 

      vector<long> mask;
      mask.resize(nslots);

      for (long j = 0; j < nslots; j++) 
        mask[j] = ((j % nrows) < (nrows - shamt));

      NewPlaintextArray emask(ea);
      encode(ea, emask, mask);

      NewPlaintextArray tmp1(*colvec[i]), tmp2(*colvec[i]);

      mul(ea, tmp1, emask);
      sub(ea, tmp2, tmp1);

      rotate(ea, tmp1, shamt);
      rotate(ea, tmp2, -(nrows-shamt));
      
      add(ea, tmp1, tmp2);
      *colvec[i] = tmp1;
    }
  }

  // combine columns...
  // optimized to avoid unnecessary constant muls
  // when hfactor == 1

  NewPlaintextArray v1(ea);

  if (hfactor == 1) {
    for (long i = 0; i < d2; i++) { 
      add(ea, v1, *colvec[i]);
    }
  }
  else {
    for (long i = 0; i < d2; i++) { 
      NewPlaintextArray const1(ea);
      vector<ZZX> vec1;
      vec1.resize(nslots);
      for (long k = 0; k < nslots; k++)
        vec1[k] = conv<ZZX>(RX(intraSlotPerm[k%nrows][i], 1) % RE::modulus());
      encode(ea, const1, vec1);

      NewPlaintextArray ctxt1(*colvec[i]);
      mul(ea, ctxt1, const1);
      
      add(ea, v1, ctxt1);
    }
  }

  v = v1;
}




template<class type>
void Step2aShuffle<type>::applyBack(Ctxt& v) const
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  REBak ebak; ebak.save(); ea.getDerived(type()).restoreContextForG();

  Tower<type> *tower = dynamic_cast<Tower<type> *>(towerBase.get());
  long nslots = ea.size();


  mat_mul(v, &Step2aShuffle<type>::iget);
  v.cleanUp();

  Mat< Vec<ZZX> > C;

  C.SetDims(d2, d2);
  for (long i = 0; i < d2; i++)
    for (long j = 0; j < d2; j++)
      C[i][j].SetLength(nrows);

  // C[i][j][k] is the j-th lin-poly coefficient
  // of the map that projects subslot i
  // onto subslot i if hfactor == 1, or
  // onto subslot 0 if hfactor != 1

  for (long k = 0; k < nrows; k++) {
    for (long i = 0; i < d2; i++) {
      long idx_in = i;
      long idx_out = (hfactor == 1) ? i : 0;

      Vec< Vec<RX> > map2;
      map2.SetLength(d2);
      map2[idx_in].SetLength(idx_out+1);
      map2[idx_in][idx_out] = 1;
      // map2 projects idx_in ontot idx_out

      Vec<RE> map1;
      map1.SetLength(d2);
      for (long j = 0; j < d2; j++)
        map1[j] = tower->convert2to1(map2[j]);

      Vec<RE> C1;
      tower->buildLinPolyCoeffs(C1, map1);

      for (long j = 0; j < d2; j++)
        C[i][j][k] = conv<ZZX>(rep(C1[j]));
    }
  }

  // mask each sub-slot

  Vec< shared_ptr<Ctxt> > frobvec; 
  frobvec.SetLength(d2);
  for (long j = 0; j < d2; j++) {
    shared_ptr<Ctxt> ptr(new Ctxt(v));
    ptr->frobeniusAutomorph(j*d1);
    frobvec[j] = ptr;
  }

  Vec< shared_ptr<Ctxt> > colvec;
  colvec.SetLength(d2);
  for (long i = 0; i < d2; i++) {
    shared_ptr<Ctxt> acc(new Ctxt(ZeroCtxtLike, v));

    for (long j = 0; j < d2; j++) {
      ZZX const1;

      vector<ZZX> vec1;
      vec1.resize(nslots);
      for (long k = 0; k < nslots; k++)
        vec1[k] = C[i][j][k % nrows];
      ea.encode(const1, vec1);

      Ctxt ctxt1(*frobvec[j]);

      ctxt1.multByConstant(const1);
      (*acc) += ctxt1;
    }

    colvec[i] = acc;
  }
  

  // rotate each subslot 

  for (long i = 0; i < d2; i++) {
    long shamt = mcMod(-cshift[i], nrows);

    if (shamt == 0) continue;

    if (nrows == nslots) {
      // simple rotation

      ea.rotate(*colvec[i], shamt);
    }
    else {
      // synthetic rotation 

      vector<long> mask;
      mask.resize(nslots);

      for (long j = 0; j < nslots; j++) 
        mask[j] = ((j % nrows) < (nrows - shamt));

      ZZX emask;
      ea.encode(emask, mask);

      Ctxt tmp1(*colvec[i]), tmp2(*colvec[i]);

      tmp1.multByConstant(emask);
      tmp2 -= tmp1;

      ea.rotate(tmp1, shamt);
      ea.rotate(tmp2, -(nrows-shamt));
      
      tmp1 += tmp2;
      *colvec[i] = tmp1;
    }
  }

  // combine columns...
  // optimized to avoid unnecessary constant muls
  // when hfactor == 1

  Ctxt v1(ZeroCtxtLike, v);

  if (hfactor == 1) {
    for (long i = 0; i < d2; i++) { 
      v1 += *colvec[i];
    }
  }
  else {
    for (long i = 0; i < d2; i++) { 
      ZZX const1;
      vector<ZZX> vec1;
      vec1.resize(nslots);
      for (long k = 0; k < nslots; k++)
        vec1[k] = conv<ZZX>(RX(intraSlotPerm[k%nrows][i], 1) % RE::modulus());
      ea.encode(const1, vec1);

      Ctxt ctxt1(*colvec[i]);
      ctxt1.multByConstant(const1);
      
      v1 += ctxt1;
    }
  }

  v = v1;
}




template<class type>
void Step2aShuffle<type>::applyFwd(Ctxt& v) const
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  REBak ebak; ebak.save(); ea.getDerived(type()).restoreContextForG();

  Tower<type> *tower = dynamic_cast<Tower<type> *>(towerBase.get());
  long nslots = ea.size();

  // cout << "starting shuffle...\n";

  // tower->print(ea, cout, v, nrows);

  // build linPolyCoeffs

  v.cleanUp();

  Mat< Vec<ZZX> > C;

  C.SetDims(d2, d2);
  for (long i = 0; i < d2; i++)
    for (long j = 0; j < d2; j++)
      C[i][j].SetLength(nrows);

  // C[i][j][k] is the j-th lin-poly coefficient
  // of the map that projects subslot intraSlotPerm[k][i]
  // onto subslot i

  for (long k = 0; k < nrows; k++) {
    for (long i = 0; i < d2; i++) {
      long idx_in = intraSlotPerm[k][i];
      long idx_out = i;

      Vec< Vec<RX> > map2;
      map2.SetLength(d2);
      map2[idx_in].SetLength(idx_out+1);
      map2[idx_in][idx_out] = 1;
      // map2 projects idx_in ontot idx_out

      Vec<RE> map1;
      map1.SetLength(d2);
      for (long j = 0; j < d2; j++)
        map1[j] = tower->convert2to1(map2[j]);

      Vec<RE> C1;
      tower->buildLinPolyCoeffs(C1, map1);

      for (long j = 0; j < d2; j++)
        C[i][j][k] = conv<ZZX>(rep(C1[j]));
    }
  }

  // FIXME: in the case where nrows != nslots, which
  // is the same as saying we are at dimension 0, we can
  // avoid the extra masking depth incurred for the 
  // synthetic rotations, by folding them into 
  // the masking/frobenius step. A similar optimization
  // applies to the applyBack routine.


  // mask each sub-slot

  Vec< shared_ptr<Ctxt> > frobvec; 
  frobvec.SetLength(d2);
  for (long j = 0; j < d2; j++) {
    shared_ptr<Ctxt> ptr(new Ctxt(v));
    ptr->frobeniusAutomorph(j*d1);
    frobvec[j] = ptr;
  }

  Vec< shared_ptr<Ctxt> > colvec;
  colvec.SetLength(d2);
  for (long i = 0; i < d2; i++) {
    shared_ptr<Ctxt> acc(new Ctxt(ZeroCtxtLike, v));

    for (long j = 0; j < d2; j++) {
      ZZX const1;

      vector<ZZX> vec1;
      vec1.resize(nslots);
      for (long k = 0; k < nslots; k++)
        vec1[k] = C[i][j][k % nrows];
      ea.encode(const1, vec1);

      Ctxt ctxt1(*frobvec[j]);

      ctxt1.multByConstant(const1);
      (*acc) += ctxt1;
    }

    colvec[i] = acc;
  }

  // for (long i = 0; i < d2; i++) {
    // cout << "column " << i << "\n";
    // tower->print(ea, cout, *colvec[i], nrows);
  // }

  // rotate each subslot 

  for (long i = 0; i < d2; i++) {
    if (cshift[i] == 0) continue;

    if (nrows == nslots) {
      // simple rotation

      ea.rotate(*colvec[i], cshift[i]);

    }
    else {
      // synthetic rotation 

      vector<long> mask;
      mask.resize(nslots);

      for (long j = 0; j < nslots; j++) 
        mask[j] = ((j % nrows) < (nrows - cshift[i]));

      ZZX emask;
      ea.encode(emask, mask);

      Ctxt tmp1(*colvec[i]), tmp2(*colvec[i]);

      tmp1.multByConstant(emask);
      tmp2 -= tmp1;

      ea.rotate(tmp1, cshift[i]);
      ea.rotate(tmp2, -(nrows-cshift[i]));
      
      tmp1 += tmp2;
      *colvec[i] = tmp1;
    }
  }

  // for (long i = 0; i < d2; i++) {
    // cout << "column " << i << "\n";
    // tower->print(ea, cout, *colvec[i], nrows);
  // }

  // conbine columns

  Ctxt v1(ZeroCtxtLike, v);
  for (long i = 0; i < d2; i++) 
    v1 += *colvec[i];


  // apply the matrix

  mat_mul(v1, &Step2aShuffle<type>::get);

  v = v1;
}




Step2aShuffleBase*
buildStep2aShuffle(const EncryptedArray& ea, 
                 shared_ptr<CubeSignature> sig,
                 const Vec<long>& reps,
                 long dim,
                 long cofactor,
                 shared_ptr<TowerBase> towerBase,
                 bool invert = false)
{
  switch (ea.getAlMod().getTag()) {
  case PA_GF2_tag: 
    return new Step2aShuffle<PA_GF2>(ea, sig, reps, dim, cofactor, towerBase, invert);

  case PA_zz_p_tag: 
    return new Step2aShuffle<PA_zz_p>(ea, sig, reps, dim, cofactor, towerBase, invert);

  default: return 0;
  }
}

/***** END Step2a stuff *****/





void init_representatives(Vec<long>& representatives, long m, long p)
{
  Vec<bool> available;
  available.SetLength(m);

  long num_available = 0;

  for (long i = 0; i < m; i++) {
    if (GCD(i, m) == 1) {
      available[i] = true;
      num_available++;
    }
    else
      available[i] = false;
  }

  representatives.SetLength(0);

  while (num_available > 0) {

    // choose next available at random
    long i;
    do {
      i = RandomBnd(m);
    } while (!available[i]);

    append(representatives, i);

    // mark all conjugates as unavailable
    long j = i;
    do {
      available[j] = false;
      num_available--;
      j = MulMod(j, p, m);
    } while (j != i);
  }
}


void init_slot_mappings(Vec<long>& slot_index, 
                        Vec<long>& slot_rotate, 
                        const Vec<long>& representatives, 
                        long m,
                        long p,
                        const FHEcontext& context)
{
   long nslots = representatives.length();

   assert(nslots == long(context.zMStar.getNSlots()));

   slot_index.SetLength(nslots);
   slot_rotate.SetLength(nslots);

   Vec<bool> used; // for debugging
   used.SetLength(nslots);
   for (long i = 0; i < nslots; i++) used[i] = false;
   
   for (long i = 0; i < nslots; i++) {
     long t = representatives[i];
     long h = 0;
     long idx;
     while ((idx = context.zMStar.indexOfRep(InvMod(t, m))) == -1) {
       t = MulMod(t, p, m);
       h++;
     }

     assert(!used[idx]);
     used[idx] = true;
     slot_index[idx] = i;
     slot_rotate[idx] = h;
   }
}

// apply p^{vec[i]} to slot i
void frobeniusAutomorph(Ctxt& ctxt, const EncryptedArray& ea, const Vec<long>& vec)
{
  long d = ea.getDegree();
  long nslots = ea.size();

  // construct masks
  Vec<ZZX> masks;
  masks.SetLength(d);

  for (long i = 0; i < d; i++) {
    vector<long> mask1_vec;
    mask1_vec.resize(nslots);
    for (long j = 0; j < nslots; j++) 
      mask1_vec[j] = (mcMod(vec[j], d) == i);

    ZZX mask1_poly;
    ea.encode(mask1_poly, mask1_vec);
    masks[i] = mask1_poly;
  }

  ctxt.cleanUp();
  Ctxt acc(ZeroCtxtLike, ctxt);
  for (long i = 0; i < d; i++) {
    if (masks[i] != 0) {
      Ctxt tmp = ctxt;
      tmp.frobeniusAutomorph(i);
      tmp.multByConstant(masks[i]);
      acc += tmp;
    }
  }

  ctxt = acc;
}

OldEvalMap::OldEvalMap(const EncryptedArray& _ea, const Vec<long>& mvec, 
                 long width, bool _invert)
  : ea(_ea), invert(_invert)
{
  const FHEcontext& context = ea.getContext();
  const PAlgebra& zMStar = context.zMStar;
  
  long p = zMStar.getP();
  long d = zMStar.getOrdP();

  // FIXME: we should check that ea was initilized with 
  // G == factors[0], but this is a slight pain to check
  // currently

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

  if (inertPrefix == nfactors-1)
    easy = true;
  else if (inertPrefix == nfactors-2)
    easy = false;
  else
    Error("OldEvalMap: case not handled: bad inertPrefix");

  Vec< Vec<long> > local_reps(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    init_representatives(local_reps[i], mvec[i], 
                         PowerMod(p % mvec[i], dprodvec[i+1], mvec[i]));




  Vec<long> crtvec(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++) 
    crtvec[i] = (m/mvec[i]) * InvMod((m/mvec[i]) % mvec[i], mvec[i]);

  Vec<long> redphivec(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    redphivec[i] = phivec[i]/dvec[i];

  CubeSignature redphisig(redphivec);


  Vec<long> global_reps(INIT_SIZE, phim/d);
  for (long i = 0; i < phim/d; i++) {
    global_reps[i] = 0;
    for (long j = 0; j < nfactors; j++) {
      long i1 = redphisig.getCoord(i, j);
      global_reps[i] = (global_reps[i] + crtvec[j]*local_reps[j][i1]) % m;
    }
  }

  Vec<long> slot_index;
  init_slot_mappings(slot_index, slot_rotate, global_reps, m, p, context);

  Vec< shared_ptr<CubeSignature> > sig_sequence;
  sig_sequence.SetLength(nfactors+1);
  sig_sequence[nfactors] = shared_ptr<CubeSignature>(new CubeSignature(phivec));

  Vec<long> reduced_phivec = phivec;

  for (long dim = nfactors-1; dim >= 0; dim--) {
    reduced_phivec[dim] /= dvec[dim];
    sig_sequence[dim] = 
      shared_ptr<CubeSignature>(new CubeSignature(reduced_phivec));
  }

  if (easy) {
    long dim = nfactors - 1;

    mat1 = shared_ptr<PlaintextBlockMatrixBaseInterface>(
      buildStep1Matrix(ea, sig_sequence[dim], local_reps[dim], dim, m/mvec[dim], invert));

    matvec.SetLength(nfactors-1);

    while (dim > 0) {
      dim--;
      matvec[dim] = shared_ptr<PlaintextMatrixBaseInterface>(
        buildStep2Matrix(ea, sig_sequence[dim], local_reps[dim], dim, m/mvec[dim], invert));
    }
  }
  else {
    long m1 = mvec[nfactors-1];
    long cofactor = m/m1;

    long d1 = dvec[nfactors-1];
    long d2 = d/d1;

    tower = shared_ptr<TowerBase>(buildTowerBase(ea, cofactor, d1, d2));
    long dim = nfactors-1;

    mat1 = shared_ptr<PlaintextBlockMatrixBaseInterface>(
      buildStep1aMatrix(ea, local_reps[dim], cofactor, d1, d2, phivec[dim], tower, invert));


    dim--;
    shuffle = shared_ptr<Step2aShuffleBase>(
      buildStep2aShuffle(ea, sig_sequence[dim], local_reps[dim], dim, m/mvec[dim], tower, invert));

    
    long phim1 = shuffle->new_order.length();

    Vec<long> no_i; // inverse function
    no_i.SetLength(phim1);
    for (long i = 0; i < phim1; i++) 
      no_i[shuffle->new_order[i]] = i;

    Vec<long> slot_index1;
    slot_index1.SetLength(nslots);
    for (long i = 0; i < nslots; i++) 
      slot_index1[i] = (slot_index[i]/phim1)*phim1 + no_i[slot_index[i] % phim1];

    slot_index = slot_index1;

    matvec.SetLength(nfactors-2);

    while (dim > 0) {
      dim--;
      matvec[dim] = shared_ptr<PlaintextMatrixBaseInterface>(
        buildStep2Matrix(ea, sig_sequence[dim], local_reps[dim], dim, m/mvec[dim], invert));
    }
  }

  if (invert) {
    Vec<long> slot_index_i; // inverse function
    slot_index_i.SetLength(nslots);
    for (long i = 0; i < nslots; i++) 
      slot_index_i[slot_index[i]] = i;

    slot_index = slot_index_i; 

    for (long i = 0; i < nslots; i++)
      slot_rotate[i] = mcMod(-slot_rotate[i], d);
  }

  Vec<GenDescriptor> gvec(INIT_SIZE, ea.dimension());
  for (long i=0; i<ea.dimension(); i++)
    gvec[i] = GenDescriptor(/*order=*/ea.sizeOfDimension(i),
                            /*good=*/ ea.nativeDimension(i), /*genIdx=*/i); 

  GeneratorTrees trees;
  long cost = trees.buildOptimalTrees(gvec, width);

  if (cost == NTL_MAX_LONG)
    Error("OldEvalMap: can't build network for given width");

  net = shared_ptr<PermNetwork>(new PermNetwork(slot_index, trees));
}

void OldEvalMap::apply(Ctxt& ctxt) const
{
  if (!invert) {
    // forward direction

    mat_mul(ea, ctxt, *mat1);
    if (!easy) shuffle->apply(ctxt);

    for (long i = matvec.length()-1; i >= 0; i--) 
      mat_mul(ea, ctxt, *matvec[i]);

    net->applyToCtxt(ctxt, ea);
    frobeniusAutomorph(ctxt, ea, slot_rotate);

  }
  else {
    frobeniusAutomorph(ctxt, ea, slot_rotate);
    net->applyToCtxt(ctxt, ea);

    for (long i = 0; i < matvec.length(); i++)
      mat_mul(ea, ctxt, *matvec[i]);

    if (!easy) shuffle->apply(ctxt);
    mat_mul(ea, ctxt, *mat1);
  }
}

//! \endcond

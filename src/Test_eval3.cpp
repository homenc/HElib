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



namespace std {} using namespace std;
namespace NTL {} using namespace NTL;


#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include "powerful.h"

#include <cassert>

#if (__cplusplus>199711L)
#include <memory>
#else
#include <tr1/memory>
using namespace tr1;
#warning "using TR1"
#endif



template<class type>
class Step1Matrix : public PlaintextBlockMatrixInterface<type> 
{
public:
  PA_INJECT(type)

private:
  const EncryptedArray& ea;
  const CubeSignature& sig;
  long dim;

  Mat< mat_R > A;

public:
  // constructor
  Step1Matrix(const EncryptedArray& _ea, 
              const CubeSignature& _sig,
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
  long i1 = sig.getCoord(i, dim);
  long j1 = sig.getCoord(j, dim);

  if (sig.addCoord(i, dim, -i1) != sig.addCoord(j, dim, -j1)) 
    return true;

  out = A[i1][j1];
  return false;
}

template<class type>
Step1Matrix<type>::Step1Matrix(const EncryptedArray& _ea, 
                               const CubeSignature& _sig,
                               const Vec<long>& reps,
                               long _dim,
                               long cofactor,
                               bool invert)
: ea(_ea), sig(_sig), dim(_dim)
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  const RX& G = ea.getDerived(type()).getG();

  assert(dim == sig.getNumDims() - 1);
  assert(sig.getSize() == ea.size());

  long sz = sig.getDim(dim);
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
                 const CubeSignature& sig,
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
  const CubeSignature& sig;
  long dim;

  Mat<RX> A;

public:
  // constructor
  Step2Matrix(const EncryptedArray& _ea, 
              const CubeSignature& _sig,
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
  long i1 = sig.getCoord(i, dim);
  long j1 = sig.getCoord(j, dim);

  if (sig.addCoord(i, dim, -i1) != sig.addCoord(j, dim, -j1)) 
    return true;

  out = A[i1][j1];
  return false;
}

template<class type>
Step2Matrix<type>::Step2Matrix(const EncryptedArray& _ea, 
                               const CubeSignature& _sig,
                               const Vec<long>& reps,
                               long _dim,
                               long cofactor,
                               bool invert)
: ea(_ea), sig(_sig), dim(_dim)
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  const RX& G = ea.getDerived(type()).getG();

  long sz = sig.getDim(dim);
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
                 const CubeSignature& sig,
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

template<class type>
class Tower {
public:
  PA_INJECT(type)

  long m1, m2, d1, d2, p, r;

  long d;
  RE zeta; // = [X^m2 mod G]

  RX H;  // = the min poly of zeta over R

  Mat<R> M2, M2i;
  // M2 is the matrix that takes us from the two-step tower
  // to the one-step tower, and M2i is its inverse.


  Tower(long _m1, long _m2, long _d1, long _d2, long _p, long _r)
    : m1(_m1), m2(_m2), d1(_d1), d2(_d2), p(_p), r(_r) 
  {
    d = RE::degree();
    assert(d == d1*d2);

    const RXModulus& G = RE::modulus();

    zeta = conv<RE>(RX(m2, 1));  // zeta = [X^m2 mod G]

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

  RE convert2to1(const Vec<RX>& v)
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

  Vec<RX> convert1to2(const RE& beta)
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

  void buildLinPolyMatrix(Mat<RE>& M)
  {
     long q = power_long(p, d1);

     M.SetDims(d2, d2);

     for (long j = 0; j < d2; j++) 
        conv(M[0][j], RX(j, 1));

     for (long i = 1; i < d2; i++)
        for (long j = 0; j < d2; j++)
           M[i][j] = power(M[i-1][j], q);
  }



  void buildLinPolyCoeffs(Vec<RE>& C_out, const Vec<RE>& L)
  {
     Mat<RE> M;
     buildLinPolyMatrix(M);

     Vec<RE> C;
     ppsolve(C, M, L, p, r);

     C_out = C;
  }

  void applyLinPoly(RE& beta, const Vec<RE>& C, const RE& alpha)
  {
     assert(d2 == C.length());
     long q = power_long(p, d1);

     RE gamma, res;

     gamma = conv<RE>(RX(1, 1));
     res = C[0]*alpha;
     for (long i = 1; i < d2; i++) {
        gamma = power(gamma, q);
        res += C[i]*conv<RE>(CompMod(rep(alpha), rep(gamma), RE::modulus()));
     }

     beta = res;
  }





};

template class Tower<PA_GF2>;
template class Tower<PA_zz_p>;


template<class type>
class Step1aMatrix : public PlaintextBlockMatrixInterface<type> 
{
public:
  PA_INJECT(type)

private:
  const EncryptedArray& ea;
  long m1, m2, d1, d2, phim1;


  copied_ptr< Tower<type> > tower;

  Mat< mat_R > A;

public:
  // constructor
  Step1aMatrix(const EncryptedArray& _ea, 
              const Vec<long>& reps,
              long _m1, long _m2, long _d1, long _d2, long _phim1,
              bool invert = false);

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

  mul(out, tmp, tower->M2);
  return false;
}

template<class type>
Step1aMatrix<type>::Step1aMatrix(const EncryptedArray& _ea, 
                               const Vec<long>& reps,
                               long _m1, long _m2, long _d1, long _d2, long _phim1,
                               bool invert)
: ea(_ea), m1(_m1), m2(_m2), d1(_d1), d2(_d2), phim1(_phim1)
{
  RBak bak; bak.save(); ea.getAlMod().restoreContext();
  REBak ebak; ebak.save(); ea.getDerived(type()).restoreContextForG();
  const RXModulus& G = RE::modulus();

  long p = ea.getAlMod().getZMStar().getP();
  long r = ea.getAlMod().getR();

  tower.set_ptr(new Tower<type>(m1, m2, d1, d2, p, r));

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
    Error("not implemented");
  }
}


PlaintextBlockMatrixBaseInterface*
buildStep1aMatrix(const EncryptedArray& ea, 
                 const Vec<long>& reps,
                 long m1, long m2, long d1, long d2, long phim1,
                 bool invert = false)
{
  switch (ea.getAlMod().getTag()) {
  case PA_GF2_tag: 
    return new Step1aMatrix<PA_GF2>(ea, reps, m1, m2, d1, d2, phim1, invert);

  case PA_zz_p_tag: 
    return new Step1aMatrix<PA_zz_p>(ea, reps, m1, m2, d1, d2, phim1, invert);

  default: return 0;
  }
}

/***** END Step1a stuff *****/




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


void alt_init_representatives(Vec<long>& rep, long m, long gen, long phim)
{
  rep.SetLength(phim);
  rep[0] = 1;
  for (long i = 1; i < phim; i++)
    rep[i] = MulMod(rep[i-1], gen, m);
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
     slot_index[i] = idx;
     slot_rotate[i] = h;
   }
}

void convertToPowerful(Vec<zz_p>& v, const zz_pX& F, const Vec<long>& mvec)
{ 
  long nfactors = mvec.length();

  long m = computeProd(mvec);
  
  Vec<long> phivec;
  phivec.SetLength(nfactors);
  for (long i = 0; i < nfactors; i++) phivec[i] = phi_N(mvec[i]);

  long phim = computeProd(phivec);

  Vec<long> divvec;
  computeDivVec(divvec, m, mvec);

  Vec<long> invvec;
  computeInvVec(invvec, divvec, mvec);

  CubeSignature shortsig(phivec);
  CubeSignature longsig(mvec);

  Vec<long> polyToCubeMap;
  Vec<long> cubeToPolyMap;
  computePowerToCubeMap(polyToCubeMap, cubeToPolyMap, m, mvec, invvec, longsig);

  Vec<long> shortToLongMap;
  computeShortToLongMap(shortToLongMap, shortsig, longsig);


  Vec<zz_pX> cycvec;
  computeCycVec(cycvec, mvec);


  ZZX PhimX = Cyclotomic(m);
  zz_pX phimX = conv<zz_pX>(PhimX);

  HyperCube<zz_p> cube(shortsig);
  HyperCube<zz_p> tmpCube(longsig);

  convertPolyToPowerful(cube, tmpCube, F, cycvec, 
                        polyToCubeMap, shortToLongMap);

  zz_pX poly1;

  convertPowerfulToPoly(poly1, cube, m, shortToLongMap, cubeToPolyMap, phimX);

  if (F == poly1)
    cout << "*********** :-)\n";
  else {
    cout << "*********** :-(\n";
    cout << F << "\n";
    cout << poly1 << "\n";
  }

  v.SetLength(phim);
  for (long i = 0; i < phim; i++) v[i] = cube[i];
}



void  TestIt(long R, long p, long r, long c, long _k, long w, 
               long L, const Vec<long>& mvec)
{
  cerr << "*** TestIt: R=" << R 
       << ", p=" << p
       << ", r=" << r
       << ", c=" << c
       << ", k=" << _k
       << ", w=" << w
       << ", L=" << L
       << ", mvec=" << mvec
       << endl;

  long nfactors = mvec.length();
  for (long i = 0; i < nfactors; i++)
    for (long j = i+1; j < nfactors; j++)
      assert(GCD(mvec[i], mvec[j]) == 1);


  long m = computeProd(mvec);
  assert(GCD(p, m) == 1); 

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

  cout << "dvec=" << dvec << "\n";

  long d = dprodvec[0];
  long nslots = phim/d;

  long inertPrefix = 0;
  for (long i = 0; i < nfactors && dvec[i] == 1; i++) {
    inertPrefix++;
  }

  cout << "inertPrefix=" << inertPrefix << "\n";

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


  FHEcontext context(m, p, r);
  buildModChain(context, L, c);
  context.zMStar.printout();
  cerr << endl;

  ZZX GG;
  GG = context.alMod.getFactorsOverZZ()[0];
  zz_p::init(context.alMod.getPPowR());
  zz_pX G = conv<zz_pX>(GG);
  zz_pE::init(G);

  Vec<zz_pE> global_points(INIT_SIZE, phim/d);
  for (long i = 0; i < phim/d; i++) 
    global_points[i] = conv<zz_pE>(zz_pX(global_reps[i], 1)); 


  zz_pX F;
  random(F, phim);

  Vec<zz_pE> global_values(INIT_SIZE, phim/d);
  for (long i = 0; i < phim/d; i++)
    global_values[i] = eval(F, global_points[i]);

  Vec<zz_p> cube;
  convertToPowerful(cube, F, mvec);

  Vec< Vec<zz_pE> > local_points(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++) {
    local_points[i].SetLength(phivec[i]/dvec[i]);
    for (long j = 0; j < phivec[i]/dvec[i]; j++)
      local_points[i][j] = conv<zz_pE>(zz_pX(local_reps[i][j]*(m/mvec[i]), 1));
  }


  Vec< Vec<zz_pE> > eval_sequence;
  eval_sequence.SetLength(nfactors+1);
  conv(eval_sequence[nfactors], cube);

  Vec< shared_ptr<CubeSignature> > sig_sequence;
  sig_sequence.SetLength(nfactors+1);
  sig_sequence[nfactors] = shared_ptr<CubeSignature>(new CubeSignature(phivec));

  Vec<long> reduced_phivec = phivec;

  for (long dim = nfactors-1; dim >= 0; dim--) {
    reduced_phivec[dim] /= dvec[dim];
    sig_sequence[dim] = 
      shared_ptr<CubeSignature>(new CubeSignature(reduced_phivec));

    shared_ptr<CubeSignature> old_sig = sig_sequence[dim+1];
    shared_ptr<CubeSignature> new_sig = sig_sequence[dim];

    

    long nslices = old_sig->getProd(0, dim); // same for both old and new
    long ncols = old_sig->getProd(dim+1);  // same for both old and new
    long old_colsz  = old_sig->getDim(dim);
    long new_colsz  = new_sig->getDim(dim);

    Vec<zz_pE> old_col(INIT_SIZE, old_colsz);
    zz_pEX old_col_as_poly;
    Vec<zz_pE> new_col(INIT_SIZE, new_colsz);

    eval_sequence[dim].SetLength(new_sig->getSize());

    for (long i = 0; i < nslices; i++) {
      for (long j = 0; j < ncols; j++) {
        // extract old column
        for (long k = 0; k < old_colsz; k++) 
          old_col[k] = eval_sequence[dim+1][i*old_colsz*ncols + j + k*ncols];

        // convert old column to a polynomial
        conv(old_col_as_poly, old_col);

        // compute new column
        for (long k = 0; k < new_colsz; k++)
          new_col[k] = eval(old_col_as_poly, local_points[dim][k]);

        // insert new column
        for (long k = 0; k < new_colsz; k++)
          eval_sequence[dim][i*new_colsz*ncols + j + k*ncols] = new_col[k];
      }
    }
  }

  if (global_values == eval_sequence[0]) 
    cout << "I win!!\n";
  else {
    cout << "I lose\n";
    cout << global_values << "\n";
    cout << eval_sequence[0] << "\n";
  }

  Vec<long> slot_index, slot_rotate;
  init_slot_mappings(slot_index, slot_rotate, global_reps, m, p, context);

  zz_pE H = conv<zz_pE>(zz_pX(p, 1));

  vector<ZZX> adjusted_values;
  adjusted_values.resize(nslots);

  for (long i = 0; i < nslots; i++) {
    zz_pE V = global_values[i];
    long h = slot_rotate[i];
    for (long j = 0; j < h; j++) 
      V = conv<zz_pE>(CompMod(rep(V), rep(H), G));
    
    adjusted_values[ slot_index[i] ] = conv<ZZX>(rep(V));
  }

  EncryptedArray ea(context, GG);

  ZZX FF1;
  ea.encode(FF1, adjusted_values);
  
  zz_pX F1 = conv<zz_pX>(FF1);

  if (F1 == F) 
    cout << "yes!!\n";
  else 
    cout << "NO!!!\n";



  for (long dim = 0; dim < inertPrefix; dim++) {
    PlaintextMatrixBaseInterface *mat = 
      buildStep2Matrix(ea, *sig_sequence[dim], local_reps[dim], dim, m/mvec[dim]);

    PlaintextMatrixBaseInterface *imat = 
      buildStep2Matrix(ea, *sig_sequence[dim], local_reps[dim], dim, m/mvec[dim], true);

    
    vector<ZZX> val1;
    val1.resize(nslots);
    for (long i = 0; i < nslots; i++) 
      val1[i] = conv<ZZX>(rep(eval_sequence[dim+1][i]));

    PlaintextArray pa1(ea);
    pa1.encode(val1);
    PlaintextArray pa1_orig(pa1);

    pa1.mat_mul(*mat);

    vector<ZZX> val2;
    val2.resize(nslots);
    for (long i = 0; i < nslots; i++) 
      val2[i] = conv<ZZX>(rep(eval_sequence[dim][i]));

    PlaintextArray pa2(ea);
    pa2.encode(val2);

    if (pa1.equals(pa2))
      cout << "dim=" << dim << " GOOD\n";
    else
      cout << "dim=" << dim << " BAD\n";

    pa1.mat_mul(*imat);
    if (pa1.equals(pa1_orig))
      cout << "dim=" << dim << " INV GOOD\n";
    else
      cout << "dim=" << dim << " INV BAD\n";

  }

  if (inertPrefix == nfactors-1) {
    cout << "easy case\n";

    long dim = nfactors-1;
    PlaintextBlockMatrixBaseInterface *mat = 
      buildStep1Matrix(ea, *sig_sequence[dim], local_reps[dim], dim, m/mvec[dim]);

    PlaintextBlockMatrixBaseInterface *imat = 
      buildStep1Matrix(ea, *sig_sequence[dim], local_reps[dim], dim, m/mvec[dim], true);


    vector<ZZX> val1;
    val1.resize(nslots);
    for (long i = 0; i < phim; i++) {
      val1[i/d] += conv<ZZX>(rep(eval_sequence[dim+1][i])) << (i % d);
    }
    PlaintextArray pa1(ea);
    pa1.encode(val1);
    PlaintextArray pa1_orig(pa1);

    pa1.mat_mul(*mat);

    vector<ZZX> val2;
    val2.resize(nslots);
    for (long i = 0; i < nslots; i++) 
      val2[i] = conv<ZZX>(rep(eval_sequence[dim][i]));

    PlaintextArray pa2(ea);
    pa2.encode(val2);

    if (pa1.equals(pa2))
      cout << "dim=" << dim << " GOOD\n";
    else
      cout << "dim=" << dim << " BAD\n";


    pa1.mat_mul(*imat);
    if (pa1.equals(pa1_orig))
      cout << "dim=" << dim << " INV GOOD\n";
    else
      cout << "dim=" << dim << " INV BAD\n";

  }
  else if (inertPrefix == nfactors-2) {
    cout << "harder case\n";

    long m1 = mvec[nfactors-1];
    long m2 = m/m1;

    long d1 = dvec[nfactors-1];
    long d2 = d/d1;

    Tower<PA_zz_p> tower(m1, m2, d1, d2, p, r);

    zz_pX g;
    random(g, d1);
    zz_pE beta = eval(g, tower.zeta) * conv<zz_pE>(zz_pX(1, 1));

    cout << g << "\n";
    cout << tower.convert1to2(beta) << "\n";

    long dim = nfactors-1;


    PlaintextBlockMatrixBaseInterface *mat =
      buildStep1aMatrix(ea, local_reps[dim], m1, m2, d1, d2, phivec[dim]);

    vector<ZZX> val1;
    val1.resize(nslots);
    for (long i = 0; i < phim; i++) {
      val1[i/d] += conv<ZZX>(rep(eval_sequence[dim+1][i])) << (i % d);
    }
    PlaintextArray pa1(ea);
    pa1.encode(val1);
    PlaintextArray pa1_orig(pa1);

    pa1.mat_mul(*mat);

    assert(eval_sequence[dim].length() == nslots*d2);

    vector<ZZX> val2;
    val2.resize(nslots);
    for (long i = 0; i < nslots; i++) {
      Vec<zz_pX> one_slot;
      one_slot.SetLength(d2);
      for (long j = 0; j < d2; j++) {
        Vec<zz_pX> v = tower.convert1to2(eval_sequence[dim][i*d2 + j]);
        for (long k = 1; k < d2; k++) assert(v[k] == 0);
        one_slot[j] = v[0];
      }
      val2[i] = conv<ZZX>(rep(tower.convert2to1(one_slot)));
    }

    PlaintextArray pa2(ea);
    pa2.encode(val2);

    if (pa1.equals(pa2))
      cout << "dim=" << dim << " GOOD\n";
    else
      cout << "dim=" << dim << " BAD\n";
  }
  else {
    cout << "case not handled\n";

  }
   

}




void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'R=1 p=2 k=80'\n\n";
  cerr << "  R is the number of rounds\n";
  cerr << "  p is the plaintext base [default=2]" << endl;
  cerr << "  r is the lifting [default=1]" << endl;
  cerr << "  d is the degree of the field extension [default==0]\n";
  cerr << "    (d == 0 => factors[0] defined the extension)\n";
  cerr << "  c is number of columns in the key-switching matrices [default=2]\n";
  cerr << "  k is the security parameter [default=80]\n";
  cerr << "  L is the # of primes in the modulus chai [default=4*R]\n";
  cerr << "  s is the minimum number of slots [default=4]\n";
  cerr << "  m defined the cyclotomic polynomial Phi_m(X)\n";
  cerr << "  seed is the PRG seed\n";
  exit(0);
}


int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["R"] = "1";
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["c"] = "2";
  argmap["k"] = "80";
  argmap["L"] = "0";
  argmap["s"] = "0";
  argmap["m1"] = "0";
  argmap["m2"] = "0";
  argmap["m3"] = "0";
  argmap["m4"] = "0";
  argmap["seed"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long R = atoi(argmap["R"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long c = atoi(argmap["c"]);
  long k = atoi(argmap["k"]);
  //  long z = atoi(argmap["z"]);
  long L = atoi(argmap["L"]);
  if (L==0) { // determine L based on R,r
    if (r==1) L = 2*R+2;
    else      L = 4*R;
  }
  long s = atoi(argmap["s"]);

  long m1 = atoi(argmap["m1"]);
  long m2 = atoi(argmap["m2"]);
  long m3 = atoi(argmap["m3"]);
  long m4 = atoi(argmap["m4"]);
  long seed = atoi(argmap["seed"]);

  long w = 64; // Hamming weight of secret key
  //  long L = z*R; // number of levels

  Vec<long> mvec;
  if (m1 != 0) append(mvec, m1);
  if (m2 != 0) append(mvec, m2);
  if (m3 != 0) append(mvec, m3);
  if (m4 != 0) append(mvec, m4);
  

  if (seed) SetSeed(conv<ZZ>(seed));

  TestIt(R, p, r, c, k, w, L, mvec);

  cerr << endl;
  printAllTimers();
  cerr << endl;

}

//   [1 1 3 8] Test_eval3_x p=2 m1=3 m2=5 m3=7 m4=17
//   [1 1 20]  Test_eval3_x p=2 m1=3 m2=11 m3=25
//   [1 1 20]  Test_eval3_x p=2 m1=3 m2=11 m3=41


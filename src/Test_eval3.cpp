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

#include <NTL/lzz_pXFactoring.h>
NTL_CLIENT

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
    cout << "special case\n";

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


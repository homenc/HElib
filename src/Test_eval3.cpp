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
class Step2Matrix : public PlaintextMatrixInterface<type> 
{
public:
  PA_INJECT(type)

private:
  const EncryptedArray& ea;
  CubeSignature sig;
  long dim;

  Mat<RX> A;

  // constructor
  Step2Matrix(const EncryptedArray& _ea, 
              const Vec<long>& dimvec,
              const Vec<long>& reps,
              long _dim,
              long cofactor);
              
};

template<class type>
Step2Matrix<type>::Step2Matrix(const EncryptedArray& _ea, 
                               const Vec<long>& dimvec,
                               const Vec<long>& reps,
                               long _dim,
                               long cofactor)
: ea(_ea), sig(dimvec), dim(_dim)
{
   RBak bak; bak.save(); ea.getContext().alMod.restoreContext();
   const RX& G = ea.getDerived(type()).getG();

   long sz = dimvec[dim];
   assert(sz == reps.length());


   Vec<RX> points;
   points.SetLength(sz);
   for (long j = 0; j < sz; j++) 
      points[j] = zz_pX(reps[j]*cofactor, 1) % G;

   A.SetDims(sz, sz);
   for (long j = 0; j < sz; j++)
      A[0][j] = 1;

   for (long i = 1; i < sz; i++)
      for (long j = 0; j < sz; j++)
         A[i][j] = (A[i-1][j] * points[j]) % G;
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

  exit(0);

/*
for i in [0..sig.getProd(0, dim))  // iterate over all dimension dim subcubes
   offset = i*sig.getProd(dim)
   for j in [0..sig.getProd(dim+1))  // iterate over all columns in the subcube
      for k = 0..sig.getDim(dim)   // iterate over the column elements
         index = offset + j + k*sig.getProd(dim+1)
*/


#if 0
    
  

  long nslots = context.zMStar.getNSlots();


  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key

  assert(d == 0); // that's the assumption here


  if (d == 0)
  else
    GG = makeIrredPoly(p, d); 

  d = deg(GG);

  cerr << "GG = " << GG << "\n";
  cerr << "generating key-switching matrices... ";
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  cerr << "done\n";

  printAllTimers();
  resetAllTimers();



  cerr << "computing masks and tables for rotation...";
  EncryptedArray ea(context, GG);
  cerr << "done\n";


  // build a Step1 matrix

  PlaintextBlockMatrixBaseInterface *step1 =
    buildStep1Matrix(ea, m1, m2);

  long gen1 = dynamic_cast<Step1MatrixSuper*>(step1)->getGen();

  Vec<long> rep1;
  alt_init_representatives(rep1, m1, gen1, phim1/d);

  // build a Step2 matrix

  PlaintextMatrixBaseInterface *step2 =
    buildStep2Matrix(ea, m1, m2);

  long gen2 = dynamic_cast<Step2MatrixSuper*>(step2)->getGen();

  Vec<long> rep2;
  alt_init_representatives(rep2, m2, gen2, phim2);


  Vec<long> representatives;
  for (long i = 0; i < rep1.length(); i++)
    for (long j = 0; j < rep2.length(); j++) {
      // chinese remaindering
      long x1 = rep1[i];
      long x2 = rep2[j];

      long x = mcMod(x1*m2*InvMod(m2, m1) + x2*m1*InvMod(m1, m2), m);

      append(representatives, x);
    }


  cout << representatives << "\n";

  Vec<long> slot_index, slot_rotate;

  init_slot_mappings(slot_index, slot_rotate, representatives, m, p, context);

  cout << slot_index << "\n";
  cout << slot_rotate << "\n";


  // evaluate F via the cube structure

  Vec<zz_pX> points1;
  points1.SetLength(phim1/d);
  for (long i = 0; i < phim1/d; i++)
    points1[i] = zz_pX(m2 * rep1[i], 1) % G;

  Vec<zz_pX> points2;
  points2.SetLength(phim2);
  for (long j = 0; j < phim2; j++)
    points2[j] = zz_pX(m1 * rep2[j], 1) % G;


  // cube represents an array with phim1 rows and phim2 colums,
  //   stored in row-major order.
  //   So, the first phim2 elements elements represent row 0,
  //   the next phim2 elements represent row 1, etc.
  //   Column i represents the polynomial f_i in the powerful 
  //   basis represenation: sum_i f_i(X_1) X_2^i, where i = 0..phim2-1
  
  Vec<zz_pX> alt_values1;
  alt_values1.SetLength(nslots);

  for (long j = 0; j < phim2; j++) {
    zz_pX tpoly;
    tpoly.SetLength(phim1);
    for (long i = 0; i < phim1; i++)
      tpoly[i] = cube[i*phim2 + j];
    tpoly.normalize();

    for (long i = 0; i < phim1/d; i++)
      alt_values1[i*phim2 + j] = CompMod(tpoly, points1[i], G);
  }

  Vec<zz_pX> alt_values2;
  alt_values2.SetLength(nslots);

  for (long i = 0; i < phim1/d; i++) {
    Vec<zz_pX> tpoly;
    tpoly.SetLength(phim2);
    for (long j = 0; j < phim2; j++)
      tpoly[j] = alt_values1[i*phim2 + j];

    for (long j = 0; j < phim2; j++)
      alt_values2[i*phim2 + j ] = EvalMod(tpoly, points2[j], G);
  }

  if (eval_values == alt_values2) 
    cout << "right on!!\n";
  else
    cout << "no way!!!\n";

  // evaluate using matrices

#if 0

  PlaintextArray pa(ea);
  vector<ZZX> Alt_values1;
  convert(Alt_values1, alt_values1);

  pa.encode(Alt_values1);
  pa.mat_mul(*step2);

  vector<ZZX> AV3;
  pa.decode(AV3);

  Vec<zz_pX> av3;
  convert(av3, AV3);

  if (eval_values == av3) 
    cout << "too cool!!\n";
  else
    cout << "nooooo!!!\n";
#endif


#if 1
  Vec<zz_pX> av2;
  av2.SetLength(nslots);

  for (long i = 0; i < nslots; i++) {
    long u = i / phim2;
    long v = i % phim2;

    // pack coeffs [ud..ud+d) of f_v into slot i
    av2[i] = 0;
    for (long j = u*d; j < u*d+d; j++) {
      av2[i] += zz_pX(j - u*d, cube[j*phim2 + v]);
    }
  }


  vector<ZZX> AV2;
  convert(AV2, av2);

  PlaintextArray pa(ea);

  pa.encode(AV2);
  PlaintextArray pa_orig = pa;

  pa.mat_mul(*step1);

  PlaintextArray pa_int = pa;

  pa.mat_mul(*step2);

  vector<ZZX> AV3;
  pa.decode(AV3);

  Vec<zz_pX> av3;
  convert(av3, AV3);

  if (eval_values == av3) 
    cout << "too cool!!\n";
  else
    cout << "nooooo!!!\n";

  

  PlaintextBlockMatrixBaseInterface *step1inv = 
    buildStep1Inverse(step1);

  PlaintextMatrixBaseInterface *step2inv =
    buildStep2Inverse(step2);

  pa.mat_mul(*step2inv);

  if (pa.equals(pa_int)) 
    cout << "step2 inverse OK\n";
  else
    cout << "step2 inverse NOT OK\n";

  pa.mat_mul(*step1inv);

  if (pa.equals(pa_orig)) 
    cout << "inverse OK\n";
  else
    cout << "inverse NOT OK\n";


#endif


  // now verify adjustments

  zz_pX H = zz_pX(p, 1) % G;

  vector<ZZX> adjusted_values;
  adjusted_values.resize(nslots);

  for (long i = 0; i < nslots; i++) {
    zz_pX V = eval_values[i];
    long h = slot_rotate[i];
    for (long j = 0; j < h; j++) 
      V = CompMod(V, H, G);
    
    adjusted_values[ slot_index[i] ] = conv<ZZX>(V);
  }

  ZZX FF1;
  ea.encode(FF1, adjusted_values);
  
  zz_pX F1 = conv<zz_pX>(FF1);

  if (F1 == F) 
    cout << "yes!!\n";
  else 
    cout << "NO!!!\n";

#endif


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

// a nice test case: Test_eval3_x p=2 m1=3 m2=5 m3=7 m4=17

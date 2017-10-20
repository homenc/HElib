/* Copyright (C) 2012-2017 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

/**
 * @file matmul.h
 * @brief some matrix / linear algenra stuff
 */
#include <cstddef>
#include <tuple>
#include "EncryptedArray.h"
#include "matmul.h"


/********************************************************************/
/****************** Linear transformation classes *******************/



class MatMul_new {
public:
  virtual ~MatMul_new() {}
  virtual const EncryptedArray& getEA() const = 0;
};

template<class type>
class MatMul_derived : public MatMul_new { 
public:
  PA_INJECT(type)

  // Should return true when the entry is a zero. 
  virtual bool get(RX& out, long i, long j) const = 0;

};

class MatMul1D {
public:
  virtual ~MatMul1D() {}
  virtual const EncryptedArray& getEA() const = 0;
  virtual long getDim() const = 0;
  virtual bool multipleTransforms() const = 0;
};

template<class type>
class MatMul1D_derived : public MatMul1D { 
public:
  PA_INJECT(type)

  // Should return true when the entry is a zero. 
  virtual bool get(RX& out, long i, long j, long k) const = 0;
};


class ConstMultiplier {
// stores a constant in either zzX or DoubleCRT format

public:

  virtual ~ConstMultiplier() {}

  virtual void mul(Ctxt& ctxt) const = 0;

  virtual ConstMultiplier *upgrade(const FHEcontext& context) const = 0;
  // Upgrade to DCRT. Returns null of no upgrade required

};

class ConstMultiplier_DoubleCRT : public ConstMultiplier {
private:

  DoubleCRT data;

public:
  ConstMultiplier_DoubleCRT(const DoubleCRT& _data) : data(_data) { }

  void mul(Ctxt& ctxt) const override {
    ctxt.multByConstant(data);
  } 

  ConstMultiplier *upgrade(const FHEcontext& context) const override {
    return nullptr;
  }

};


class ConstMultiplier_zzX : public ConstMultiplier {
private:

  zzX data;

public:

  ConstMultiplier_zzX(const zzX& _data) : data(_data) { }

  void mul(Ctxt& ctxt) const override {
    ctxt.multByConstant(data);
  } 

  ConstMultiplier *upgrade(const FHEcontext& context) const override {
    return new ConstMultiplier_DoubleCRT(DoubleCRT(data, context));
  }

};



class ConstMultiplierCache {
public:
  vector<shared_ptr<ConstMultiplier>> multiplier;

  void upgrade(const FHEcontext& context) {
    for (auto&& ptr: multiplier) 
      if (ptr) 
        if (auto newptr = ptr->upgrade(context)) 
          ptr = shared_ptr<ConstMultiplier>(newptr); 
  }

};

class MatMul1DExec {
public:

  const EncryptedArray& ea;
  long dim;
  long D;
  bool native;
  long g;

  ConstMultiplierCache cache;
  ConstMultiplierCache cache1; // only for non-native dimension


  MatMul1DExec(const MatMul1D& mat);
  void mul(Ctxt& ctxt);
  void upgrade() { 
    cache.upgrade(ea.getContext()); 
    cache1.upgrade(ea.getContext()); 
  }
};


static inline long dimSz(const EncryptedArray& ea, long dim)
{
   return (dim==ea.dimension())? 1 : ea.sizeOfDimension(dim);
}

static inline long dimSz(const EncryptedArrayBase& ea, long dim)
{
   return (dim==ea.dimension())? 1 : ea.sizeOfDimension(dim);
}

static inline long dimNative(const EncryptedArray& ea, long dim)
{
   return (dim==ea.dimension())? true : ea.nativeDimension(dim);
}

static inline long dimNative(const EncryptedArrayBase& ea, long dim)
{
   return (dim==ea.dimension())? true : ea.nativeDimension(dim);
}


template<class type>
struct MatMul1DExec_construct {
  PA_INJECT(type)

  static
  void processDiagonal1(RX& poly, long i, long rotAmt,
                        const EncryptedArrayDerived<type>& ea,
                        const MatMul1D_derived<type>& mat)
  {
    long dim = mat.getDim();
    long D = dimSz(ea, dim);

    vector<RX> tmpDiag(D);
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry
    RX entry;

    // Process the entries in this diagonal one at a time
    for (long j = 0; j < D; j++) { // process entry j
      long rotJ = (j+rotAmt) % D;  // need to rotate constant by rotAmt
      bool zEntry = mat.get(entry, mcMod(rotJ-i, D), rotJ, 0); 
        // entry [j-i mod D, j]

      assert(zEntry || deg(entry) < ea.getDegree());
      // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry = true;// zero is an empty entry too

      if (!zEntry) {   // not a zero entry
        zDiag = false; // mark diagonal as non-empty

        // clear entries between last nonzero entry and this one
        for (long jj = nzLast+1; jj < j; jj++) clear(tmpDiag[jj]);
        tmpDiag[j] = entry;
        nzLast = j;
      }
    }    
    if (zDiag) {
      clear(poly);
    } 
    else {

      // clear trailing zero entries
      for (long jj = nzLast+1; jj < D; jj++) clear(tmpDiag[jj]);
      
      vector<RX> diag(ea.size());
      if (D==1) 
	diag.assign(ea.size(), tmpDiag[0]); // dimension of size one
      else {
	for (long j = 0; j < ea.size(); j++)
	  diag[j] = tmpDiag[ ea.coordinate(dim,j) ];
	  // rearrange the indexes based on the current dimension
      }

      ea.encode(poly, diag);
    }
  }


  static
  void processDiagonal2(RX& poly, long idx, long rotAmt,
                        const EncryptedArrayDerived<type>& ea,
                        const MatMul1D_derived<type>& mat)
  {
    long dim = mat.getDim();
    long D = dimSz(ea, dim);

    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry
    RX entry;

    long n = ea.size();

    // Process the entries in this diagonal one at a time
    long blockIdx, innerIdx;
    vector<RX> diag(n);
    for (long j=0; j < n; j++) {
      if (D==1) {
	blockIdx=j; 
        innerIdx = 0;
      } 
      else {
	std::tie(blockIdx, innerIdx) // std::pair<long,long> idxes
	  = ea.getContext().zMStar.breakIndexByDim(j, dim);
	//	blockIdx = idxes.first;  // which transformation
	//	innerIdx = idxes.second; // index along dimension dim
        innerIdx = (innerIdx+rotAmt) % D;  // need to rotate constant by rotAmt
      }
      // process entry j
      bool zEntry=mat.get(entry, mcMod(innerIdx-idx,D), innerIdx, blockIdx);
      // entry [i,j-i mod D] in the block corresponding to blockIdx
      // get(...) returns true if the entry is empty, false otherwise

      // If non-zero, make sure the degree is not too large
      assert(zEntry || deg(entry) < ea.getDegree());

      if (!zEntry && IsZero(entry)) zEntry = true; // zero is an empty entry too

      if (!zEntry) {   // not a zero entry
	zDiag = false; // mark diagonal as non-empty

	// clear entries between last nonzero entry and this one
	for (long jj = nzLast+1; jj < j; jj++) clear(diag[jj]);
	nzLast = j;
	diag[j] = entry;
      }
    }    
    if (zDiag) {
      clear(poly);
    }
    else {

      // clear trailing zero entries
      for (long jj = nzLast+1; jj < ea.size(); jj++) clear(diag[jj]);

      ea.encode(poly, diag);
    }
  }

  // Get the i'th diagonal, encoded as a single constant. 
  static
  void processDiagonal(RX& poly, long i, long rotAmt,
                        const EncryptedArrayDerived<type>& ea,
                        const MatMul1D_derived<type>& mat)
  {
    if (mat.multipleTransforms())
      processDiagonal2(poly, i, rotAmt, ea, mat);
    else
      processDiagonal1(poly, i, rotAmt, ea, mat);
  }

  static
  void apply(const EncryptedArrayDerived<type>& ea,
             const MatMul1D& mat_basetype,
             vector<shared_ptr<ConstMultiplier>>& vec,
             vector<shared_ptr<ConstMultiplier>>& vec1,
             long g)
  {

    const MatMul1D_derived<type>& mat =
      dynamic_cast< const MatMul1D_derived<type>& >(mat_basetype);

    long dim = mat.getDim();
    long D = dimSz(ea, dim);
    bool native = dimNative(ea, dim);

    RBak bak; bak.save(); ea.getTab().restoreContext();

    vec.resize(D);
    if (!native) vec1.resize(D);

    for (long i = 0; i < D; i++) {
      // i == j + g*k
      long j = i % g;
      long k = i / g;

      // long rotAmt = g*k;
      // This assumes we process baby steps first, then giant steps.
      // For the reverse, set rotAmt = j

      RX poly;
      processDiagonal(poly, i, 0, ea, mat);

      // FIXME: not quite right. Need to multiply by mask constants
      // in non-native dimension

      if (IsZero(poly)) {
        vec[i] = nullptr; 
        if (!native) vec1[i] = nullptr;
      }
      else {
        RX poly1;
        plaintextAutomorph(poly1, poly, dim, -g*k, ea); 
        vec[i] = shared_ptr<ConstMultiplier>(new ConstMultiplier_zzX(convert<zzX>(poly1)));

        if (!native) {
          RX poly2;
          plaintextAutomorph(poly1, poly, dim, D, ea); 
          vec1[i] = shared_ptr<ConstMultiplier>(new ConstMultiplier_zzX(convert<zzX>(poly2)));
        }
      }
    }
  }






};



MatMul1DExec::MatMul1DExec(const MatMul1D& mat)
  : ea(mat.getEA())
{
    dim = mat.getDim();
    assert(dim >= 0 && dim <= ea.dimension());
    D = dimSz(ea, dim);
    native = dimNative(ea, dim);
    g = KSGiantStepSize(D);

    ea.dispatch<MatMul1DExec_construct>(mat, Fwd(cache.multiplier), g);
}


/***************************************************************************

BS/GS logic:

  \sum_{i=0}^{D-1} const_i rot^i(v)
    = \sum_k \sum_j const_{j+g*k} rot^{j+g*k}(v)
    = \sum_k rot^{g*k}[ \sum_j rot^{-g*k}(const_{j+g*k}) rot^j(v) ]

So we first compute baby_steps[j] = rot^j(v) for j in [0..g).
Then for each k in [0..ceil(D/g)), we compute 
   giant_steps[k] = \rot^{g*k}[ rot^{-g*k}(const_{j+g*k}) baby_steps[j] ] 
Then we add up all the giant_steps.

In bad dimesnions:

We need to compute
\[
  \sum_{j,k} c_{j+gk} r^{j+gk}(x)
\]
where $r^i$ denotes rotation by $i$.
In bad dimensions, we have
\[
 r^i(x) = d_i \rho^i(x) + e_i \rho^{i-D}(x)
\]
for constants $d_i$ and $e_i$.
Here, d_i is maskTable[i][amt] and e_i = 1-d_i

So putting it all together
\[
  \sum_{j,k} c_{j+gk} r^{j+gk}(x)
= \sum_{j,k} d'_{j+gk} \rho^{j+gk}(x) + e'_{j+gk} \rho^{j+gk-D}(x) 
     \text{where $d'_i=c_i d_i$ and $e'_i = c_i e_i$}


=               \sum_k \rho^{gk}[ \sum_j d''_{j+gk} \rho^j(x) ]
   + \rho^{-D}[ \sum_k \rho^{gk}[ \sum_j e''_{j+gk} \rho^j(x) ] ]
      \text{where $d''_{j+gk} = \rho^{-gk}(d'_{j+gk})$ and
                  $e''_{j+gk} = \rho^{D-gk}(d'_{j+gk})$}
 
\]

***************************************************************************/


void
MatMul1DExec::mul(Ctxt& ctxt)
{
   ctxt.cleanUp();

   long nintervals = divc(D, g);

   vector<Ctxt> baby_steps(g, ctxt);

   // FIXME: use parallel for loop
   for (long j = 1; j < g; j++) {
     ea.rotate1D(baby_steps[j], dim, j);
     baby_steps[j].cleanUp();
   }

   Ctxt acc(ZeroCtxtLike, ctxt);

   // FIXME: use parallel for loop
   for (long k = 0; k < nintervals; k++) {
      Ctxt acc1(ZeroCtxtLike, ctxt);

      for (long j = 0; j < g; j++) {
         long i = j + g*k;
         if (i >= D) break;
         if (cache.multiplier[i]) {
            Ctxt tmp(baby_steps[j]);
            cache.multiplier[i]->mul(tmp);
            acc1 += tmp;
         }
      }

      if (k > 0) ea.rotate1D(acc1, dim, g*k);
      acc += acc1;
   }

   ctxt = acc;
}

template<class type> class RandomMatrix_new : public  MatMul1D_derived<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< RX > > data;
  const EncryptedArray& ea;
  long dim;

public:
  virtual ~RandomMatrix_new() {}
  RandomMatrix_new(const EncryptedArray& _ea, long _dim): 
    ea(_ea), dim(_dim)
  {
    RBak bak; bak.save(); ea.getAlMod().restoreContext();
    long n = ea.size();
    long d = ea.getDegree();
    long D = ea.sizeOfDimension(dim);

    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(D);
    for (long i = 0; i < D; i++) {
      data[i].resize(D);
      for (long j = 0; j < D; j++) {
        random(data[i][j], d);
      }
    }
  }

  const EncryptedArray& getEA() const override { return ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }

  bool get(RX& out, long i, long j, long k) const override {
    long D = ea.sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};


static MatMul1D*
buildRandomMatrix_new(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMatrix_new<PA_GF2>(ea, dim);
    }
    case PA_zz_p_tag: {
      return new RandomMatrix_new<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}


template<class type> class RandomMatrix : public  MatMul<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< RX > > data;
  long dim;

public:
  virtual ~RandomMatrix() {}
  RandomMatrix(const EncryptedArray& _ea, long _dim, long g): 
    MatMul<type>(_ea,g), dim(_dim)
  {
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long n = _ea.size();
    long d = _ea.getDegree();
    long D = _ea.sizeOfDimension(dim);

    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(D);
    for (long i = 0; i < D; i++) {
      data[i].resize(D);
      for (long j = 0; j < D; j++) {
        random(data[i][j], d);
      }
    }
  }

  virtual bool get(RX& out, long i, long j) const {
    long D = this->getEA().sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};


static MatMulBase*
buildRandomMatrix(const EncryptedArray& ea, long dim, long giantStep)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMatrix<PA_GF2>(ea, dim, giantStep);
    }
    case PA_zz_p_tag: {
      return new RandomMatrix<PA_zz_p>(ea, dim, giantStep);
    }
    default: return 0;
  }
}


void  TestIt(FHEcontext& context, long g, long dim, bool verbose)
{
  if (verbose) {
    context.zMStar.printout();
    cout << endl;
  }

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64); // A Hamming-weight-w secret key

  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  EncryptedArray ea(context, context.alMod);


  // choose a random plaintext square matrix
  std::unique_ptr< MatMulBase > ptr(buildRandomMatrix(ea,dim,g));
  std::unique_ptr< MatMul1D > ptr_new(buildRandomMatrix_new(ea,dim));

  MatMul1DExec mat_exec(*ptr_new);
  mat_exec.upgrade();

  // choose a random plaintext vector
  NewPlaintextArray v(ea);
  random(ea, v);

  // encrypt the random vector
  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, v);
  Ctxt ctxt2 = ctxt;

  resetAllTimers();

  { FHE_NTIMER_START(AAA_matmul1D);
  //matMul1D(ctxt, *ptr, dim);               // then use it
  mat_exec.mul(ctxt);
  }

  printAllTimers();

  matMul1D(v, *ptr, dim);     // multiply the plaintext vector

  NewPlaintextArray v1(ea);
  ea.decrypt(ctxt, secretKey, v1); // decrypt the ciphertext vector

  if (equals(ea, v, v1))        // check that we've got the right answer
    cout << "Nice!!\n";
  else
    cout << "Grrr@*\n";


}


int main(int argc, char *argv[]) 
{
  ArgMapping amap;

  long m=2047;
  amap.arg("m", m, "defines the cyclotomic polynomial Phi_m(X)");
  long p=2;
  amap.arg("p", p, "plaintext base");
  long r=1;
  amap.arg("r", r,  "lifting");
  long g=2;
  amap.arg("g", g,  "giant-step parameter");
  long L=4;
  amap.arg("L", L, "# of levels in the modulus chain");
  long dim=0;
  amap.arg("dim", dim, "dimension along which to multiply");
  long verbose=0;
  amap.arg("verbose", verbose, "print timing and other info");

  NTL::Vec<long> gens;
  amap.arg("gens", gens, "use specified vector of generators", NULL);
  amap.note("e.g., gens='[562 1871 751]'");
  NTL::Vec<long> ords;
  amap.arg("ords", ords, "use specified vector of orders", NULL);
  amap.note("e.g., ords='[4 2 -4]', negative means 'bad'");

  amap.parse(argc, argv);

  cout << "*** matmul1D: m=" << m
       << ", p=" << p
       << ", r=" << r
       << ", L=" << L
       << ", g=" << g
       << ", dim=" << dim
       // << ", gens=" << gens
       // << ", ords=" << ords
       << endl;

  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);

  setTimersOn();

  FHEcontext context(m, p, r, gens1, ords1);
  buildModChain(context, L, /*c=*/3);

  TestIt(context, g, dim, verbose);
  cout << endl;
  if (0 && verbose) {
    printAllTimers();
    cout << endl;
  }
}

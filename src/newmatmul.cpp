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
#include <NTL/BasicThreadPool.h>


/********************************************************************/
/****************** Auxiliary stuff: should go elsewhere   **********/



// FIXME: this is copied verbatim comes from Ctxt.cpp
// Compute the number of digits that we need and the esitmated
// added noise from switching this ciphertext part.
static std::pair<long, NTL::xdouble>
computeKSNoise(const CtxtPart& p, const FHEPubKey& pubKey, long pSpace)
{
  const FHEcontext& context = p.getContext();
  long nDigits = 0;
  xdouble addedNoise = to_xdouble(0.0);
  double sizeLeft = context.logOfProduct(p.getIndexSet());
  for (size_t i=0; i<context.digits.size() && sizeLeft>0.0; i++) {    
    nDigits++;

    double digitSize = context.logOfProduct(context.digits[i]);
    if (sizeLeft<digitSize) digitSize=sizeLeft;// need only part of this digit

    // Added noise due to this digit is phi(m) *sigma^2 *pSpace^2 *|Di|^2/4, 
    // where |Di| is the magnitude of the digit

    // WARNING: the following line is written just so to prevent overflow
    addedNoise += to_xdouble(context.zMStar.getPhiM()) * pSpace*pSpace
      * xexp(2*digitSize) * context.stdev*context.stdev / 4.0;

    sizeLeft -= digitSize;
  }

  // Sanity-check: make sure that the added noise is not more than the special
  // primes can handle: After dividing the added noise by the product of all
  // the special primes, it should be smaller than the added noise term due
  // to modulus switching, i.e., keyWeight * phi(m) * pSpace^2 / 12

  long keyWeight = pubKey.getSKeyWeight(p.skHandle.getSecretKeyID());
  double phim = context.zMStar.getPhiM();
  double logModSwitchNoise = log((double)keyWeight) 
    +2*log((double)pSpace) +log(phim) -log(12.0);
  double logKeySwitchNoise = log(addedNoise) 
    -2*context.logOfProduct(context.specialPrimes);
  assert(logKeySwitchNoise < logModSwitchNoise);

  return std::pair<long, NTL::xdouble>(nDigits,addedNoise);
}

class BasicAutomorphPrecon {
  Ctxt ctxt;
  NTL::xdouble noise;
  std::vector<DoubleCRT> polyDigits;

public:
  BasicAutomorphPrecon(const Ctxt& _ctxt) : ctxt(_ctxt)
  {
    FHE_TIMER_START;
    ctxt.cleanUp();

    const FHEcontext& context = ctxt.getContext();
    const FHEPubKey& pubKey = ctxt.getPubKey();
    long keyID = ctxt.getKeyID();

    // The call to cleanUp() should ensure that these assertions pass.

    assert(ctxt.parts.size() == 2); 
    assert(ctxt.parts[0].skHandle.isOne());
    assert(ctxt.parts[1].skHandle.isBase(keyID));
    assert(ctxt.getPrimeSet().disjointFrom(context.specialPrimes));
    

    // Compute the number of digits that we need and the esitmated
    // added noise from switching this ciphertext.
    long nDigits;
    std::tie(nDigits, noise)
      = computeKSNoise(ctxt.parts[1], pubKey, pubKey.keySWlist().at(0).ptxtSpace);

    double logProd = context.logOfProduct(context.specialPrimes);
    noise += ctxt.getNoiseVar() * xexp(2*logProd);

    // Break the ciphertext part into digits, if needed, and scale up these
    // digits using the special primes.

    ctxt.parts[1].breakIntoDigits(polyDigits, nDigits);
  }

  
  shared_ptr<Ctxt>
  automorph(long k) const
  {
    FHE_TIMER_START;

    if (k == 1) return make_shared<Ctxt>(ctxt);

    const FHEcontext& context = ctxt.getContext();
    const FHEPubKey& pubKey = ctxt.getPubKey();
    long keyID = ctxt.getKeyID();

    assert(pubKey.haveKeySWmatrix(1,k,keyID,keyID));
    const KeySwitch& W = pubKey.getKeySWmatrix(1,k,keyID,keyID);

    shared_ptr<Ctxt> result = make_shared<Ctxt>(Ctxt(ZeroCtxtLike, ctxt));

    result->noiseVar = noise; // noise estimate


    // Add in the constant part
    CtxtPart tmpPart = ctxt.parts[0];
    tmpPart.automorph(k);
    tmpPart.addPrimesAndScale(context.specialPrimes);
    result->addPart(tmpPart, /*matchPrimeSet=*/true);

    // "rotate" the digits before key-switching them
    vector<DoubleCRT> tmpDigits = polyDigits;
    for (auto&& tmp: tmpDigits) // rotate each of the digits
      tmp.automorph(k);

    result->keySwitchDigits(W, tmpDigits);

    return result;
  }
};


class GeneralAutomorphPrecon {
public:
  virtual ~GeneralAutomorphPrecon() {}

  virtual shared_ptr<Ctxt> operator()(long i) const = 0;

};

class GeneralAutomorphPrecon_UNKNOWN : public GeneralAutomorphPrecon {
private:
  Ctxt ctxt;
  long dim;
  const PAlgebra& zMStar;

public:
  GeneralAutomorphPrecon_UNKNOWN(const Ctxt& _ctxt, long _dim) :
    ctxt(_ctxt), dim(_dim), zMStar(_ctxt.getContext().zMStar)
  {
    ctxt.cleanUp();
  }

  shared_ptr<Ctxt> operator()(long i) const override
  {
    shared_ptr<Ctxt> result = make_shared<Ctxt>(ctxt);

    // guard against i == 0, as dim may be #gens
    if (i != 0) result->smartAutomorph(zMStar.genToPow(dim, i));

    return result;
  }
};

class GeneralAutomorphPrecon_FULL : public GeneralAutomorphPrecon {
private:
  BasicAutomorphPrecon precon;
  long dim;
  const PAlgebra& zMStar;

public:
  GeneralAutomorphPrecon_FULL(const Ctxt& _ctxt, long _dim) :
    precon(_ctxt), dim(_dim), zMStar(_ctxt.getContext().zMStar)
  { }

  shared_ptr<Ctxt> operator()(long i) const override
  {
    return precon.automorph(zMStar.genToPow(dim, i));
  }

};

class GeneralAutomorphPrecon_BSGS : public GeneralAutomorphPrecon {
private:
  long dim;
  const PAlgebra& zMStar;

  long D;
  long g;
  long nintervals;
  vector<shared_ptr<BasicAutomorphPrecon>> precon;

public:
  GeneralAutomorphPrecon_BSGS(const Ctxt& _ctxt, long _dim) :
    dim(_dim), zMStar(_ctxt.getContext().zMStar)
  { 
    D = (dim == -1) ? zMStar.getOrdP() : zMStar.OrderOf(dim);
    g = KSGiantStepSize(D);
    nintervals = divc(D, g);

    BasicAutomorphPrecon precon0(_ctxt);
    precon.resize(nintervals);

    // parallel for k in [0..nintervals)
    NTL_EXEC_RANGE(nintervals, first, last)
      for (long k = first; k < last; k++) {
	shared_ptr<Ctxt> p = precon0.automorph(zMStar.genToPow(dim, g*k));
	precon[k] = make_shared<BasicAutomorphPrecon>(*p);
      }
    NTL_EXEC_RANGE_END
  }

  shared_ptr<Ctxt> operator()(long i) const override
  {
    assert(i >= 0 && i < D);
    long j = i % g;
    long k = i / g;
    // i == j + g*k
    return precon[k]->automorph(zMStar.genToPow(dim, j));
  }

};

shared_ptr<GeneralAutomorphPrecon>
buildGeneralAutomorphPrecon(const Ctxt& ctxt, long dim)
{
  // allow dim == -1 (Frobenius)
  // allow dim == #gens (the dummy generator of order 1)
  assert(dim >= -1 && dim <= long(ctxt.getContext().zMStar.numOfGens()));

  switch (ctxt.getPubKey().getKSStrategy(dim)) {
    case FHE_KSS_BSGS:
      return make_shared<GeneralAutomorphPrecon_BSGS>(ctxt, dim);

    case FHE_KSS_FULL:
      return make_shared<GeneralAutomorphPrecon_FULL>(ctxt, dim);
      
    default:
      return make_shared<GeneralAutomorphPrecon_UNKNOWN>(ctxt, dim);
  }
}


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

class BlockMatMul1D {
public:
  virtual ~BlockMatMul1D() {}
  virtual const EncryptedArray& getEA() const = 0;
  virtual long getDim() const = 0;
  virtual bool multipleTransforms() const = 0;
};

template<class type>
class BlockMatMul1D_derived : public BlockMatMul1D { 
public:
  PA_INJECT(type)

  // Should return true when the entry is a zero. 
  virtual bool get(mat_R& out, long i, long j, long k) const = 0;
};


class ConstMultiplier {
// stores a constant in either zzX or DoubleCRT format

public:

  virtual ~ConstMultiplier() {}

  virtual void mul(Ctxt& ctxt) const = 0;

  virtual shared_ptr<ConstMultiplier> upgrade(const FHEcontext& context) const = 0;
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

  shared_ptr<ConstMultiplier> upgrade(const FHEcontext& context) const override {
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

  shared_ptr<ConstMultiplier> upgrade(const FHEcontext& context) const override {
    return make_shared<ConstMultiplier_DoubleCRT>(DoubleCRT(data, context));
  }

};

shared_ptr<ConstMultiplier> 
build_ConstMultiplier(const GF2X& poly)
{
   return make_shared<ConstMultiplier_zzX>(convert<zzX>(poly));
}

shared_ptr<ConstMultiplier> 
build_ConstMultiplier(const zz_pX& poly)
{
   return make_shared<ConstMultiplier_zzX>(convert<zzX>(poly));
}



class ConstMultiplierCache {
public:
  vector<shared_ptr<ConstMultiplier>> multiplier;

  void upgrade(const FHEcontext& context) {
    FHE_TIMER_START;

    long n = multiplier.size();
    NTL_EXEC_RANGE(n, first, last)
    for (long i: range(first, last)) {
      if (multiplier[i]) 
        if (auto newptr = multiplier[i]->upgrade(context)) 
          multiplier[i] = shared_ptr<ConstMultiplier>(newptr); 
    }
    NTL_EXEC_RANGE_END
  }

};

class MatMul1DExec {
public:

  const EncryptedArray& ea;
  bool minimal;

  long dim;
  long D;
  bool native;
  long g;

  ConstMultiplierCache cache;
  ConstMultiplierCache cache1; // only for non-native dimension


  // If minimal, then it is assumed minimal KS matrices will
  // be present (one for the generator g, and one for g^{-D} 
  // for non-native dimensions).  With this flag set, all BS/GS
  // and parallel strategies area voided.
  explicit
  MatMul1DExec(const MatMul1D& mat, bool minimal=false);

  void mul(Ctxt& ctxt);

  void upgrade() { 
    cache.upgrade(ea.getContext()); 
    cache1.upgrade(ea.getContext()); 
  }
};

class BlockMatMul1DExec {
public:

  const EncryptedArray& ea;

  long dim;
  long D;
  long d;
  bool native;
  long strategy;

  ConstMultiplierCache cache;
  ConstMultiplierCache cache1; // only for non-native dimension


  // If minimal, then it is assumed minimal KS matrices will
  // be present (one for the generator g, and one for g^{-D} 
  // for non-native dimensions).  With this flag set, all BS/GS
  // and parallel strategies area voided.
  explicit
  BlockMatMul1DExec(const BlockMatMul1D& mat);

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

    if (native) {

      vec.resize(D);

      for (long i = 0; i < D; i++) {
	// i == j + g*k
        long j, k;
      
        if (g) {
          j = i % g;
          k = i / g;
        }
        else {
          j = i;
          k = 1;
        }

	RX poly;
	processDiagonal(poly, i, 0, ea, mat);

        // vec[i] = rho_dim^{-g*k}(poly)

	if (IsZero(poly)) {
	  vec[i] = nullptr; 
          continue;
	}

	plaintextAutomorph(poly, poly, dim, -g*k, ea); 
	vec[i] = build_ConstMultiplier(poly);
      }
    }
    else {
      vec.resize(D);
      vec1.resize(D);

      for (long i = 0; i < D; i++) {
	// i == j + g*k
        long j, k;
      
        if (g) {
          j = i % g;
          k = i / g;
        }
        else {
          j = i;
          k = 1;
        }

	RX poly;
	processDiagonal(poly, i, 0, ea, mat);

        if (IsZero(poly)) {
          vec[i] = nullptr;
          vec1[i] = nullptr;
          continue;
        }

        const RX& mask = ea.getTab().getMaskTable()[dim][i];
        const RXModulus& PhimXMod = ea.getTab().getPhimXMod();

        RX poly1, poly2;
        MulMod(poly1, poly, mask, PhimXMod);
        sub(poly2, poly, poly1);

        // poly1 = poly w/ first i slots zeroed out
        // poly2 = poly w/ last D-i slots zeroed out

        // vec[i] = rho_dim^{-g*k}(poly1)

        if (IsZero(poly1)) {
          vec[i] = nullptr;
        }
        else {
	  plaintextAutomorph(poly1, poly1, dim, -g*k, ea); 
	  vec[i] = build_ConstMultiplier(poly1);
	}

        // vec1[i] = rho_dim^{D-g*k}(poly2)

        if (IsZero(poly2)) {
          vec1[i] = nullptr;
        }
        else {
	  plaintextAutomorph(poly2, poly2, dim, D-g*k, ea); 
	  vec1[i] = build_ConstMultiplier(poly2); 
	}
      }
    }
  }






};


#define FHE_BSGS_MUL_THRESH FHE_KEYSWITCH_THRESH
// uses a BSGS multiplication strategy if sizeof(dim) > FHE_BSGS_MUL_THRESH;
// otherwise uses the old strategy (but potentially with hoisting)

// For performace purposes, should not exceed FHE_KEYSWITCH_THRESH
// For testing purposes: 
//    set to 1 to always use BSGS
//    set to infty to never use BSGS            



MatMul1DExec::MatMul1DExec(const MatMul1D& mat, bool _minimal)
  : ea(mat.getEA()), minimal(_minimal)
{
    FHE_NTIMER_START(MatMul1DExec);

    dim = mat.getDim();
    assert(dim >= 0 && dim <= ea.dimension());
    D = dimSz(ea, dim);
    native = dimNative(ea, dim);

    // FIXME: performance tune
    if (D <= FHE_BSGS_MUL_THRESH || minimal)
       g = 0; // do not use BSGS
    else
       g = KSGiantStepSize(D); // use BSGS

    ea.dispatch<MatMul1DExec_construct>(mat, Fwd(cache.multiplier), 
                                        Fwd(cache1.multiplier), g);
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

void GenBabySteps(vector<shared_ptr<Ctxt>>& v, const Ctxt& ctxt, long dim, 
                  bool clean)
{
  long n = v.size();
  assert(n > 0);

  if (n == 1) {
    v[0] = make_shared<Ctxt>(ctxt);
    if (clean) v[0]->cleanUp();
    return;
  }

  const PAlgebra& zMStar = ctxt.getContext().zMStar;

  if (ctxt.getPubKey().getKSStrategy(dim) != FHE_KSS_UNKNOWN) {
    BasicAutomorphPrecon precon(ctxt);

    NTL_EXEC_RANGE(n, first, last)
      for (long j = first; j < last; j++) {
	 v[j] = precon.automorph(zMStar.genToPow(dim, j));
	 if (clean) v[j]->cleanUp();
      }
    NTL_EXEC_RANGE_END
  }
  else {
    Ctxt ctxt0(ctxt);
    ctxt0.cleanUp();
 
    NTL_EXEC_RANGE(n, first, last)
      for (long j = first; j < last; j++) {
	 v[j] = make_shared<Ctxt>(ctxt0);
	 v[j]->smartAutomorph(zMStar.genToPow(dim, j));
	 if (clean) v[j]->cleanUp();
      }
    NTL_EXEC_RANGE_END
  }
  
}

void
MatMul1DExec::mul(Ctxt& ctxt)
{
   assert(&ea.getContext() == &ctxt.getContext());
   const PAlgebra& zMStar = ea.getContext().zMStar;

   ctxt.cleanUp();

   if (g != 0) {
      // baby-step / giant-step

      if (native) {
	 long nintervals = divc(D, g);
	 vector<shared_ptr<Ctxt>> baby_steps(g);
	 GenBabySteps(baby_steps, ctxt, dim, true);

	 PartitionInfo pinfo(nintervals);
	 long cnt = pinfo.NumIntervals();

	 vector<Ctxt> acc(cnt, Ctxt(ZeroCtxtLike, ctxt));

	 // parallel for loop: k in [0..nintervals)
	 NTL_EXEC_INDEX(cnt, index)
	    long first, last;
	    pinfo.interval(first, last, index);

	    for (long k = first; k < last; k++) {
	       Ctxt acc_inner(ZeroCtxtLike, ctxt);

	       for (long j = 0; j < g; j++) {
		  long i = j + g*k;
		  if (i >= D) break;
		  if (cache.multiplier[i]) {
		     Ctxt tmp(*baby_steps[j]);
		     cache.multiplier[i]->mul(tmp);
		     acc_inner += tmp;
		  }
	       }

	       if (k > 0) acc_inner.smartAutomorph(zMStar.genToPow(dim, g*k));
	       acc[index] += acc_inner;
	    }
	 NTL_EXEC_INDEX_END

	 ctxt = acc[0];
	 for (long i = 1; i < cnt; i++)
	    ctxt += acc[i];
      }
      else {
	 long nintervals = divc(D, g);
	 vector<shared_ptr<Ctxt>> baby_steps(g);
	 GenBabySteps(baby_steps, ctxt, dim, true);

	 PartitionInfo pinfo(nintervals);
	 long cnt = pinfo.NumIntervals();

	 vector<Ctxt> acc(cnt, Ctxt(ZeroCtxtLike, ctxt));
	 vector<Ctxt> acc1(cnt, Ctxt(ZeroCtxtLike, ctxt));

	 // parallel for loop: k in [0..nintervals)
	 NTL_EXEC_INDEX(cnt, index)

	    long first, last;
	    pinfo.interval(first, last, index);

	    for (long k = first; k < last; k++) {
	       Ctxt acc_inner(ZeroCtxtLike, ctxt);
	       Ctxt acc_inner1(ZeroCtxtLike, ctxt);

	       for (long j = 0; j < g; j++) {
		  long i = j + g*k;
		  if (i >= D) break;
		  if (cache.multiplier[i]) {
		     Ctxt tmp(*baby_steps[j]);
		     cache.multiplier[i]->mul(tmp);
		     acc_inner += tmp;
		  }
		  if (cache1.multiplier[i]) {
		     Ctxt tmp(*baby_steps[j]);
		     cache1.multiplier[i]->mul(tmp);
		     acc_inner1 += tmp;
		  }
	       }

	       if (k > 0) {
		  acc_inner.smartAutomorph(zMStar.genToPow(dim, g*k));
		  acc_inner1.smartAutomorph(zMStar.genToPow(dim, g*k));
	       }

	       acc[index] += acc_inner;
	       acc1[index] += acc_inner1;
	    }

	 NTL_EXEC_INDEX_END

	 for (long i = 1; i < cnt; i++) acc[0] += acc[i];
	 for (long i = 1; i < cnt; i++) acc1[0] += acc1[i];

	 acc1[0].smartAutomorph(zMStar.genToPow(dim, -D));
	 acc[0] += acc1[0];
	 ctxt = acc[0];
      }
   }
   else if (!minimal) {
      if (native) {
         shared_ptr<GeneralAutomorphPrecon> precon =
            buildGeneralAutomorphPrecon(ctxt, dim);

	 PartitionInfo pinfo(D);
	 long cnt = pinfo.NumIntervals();

	 vector<Ctxt> acc(cnt, Ctxt(ZeroCtxtLike, ctxt));

	 // parallel for loop: i in [0..D)
	 NTL_EXEC_INDEX(cnt, index)
	    long first, last;
	    pinfo.interval(first, last, index);

	    for (long i = first; i < last; i++) {
	       if (cache.multiplier[i]) {
		  shared_ptr<Ctxt> tmp = (*precon)(i);
		  cache.multiplier[i]->mul(*tmp);
		  acc[index] += *tmp;
	       }
	    }
	 NTL_EXEC_INDEX_END

	 ctxt = acc[0];
	 for (long i = 1; i < cnt; i++)
	    ctxt += acc[i];
      }
      else {
         shared_ptr<GeneralAutomorphPrecon> precon =
            buildGeneralAutomorphPrecon(ctxt, dim);

	 PartitionInfo pinfo(D);
	 long cnt = pinfo.NumIntervals();

	 vector<Ctxt> acc(cnt, Ctxt(ZeroCtxtLike, ctxt));
	 vector<Ctxt> acc1(cnt, Ctxt(ZeroCtxtLike, ctxt));

	 // parallel for loop: i in [0..D)
	 NTL_EXEC_INDEX(cnt, index)
	    long first, last;
	    pinfo.interval(first, last, index);

	    for (long i = first; i < last; i++) {
	       if (cache.multiplier[i] || cache1.multiplier[i]) {
		  shared_ptr<Ctxt> tmp = (*precon)(i);
                  shared_ptr<Ctxt> tmp1 = make_shared<Ctxt>(*tmp);
		  if (cache.multiplier[i]) {
                     cache.multiplier[i]->mul(*tmp);
		     acc[index] += *tmp;
                  }
		  if (cache1.multiplier[i]) {
                     cache1.multiplier[i]->mul(*tmp1);
		     acc1[index] += *tmp1;
                  }
	       }
	    }
	 NTL_EXEC_INDEX_END

	 for (long i = 1; i < cnt; i++) acc[0] += acc[i];
	 for (long i = 1; i < cnt; i++) acc1[0] += acc1[i];

	 acc1[0].smartAutomorph(zMStar.genToPow(dim, -D));
	 acc[0] += acc1[0];
	 ctxt = acc[0];
      }
   }
   else /* minimal */ {
      if (native) {
	 Ctxt acc(ZeroCtxtLike, ctxt);
         Ctxt sh_ctxt(ctxt);
 
         for (long i = 0; i < D; i++) {
	    if (i > 0) sh_ctxt.smartAutomorph(zMStar.genToPow(dim, 1));

	    if (cache.multiplier[i]) {
	       Ctxt tmp(sh_ctxt);
	       cache.multiplier[i]->mul(tmp);
	       acc += tmp;
	    }
         }

	 ctxt = acc;
      }
      else {
	 Ctxt acc(ZeroCtxtLike, ctxt);
	 Ctxt acc1(ZeroCtxtLike, ctxt);
         Ctxt sh_ctxt(ctxt);
 
         for (long i = 0; i < D; i++) {
	    if (i > 0) sh_ctxt.smartAutomorph(zMStar.genToPow(dim, 1));

	    if (cache.multiplier[i]) {
	       Ctxt tmp(sh_ctxt);
	       cache.multiplier[i]->mul(tmp);
	       acc += tmp;
	    }
	    if (cache1.multiplier[i]) {
	       Ctxt tmp(sh_ctxt);
	       cache1.multiplier[i]->mul(tmp);
	       acc1 += tmp;
	    }
         }

	 acc1.smartAutomorph(zMStar.genToPow(dim, -D));
	 acc += acc1;
	 ctxt = acc;
      }
   }
}



// ========================== BlockMatMul1D stuff =====================

template<class type>
struct BlockMatMul1DExec_construct {
  PA_INJECT(type)

  // return true if zero
  static
  bool processDiagonal1(vector<RX>& poly, long i, 
                        const EncryptedArrayDerived<type>& ea,
                        const BlockMatMul1D_derived<type>& mat)
  {
    long dim = mat.getDim();
    long D = dimSz(ea, dim);
    long nslots = ea.size();
    long d = ea.getDegree();

    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry

    mat_R entry(INIT_SIZE, d, d);
    std::vector<RX> entry1(d);
    std::vector< std::vector<RX> > tmpDiag(D);

    vector<vector<RX>> diag(nslots);

    // Process the entries in this diagonal one at a time
    for (long j: range(D)) { // process entry j
      bool zEntry = mat.get(entry, mcMod(j-i, D), j, 0); // entry [j-i mod D, j]
      // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry = true;// zero is an empty entry too
      assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));

      if (!zEntry) {   // not a zero entry
        zDiag = false; // mark diagonal as non-empty

	for (long jj: range(nzLast+1, j)) {// clear from last nonzero entry
          tmpDiag[jj].assign(d, RX());
        }
        nzLast = j; // current entry is the last nonzero one

        // recode entry as a vector of polynomials
        for (long k: range(d)) conv(entry1[k], entry[k]);

        // compute the linearlized polynomial coefficients
	ea.buildLinPolyCoeffs(tmpDiag[j], entry1);
      }
    }
    if (zDiag) return true; // zero diagonal, nothing to do

    // clear trailing zero entries
    for (long jj: range(nzLast+1, D)) {
      tmpDiag[jj].assign(d, RX());
    }

    if (D==1) 
       diag.assign(nslots, tmpDiag[0]); // dimension of size one
    else {
      for (long j: range(nslots))
        diag[j] = tmpDiag[ ea.coordinate(dim,j) ];
           // rearrange the indexes based on the current dimension
    }

    // transpose and encode diag to form polys

    vector<RX> slots(nslots);
    poly.resize(d);
    for (long i: range(d)) {
      for (long j: range(nslots)) slots[j] = diag[j][i];
      ea.encode(poly[i], slots);
    }

    return false; // a nonzero diagonal


  }

  // return true if zero
  static
  bool processDiagonal2(vector<RX>& poly, long idx,
                        const EncryptedArrayDerived<type>& ea,
                        const BlockMatMul1D_derived<type>& mat)
  {
    long dim = mat.getDim();
    long D = dimSz(ea, dim);
    long nslots = ea.size();
    long d = ea.getDegree();

    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry

    mat_R entry(INIT_SIZE, d, d);
    std::vector<RX> entry1(d);

    vector<vector<RX>> diag(nslots);

    // Get the slots in this diagonal one at a time
    long blockIdx, rowIdx, colIdx;
    for (long j: range(nslots)) { // process entry j
      if (dim == ea.dimension()) { // "special" last dimenssion of size 1
	rowIdx = colIdx = 0; blockIdx=j;
      } 
      else {
        std::tie(blockIdx, colIdx)
	  = ea.getContext().zMStar.breakIndexByDim(j, dim);
	rowIdx = mcMod(colIdx-idx,D);
      }
      bool zEntry = mat.get(entry,rowIdx,colIdx,blockIdx);
      // entry [i,j-i mod D] in the block corresponding to blockIdx
      // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too
      assert(zEntry ||
             (entry.NumRows() == d && entry.NumCols() == d));

      if (!zEntry) {    // non-empty entry
	zDiag = false;  // mark diagonal as non-empty

	for (long jj: range(nzLast+1, j)) // clear from last nonzero entry
          diag[jj].assign(d, RX());

	nzLast = j; // current entry is the last nonzero one

	// recode entry as a vector of polynomials
	for (long k: range(d)) conv(entry1[k], entry[k]);

        // compute the linearlized polynomial coefficients
	ea.buildLinPolyCoeffs(diag[j], entry1);
      }
    }
    if (zDiag) return true; // zero diagonal, nothing to do

    // clear trailing zero entries
    for (long jj: range(nzLast+1, nslots))
      diag[jj].assign(d, RX());

    // transpose and encode diag to form polys

    vector<RX> slots(nslots);
    poly.resize(d);
    for (long i: range(d)) {
      for (long j: range(nslots)) slots[j] = diag[j][i];
      ea.encode(poly[i], slots);
    }

    return false; // a nonzero diagonal
  }

  // return true if zero
  static
  bool processDiagonal(vector<RX>& poly, long i,
                        const EncryptedArrayDerived<type>& ea,
                        const BlockMatMul1D_derived<type>& mat)
  {
    if (mat.multipleTransforms())
      return processDiagonal2(poly, i, ea, mat);
    else
      return processDiagonal1(poly, i,  ea, mat);

  }


  // Basic logic:
  // We are computing \sum_i \sum_j c_{ij} \sigma^j rot^i(v),
  // where \sigma = frobenius, rot = rotation 
  // For good dimensions, rot = \rho (the basic automorphism),
  // so we need to compute
  //     \sum_i \sum_j c_{ij} \sigma^j \rho^i(v)
  //  =  \sum_j \sigma^j[ \sigma^{-j}(c_{ij}) \rho^i(v) ]

  // For bad dimensions, we have 
  //   rot^i(v) = d_i \rho^i(v) + e_i \rho^{i-D}(v)
  // and so we need to compute
  //   \sum_i \sum_j c_{ij} \sigma^j (d_i \rho^i(v) + e_i \rho^{i-D}(v))
  // =      \sum_j \sigma_j[  \sigma^{-j}(c_{ij}) d_i \rho^i(v) ] +
  //   \rho^{-D}[ \sum_j \sigma_j[ \rho^{D}{\sigma^{-j}(c_{ij}) e_i} \rho^i(v) ] ]


  // strategy == +1 : factor \sigma
  // strategy == -1 : factor \rho
  // strategy ==  0 : no factoring

   

  static
  void apply(const EncryptedArrayDerived<type>& ea,
             const BlockMatMul1D& mat_basetype,
             vector<shared_ptr<ConstMultiplier>>& vec,
             vector<shared_ptr<ConstMultiplier>>& vec1,
             long strategy)
  {
    const BlockMatMul1D_derived<type>& mat =
      dynamic_cast< const BlockMatMul1D_derived<type>& >(mat_basetype);

    long dim = mat.getDim();
    long D = dimSz(ea, dim);
    long d = ea.getDegree();
    bool native = dimNative(ea, dim);

    RBak bak; bak.save(); ea.getTab().restoreContext();

    vector<RX> poly;

    switch (strategy) {
    case +1: // factor \sigma

      if (native) {
        vec.resize(D*d);
        for (long i: range(D)) {
          bool zero = processDiagonal(poly, i, ea, mat);
          if (zero) {
            for (long j: range(d)) vec[i*d+j] = nullptr;
          }
          else {
	    for (long j: range(d)) {
	      plaintextAutomorph(poly[j], poly[j], -1, -j, ea);
	      vec[i*d+j] = build_ConstMultiplier(poly[j]);
	    }
          }
        }
      }
      else {
        vec.resize(D*d);
        vec1.resize(D*d);
        for (long i: range(D)) {
          bool zero = processDiagonal(poly, i, ea, mat);
          if (zero) {
            for (long j: range(d)) {
              vec [i*d+j] = nullptr;
              vec1[i*d+j] = nullptr;
            }
          }
          else {
	    const RX& mask = ea.getTab().getMaskTable()[dim][i];
	    const RXModulus& F = ea.getTab().getPhimXMod();

            for (long j: range(d)) {
              plaintextAutomorph(poly[j], poly[j], -1, -j, ea);

              RX poly1;
              MulMod(poly1, poly[j], mask, F); // poly[j] w/ first i slots zeroed out
              if (IsZero(poly1)) 
                vec[i*d+j] = nullptr;
              else 
	        vec[i*d+j] = build_ConstMultiplier(poly1);

              sub(poly1, poly[j], poly1); // poly[j] w/ last D-i slots zeroed out
              if (IsZero(poly1)) {
                vec1[i*d+j] = nullptr;
              }
              else {
                plaintextAutomorph(poly1, poly1, dim, D, ea);
	        vec1[i*d+j] = build_ConstMultiplier(poly1);
              }
            }
          }
        }
      }
      break;

    case -1: // factor \rho

      if (native) {
        vec.resize(D*d);
        for (long i: range(D)) {
          bool zero = processDiagonal(poly, i, ea, mat);
          if (zero) {
            for (long j: range(d)) vec[i+j*D] = nullptr;
          }
          else {
	    for (long j: range(d)) {
	      plaintextAutomorph(poly[j], poly[j], dim, -i, ea);
	      vec[i+j*D] = build_ConstMultiplier(poly[j]);
	    }
          }
        }
      }
      else {
        vec.resize(D*d);
        vec1.resize(D*d);
        for (long i: range(D)) {
          bool zero = processDiagonal(poly, i, ea, mat);
          if (zero) {
            for (long j: range(d)) {
              vec [i+j*D] = nullptr;
              vec1[i+j*D] = nullptr;
            }
          }
          else {
	    const RX& mask = ea.getTab().getMaskTable()[dim][i];
	    const RXModulus& F = ea.getTab().getPhimXMod();

            for (long j: range(d)) {
              RX poly1, poly2;
              MulMod(poly1, poly[j], mask, F); // poly[j] w/ first i slots zeroed out
              sub(poly2, poly[j], poly1);      // poly[j] w/ last D-i slots zeroed out

              if (IsZero(poly1)) {
                vec[i+j*D] = nullptr;
              }
              else {
                plaintextAutomorph(poly1, poly1, dim, -i, ea);
	        vec[i+j*D] = build_ConstMultiplier(poly1);
              }

              if (IsZero(poly2)) {
                vec1[i+j*D] = nullptr;
              }
              else {
                plaintextAutomorph(poly2, poly2, dim, D-i, ea);
	        vec1[i+j*D] = build_ConstMultiplier(poly2);
              }
            }
          }
        }
      }

      break;

    case 0:
      Error("not implemented");

    default:
      Error("unknown strategy");
    }
      
  }

};


BlockMatMul1DExec::BlockMatMul1DExec(const BlockMatMul1D& mat)
  : ea(mat.getEA())
{
    FHE_TIMER_START;

    dim = mat.getDim();
    assert(dim >= 0 && dim <= ea.dimension());
    D = dimSz(ea, dim);
    d = ea.getDegree();
    native = dimNative(ea, dim);
    
    if (D >= d) 
      strategy = +1;
    else
      strategy = -1;

    ea.dispatch<BlockMatMul1DExec_construct>(mat, Fwd(cache.multiplier), 
                                        Fwd(cache1.multiplier), strategy);
}


void
BlockMatMul1DExec::mul(Ctxt& ctxt)
{
   assert(&ea.getContext() == &ctxt.getContext());
   const PAlgebra& zMStar = ea.getContext().zMStar;

   ctxt.cleanUp();

   long d0, d1;
   long dim0, dim1;

   if (strategy == +1) {
      d0 = D;
      dim0 = dim;
      d1 = d;
      dim1 = -1;
   }
   else {
      d1 = D;
      dim1 = dim;
      d0 = d;
      dim0 = -1;
   }

   const long par_buf_max = 50;

   if (native) {
      vector<Ctxt> acc(d1, Ctxt(ZeroCtxtLike, ctxt));

      shared_ptr<GeneralAutomorphPrecon> precon =
	       buildGeneralAutomorphPrecon(ctxt, dim0);

#if 0
      // This is the original code

      for (long i: range(d0)) {
	 shared_ptr<Ctxt> tmp = (*precon)(i);
	 for (long j: range(d1)) {
	    if (cache.multiplier[i*d1+j]) {
	       Ctxt tmp1(*tmp);
	       cache.multiplier[i*d1+j]->mul(tmp1);
	       acc[j] += tmp1;
	    }
	 }
      }

      Ctxt sum(ZeroCtxtLike, ctxt);
      for (long j: range(d1)) {
	 if (j > 0) acc[j].smartAutomorph(zMStar.genToPow(dim1, j));
	 sum += acc[j];
      }

      ctxt = sum;
#endif

      long par_buf_sz = 1;
      if (AvailableThreads() > 1) 
         par_buf_sz = min(d0, par_buf_max);

      vector<shared_ptr<Ctxt>> par_buf(par_buf_sz);

      for (long first_i = 0; first_i < d0; first_i += par_buf_sz) {
         long last_i = min(first_i + par_buf_sz, d0);

         // for i in [first_i..last_i), generate automorphosm i and store
         // in par_buf[i-first_i]

         NTL_EXEC_RANGE(last_i-first_i, first, last) 
  
            for (long idx: range(first, last)) {
              long i = idx + first_i;
              par_buf[idx] = (*precon)(i);
            }

         NTL_EXEC_RANGE_END

         NTL_EXEC_RANGE(d1, first, last)

            for (long j: range(first, last)) {
               for (long i: range(first_i, last_i)) {
                  // acc[j] += cache.multiplier[i*d1+j]*par_buf[i-first_i] 

                  if (cache.multiplier[i*d1+j]) {
                     Ctxt tmp1(*par_buf[i-first_i]);
                     cache.multiplier[i*d1+j]->mul(tmp1);
                     acc[j] += tmp1;
                  }
               }
            }

         NTL_EXEC_RANGE_END
      }

      par_buf.resize(0); // free-up space

      PartitionInfo pinfo(d1);
      long cnt = pinfo.NumIntervals();

      vector<Ctxt> sum(cnt, Ctxt(ZeroCtxtLike, ctxt));

      // for j in [0..d1)
      NTL_EXEC_INDEX(cnt, index)
         long first, last;
         pinfo.interval(first, last, index);
         for (long j: range(first, last)) {
	    if (j > 0) acc[j].smartAutomorph(zMStar.genToPow(dim1, j));
	    sum[index] += acc[j];
         }
      NTL_EXEC_INDEX_END

      ctxt = sum[0];
      for (long i: range(1, cnt)) ctxt += sum[i];
   }
   else {
      vector<Ctxt> acc(d1, Ctxt(ZeroCtxtLike, ctxt));
      vector<Ctxt> acc1(d1, Ctxt(ZeroCtxtLike, ctxt));

      shared_ptr<GeneralAutomorphPrecon> precon =
	       buildGeneralAutomorphPrecon(ctxt, dim0);

      for (long i: range(d0)) {
	 shared_ptr<Ctxt> tmp = (*precon)(i);
	 for (long j: range(d1)) {
	    if (cache.multiplier[i*d1+j]) {
	       Ctxt tmp1(*tmp);
	       cache.multiplier[i*d1+j]->mul(tmp1);
	       acc[j] += tmp1;
	    }
	    if (cache1.multiplier[i*d1+j]) {
	       Ctxt tmp1(*tmp);
	       cache1.multiplier[i*d1+j]->mul(tmp1);
	       acc1[j] += tmp1;
	    }
	 }
      }

      Ctxt sum(ZeroCtxtLike, ctxt);
      Ctxt sum1(ZeroCtxtLike, ctxt);
      for (long j: range(d1)) {
	 if (j > 0) {
            acc[j].smartAutomorph(zMStar.genToPow(dim1, j));
            acc1[j].smartAutomorph(zMStar.genToPow(dim1, j));
         }
	 sum += acc[j];
	 sum1 += acc1[j];
      }

      sum1.smartAutomorph(zMStar.genToPow(dim, -D));
      sum += sum1;
      ctxt = sum;
   }
}


// ====================================================================

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

template<class type> class RandomMultiMatrix_new : public  MatMul1D_derived<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< vector< RX > > > data;
  const EncryptedArray& ea;
  long dim;

public:
  virtual ~RandomMultiMatrix_new() {}
  RandomMultiMatrix_new(const EncryptedArray& _ea, long _dim): 
    ea(_ea), dim(_dim)
  {
    RBak bak; bak.save(); ea.getAlMod().restoreContext();
    long n = ea.size();
    long d = ea.getDegree();
    long D = ea.sizeOfDimension(dim);

    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(n/D);
    for (long k = 0; k < n/D; k++) {
      data[k].resize(D);
      for (long i = 0; i < D; i++) {
	data[k][i].resize(D);
	for (long j = 0; j < D; j++) {
	  random(data[k][i][j], d);
	}
      }
    }
  }

  const EncryptedArray& getEA() const override { return ea; }
  bool multipleTransforms() const override { return true; }
  long getDim() const override { return dim; }

  bool get(RX& out, long i, long j, long k) const override {
    long n = ea.size();
    long D = ea.sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    assert(k >= 0 && k < n/D);
    if (IsZero(data[k][i][j])) return true;
    out = data[k][i][j];
    return false;
  }
};


static MatMul1D*
buildRandomMultiMatrix_new(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMultiMatrix_new<PA_GF2>(ea, dim);
    }
    case PA_zz_p_tag: {
      return new RandomMultiMatrix_new<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}

template<class type> class RandomMultiMatrix : public  MatMul<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< vector< RX > > > data;
  long dim;

public:
  virtual ~RandomMultiMatrix() {}
  RandomMultiMatrix(const EncryptedArray& _ea, long _dim, long g)
    : MatMul<type>(_ea, g), dim(_dim)
  {
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long n = _ea.size();
    long d = _ea.getDegree();
    long D = _ea.sizeOfDimension(dim);

    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(n/D);
    for (long k = 0; k < n/D; k++) {
      data[k].resize(D);
      for (long i = 0; i < D; i++) {
	data[k][i].resize(D);
	for (long j = 0; j < D; j++) {
	  random(data[k][i][j], d);
	}
      }
    }
  }

  virtual bool multiGet(RX& out, long i, long j, long k) const
  {
    long D = this->getEA().sizeOfDimension(dim);
    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    assert(k >= 0 && k < this->getEA().size()/D);
    if (IsZero(data[k][i][j])) return true;
    out = data[k][i][j];
    return false;
  }
};

static MatMulBase*
buildRandomMultiMatrix(const EncryptedArray& ea, long dim, long giantStep)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMultiMatrix<PA_GF2>(ea, dim, giantStep);
    }
    case PA_zz_p_tag: {
      return new RandomMultiMatrix<PA_zz_p>(ea, dim, giantStep);
    }
    default: return 0;
  }
}



//********************************

template<class type> 
class RandomBlockMatrix : public BlockMatMul<type> {
  PA_INJECT(type) 

  vector< vector< mat_R > > data;
  long dim;

public:
  ~RandomBlockMatrix() { /*cout << "destructor: random matrix\n";*/ }

  RandomBlockMatrix(const EncryptedArray& _ea, long _dim):
    BlockMatMul<type>(_ea), dim(_dim)
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
        data[i][j].SetDims(d, d);
        for (long u = 0; u < d; u++)
          for (long v = 0; v < d; v++) 
            random(data[i][j][u][v]);
      }
    }
  }

  virtual bool get(mat_R& out, long i, long j) const
  {
    const EncryptedArray& ea = this->getEA();
    long D = ea.sizeOfDimension(dim);
    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};

static MatMulBase*
buildRandomBlockMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomBlockMatrix<PA_GF2>(ea,dim);
    }
    case PA_zz_p_tag: {
      return new RandomBlockMatrix<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}

template<class type> 
class RandomBlockMatrix_new : public BlockMatMul1D_derived<type> {
  PA_INJECT(type) 

  const EncryptedArray& ea;
  long dim;

  vector< vector< mat_R > > data;

public:

  RandomBlockMatrix_new(const EncryptedArray& _ea, long _dim):
    ea(_ea), dim(_dim)
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
        data[i][j].SetDims(d, d);
        for (long u = 0; u < d; u++)
          for (long v = 0; v < d; v++) 
            random(data[i][j][u][v]);
      }
    }
  }

  bool get(mat_R& out, long i, long j, long k) const override
  {
    long D = ea.sizeOfDimension(dim);
    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return ea; }
  long getDim() const override { return dim; }
  bool multipleTransforms() const override { return false; }
};

static BlockMatMul1D*
buildRandomBlockMatrix_new(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomBlockMatrix_new<PA_GF2>(ea,dim);
    }
    case PA_zz_p_tag: {
      return new RandomBlockMatrix_new<PA_zz_p>(ea, dim);
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

  // encrypted array with "full slots"
  EncryptedArray ea(context, context.alMod);


  // choose a random plaintext square matrix
  //std::unique_ptr< MatMulBase > ptr(buildRandomMultiMatrix(ea,dim,g));
  //std::unique_ptr< MatMul1D > ptr_new(buildRandomMultiMatrix_new(ea,dim));

  //std::unique_ptr< MatMulBase > ptr(buildRandomMatrix(ea,dim,g));
  //std::unique_ptr< MatMul1D > ptr_new(buildRandomMatrix_new(ea,dim));

  std::unique_ptr< MatMulBase > ptr(buildRandomBlockMatrix(ea,dim));
  std::unique_ptr< BlockMatMul1D > ptr_new(buildRandomBlockMatrix_new(ea,dim));

  resetAllTimers();
  BlockMatMul1DExec mat_exec(*ptr_new);
  mat_exec.upgrade();
  printAllTimers();

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

  //matMulti1D(v, *ptr, dim);     // multiply the plaintext vector
  blockMatMul1D(v, *ptr, dim);     // multiply the plaintext vector

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
  long nt=1;
  amap.arg("nt", nt, "# threads");

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
       << ", nt=" << nt
       // << ", gens=" << gens
       // << ", ords=" << ords
       << endl;

  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);

  SetNumThreads(nt);

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

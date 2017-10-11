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


/********************************************************************/
/****************** Linear transformation classes *******************/



class MatMul {
public:
  virtual ~MatMul() {}
  virtual const EncryptedArray& getEA() const = 0;
};

template<class type>
class MatMul_derived : public MatMul { 
public:
  PA_INJECT(type)

  // Should return true when the entry is a zero. 
  virtual bool get(RX& out, long i, long j) const = 0;

};

class MatMul1D {
public:
  virtual ~MatMul1D() {}
  virtual const EncryptedArray& getEA() const = 0;
  virtual bool multipleTransforms() const = 0;
  virtual long getDim() const = 0;
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

  virtual void mul(Ctxt& ctxt) = 0;

};

class ConstMultiplier_zzX : public ConstMultiplier {
private:

  zzX data;

public:

  ConstMultiplier_zzX(const zzX& _data) : data(_data) { }

  virtual void mul(Ctxt& ctxt) {
    ctxt.multByConstant(data);
  } 

};

class ConstMultiplier_DoubleCRT : public ConstMultiplier {
private:

  DoubleCRT data;

public:
  ConstMultiplier_DoubleCRT(const DoubleCRT& _data) : data(_data) { }

  virtual void mul(Ctxt& ctxt) {
    ctxt.multByConstant(data);
  } 

};


class ConstMultiplierCache {
public:

  vector<shared_ptr<ConstMultiplier>> multiplier;

};

class MatMul1DExec {
public:

  const EncryptedArray& ea;
  long dim;
  long D;

  long g;


  ConstMultiplierCache cache;


  MatMul1DExec(const MatMul1D& data);
  void apply(Ctxt& ctxt);
};


static inline long dimSz(const EncryptedArray& ea, long dim)
{
   return (dim==ea.dimension())? 1 : ea.sizeOfDimension(dim);
}

static inline long dimSz(const EncryptedArrayBase& ea, long dim)
{
   return (dim==ea.dimension())? 1 : ea.sizeOfDimension(dim);
}


template<class type>
struct MatMul1DExec_construct {
  PA_INJECT(type)

  static
  void processDiagonal1(zzX& poly, long i, long rotAmt,
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
  void processDiagonal2(zzX& poly, long idx, long rotAmt,
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
  void processDiagonal(zzX& poly, long i, long rotAmt,
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
             long g)
  {

    const MatMul1D_derived<type>& mat =
      dynamic_cast< const MatMul1D_derived<type>& >(mat_basetype);

    long dim = mat.getDim();
    long D = dimSz(ea, dim);

    RBak bak; bak.save(); ea.getTab().restoreContext();

    vec.resize(D);

    for (long i = 0; i < D; i++) {
      // i == j + g*k
      long j = i % g;
      long k = i / g;

      long rotAmt = g*k;
      // This assumes we process baby steps first, then giant steps.
      // For the reverse, set rotAmt = j

      zzX poly;
      processDiagonal(poly, i, rotAmt, ea, mat);

      if (IsZero(poly))
         vec[i] = nullptr;
      else
         vec[i] = shared_ptr<ConstMultiplier>(new ConstMultiplier_zzX(poly));
    }

  }



    // Process the diagonals in baby-step/giant-step ordering.
    //   sum_{i=0}^{D-1} const_i rot^i(X)
    //   = \sum_{i=0}^{g-1} \sum_{j=0}^{ngs-1} const_{i+g*j} rot^{i+g*j}(X)
    //   = \sum_i rot^i(sum_j rot^{-i}(const_{i+g*j}) rot^{g*j}(X))
    //
    // so for i=0..g-1 we let
    //    Y_i = sum_j rot^{-i}(const_{i+g*j}) rot^{g*j}(X)
    // then compute \sum_{i=0}^{g-1} rot^i(Y_i).
    //
    // Computing the Y_i's, we initialize an accumulator for each Y_i,
    // then compute the rotations X_j = rot^{g*j}(X), j=0,...,ngs-1.
    // Each X_j is multiplied by all the constants rot^{-i}(const_{i+g*j}),
    // i=0,...,g-1, and the i'th product is added to the accumulator for
    // the corresponding Y_i.



};



MatMul1DExec::MatMul1DExec(const MatMul1D& mat)
  : ea(mat.getEA())
{
    dim = mat.getDim();
    assert(dim >= 0 && dim <= ea.dimension());
    D = dimSz(ea, dim);

    g = SqrRoot(D);
    if (g*g < D) g++;  // g = ceiling(sqrt(D))

    ea.dispatch<MatMul1DExec_construct>(mat, Fwd(cache.multiplier), g);
}

void
MatMul1DExec::apply(Ctxt& ctxt)
{
   ctxt.cleanUp();

   long ngs = divc(D, g);

   vector<Ctxt> baby_steps(g, ctxt);

   // FIXME: use parallel for loop
   for (long j = 1; j < g; j++) {
     ea.rotate1D(baby_steps[j], dim, j);
     baby_steps[j].cleanUp();
   }

   vector<Ctxt> giant_steps(ngs, Ctxt(ZeroCtxtLike, ctxt));

   // FIXME: use parallel for loop
   for (long k = 0; k < ngs; k++) {
      for (long j = 0; j < g; j++) {
         long i = j + g*k;
         if (i >= D) break;
         if (cache.multiplier[i]) {
            Ctxt tmp(baby_steps[j]);
            cache.multiplier[i]->mul(tmp);
            giant_steps[k] += tmp;
         }
      }

      if (k > 0) ea.rotate1D(giant_steps[k], dim, g*k);
   }

   Ctxt acc(ZeroCtxtLike, ctxt);
   for (long k = 0; k < ngs; k++) {
      acc += giant_steps[k];
   }

   ctxt = acc;
}


int main() { }

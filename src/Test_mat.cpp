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

/* Test_matmul.cpp - Testing the functionality of multiplying an encrypted
 * vector by a plaintext matrix, either over the extension- or the
 * base-field/ring.
 */

#include <cassert>
#include <NTL/lzz_pXFactoring.h>
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include "matmul.h"


//*******************************************


template<class type>
class mat_multi1D_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const PlaintextMatrixBaseInterface& mat, long dim)
  {
    PA_BOILER

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    const PlaintextMultiMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMultiMatrixInterface<type>& >( mat );

    vector< vector<RX> > data1;
    data1.resize(n/D);
    for (long k = 0; k < n/D; k++)
      data1[k].resize(D);

    for (long i = 0; i < n; i++) {
      auto p = zMStar.breakIndexByDim(i, dim);
      long k = p.first; 
      long j = p.second;
      data1[k][j] = data[i];
    }

    for (long k = 0; k < n/D; k++) {
      for (long j = 0; j < D; j++) {
	RX acc, val, tmp; 
	acc = 0;
	for (long i = 0; i < D; i++) {
	  if (!mat1.get(val, i, j, k)) {
	    NTL::mul(tmp, data1[k][i], val);
	    NTL::add(acc, acc, tmp);
	  }
	}
	rem(acc, acc, G);
	data[zMStar.assembleIndexByDim(make_pair(k,j), dim)] = acc;
      }
    }
  }

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const PlaintextBlockMatrixBaseInterface& mat, long dim)
  {
    PA_BOILER

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    const PlaintextMultiBlockMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMultiBlockMatrixInterface<type>& >( mat );

    vector< vector<RX> > data1;
    data1.resize(n/D);
    for (long k = 0; k < n/D; k++)
      data1[k].resize(D);

    for (long i = 0; i < n; i++) {
      auto p = zMStar.breakIndexByDim(i, dim);
      long k = p.first; 
      long j = p.second;
      data1[k][j] = data[i];
    }

    for (long k = 0; k < n/D; k++) {
      for (long j = 0; j < D; j++) {
	vec_R acc, tmp, tmp1;
	mat_R val;
	acc.SetLength(d);
	for (long i = 0; i < D; i++) {
	  if (!mat1.get(val, i, j, k)) {
            VectorCopy(tmp1, data1[k][i], d);
            mul(tmp, tmp1, val);
            add(acc, acc, tmp);
	  }
	}
        conv(data[zMStar.assembleIndexByDim(make_pair(k,j), dim)], acc);
      }
    }
  }

}; 


void mat_multi1D(const EncryptedArray& ea, NewPlaintextArray& pa, 
  const PlaintextMatrixBaseInterface& mat, long dim)
{
  ea.dispatch<mat_multi1D_pa_impl>(Fwd(pa), mat, dim); 
}


void mat_multi1D(const EncryptedArray& ea, NewPlaintextArray& pa, 
  const PlaintextBlockMatrixBaseInterface& mat, long dim)
{
  ea.dispatch<mat_multi1D_pa_impl>(Fwd(pa), mat, dim); 
}




template<class type> 
class RandomMultiMatrix : public  PlaintextMultiMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

  vector< vector< vector< RX > > > data;
  long dim;

public:
  ~RandomMultiMatrix() { /*cout << "destructor: random matrix\n";*/ }

  RandomMultiMatrix(const EncryptedArray& _ea, long _dim) : ea(_ea), dim(_dim) { 
    long n = ea.size();
    long d = ea.getDegree();

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();

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

  const vector< vector< vector< RX > > >& getData() const { return data; }
  long getDim() const { return dim; }
  virtual const EncryptedArray& getEA() const { return ea; }

  virtual long size() const {
    long n = ea.size();
    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);
    return n/D;
  }

  virtual bool get(RX& out, long i, long j, long k) const {
    long n = ea.size();
    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    assert(k >= 0 && k < n/D);
    if (IsZero(data[k][i][j])) return true;
    out = data[k][i][j];
    return false;
  }
};

PlaintextMatrixBaseInterface *
buildRandomMultiMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new RandomMultiMatrix<PA_GF2>(ea, dim);
    }

    case PA_zz_p_tag: {
      return new RandomMultiMatrix<PA_zz_p>(ea, dim);
    }

    default: return 0;
  }
}

template<class type> class RandomMultiMatrix2 : public  MatMul<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< vector< RX > > > data;
  long dim;

public:
  ~RandomMultiMatrix2() { /*cout << "destructor: random matrix\n";*/ }

  RandomMultiMatrix2(const EncryptedArray& _ea, long _dim,
		     const vector< vector< vector< RX > > >& _data)
    : MatMul<type>(_ea), dim(_dim), data(_data) {} 

  virtual bool multiGet(RX& out, long i, long j, long k) const {
    long D = this->getEA().sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    assert(k >= 0 && k < this->getEA().size()/D);
    if (IsZero(data[k][i][j])) return true;
    out = data[k][i][j];
    return false;
  }
};

MatMulBase *buildRandomMultiMatrix2(const PlaintextMatrixBaseInterface& rMat)
{
  switch (rMat.getEA().getTag()) {
    case PA_GF2_tag: {
      const RandomMultiMatrix<PA_GF2>& rMat2 = 
	dynamic_cast< const RandomMultiMatrix<PA_GF2>& >( rMat );
      return new RandomMultiMatrix2<PA_GF2>(rMat2.getEA(),
				       rMat2.getDim(), rMat2.getData());
    }
    case PA_zz_p_tag: {
      const RandomMultiMatrix<PA_zz_p>& rMat2 = 
	dynamic_cast< const RandomMultiMatrix<PA_zz_p>& >( rMat );
      return new RandomMultiMatrix2<PA_zz_p>(rMat2.getEA(),
					rMat2.getDim(), rMat2.getData());
    }
    default: return 0;
  }
}



template<class type> 
class RandomMultiBlockMatrix : public  PlaintextMultiBlockMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

  vector< vector< vector< mat_R > > > data;
  long dim;

public:
  ~RandomMultiBlockMatrix() { /*cout << "destructor: random matrix\n";*/ }

  RandomMultiBlockMatrix(const EncryptedArray& _ea, long _dim) : ea(_ea), dim(_dim) { 
    long n = ea.size();
    long d = ea.getDegree();

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();

    data.resize(n/D);
    for (long k = 0; k < n/D; k++) {
      data[k].resize(D);
      for (long i = 0; i < D; i++) {
	data[k][i].resize(D);
	for (long j = 0; j < D; j++) {
          data[k][i][j].SetDims(d, d);
	  for (long u = 0; u < d; u++)
	    for (long v = 0; v < d; v++) 
	      random(data[k][i][j][u][v]);
	}
      }
    }
  }

  virtual const EncryptedArray& getEA() const { return ea; }
  const long getDim() const { return dim; }
  const vector< vector< vector< mat_R > > >& getData() const { return data; }

  virtual long size() const {
    long n = ea.size();
    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);
    return n/D;
  }

  virtual bool get(mat_R& out, long i, long j, long k) const {
    long n = ea.size();
    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    assert(k >= 0 && k < n/D);
    if (IsZero(data[k][i][j])) return true;
    out = data[k][i][j];
    return false;
  }
};

PlaintextBlockMatrixBaseInterface *
buildRandomMultiBlockMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new RandomMultiBlockMatrix<PA_GF2>(ea, dim);
    }

    case PA_zz_p_tag: {
      return new RandomMultiBlockMatrix<PA_zz_p>(ea, dim);
    }

    default: return 0;
  }
}

template<class type> 
class RandomMultiBlockMatrix2 : public BlockMatMul<type> {
  PA_INJECT(type) 

  vector< vector< vector< mat_R > > > data;
  long dim;

public:
  ~RandomMultiBlockMatrix2() { /*cout << "destructor: random matrix\n";*/ }

  RandomMultiBlockMatrix2(const EncryptedArray& _ea, long _dim,
                          const vector< vector< vector< mat_R > > >& _data):
    BlockMatMul<type>(_ea), dim(_dim), data(_data)
  {}

  virtual long size() const { // how many transformations
    const EncryptedArray& ea = this->getEA();
    long n = ea.size();
    long D = ea.sizeOfDimension(dim);
    return n/D;
  }

  virtual bool multiGet(mat_R& out, long i, long j, long k) const {
    const EncryptedArray& ea = this->getEA();
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

MatMulBase *
buildRandomMultiBlockMatrix2(const PlaintextBlockMatrixBaseInterface& rMat)
{
  switch (rMat.getEA().getTag()) {
    case PA_GF2_tag: {
      const RandomMultiBlockMatrix<PA_GF2>& rMat2 = 
	dynamic_cast< const RandomMultiBlockMatrix<PA_GF2>& >(rMat);
      return new RandomMultiBlockMatrix2<PA_GF2>(rMat2.getEA(),
					rMat2.getDim(), rMat2.getData());
    }
    case PA_zz_p_tag: {
      const RandomMultiBlockMatrix<PA_zz_p>& rMat2 = 
	dynamic_cast< const RandomMultiBlockMatrix<PA_zz_p>& >( rMat );
      return new RandomMultiBlockMatrix2<PA_zz_p>(rMat2.getEA(),
					rMat2.getDim(), rMat2.getData());
    }
    default: return 0;
  }
}


//*******************************************


template<class type>
class mat_mul1D_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const PlaintextMatrixBaseInterface& mat, long dim)
  {
    PA_BOILER

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    vector< vector<RX> > data1;
    data1.resize(n/D);
    for (long k = 0; k < n/D; k++)
      data1[k].resize(D);

    for (long i = 0; i < n; i++) {
      auto p = zMStar.breakIndexByDim(i, dim);
      long k = p.first; 
      long j = p.second;
      data1[k][j] = data[i];
    }

    for (long k = 0; k < n/D; k++) {
      for (long j = 0; j < D; j++) {
	RX acc, val, tmp; 
	acc = 0;
	for (long i = 0; i < D; i++) {
	  if (!mat1.get(val, i, j)) {
	    NTL::mul(tmp, data1[k][i], val);
	    NTL::add(acc, acc, tmp);
	  }
	}
	rem(acc, acc, G);
	data[zMStar.assembleIndexByDim(make_pair(k,j), dim)] = acc;
      }
    }
  }



  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const PlaintextBlockMatrixBaseInterface& mat, long dim)
  {
    PA_BOILER

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    const PlaintextBlockMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

    vector< vector<RX> > data1;
    data1.resize(n/D);
    for (long k = 0; k < n/D; k++)
      data1[k].resize(D);

    for (long i = 0; i < n; i++) {
      auto p = zMStar.breakIndexByDim(i, dim);
      long k = p.first; 
      long j = p.second;
      data1[k][j] = data[i];
    }

    for (long k = 0; k < n/D; k++) {
      for (long j = 0; j < D; j++) {
	vec_R acc, tmp, tmp1;
	mat_R val;
	acc.SetLength(d);
	for (long i = 0; i < D; i++) {
	  if (!mat1.get(val, i, j)) {
            VectorCopy(tmp1, data1[k][i], d);
            mul(tmp, tmp1, val);
            add(acc, acc, tmp);
	  }
	}
        conv(data[zMStar.assembleIndexByDim(make_pair(k,j), dim)], acc);
      }
    }
  }

}; 


void mat_mul1D(const EncryptedArray& ea, NewPlaintextArray& pa, 
  const PlaintextMatrixBaseInterface& mat, long dim)
{
  ea.dispatch<mat_mul1D_pa_impl>(Fwd(pa), mat, dim); 
}

void mat_mul1D(const EncryptedArray& ea, NewPlaintextArray& pa, 
  const PlaintextBlockMatrixBaseInterface& mat, long dim)
{
  ea.dispatch<mat_mul1D_pa_impl>(Fwd(pa), mat, dim); 
}









template<class type> 
class RandomMatrix : public  PlaintextMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

  vector< vector< RX > > data;
  long dim;

public:
  ~RandomMatrix() { /*cout << "destructor: random matrix\n";*/ }

  RandomMatrix(const EncryptedArray& _ea, long _dim) : ea(_ea), dim(_dim) { 
    long n = ea.size();
    long d = ea.getDegree();

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();

    data.resize(D);
    for (long i = 0; i < D; i++) {
      data[i].resize(D);
      for (long j = 0; j < D; j++) {
        random(data[i][j], d);
      }
    }
  }
  const vector< vector< RX > >& getData() const { return data; }
  long getDim() const { return dim; }

  virtual const EncryptedArray& getEA() const { return ea; }

  virtual bool get(RX& out, long i, long j) const {
    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};

PlaintextMatrixBaseInterface *
buildRandomMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new RandomMatrix<PA_GF2>(ea, dim);
    }

    case PA_zz_p_tag: {
      return new RandomMatrix<PA_zz_p>(ea, dim);
    }

    default: return 0;
  }
}

template<class type> class RandomMatrix2 : public  MatMul<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< RX > > data;
  long dim;

public:
  ~RandomMatrix2() { /*cout << "destructor: random matrix\n";*/ }

  RandomMatrix2(const EncryptedArray& _ea, long _dim,
		const vector< vector< RX > >& _data)
    : MatMul<type>(_ea), dim(_dim), data(_data) {} 

  virtual bool get(RX& out, long i, long j) const {
    long D = this->getEA().sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};


MatMulBase *buildRandomMatrix2(const PlaintextMatrixBaseInterface& rMat)
{
  switch (rMat.getEA().getTag()) {
    case PA_GF2_tag: {
      const RandomMatrix<PA_GF2>& rMat2 = 
	dynamic_cast< const RandomMatrix<PA_GF2>& >( rMat );
      return new RandomMatrix2<PA_GF2>(rMat2.getEA(),
				       rMat2.getDim(), rMat2.getData());
    }
    case PA_zz_p_tag: {
      const RandomMatrix<PA_zz_p>& rMat2 = 
	dynamic_cast< const RandomMatrix<PA_zz_p>& >( rMat );
      return new RandomMatrix2<PA_zz_p>(rMat2.getEA(),
					rMat2.getDim(), rMat2.getData());
    }
    default: return 0;
  }
}

template<class type> 
class RandomBlockMatrix : public  PlaintextBlockMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

  vector< vector< mat_R > > data;
  long dim;

public:
  ~RandomBlockMatrix() { /*cout << "destructor: random block matrix\n";*/ }
  RandomBlockMatrix(const EncryptedArray& _ea, long _dim) : ea(_ea), dim(_dim) { 
    long n = ea.size();
    long d = ea.getDegree();

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();

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

  virtual const EncryptedArray& getEA() const { return ea; }
  long getDim() const { return dim; }
  const vector< vector< mat_R > >& getData() const { return data; }

  virtual bool get(mat_R& out, long i, long j) const {
    const PAlgebra& zMStar = ea.getContext().zMStar;
    long D = zMStar.OrderOf(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};

PlaintextBlockMatrixBaseInterface *
buildRandomBlockMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new RandomBlockMatrix<PA_GF2>(ea, dim);
    }

    case PA_zz_p_tag: {
      return new RandomBlockMatrix<PA_zz_p>(ea, dim);
    }

    default: return 0;
  }
}

template<class type> 
class RandomBlockMatrix2 : public BlockMatMul<type> {
  PA_INJECT(type) 

  vector< vector< mat_R > > data;
  long dim;

public:
  ~RandomBlockMatrix2() { /*cout << "destructor: random matrix\n";*/ }

  RandomBlockMatrix2(const EncryptedArray& _ea, long _dim,
                          const vector< vector< mat_R > >& _data):
    BlockMatMul<type>(_ea), dim(_dim), data(_data)
  {}

  virtual bool get(mat_R& out, long i, long j) const {
    const EncryptedArray& ea = this->getEA();
    long D = ea.sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};

MatMulBase *
buildRandomBlockMatrix2(const PlaintextBlockMatrixBaseInterface& rMat)
{
  switch (rMat.getEA().getTag()) {
    case PA_GF2_tag: {
      const RandomBlockMatrix<PA_GF2>& rMat2 = 
	dynamic_cast< const RandomBlockMatrix<PA_GF2>& >(rMat);
      return new RandomBlockMatrix2<PA_GF2>(rMat2.getEA(),
					rMat2.getDim(), rMat2.getData());
    }
    case PA_zz_p_tag: {
      const RandomBlockMatrix<PA_zz_p>& rMat2 = 
	dynamic_cast< const RandomBlockMatrix<PA_zz_p>& >( rMat );
      return new RandomBlockMatrix2<PA_zz_p>(rMat2.getEA(),
					rMat2.getDim(), rMat2.getData());
    }
    default: return 0;
  }
}









void  TestIt(long m, long p, long r, long d, long L, long dim)
{
  cout << "*** TestIt: m=" << m
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << ", L=" << L
       << ", dim=" << dim
       << endl;

  FHEcontext context(m, p, r);
  buildModChain(context, L, /*c=*/3);

  context.zMStar.printout();
  cout << endl;

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64); // A Hamming-weight-w secret key

  ZZX G;
  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d); 

  cout << "G = " << G << "\n";
  cout << "generating key-switching matrices... ";
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  cout << "done\n";

  cout << "computing masks and tables for rotation...";
  EncryptedArray ea(context, G);
  cout << "done\n";

  // Test a "normal" matrix over the extension field
  {
    // choose a random plaintext square matrix
    shared_ptr<PlaintextMatrixBaseInterface> ptr(buildRandomMatrix(ea, dim));
    std::unique_ptr< MatMulBase > ptr2(buildRandomMatrix2(*ptr));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying 1D with MatMulBase... " << std::flush;
    // mat_mul1D(ea, v, *ptr, dim);     // multiply the plaintext vector
    // mat_mul1D(ea, ctxt2, *ptr, dim); // multiply the ciphertext vector
    matMul1D(v, *ptr2, dim);
    matMul1D(ctxt2, *ptr2, dim, cachezzX);// multiply ciphertext and build cache
    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";

    cout << " Multiplying 1D with MatMulBase+dcrt cache... " << std::flush;
    ctxt2 = ctxt;
    matMul1D(ctxt2, *ptr2, dim, cacheDCRT); // upgrade cache and use in multiplication

    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  {
    // choose a random plaintext square matrix
    shared_ptr<PlaintextMatrixBaseInterface> ptr(buildRandomMatrix(ea, dim));
    std::unique_ptr< MatMulBase > ptr2(buildRandomMatrix2(*ptr));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying 1D with MatMulBase+zzx cache... " << std::flush;
    buildCache4MatMul1D(*ptr2, dim, cachezzX);// build the cache
    matMul1D(ctxt, *ptr2, dim);               // then use it
    matMul1D(v, *ptr2, dim);     // multiply the plaintext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }

  // Test a "multi" matrix over the extension field
  {
    // choose a random plaintext square matrix
    shared_ptr<PlaintextMatrixBaseInterface>
      ptr(buildRandomMultiMatrix(ea,dim));
    std::unique_ptr< MatMulBase > ptr2(buildRandomMultiMatrix2(*ptr));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying multi 1D with MatMulBase... " << std::flush;
    // mat_multi1D(ea, v, *ptr, dim);     // multiply the plaintext vector
    // mat_multi1D(ctxt2, ea, dim, *ptr); // multiply the ciphertext vector
    matMulti1D(v, *ptr2, dim);
    matMulti1D(ctxt2, *ptr2, dim, cachezzX);// multiply ciphertext and build cache
    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";

    cout <<" Multiplying multi 1D with MatMulBase+dcrt cache... "<< std::flush;
    ctxt2 = ctxt;
    matMulti1D(ctxt2, *ptr2, dim, cacheDCRT); // upgrade cache and use in multiplication

    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))    // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  {
    // choose a random plaintext square matrix
    shared_ptr<PlaintextMatrixBaseInterface>
      ptr(buildRandomMultiMatrix(ea,dim));
    std::unique_ptr< MatMulBase > ptr2(buildRandomMultiMatrix2(*ptr));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying multi 1D with MatMulBase+zzx cache... " << std::flush;
    buildCache4MatMulti1D(*ptr2, dim, cachezzX);// build the cache
    matMulti1D(ctxt, *ptr2, dim);               // then use it
    matMulti1D(v, *ptr2, dim);     // multiply the plaintext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }

  // Test a "block matrix" over the base field
  {
    // choose a random plaintext square matrix
    shared_ptr<PlaintextBlockMatrixBaseInterface>
      ptr(buildRandomBlockMatrix(ea, dim));
    shared_ptr<MatMulBase> ptr2(buildRandomBlockMatrix2(*ptr));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << endl << " Multiplying with BlockMatMul 1D... " 
	 << std::flush;
    // mat_mul1D(ea, ctxt2, *ptr, dim);  // multiply the ciphertext vector
    // mat_mul1D(ea, v, *ptr, dim);      // multiply the plaintext vector
    blockMatMul1D(ctxt2, *ptr2, dim, cachezzX);
                                  // multiply ciphertext and build cache
    blockMatMul1D(v, *ptr2, dim); // multiply the plaintext vector
    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr...\n";

    cout << " Multiplying with BlockMatMul 1D+dcrt cache... " << std::flush;
    ctxt2 = ctxt;
    blockMatMul1D(ctxt2, *ptr2, dim, cacheDCRT);
                                  // upgrade cache and use in multiplication
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  {
    // choose a random plaintext square matrix
    shared_ptr<PlaintextBlockMatrixBaseInterface>
      ptr(buildRandomBlockMatrix(ea, dim));
    shared_ptr<MatMulBase> ptr2(buildRandomBlockMatrix2(*ptr));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying with BlockMatMul 1D+zzx cache... " << std::flush;
    buildCache4BlockMatMul1D(*ptr2, dim, cachezzX);// build the cache
    blockMatMul1D(ctxt2, *ptr2, dim);              // then use it
    blockMatMul1D(v, *ptr2, dim); // multiply the plaintext vector
    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr...\n";
  }
  // Test multiple "block" matrices over the base field
  {
    // choose a random plaintext square matrix
    shared_ptr<PlaintextBlockMatrixBaseInterface> ptr(buildRandomMultiBlockMatrix(ea, dim));
    shared_ptr<MatMulBase> ptr2(buildRandomMultiBlockMatrix2(*ptr));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying with multi BlockMatMul... " << std::flush;
    // mat_multi1D(ea, v, *ptr, dim);         // multiply the plaintext vector
    // mat_multi1D_block(ctxt2, ea,dim,*ptr); // multiply the ciphertext vector
    blockMatMulti1D(v, *ptr2, dim);     // multiply the plaintext vector
    blockMatMulti1D(ctxt2, *ptr2, dim); // multiply the ciphertext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";

  }
}


void usage(char *prog) 
{
  cout << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cout << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cout << "  e.g, 'm=2047 p=2 L=4'\n\n";
  cout << "  m defines the cyclotomic polynomial Phi_m(X)\n";
  cout << "  p is the plaintext base [default=2]" << endl;
  cout << "  r is the lifting [default=1]" << endl;
  cout << "  d is the degree of the field extension [default==1]\n";
  cout << "    (d == 0 => factors[0] defined the extension)\n";
  cout << "  L is the # of primes in the modulus chain [default=4]\n";
  exit(0);
}

/* Testing the functionality of multiplying an encrypted vector by a plaintext
 * matrix, either over the extension- or the base-field/ring.
 */
int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["m"] = "2047";
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["d"] = "0";
  argmap["L"] = "4";
  argmap["dim"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long m = atoi(argmap["m"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long d = atoi(argmap["d"]);
  long L = atoi(argmap["L"]);
  long dim = atoi(argmap["dim"]);

  //  setTimersOn();
  setTimersOn();
  TestIt(m, p, r, d, L, dim);
  cout << endl;
  //printAllTimers();
  //cout << endl;

}

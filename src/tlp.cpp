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
#include "matrix.h"
#include <NTL/ZZ.h>
#include <NTL/lzz_pXFactoring.h>
#include <cassert>
#include <cstdio>




template<class type> 
class SingleBlockMatrix : public  PlaintextBlockMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

  Mat<R> data;

public:
  SingleBlockMatrix(const EncryptedArray& _ea, const vector<ZZX>& vec) : ea(_ea) { 
    long d = ea.getDegree();

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();

    data.SetDims(d, d);
    for (long i = 0; i < d; i++) 
      for (long j = 0; j < d; j++) 
         conv(data[i][j], coeff(vec[i], j));
  }

  virtual const EncryptedArray& getEA() const {
    return ea;
  }

  virtual bool get(Mat<R>& out, long i, long j) const {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    if (i != j) return true;
    out = data;
    return false;
  }

};



PlaintextBlockMatrixBaseInterface *buildSingleBlockMatrix(const EncryptedArray& ea,
                                                          const vector<ZZX>& vec)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new SingleBlockMatrix<PA_GF2>(ea, vec);
    }

    case PA_zz_p_tag: {
      return new SingleBlockMatrix<PA_zz_p>(ea, vec);
    }

    default: return 0;
  }
}





template<class type> 
class MultiBlockMatrix : public  PlaintextBlockMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

  Vec< Mat<R> > data;

public:
  MultiBlockMatrix(const EncryptedArray& _ea, const vector< vector<ZZX> >& vec) : ea(_ea) { 
    long n = ea.size();
    long d = ea.getDegree();

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();

    data.SetLength(n);
    for (long k = 0; k < n; k++) {
      data[k].SetDims(d, d);
      for (long i = 0; i < d; i++) 
        for (long j = 0; j < d; j++) 
           conv(data[k][i][j], coeff(vec[k][i], j));
    }
  }

  virtual const EncryptedArray& getEA() const {
    return ea;
  }

  virtual bool get(Mat<R>& out, long i, long j) const {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    if (i != j) return true;
    out = data[i];
    return false;
  }
};

PlaintextBlockMatrixBaseInterface *buildMultiBlockMatrix(const EncryptedArray& ea,
                                                          const vector< vector<ZZX> >& vec)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new MultiBlockMatrix<PA_GF2>(ea, vec);
    }

    case PA_zz_p_tag: {
      return new MultiBlockMatrix<PA_zz_p>(ea, vec);
    }

    default: return 0;
  }
}



void  TestIt(long m, long p, long r, long d)
{
  cout << "\n\n******** TestIt" << (isDryRun()? "(dry run):" : ":")
       << " m=" << m 
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << endl;

  FHEcontext context(m, p, r);
  buildModChain(context, /*L=*/10, /*c=*/2);
  context.zMStar.printout();
  cout << endl;

  cout << "generating keys and key-switching matrices... " << std::flush;
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64);// A Hamming-weight-w secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  cout << "done\n";

  ZZX G;
  if (d == 0) {
    G = context.alMod.getFactorsOverZZ()[0];
    d = deg(G);
  }
  else
    G = makeIrredPoly(p, d);
  cout << "G = " << G << "\n";
  cout << "computing masks and tables for rotation... " << std::flush;
  EncryptedArray ea(context, G);
  cout << "done\n";

  long nslots = ea.size();

  PlaintextArray p0(ea);
  PlaintextArray pp0(ea);

  Ctxt c0(publicKey), cc0(publicKey);

  cout << "\nTest #1: Apply the same linear transformation to all slots\n";
  {
  vector<ZZX> LM(d); // LM selects even coefficients
  for (long j = 0; j < d; j++) 
    if (j % 2 == 0) LM[j] = ZZX(j, 1);

  // "building" the linearized-polynomial coefficients
  vector<ZZX> C;
  ea.buildLinPolyCoeffs(C, LM);

  p0.random();  
  ea.encrypt(c0, publicKey, p0);
  ea.encrypt(cc0, publicKey, p0);

  applyLinPoly1(ea, c0, C);
  ea.decrypt(c0, secretKey, p0);

  shared_ptr<PlaintextBlockMatrixBaseInterface> mat(buildSingleBlockMatrix(ea, LM));
  free_mat_mul(ea, cc0, *mat);
  ea.decrypt(cc0, secretKey, pp0);

  if (pp0.equals(p0))
    cout << "GOOD\n";
  else
    cout << "BAD\n";
  }


  cout << "\nTest #2: Apply different transformations to the different slots\n";
  {
  vector< vector<ZZX> > LM(nslots); 
  // LM[i] rotates the coefficients in the i'th slot by (i % d)
  for (long i = 0; i < nslots; i++) {
    LM[i].resize(d);
    for (long j = 0; j < d; j++)  {
      long jj = (i+j) % d;
      LM[i][j] = ZZX(jj, 1);
    }
  }

  // "building" the linearized-polynomial coefficients
  vector< vector<ZZX> > C(nslots);
  for (long i = 0; i < nslots; i++)
    ea.buildLinPolyCoeffs(C[i], LM[i]);

  p0.random();
  ea.encrypt(c0, publicKey, p0);
  ea.encrypt(cc0, publicKey, p0);

  applyLinPolyMany(ea, c0, C); // apply the linearized polynomials
  ea.decrypt(c0, secretKey, p0);

  shared_ptr<PlaintextBlockMatrixBaseInterface> mat(buildMultiBlockMatrix(ea, LM));
  free_mat_mul(ea, cc0, *mat);
  ea.decrypt(cc0, secretKey, pp0);

  if (pp0.equals(p0))
    cout << "GOOD\n";
  else
    cout << "BAD\n";
  }

  cout << "\nTest #3: Testing low-level (cached) implementation\n";
  {
  vector< vector<ZZX> > LM(nslots); 
  // LM[i] adds coefficients (i % d) and (i+1 % d) in the i'th slot
  for (long i = 0; i < nslots; i++) {
    LM[i].resize(d);
    for (long j = 0; j < d; j++)  {
      if ( j == (i % d) || j == ((i+1)%d) )
	LM[i][j] = conv<ZZX>(1L);
    }
  }

  // "building" the linearized-polynomial coefficients
  vector< vector<ZZX> > C(nslots);
  for (long i = 0; i < nslots; i++)
    ea.buildLinPolyCoeffs(C[i], LM[i]);

  // "encoding" the linearized-polynomial coefficients
  vector<ZZX> encodedC(d);
  for (long j = 0; j < d; j++) {
    vector<ZZX> v(nslots);
    for (long i = 0; i < nslots; i++) v[i] = C[i][j];
    ea.encode(encodedC[j], v);
  }

  p0.random();  
  ea.encrypt(c0, publicKey, p0);
  ea.encrypt(cc0, publicKey, p0);

  applyLinPolyLL(c0, encodedC, ea.getDegree()); // apply linearized polynomials
  ea.decrypt(c0, secretKey, p0);

  shared_ptr<PlaintextBlockMatrixBaseInterface> mat(buildMultiBlockMatrix(ea, LM));
  free_mat_mul(ea, cc0, *mat);
  ea.decrypt(cc0, secretKey, pp0);

  if (pp0.equals(p0))
    cout << "GOOD\n";
  else
    cout << "BAD\n";
  }
}


int main(int argc, char *argv[]) 
{
  ArgMapping amap;

  bool dry = false;
  amap.arg("dry", dry, "dry=1 for a dry-run");

  long m=91;
  amap.arg("m", m, "use specified value as modulus");

  long p=2;
  amap.arg("p", p, "plaintext base");

  long r=1;
  amap.arg("r", r,  "lifting");

  long d=0;
  amap.arg("d", d, "degree of the field extension");
  amap.note("d == 0 => factors[0] defines extension");

  amap.parse(argc, argv);

  long repeat = 2;
  setTimersOn();
  setDryRun(dry);
  for (long repeat_cnt = 0; repeat_cnt < repeat; repeat_cnt++) {
    TestIt(m, p, r, d);
  }

}

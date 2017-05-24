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
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include "matmul.h"
#include <NTL/ZZ.h>
#include <NTL/lzz_pXFactoring.h>
#include <cassert>
#include <cstdio>

static bool noPrint = false;

static MatMulBase*
buildSingleBlockMatrix(const EncryptedArray& ea, const vector<ZZX>& vec);


template<class type> class SingleBlockMatrix : public BlockMatMul<type> {
  PA_INJECT(type) 

  Mat<R> data;

public:
  SingleBlockMatrix(const EncryptedArray& _ea, const vector<ZZX>& vec) :
    BlockMatMul<type>(_ea)
  { 
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long d = _ea.getDegree();

    data.SetDims(d, d);
    for (long i = 0; i < d; i++) 
      for (long j = 0; j < d; j++) 
         conv(data[i][j], coeff(vec[i], j));
  }

  virtual bool get(Mat<R>& out, long i, long j) const
  {
    assert(i >= 0 && i < this->getEA().size());
    assert(j >= 0 && j < this->getEA().size());
    if (i != j) return true;
    out = data;
    return false;
  }
};

static MatMulBase*
buildSingleBlockMatrix(const EncryptedArray& ea, const vector<ZZX>& vec)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new SingleBlockMatrix<PA_GF2>(ea, vec);
    }
    case PA_zz_p_tag: {
      return new SingleBlockMatrix<PA_zz_p>(ea, vec);
    }
    default: return nullptr;
  }
}


template<class type> class MultiBlockMatrix : public BlockMatMul<type> {
  PA_INJECT(type) 

  Vec< Mat<R> > data;

public:
  MultiBlockMatrix(const EncryptedArray& _ea, const vector<vector<ZZX> >& vec):
    BlockMatMul<type>(_ea)
  { 
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long n = _ea.size();
    long d = _ea.getDegree();

    data.SetLength(n);
    for (long k = 0; k < n; k++) {
      data[k].SetDims(d, d);
      for (long i = 0; i < d; i++) 
        for (long j = 0; j < d; j++) 
           conv(data[k][i][j], coeff(vec[k][i], j));
    }
  }

  virtual bool get(Mat<R>& out, long i, long j) const {
    assert(i >= 0 && i < this->getEA().size());
    assert(j >= 0 && j < this->getEA().size());
    if (i != j) return true;
    out = data[i];
    return false;
  }
};

static MatMulBase* buildMultiBlockMatrix(const EncryptedArray& ea,
					 const vector< vector<ZZX> >& vec)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new MultiBlockMatrix<PA_GF2>(ea, vec);
    }
    case PA_zz_p_tag: {
      return new MultiBlockMatrix<PA_zz_p>(ea, vec);
    }
    default: return nullptr;
  }
}



void  TestIt(long m, long p, long r, long d)
{
  if (!noPrint)
    cout << "\n\n******** TestIt" << (isDryRun()? "(dry run):" : ":")
       << " m=" << m 
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << endl;

  FHEcontext context(m, p, r);
  buildModChain(context, /*L=*/3, /*c=*/2);

  ZZX G;
  if (d == 0) {
    G = context.alMod.getFactorsOverZZ()[0];
    d = deg(G);
  }
  else
    G = makeIrredPoly(p, d);

  if (!noPrint) {
    context.zMStar.printout();
    cout << endl;
    cout << "G = " << G << "\n";

    cout << "generating keys and key-switching matrices... " << std::flush;
  }
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64);// A Hamming-weight-w secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  if (!noPrint) {
    cout << "done\n";
    cout << "computing masks and tables for rotation... " << std::flush;
  }
  EncryptedArray ea(context, G);
  if (!noPrint)
    cout << "done\n";

  long nslots = ea.size();

  NewPlaintextArray p0(ea);
  NewPlaintextArray pp0(ea);

  Ctxt c0(publicKey);

  if (!noPrint)
    cout << "\nTest #1: Apply the same linear transformation to all slots\n";
  {
  vector<ZZX> LM(d); // LM selects even coefficients
  for (long j = 0; j < d; j++) 
    if (j % 2 == 0) LM[j] = ZZX(j, 1);

  // "building" the linearized-polynomial coefficients
  vector<ZZX> C;
  ea.buildLinPolyCoeffs(C, LM);

  random(ea, p0);  
  ea.encrypt(c0, publicKey, p0);
  applyLinPoly1(ea, c0, C);
  ea.decrypt(c0, secretKey, pp0);

  shared_ptr<MatMulBase> mat(buildSingleBlockMatrix(ea, LM));
  NewPlaintextArray p1(p0);
  blockMatMul(p1, *mat);
  if (equals(ea, pp0, p1))
    cout << "GOOD\n";
  else
    cout << "BAD\n";
  }


  if (!noPrint)
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

  random(ea, p0);
  ea.encrypt(c0, publicKey, p0);
  applyLinPolyMany(ea, c0, C); // apply the linearized polynomials
  ea.decrypt(c0, secretKey, pp0);

  shared_ptr<MatMulBase> mat(buildMultiBlockMatrix(ea, LM));
  NewPlaintextArray p1(p0);
  blockMatMul(p1, *mat);
  if (equals(ea, pp0, p1))
    cout << "GOOD\n";
  else
    cout << "BAD\n";
  }

  if (!noPrint)
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

  random(ea, p0);  
  ea.encrypt(c0, publicKey, p0);
  applyLinPolyLL(c0, encodedC, ea.getDegree()); // apply linearized polynomials
  ea.decrypt(c0, secretKey, pp0);

  shared_ptr<MatMulBase> mat(buildMultiBlockMatrix(ea, LM));
  NewPlaintextArray p1(p0);
  blockMatMul(p1, *mat);
  if (equals(ea, pp0, p1))
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

  amap.arg("noPrint", noPrint, "suppress printouts");

  amap.parse(argc, argv);

  long repeat = 2;
  setTimersOn();
  setDryRun(dry);
  for (long repeat_cnt = 0; repeat_cnt < repeat; repeat_cnt++) {
    TestIt(m, p, r, d);
  }

}

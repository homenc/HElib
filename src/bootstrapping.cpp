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
#include "NTL/ZZ_pE.h"
#include "NumbTh.h"    // defines argmax(...)
#include "PAlgebra.h"
#include "EncryptedArray.h"

// NTL_CLIENT

// Computing an array containing the powers X^{e*i} mod G
template<class type> class PolyPowers {
public:
  PA_INJECT(type)

  static long computePolyPowers(vector<ZZX>& v,
				const PAlgebraMod& alMod, long e, long size)
  {
    const PAlgebraModDerived<type>& tab = alMod.getDerived(type());
    zz_pBak bak; bak.save(); tab.restoreContext();
    RXModulus G = tab.getFactors()[0];
    RX alpha; SetX(alpha); PowerMod(alpha,alpha,e,G); // alpha = X^e mod G
    v.resize(size);
    conv(v[0], 1);
    conv(v[1], alpha);
    RX alpha2i = alpha;
    for (long i=2; i<size; i++) {
      MulMod(alpha2i, alpha2i, alpha, G); // alpha2i = alpha^i mod G
      conv(v[i], alpha2i);
    }
    return 0;
  }
};
long
buildPolyPowers(vector<ZZX>& v, const PAlgebraMod& alMod, long e, long size)
{
  switch (alMod.getTag()) {
  case PA_GF2_tag:
    return PolyPowers<PA_GF2>::computePolyPowers(v, alMod, e, size);

  case PA_zz_p_tag:
    return PolyPowers<PA_zz_p>::computePolyPowers(v, alMod, e, size);
  
  default:
    v.clear(); 
    return -1;
  }
}


void getStepOneMatrix(vector< vector<ZZX> >&matrix,
                      long m1, long m2, const EncryptedArray& ea)
{
  // Find generator sets for Zm1* /<p>
  vector<long> gens, ords;
  long p = ea.getContext().zMStar.getP();
  long d = findGenerators(gens, ords, m1, p);
  if (ea.getDegree() != d) { // verify that p has the same order mod m,m1
     cerr << "  Cannod handle the case where d="<<ea.getDegree()
	  << "!=d1="<<d<<endl;
     exit(0);
  }
  assert(gens.size()==1); // Zm1* /<p> is cyclic

  long n = abs(ords[0]); // the order of Zm1* /<p>
  vector<long> T(n);     // representative set for Zm1* /<p>
  for (long i=0; i<(long)T.size(); i++) 
    T[i] = PowerMod(gens[0], i, m1);

  vector< ZZX > eta1Powers;      // eta1 = X^m2 mod G is an m1'th root of unity
  buildPolyPowers(eta1Powers, ea.getContext().alMod, m2, m1); // powers of eta1

  // Now prepare the matrix to use for Step 1
  //
  // We return an nSlots-by-nSlots matrix consisting of n^2 submatrices B_{i,j}
  // of dimension phi(m2)-by-phi(m2), where each B_{i,j} is fixed along the
  // main diagonal and zero elsewhere (so we have total n nonzero diagonals).
  // The main diagonal of each B_{i,j} consists of phi(m2) repetitions of the
  // same A_{i,j}. For example with phi(m2)=3 and n=2 we have the 6x6 matrix
  //
  //   ( A_{0,0}              A_{0,1}               )
  //   (        A_{0,0}              A_{0,1}        )
  //   (               A_{0,0}              A_{0,1} )
  //   ( A_{1,0}              A_{1,1}               )
  //   (        A_{1,0}              A_{1,1}        )
  //   (               A_{1,0}              A_{1,1} )
  //
  // Each A_{i,j} is itself a d-by-d block over Z_p with columns
  //
  //   A_{i,j}= ( ti^{jd} | ti^{jd+1} | ... | ti^{(j+1)d-1} )
  //
  // where ti = eta1^{T[i]}\in GF(p^d), represented as a vector over Z_p.

  // We first build a linearized polynomial for each A_{i,j}

  vector<ZZX> A[n][n];
  for (long i=0; i<n; i++) { // Go over the rows A[i]
    long ti = T[i];

    for (long j=0; j<n; j++) { // Go over columns and build A[i][j]
      long jd = MulMod(j,d,m1);

      vector<ZZX> L(d); // L[k] is the k'th column of A[i][j]
      for (long k=0; k<d; k++) {
	long idx = MulMod(ti, AddMod(jd, k, m1), m1);
	L[k] = eta1Powers[idx]; // L[k] = eta1^{ti*(jd+k)}
      }
      ea.buildLinPolyCoeffs(A[i][j], L);
    }
  }

  // Next we encode the diagonals of the matrix: In the i'th diagonal, 
  // slots j*phi(m2),...,(j+1)*phi(m2)-1 encode the matrix A[j][i+j mod n]

  long nSlots = ea.size();
  long phi_m2 = nSlots/n;    // only works when d1==d
  assert(phi_N(m2)==phi_m2); // sanity check

  matrix.resize(nSlots);
  for (long i=0; i<n; i++) { // prepare the i'th nonempty diagonal
    vector<ZZX>& diag = matrix[i*phi_m2];
    diag.resize(d); // each diagonal has d constants

    for (long k=0; k<d; k++) { // prepare the k'th constant in the diagonal
      vector<ZZX> slots(nSlots); // each conatant has nSlots slots

      // The slots are grouped in n batches of phi(m2) slots each
      for (long j=0; j<n; j++) { // set the j'th batch
	long i_plus_j = AddMod(i,j,n);
        ZZX& slotContent = A[j][i_plus_j][k];
	for (long b=0; b<phi_m2; b++) // All slots in this batch are the same
	  slots[j*phi_m2 +b] = slotContent;
      }
      ea.encode(diag[k], slots);
    }
  }
}

template<class type> 
class Step2Matrix : public  PlaintextMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;
  long m2;
  vector<long> T;         // representative set for Zm2*/ <p^d1>
  vector<ZZX> eta2Powers; // powers of eta2 = X^m1 mod G

public:
  ~Step2Matrix() { cerr << "destructor: Step-2 matrix\n"; }

  Step2Matrix(const EncryptedArray& _ea, long _m1, long _m2): ea(_ea)
  { 
    m2 = _m2;
    long d = ea.getDegree();
    long p = ea.getContext().zMStar.getP();
    long p2d = PowerMod(p % _m2, d, _m2);
    vector<long> gens, ords;     // Find generator sets for Zm2* /<p^d>
    assert(findGenerators(gens, ords, _m2, p2d) == 1); // p^d==1 mod m2
    assert(gens.size()==1);      // Zm2* = Zm2* /<p^d> is cyclic

    T.resize(abs(ords[0]));     // representative set for Zm2* / <p^d>
    T[0] = 1;
    for (long i=1; i<(long)T.size(); i++)
      T[i] = MulMod(T[i-1], gens[0], _m2);

    PolyPowers<type>::computePolyPowers(eta2Powers,
					ea.getContext().alMod, _m1, _m2);
  }

  virtual const EncryptedArray& getEA() const { return ea; }

  virtual void get(RX& out, long i, long j) const {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    if ((i/T.size()) != (j/T.size())) // return zero outside the main blocks
      clear(out);
    else {
      // find the position of (i,j) inside their block
      i %= T.size();
      j %= T.size();

      // return out = eta2^{T[j]*i mod m2}
      long exp = MulMod(i, T[j], m2);
      conv(out, eta2Powers[exp]);
    }
  }
};
PlaintextMatrixBaseInterface*
buildStep2Matrix(const EncryptedArray& ea, long m1, long m2)
{
  switch (ea.getContext().alMod.getTag()) {
  case PA_GF2_tag: 
    return new Step2Matrix<PA_GF2>(ea, m1, m2);

  case PA_zz_p_tag: 
    return new Step2Matrix<PA_zz_p>(ea, m1, m2);

  default: return 0;
  }
}


// An elementary multiplication by a Z_p plaintext matrix.
// matrix is a list of diagonals, each diagonal is either empty or
// of the format needed for applyLinPolyLL (in EncryptedArray.cpp)
void elementaryMatrixMultiply(Ctxt& res, const vector< vector<ZZX> >& matrix,
			      const Ctxt& in, const EncryptedArray& ea)
{
  res.clear(); // initialize to zero
  for (long i=0; i<(long)matrix.size(); i++) {
    if (matrix[i].size()==0) continue; // a zero diagonal

    // non-zero diagonal: rotate by i and apply the linear transfromation
    Ctxt tmp(in);
    ea.rotate(tmp, i);                  // rotate by i
    applyLinPolyLL(ea, tmp, matrix[i]); // multiply by diagonal
    res += tmp;                         // add to result
  }
}



void testIt(long m1, long m2, long p, long r)
{
  long m=m1*m2;
  FHEcontext context(m, p, r);
  EncryptedArray ea(context, context.alMod.getFactorsOverZZ()[0]);

  // Compute the matrix for step 1
  vector< vector<ZZX> > step1Matrix;
  getStepOneMatrix(step1Matrix, m1, m2, ea);

  // Compute the matrix for step 2
  PlaintextMatrixBaseInterface* step2Matrix_pt = buildStep2Matrix(ea, m1, m2);

  cout << "okay";
}

void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'm1=257 m2=17 p=2 r=1'\n\n";
  cerr << "  m = m1*m2 is the cyclotomic ring [defaults=4369]\n";
  cerr << "  p is the plaintext base [default=2]\n";
  cerr << "  r is the lifting [default=1]\n";
  exit(0);
}

int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["m1"] = "257";
  argmap["m2"] = "17";
  argmap["p"] = "2";
  argmap["r"] = "1";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long m1 = atoi(argmap["m1"]);
  long m2 = atoi(argmap["m2"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);

  assert(GCD(m1,m2)==1);
  cout << "m1="<<m1<<", m2="<<m2<<", p="<<p<<", r="<<r<<endl;
  testIt(m1,m2,p,r);

  return 0;
}

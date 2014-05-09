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

/***************************************************************/
/*   Computing an array containing the powers X^{e*i} mod G    */
/***************************************************************/
template<class type> class PolyPowers {
public:
  PA_INJECT(type)

  static long computePolyPowers(vector<RX>& v,
				const EncryptedArrayDerived<type>& ea,
				long e, long size)
  {
    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();
    RXModulus G = ea.getG();
    RX alpha; SetX(alpha); PowerMod(alpha,alpha,e,G); // alpha = X^e mod G
    v.resize(size);
    conv(v[0], 1);
    v[1] = alpha;
    for (long i=2; i<size; i++) {
      MulMod(v[i], v[i-1], alpha, G); // v[i] = alpha^i mod G
    }
    return 0;
  }
};

/*************************************************************/
/****      Building the "Step-1 matrix" over Z_{p^r}      ****/
/*************************************************************/
template<class type> 
class Step1Matrix : public  PlaintextBlockMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;
  // long ord1; // order of quotient group Z_m1* / <p>    = phi(m1)/d
  long ord2; // order of quotient group Z_m2* / <p^d> [[ = phi(m2) for now ]]

  vector< vector< mat_R > > A; // basic ord1 x ord1 block of this matrix

  // return the inverse of the matrix encoded in the A data member
  void invert(vector< vector< mat_R > >& Ainv) const;

public:
  Step1Matrix(const EncryptedArray& ea, long m1, long m2); // constructor

  // copy/inverse constructor
  Step1Matrix(const Step1Matrix& other, bool invFlag=false) : ea(other.ea) {
    ord2 = other.ord2;
    if (invFlag) other.invert(A); // Invert the matrix of other
    else         A = other.A;     // Copy the matrix of other
  }

  virtual const EncryptedArray& getEA() const { return ea; }

  // The matrix represented by this object has dimension nSlots x nSlots,
  // and it consists of ord1^2 submatrices B_{i,j} of dimension ord2 x ord2,
  // where each B_{i,j} is fixed along the main diagonal and zero elsewhere
  // (so we have total ord1 nonzero diagonals).
  // The main diagonal of each B_{i,j} consists of ord2 repetitions of the
  // same A_{i,j}. For example with ord1=2 and ord2=3 we have the 6x6 matrix
  //
  //   ( A_{0,0}              A_{0,1}               )
  //   (        A_{0,0}              A_{0,1}        )
  //   (               A_{0,0}              A_{0,1} )
  //   ( A_{1,0}              A_{1,1}               )
  //   (        A_{1,0}              A_{1,1}        )
  //   (               A_{1,0}              A_{1,1} )
  //
  // Each A_{i,j} is itself a d-by-d block over Z_{p^r}

  virtual bool get(mat_R& out, long i, long j) const { // callback function
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    if ((i % ord2)!=(j % ord2)) // Not on the diagonal of the B_{i,j}'s
      return true;              // return true for an empty entry

    out = A[i/ord2][j/ord2]; // The block indexes are i,j div ord2
    return false;
  }
};
PlaintextBlockMatrixBaseInterface *
buildStep1Matrix(const EncryptedArray& ea, long m1, long m2)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new Step1Matrix<PA_GF2>(ea, m1, m2);
    }
    case PA_zz_p_tag: {
      return new Step1Matrix<PA_zz_p>(ea, m1, m2);
    }
    default: return 0;
  }
}
PlaintextBlockMatrixBaseInterface *
buildStep1Inverse(const PlaintextBlockMatrixBaseInterface* other)
{
  switch (other->getEA().getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      const Step1Matrix<PA_GF2>* ptr 
	= dynamic_cast< const Step1Matrix<PA_GF2>* >(other);
      return new Step1Matrix<PA_GF2>(*ptr, /*invert=*/true);
    }
    case PA_zz_p_tag: {
      const Step1Matrix<PA_zz_p>* ptr
	= dynamic_cast< const Step1Matrix<PA_zz_p>* >(other);
      return new Step1Matrix<PA_zz_p>(*ptr, /*invert=*/true);
    }
    default: return 0;
  }
}

template<class type>
Step1Matrix<type>::Step1Matrix(const EncryptedArray& _ea, long m1, long m2)
  : ea(_ea)
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

  long nSlots = ea.size();
  long ord1 = abs(ords[0]); // the order of Zm1* /<p>
  ord2 = nSlots / ord1; // the order of Zm2* /<p^d> [[ = phi(m2) for now ]]

  vector<long> T(ord1); // representative set for Zm1* /<p>
  T[0] = 1;
  T[1] = gens[0];
  for (long i=2; i<ord1; i++) 
    T[i] = MulMod(T[i-1], gens[0], m1); // T[i] = g^i mod m1

  // prepare a table of eta1^i for i=0,1,...,m2-1
  vector< RX > eta1Powers; // eta1 = X^m2 mod G is an m1'th root of unity
  PolyPowers<type>::computePolyPowers(eta1Powers,       // powers of eta1
				      ea.getDerived(type()),m2,m1);

  const PAlgebraModDerived<type>& tab=ea.getContext().alMod.getDerived(type());
  RBak bak; bak.save();
  tab.restoreContext();
  A.resize(ord1);
  for (long i=0; i<ord1; i++) { // Go over the rows A_i
    A[i].resize(ord1);
    long di= MulMod(d,i,m1);    // d*i mod m1

    for (long j=0; j<ord1; j++) { // Go over columns and build A_{i,j}
      A[i][j].SetDims(d,d);
      long tj = T[j];

      for (long k=0; k<d; k++) { // Go over the d rows of A_{i,j}
	// set the k'th row to eta1^{tj*(di+k)}
	long idx = MulMod(tj, AddMod(di, k, m1), m1);  // tj*(di+k)
	VectorCopy(A[i][j][k], eta1Powers[idx], d); // eta1^{tj*(di+k)}
      }
    }
  }
};

// returns the inverse of the matrix encoded in the A data member
template <class type>
void Step1Matrix<type>::invert(vector< vector< mat_R > >& Ainv) const
{
  const PAlgebraModDerived<type>& tab=ea.getContext().alMod.getDerived(type());
  zz_pBak bak; bak.save(); tab.restoreContext();

  // Prepare a single matrix with all the A_{i,j}'s
  long d = A[0][0].NumRows(); // dimension of small A_{i,j}'s
  long dim = A.size() * d;    // dimension of big matrix
  mat_R bigA(INIT_SIZE, dim, dim);

  // tile bigA with all the A_{i,j}'s
  for (long i=0; i<dim; i++) for (long j=0; j<dim; j++)
      bigA[i][j] = A[i/d][j/d][i%d][j%d];

  // invert the big matrix
  ppInvert(bigA, bigA, tab.getZMStar().getP(), tab.getR());

  // Copy the inverted matrix back into the small matrices in Ainv

  // begin by allocating space
  Ainv.resize(A.size());
  for (long i=0; i<(long)A.size(); i++) {
    Ainv[i].resize(A.size());
    for (long j=0; j<(long)A.size(); j++)
      Ainv[i][j].SetDims(d,d);
  }
  // Then copy the data
  for (long i=0; i<dim; i++) for (long j=0; j<dim; j++)
      Ainv[i/d][j/d][i%d][j%d] = bigA[i][j];
}


/*******************************************************************/
/** Building the "Step-2 matrix" over GF(p^d) (lifted to mod p^r) **/
/*******************************************************************/
template<class type> 
class Step2Matrix : public  PlaintextMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;
  long ord1; // order of quotient group Z_m1* / <p>    = phi(m1)/d
  // long ord2; // order of quotient group Z_m2* / <p^d> [[=phi(m2) for now ]]

  mat_RE A;

  // return the inverse of the matrix [get(i,j)]_{i,j}
  void invert(mat_RE& Minv) const;

public:
  Step2Matrix(const EncryptedArray& _ea, long m1, long m2); // constructor

  // copy/inverse constructor
  Step2Matrix(const Step2Matrix& other, bool invFlag=false) : ea(other.ea) {
    ord1 = other.ord1;
    if (invFlag) other.invert(A); // Invert the matrix of other
    else         A = other.A;     // Copy the matrix of other
  }

  virtual const EncryptedArray& getEA() const { return ea; }

  // The matrix represented by this object has dimension nSlots x nSlots over
  // GF(p^d), which is a block matrix with ord2 x ord2 blocks on the main
  // diagonal and zero elsewhere. Each block is a Vandermonde matrix with
  // all the powers 0,1,...ord2 of eta2^i for all i \in Zm2*. For example
  // for m2=5 we have ord2=phi(m2)=4 and eta2 of order 5. Id we have nSlots=8
  // then we get the followin g8x8 matrix:
  //
  //   ( 1 eta2   eta2^2 eta2^3                        )
  //   ( 1 eta2^2 eta2^4 eta2^3                        )
  //   ( 1 eta2^4 eta2^3 eta2^2                        )
  //   ( 1 eta2^3 eta2   eta2^4                        )
  //   (                        1 eta2   eta2^2 eta2^3 )
  //   (                        1 eta2^2 eta2^4 eta2^3 )
  //   (                        1 eta2^4 eta2^3 eta2^2 )
  //   (                        1 eta2^3 eta2   eta2^4)
  //
  // where each block is V(eta, eta^2, eta^4, eta^3). Note that the order
  // or the rows in the block is determined by the ordering of Zm2* which
  // is determined by some generator. 
  //
  // In the example above the generator was 2 so we get th eorder 1,2,4,3.

  virtual bool get(RX& out, long i, long j) const {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    if ((i/ord1) != (j/ord1)) // zero outside the main blocks
      return true; // return true for an empty entry

    // Return the position of (i,j) inside their block
    conv(out, A[i % ord1][j % ord1]);
    return false; // false for a non-empty entry
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
PlaintextMatrixBaseInterface *
buildStep2Inverse(const PlaintextMatrixBaseInterface* other)
{
  switch (other->getEA().getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      const Step2Matrix<PA_GF2>* ptr 
	= dynamic_cast< const Step2Matrix<PA_GF2>* >(other);
      return new Step2Matrix<PA_GF2>(*ptr, /*invert=*/true);
    }
    case PA_zz_p_tag: {
      const Step2Matrix<PA_zz_p>* ptr 
	= dynamic_cast< const Step2Matrix<PA_zz_p>* >(other);
      return new Step2Matrix<PA_zz_p>(*ptr, /*invert=*/true);
    }
    default: return 0;
  }
}

template<class type>
Step2Matrix<type>::Step2Matrix(const EncryptedArray& _ea, long m1, long m2)
  : ea(_ea)
{ 
  // Find generator sets for Zm2* /<p^d>
  long d = ea.getDegree();
  long p = ea.getContext().zMStar.getP();
  long p2d = PowerMod(p % m2, d, m2);
  vector<long> gens, ords;     // Find generator sets for Zm2* /<p^d>
  assert(findGenerators(gens, ords, m2, p2d) == 1); // p^d==1 mod m2
  assert(gens.size()==1);      // Zm2* = Zm2* /<p^d> is cyclic

  long nSlots = ea.size();
  long ord2 = abs(ords[0]); // the order of Zm2* /<p^d> = Zm2*
  ord1 = nSlots / ord2;     // the order of Zm1* /<p>

  vector<long> T(ord2);     // representative set for Zm2*/ <p^d>
  T[0] = 1;
  T[1] = gens[0];
  for (long i=2; i<ord2; i++)
    T[i] = MulMod(T[i-1], gens[0], m2);

  // prepare a table of eta2^i for i=0,1,...,m2-1
  vector<RX> eta2Powers; // powers of eta2 = X^m1 mod G
  PolyPowers<type>::computePolyPowers(eta2Powers,
				      ea.getDerived(type()), m1, m2);

  const PAlgebraModDerived<type>& tab=ea.getContext().alMod.getDerived(type());
  RBak bak; bak.save();
  REBak bakE; bakE.save();
  tab.restoreContext();
  RE::init(ea.getDerived(type()).getG());
  A.SetDims(ord2, ord2);
  for (long i=0; i<ord2; i++) for (long j=0; j<ord2; j++) {
      long exp = MulMod(i, T[j], m2);
      convert(A[i][j], eta2Powers[exp]);
    }
}

// returns the inverse of the matrix encoded in the A data member
template<class type> void Step2Matrix<type>::invert(mat_RE& Ainv) const
{
  const PAlgebraModDerived<type>& tab=ea.getContext().alMod.getDerived(type());
  RBak bak; bak.save();
  REBak bakE; bakE.save();
  tab.restoreContext();
  RE::init(ea.getDerived(type()).getG());

  ppInvert(Ainv, A, tab.getZMStar().getP(), tab.getR()); // invert mod p^r
}





void testIt(long m1, long m2, long p, long r)
{
  long m=m1*m2;
  FHEcontext context(m, p, r);
  EncryptedArray ea(context, context.alMod.getFactorsOverZZ()[0]);

#if 0
  if (r>1) {
    const PAlgebraModDerived<PA_zz_p>& tab = context.alMod.getDerived(PA_zz_p());
    zz_pBak bak; bak.save(); tab.restoreContext();
    zz_pEBak ebak; ebak.save(); zz_pE::init(tab.getFactors()[0]);

    long n=0;
    //    mat_zz_pE A(INIT_SIZE, n, n);
    //    mat_zz_pE Ainv;
    mat_zz_p A(INIT_SIZE, n, n);
    mat_zz_p Ainv;
    for (long i=0; i<n; i++) for (long j=0; j<n; j++) random(A[i][j]);
    ppInvert(Ainv, A, context.zMStar.getP(), tab.getR());
    cout << "okay\n";
    exit(0);
  }
#endif

  // Compute the matrix for step 1 and its inverse
  cerr << " + Computing step-1 matrix..." << std::flush;
  PlaintextBlockMatrixBaseInterface* step1Matrix_pt
    = buildStep1Matrix(ea, m1, m2);
  cerr << "done\n + Inverting step-1 matrix..." << std::flush;
  PlaintextBlockMatrixBaseInterface* step1Inverse_pt
    = buildStep1Inverse(step1Matrix_pt);

  // Compute the matrix for step 2
  cerr << "done\n + Computing step-2 matrix..." << std::flush;
  PlaintextMatrixBaseInterface* step2Matrix_pt
    = buildStep2Matrix(ea, m1, m2);
  cerr << "done\n + Inverting step-2 matrix..." << std::flush;
  PlaintextMatrixBaseInterface* step2Inverse_pt
    = buildStep2Inverse(step2Matrix_pt);

  cerr << "done";
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

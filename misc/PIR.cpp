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
 * @file PIR.cpp
 * @brief A simple PIR implementation
 */
#include <cstdint>
#include <NTL/GF2X.h>
#include <NTL/BasicThreadPool.h>
#include "intraSlot.h"
#include "matmul.h"

#if (defined(__unix__) || defined(__unix) || defined(unix))
#include <sys/time.h>
#include <sys/resource.h>
#endif

NTL_CLIENT

// Simplistic un-optimized encoding/decoding routines
static void encodeBits(GF2X& poly,
                       const vector<uint8_t>& data, long idx, long d)
{
  FHE_TIMER_START;
  for (long i=0; i<d; i++, idx++) {
    // compute byte index from bit index
    long byteIdx = idx / 8;
    if (byteIdx < lsize(data)) {
      long coef = (data[byteIdx] >> (idx % 8))&1;
      NTL::SetCoeff(poly, i, coef);
    }
  }
}
static void decodeBits(vector<uint8_t>& data, long idx,
                       const GF2X& poly, long d)
{
  FHE_TIMER_START;
  for (long i=0; i<d; i++, idx++) {
    // compute byte index from bit index
    long byteIdx = idx / 8;
    if (byteIdx < lsize(data)) {
      uint8_t coef = NTL::conv<uint>(NTL::coeff(poly,i)) & 1;
      data[byteIdx] |= (coef << (idx % 8));
    }
  }
}

static void decodeDBentry(vector<uint8_t>& entry,
                          const vector<GF2X>& encoded, long d)
{
  long dbSize = divc(lsize(encoded)*d, 8); // size in bytes

}


// Encoding a "database" in a matrix, the number of entries in the
// database is upto the size of the 1st dimension, i.e. sizeOfDim(0).
// The number of bits per entry is upto phi(m)/sizeOfDim(0).
class DBMatrix : public MatMul1D_derived<PA_GF2> {
public:
  PA_INJECT(PA_GF2)

private:
  vector< vector< vector< RX > > > data;
  const EncryptedArray& ea;

public:
  virtual ~DBMatrix() {}
  DBMatrix(const EncryptedArray& _ea, const vector< vector<uint8_t> >&db):
    ea(_ea)
  {
    RBak bak; bak.save(); ea.getAlMod().restoreContext();
    long n = ea.size();
    long d = ea.getDegree();
    long D = ea.sizeOfDimension(0);

    // vector< vector<uint8_t> > db2(lsize(db)); // sanity check
    // for (long i=0; i<lsize(db); i++)
    //   db2[i].resize(db[i].size(), 0);

    // Allocate space for the matrices
    data.resize(n/D);
    for (long k=0; k < n/D; k++) {
      data[k].resize(D);
      for (long i=0; i<D; i++)
	data[k][i].resize(D, RX()); // initialize to D empty polynomials
    }

    /* Each database entry is encoded in upto n/D partial rows, each partial
     * row has D cells of d bits each. In the picture here D=3, n/D=2.
     *
     *  (X1 X2 X3 | X4 X5 X6)  // database entry X
     *  (Y1 Y2 Y3 | Y4 Y5 Y6)  // database entry Y
     *  (Z1 Z2 Z3 | Z4 Z5 Z6)  // database entry Z
     *
     *  For example, the database entry Z is encoded in the two row vectors
     *  (Z1 Z2 Z3), (Z4 Z5 Z6). Each Zi is a degree-d binary polynomial.
     */

    for (long k=0; k<n/D; k++) { // go over row vectors

      // Encode the next d*D bits of each entry db[i]
      for (long i=0; i<D; i++) {
	for (long j=0; j<D; j++) { // encode next d bit of db[i]
          long dataIdx= (k*D+j)*d; // how many bits were already encoded
          if (dataIdx < 8*lsize(db[i])) { // more bits to encode
            encodeBits(data[k][i][j], db[i], dataIdx, d);
            // decodeBits(db2[i], dataIdx, data[k][i][j], d);
          }
	}
      }
    }
    // for (long i=0; i<lsize(db); i++)
    //   for (long j=0; j<lsize(db[i]); j++)
    //     assert(db[i][j] == db2[i][j]);
  }

  const EncryptedArray& getEA() const override { return ea; }
  bool multipleTransforms() const override { return true; }
  long getDim() const override { return 0; }

  bool get(RX& out, long i, long j, long k) const override {
    long n = ea.size();
    long D = ea.sizeOfDimension(0);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    assert(k >= 0 && k < n/D);
    if (IsZero(data[k][i][j])) return true;
    out = data[k][i][j];
    return false;
  }
};


static MatMul1D*
buildDBMatrix(const EncryptedArray& ea, const vector< vector<uint8_t> >&db)
{
  FHE_TIMER_START;
  assert (ea.getTag()==PA_GF2_tag);
  return new DBMatrix(ea, db);
}

static vector< vector<uint8_t> >*
buildRandomDB(const EncryptedArray& ea)
{
  FHE_TIMER_START;
  assert (ea.getTag()==PA_GF2_tag);
  vector< vector<uint8_t> >* v = new(vector< vector<uint8_t> >);
  v->resize(ea.sizeOfDimension(0));
  long entrySize = ea.size()/(ea.sizeOfDimension(0)*8);
  for (long i=0; i<ea.sizeOfDimension(0); i++) {
    (*v)[i].resize(entrySize);
    for (long j=0; j<entrySize; j++)
      (*v)[i][j] = NTL::RandomBits_long(8);
  }
  return v;
}

bool DoTest(const EncryptedArray& ea,
            const SecKey& secretKey, bool minimal, bool verbose)
{
  // choose a random database that fits in one ciphertext
  std::unique_ptr<vector< vector<uint8_t> > > db(buildRandomDB(ea));

  // Encode the database as a matrix
  FHE_NTIMER_START(PIR_EncodeDBasMatrix);
  std::unique_ptr< MatMul1D > mat(buildDBMatrix(ea, *db));
  MatMul1D::ExecType mat_exec(*mat, minimal);
  mat_exec.upgrade();
  FHE_NTIMER_STOP(PIR_EncodeDBasMatrix);

  // Choose a random index into the database
  long idx = NTL::RandomBnd(lsize(*db));

  // Encrypt index as a unit vector
  FHE_NTIMER_START(PIR_EncryptSelectionAsVector);
  Ctxt ctxt(secretKey);
  {vector<long> slots(ea.size(), 0);
  for (long i=0; i<ea.size(); i++) {
    if (ea.coordinate(0,i) == idx)
      slots[i]=1;
  }
  ea.encrypt(ctxt, secretKey, slots);
  }
  FHE_NTIMER_STOP(PIR_EncryptSelectionAsVector);

  // Do the matrix-vector multiply
  mat_exec.mul(ctxt);

  // Decrypt and check the result
  FHE_NTIMER_START(PIR_Decrypt_Result);

  vector<ZZX> slots;
  ea.decrypt(ctxt, secretKey, slots); // decrypt the ciphertext vector

  // maximum size of an entry in bits
  long d = ea.getDegree();
  long nOverD = ea.size()/ea.sizeOfDimension(0);
  long entrySize = (ea.size()*d) / ea.sizeOfDimension(0);
  vector<uint8_t> dbEntry(divc(entrySize,8), 0);

  // Go over the slots and put each one in its right place
  for (long i=0; i<ea.size(); i++) {
    // represent i as (ii,jj) with ii index along dim 0 and jj the rest
    std::pair<long,long> pp = ea.getAlMod().getZMStar().breakIndexByDim(i,0);
    long idxx = pp.first*ea.sizeOfDimension(0) + pp.second;
    decodeBits(dbEntry, idxx*d, conv<GF2X>(slots[i]), d);
  }
  FHE_NTIMER_STOP(PIR_Decrypt_Result);

  // check that we've got the right answer
  for (long i=0; i<std::min(lsize(dbEntry), lsize((*db)[idx])); i++)
    if (dbEntry[i] != (*db)[idx][i]) {
      cout << "Grrr@*, entry["<<i<<"]="<<((int)dbEntry[i])<<"!=db["<<idx<<"]["
           <<i<<"]="<<((int)(*db)[idx][i])<<endl;
      return false;
    }
  return true;
}


int ks_strategy = 0;
// 0 == default
// 1 == full
// 2 == BSGS
// 3 == minimal


void  TestIt(Context& context, bool verbose)
{
  if (verbose) {
    context.zMStar.printout();
    std::cout <<"  #threads="<<NTL::AvailableThreads()
              <<", security=" << context.securityLevel()<<endl<< endl;
  }

  SecKey secretKey(context);
  const PubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64); // A Hamming-weight-w secret key

  bool minimal = (ks_strategy == 3);

  // we call addSomeFrbMatrices for all strategies except minimal

  switch (ks_strategy) {
  case 0:
    addSome1DMatrices(secretKey);
    addSomeFrbMatrices(secretKey);
    break;
  case 1:
    add1DMatrices(secretKey);
    addSomeFrbMatrices(secretKey);
    break;
  case 2:
    addBSGS1DMatrices(secretKey);
    addSomeFrbMatrices(secretKey);
    break;
  case 3:
    addMinimal1DMatrices(secretKey);
    addMinimalFrbMatrices(secretKey);
    break;

   default:
     Error("bad ks_strategy");
   }
  EncryptedArray ea(context, context.alMod);
  bool okSoFar=true;

  cout << " * "<<ea.sizeOfDimension(0)<<"-entry DB, |entry|="
       << (ea.size()*ea.getDegree()) << " bits. ";
  for (long i=0; i<5; i++)
    if (!DoTest(ea, secretKey, minimal, verbose)) {
      okSoFar = false;
      break;
    }
  if (okSoFar)
    cout << "Nice!!\n\n";

  if (verbose) {
    printAllTimers(cout);
#if (defined(__unix__) || defined(__unix) || defined(unix))
      struct rusage rusage;
      getrusage( RUSAGE_SELF, &rusage );
      cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif
  }

}


int main(int argc, char *argv[])
{
  ArgMapping amap;

  long m=2047;
  amap.arg("m", m, "defines the cyclotomic polynomial Phi_m(X)");
  long p=2;
  long L=3;
  amap.arg("L", L, "# of levels in the modulus chain");
  long verbose=0;
  amap.arg("verbose", verbose, "print timing and other info");
  long nt=1;
  amap.arg("nt", nt, "# threads");

  amap.arg("force_bsgs", fhe_test_force_bsgs,
           "1 to force on, -1 to force off");
  amap.arg("force_hoist", fhe_test_force_hoist,
           "-1 to force off");
  amap.arg("ks_strategy", ks_strategy,
           "0: default, 1:full, 2:bsgs, 3:minimal");

  NTL::Vec<long> gens;
  amap.arg("gens", gens, "use specified vector of generators", nullptr);
  amap.note("e.g., gens='[562 1871 751]'");
  NTL::Vec<long> ords;
  amap.arg("ords", ords, "use specified vector of orders", nullptr);
  amap.note("e.g., ords='[4 2 -4]', negative means 'bad'");

  amap.parse(argc, argv);

  if (verbose) {
    cout << "*** matmul1D: m=" << m
	 << ", p^r=2"
	 << ", L=" << L
	 << ", nt=" << nt
	 << ", force_bsgs=" << fhe_test_force_bsgs
	 << ", force_hoist=" << fhe_test_force_hoist
	 << ", ks_strategy=" << ks_strategy
	 << endl;
   }

  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);

  if (nt > 1) SetNumThreads(nt);

  setTimersOn();

  Context context(m, 2, 1, gens1, ords1);
  buildModChain(context, L, /*c=*/3);

  TestIt(context, verbose);
}

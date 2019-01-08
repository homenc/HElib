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
#include <NTL/ZZ.h>
#include <algorithm>
#include <complex>

#include "norms.h"
#include "EncryptedArray.h"
#include "FHE.h"
#include "debugging.h"

NTL_CLIENT

bool verbose=false;

// Compute the L-infinity distance between two vectors
double calcMaxDiff(const vector<cx_double>& v1, 
                   const vector<cx_double>& v2){

  if(lsize(v1)!=lsize(v2))
    NTL::Error("Vector sizes differ.\nFAILED\n");

  double maxDiff = 0.0;  
  for (long i=0; i<lsize(v1); i++) {
    double diffAbs = std::abs(v1[i]-v2[i]);
    if (diffAbs > maxDiff)
      maxDiff = diffAbs;
  }

  return maxDiff;
}

inline bool cx_equals(const vector<cx_double>& v1, 
                      const vector<cx_double>& v2, 
                      double epsilon)
{
  return (calcMaxDiff(v1,v2) < epsilon);
}



void testBasicArith(const FHEPubKey& publicKey, 
                    const FHESecKey& secretKey, 
                    const EncryptedArrayCx& ea, double epsilon);
void testComplexArith(const FHEPubKey& publicKey, 
                      const FHESecKey& secretKey, 
                      const EncryptedArrayCx& ea, double epsilon);
void testRotsNShifts(const FHEPubKey& publicKey, 
                     const FHESecKey& secretKey, 
                     const EncryptedArrayCx& ea, double epsilon);


int main(int argc, char *argv[]) 
{

  // Commandline setup

  ArgMapping amap;

  long m=16;
  long r=8;
  long L=150;
  double epsilon=0.01; // Accepted accuracy

  amap.arg("m", m, "Cyclotomic index");
  amap.note("e.g., m=1024, m=2047");
  amap.arg("r", r, "Bits of precision");
  amap.arg("L", L, "Number of levels");
  amap.arg("ep", epsilon, "Accepted accuracy");
  amap.arg("verbose", verbose, "more printouts");

  amap.parse(argc, argv);

  try{

    // FHE setup keys, context, SKMs, etc

    FHEcontext context(m, /*p=*/-1, r);
    buildModChain(context, L, /*c=*/2);

    FHESecKey secretKey(context);
    secretKey.GenSecKey(); // A +-1/0 secret key
    addSome1DMatrices(secretKey); // compute key-switching matrices

    const FHEPubKey publicKey = secretKey;
    const EncryptedArrayCx& ea = context.ea->getCx();

    if (verbose) {
      ea.getPAlgebra().printout();
      cout << "r = " << context.alMod.getR() << endl;
      cout << "ctxtPrimes="<<context.ctxtPrimes
           << ", specialPrimes="<<context.specialPrimes<<endl<<endl;
    }

    // Run the tests.
    testBasicArith(publicKey, secretKey, ea, epsilon);
    testComplexArith(publicKey, secretKey, ea, epsilon);
    testRotsNShifts(publicKey, secretKey, ea, epsilon);

  } 
  catch (exception& e) {
    cerr << e.what() << endl;
    cerr << "***Major FAIL***" << endl;  
  }

  return 0;
}


void testBasicArith(const FHEPubKey& publicKey,
                    const FHESecKey& secretKey,
                    const EncryptedArrayCx& ea, double epsilon)
{
  if (verbose)  cout << "Test Arithmetic ";
  // Test objects

  Ctxt c1(publicKey), c2(publicKey), c3(publicKey);
  
  vector<cx_double> vd;
  vector<cx_double> vd1, vd2, vd3;
  ea.random(vd1);
  ea.random(vd2);

  // test encoding of shorter vectors
  vd1.resize(vd1.size()-2);
  ea.encrypt(c1, publicKey, vd1);
  vd1.resize(vd1.size()+2, 0.0);

  ea.encrypt(c2, publicKey, vd2);


  // Test - Multiplication  
  c1 *= c2;
  for (long i=0; i<lsize(vd1); i++) vd1[i] *= vd2[i];

  ZZX poly;
  ea.random(vd3);
  ea.encode(poly,vd3);
  c1.addConstant(poly); // vd1*vd2 + vd3
  for (long i=0; i<lsize(vd1); i++) vd1[i] += vd3[i];

  // Test encoding, encryption of a single number
  double xx = NTL::RandomLen_long(16)/double(1L<<16); // random in [0,1]
  ea.encryptOneNum(c2, publicKey, xx);
  c1 += c2;
  for (auto& x : vd1) x += xx;

  // Test - Multiply by a mask
  vector<long> mask(lsize(vd1), 1);
  for (long i=0; i*(i+1)<lsize(mask); i++) {
    mask[i*i] = 0;
    mask[i*(i+1)] = -1;
  }

  ea.encode(poly,mask);
  c1.multByConstant(poly); // mask*(vd1*vd2 + vd3)
  for (long i=0; i<lsize(vd1); i++) vd1[i] *= mask[i];

  // Test - Addition
  ea.random(vd3);
  ea.encrypt(c3, publicKey, vd3);
  c1 += c3;
  for (long i=0; i<lsize(vd1); i++) vd1[i] += vd3[i];

  c1.negate();
  c1.addConstant(to_ZZ(1));
  for (long i=0; i<lsize(vd1); i++) vd1[i] = 1.0 - vd1[i];

  // Diff between approxNums HE scheme and plaintext floating  
  ea.decrypt(c1, secretKey, vd);
#ifdef DEBUG_PRINTOUT
  printVec(cout<<"res=", vd, 10)<<endl;
  printVec(cout<<"vec=", vd1, 10)<<endl;
#endif
  if (verbose)
    cout << "(max |res-vec|_{infty}="<< calcMaxDiff(vd, vd1) << "): ";

  cx_equals(vd, vd1, epsilon)?
    cout << "GOOD\n":
    cout << "BAD\n";
}


void testComplexArith(const FHEPubKey& publicKey,
                      const FHESecKey& secretKey,
                      const EncryptedArrayCx& ea, double epsilon)
{

  // Test complex conjugate
  Ctxt c1(publicKey), c2(publicKey);

  vector<cx_double> vd;
  vector<cx_double> vd1, vd2;
  ea.random(vd1);
  ea.random(vd2);
   
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);

  if (verbose)
    cout << "Test Conjugate: ";
  for_each(vd1.begin(), vd1.end(), [](cx_double& d){d=std::conj(d);});
  c1.complexConj();  
  ea.decrypt(c1, secretKey, vd);
#ifdef DEBUG_PRINTOUT
  printVec(cout<<"vd1=", vd1, 10)<<endl;
  printVec(cout<<"res=", vd, 10)<<endl;
#endif
  cx_equals(vd, vd1, epsilon)?
    cout << "GOOD\n":
    cout << "BAD\n";

  // Test that real and imaginary parts are actually extracted.
  Ctxt realCtxt(c2), imCtxt(c2);
  vector<cx_double> realParts(vd2), real_dec;
  vector<cx_double> imParts(vd2), im_dec;

  if (verbose)
    cout << "Test Real and Im parts: ";
  for_each(realParts.begin(), realParts.end(), [](cx_double& d){d=std::real(d);});
  for_each(imParts.begin(), imParts.end(), [](cx_double& d){d=std::imag(d);});

  ea.extractRealPart(realCtxt);
  ea.decrypt(realCtxt, secretKey, real_dec);

  ea.extractImPart(imCtxt);
  ea.decrypt(imCtxt, secretKey, im_dec);

#ifdef DEBUG_PRINTOUT
  printVec(cout<<"vd2=", vd2, 10)<<endl;
  printVec(cout<<"real=", realParts, 10)<<endl;
  printVec(cout<<"res=", real_dec, 10)<<endl;
  printVec(cout<<"im=", imParts, 10)<<endl;
  printVec(cout<<"res=", im_dec, 10)<<endl;
#endif
  cx_equals(realParts,real_dec,epsilon) && cx_equals(imParts,im_dec,epsilon)?
    cout << "GOOD\n":
    cout << "BAD\n";
}

void testRotsNShifts(const FHEPubKey& publicKey, 
                     const FHESecKey& secretKey,
                     const EncryptedArrayCx& ea, double epsilon)
{

  std::srand(std::time(0)); // set seed, current time.
  int nplaces = rand() % static_cast<int>(ea.size()/2.0) + 1;

  if (verbose)
    cout << "Test Rotation of " << nplaces << ": ";  

  Ctxt c1(publicKey);
  vector<cx_double> vd1;
  vector<cx_double> vd_dec;
  ea.random(vd1);
  ea.encrypt(c1, publicKey, vd1);

#ifdef DEBUG_PRINTOUT
  printVec(cout<< "vd1=", vd1, 10)<<endl;
#endif
  std::rotate(vd1.begin(), vd1.end()-nplaces, vd1.end());
  ea.rotate(c1, nplaces);
  ea.decrypt(c1, secretKey, vd_dec);
#ifdef DEBUG_PRINTOUT
  printVec(cout<< "vd1(rot)=", vd1, 10)<<endl;
  printVec(cout<<"res: ", vd_dec, 10)<<endl;
#endif

  cx_equals(vd1, vd_dec, epsilon)?
    cout << "GOOD\n":
    cout << "BAD\n";
}

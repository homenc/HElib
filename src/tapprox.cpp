/* Copyright (C) 2012-2020 IBM Corp.
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

#include <helib/norms.h>
#include <helib/helib.h>
#include <helib/debugging.h>
#include <helib/ArgMap.h>

NTL_CLIENT
using namespace helib;

bool verbose = false;

bool reset = false;

// Compute the L-infinity distance between two vectors
double calcMaxDiff(const vector<cx_double>& v1, const vector<cx_double>& v2)
{

  if (lsize(v1) != lsize(v2))
    NTL::Error("Vector sizes differ.\nFAILED\n");

  double maxDiff = 0.0;
  for (long i = 0; i < lsize(v1); i++) {
    double diffAbs = std::abs(v1[i] - v2[i]);
    if (diffAbs > maxDiff)
      maxDiff = diffAbs;
  }

  return maxDiff;
}
// Compute the max relative difference between two vectors
double calcMaxRelDiff(const vector<cx_double>& v1, const vector<cx_double>& v2)
{
  if (lsize(v1) != lsize(v2))
    NTL::Error("Vector sizes differ.\nFAILED\n");

  // Compute the largest-magnitude value in the vector
  double maxAbs = 0.0;
  for (auto& x : v1) {
    if (std::abs(x) > maxAbs)
      maxAbs = std::abs(x);
  }
  if (maxAbs < 1e-10)
    maxAbs = 1e-10;

  double maxDiff = 0.0;
  for (long i = 0; i < lsize(v1); i++) {
    double relDiff = std::abs(v1[i] - v2[i]) / maxAbs;
    if (relDiff > maxDiff)
      maxDiff = relDiff;
  }

  return maxDiff;
}

inline bool cx_equals(const vector<cx_double>& v1,
                      const vector<cx_double>& v2,
                      double epsilon)
{
  return (calcMaxRelDiff(v1, v2) < epsilon);
}

void testBasicArith(const PubKey& publicKey,
                    const SecKey& secretKey,
                    const EncryptedArrayCx& ea,
                    double epsilon);
void testComplexArith(const PubKey& publicKey,
                      const SecKey& secretKey,
                      const EncryptedArrayCx& ea,
                      double epsilon);
void testRotsNShifts(const PubKey& publicKey,
                     const SecKey& secretKey,
                     const EncryptedArrayCx& ea,
                     double epsilon);

void debugCompare(const EncryptedArrayCx& ea,
                  const SecKey& sk,
                  vector<cx_double>& p,
                  const Ctxt& c,
                  double epsilon)
{
  double maxAbs = 0.0;
  for (auto& x : p) {
    if (std::abs(x) > maxAbs)
      maxAbs = std::abs(x);
  }
  vector<cx_double> pp;
  ea.decrypt(c, sk, pp);

  double err = calcMaxDiff(p, pp);
  double err_bnd = NTL::conv<double>(c.getNoiseBound() / c.getRatFactor());
  double log2_ratio = std::log(err / err_bnd) / std::log(2.0);
  std::cout << "    "
            //<< " relative-error="<<calcMaxRelDiff(p,pp)
            << " err=" << err << " ptxtMag=" << c.getPtxtMag()
            << " maxAbs=" << maxAbs << " ratFactor="
            << c.getRatFactor()
            //<< " noiseBound=" << c.getNoiseBound()
            << " log2(err/bound)=" << log2_ratio
            << " log2(ptxt/bound)=" << (log(maxAbs / c.getPtxtMag()) / log(2.0))
            << endl;
  //  if (!cx_equals(pp, p, epsilon)) {
  //    std::cout << "oops:\n"; std::cout << p << "\n";
  //    std::cout << pp << "\n";
  //    exit(0);
  //  }
}

void negateVec(vector<cx_double>& p1)
{
  for (auto& x : p1)
    x = -x;
}
void add(vector<cx_double>& to, const vector<cx_double>& from)
{
  if (to.size() < from.size())
    to.resize(from.size(), 0);
  for (long i = 0; i < from.size(); i++)
    to[i] += from[i];
}
void sub(vector<cx_double>& to, const vector<cx_double>& from)
{
  if (to.size() < from.size())
    to.resize(from.size(), 0);
  for (long i = 0; i < from.size(); i++)
    to[i] -= from[i];
}
void mul(vector<cx_double>& to, const vector<cx_double>& from)
{
  if (to.size() < from.size())
    to.resize(from.size(), 0);
  for (long i = 0; i < from.size(); i++)
    to[i] *= from[i];
}
void rotate(vector<cx_double>& p, long amt)
{
  long sz = p.size();
  vector<cx_double> tmp(sz);
  for (long i = 0; i < sz; i++)
    tmp[((i + amt) % sz + sz) % sz] = p[i];
  p = tmp;
}

void resetPtxtMag(Ctxt& c, const vector<cx_double>& p)
{
  double maxAbs = 0.0;
  for (auto& x : p) {
    if (std::abs(x) > maxAbs)
      maxAbs = std::abs(x);
  }

  if (maxAbs < 1.0)
    maxAbs = 1.0;
  else
    maxAbs = std::pow(
        2,
        std::ceil(std::log(maxAbs) / std::log(2))); // next power of two

  c.setPtxtMag(NTL::xdouble(maxAbs));
}

/************** Each round consists of the following:
tmp1 = rotate(c0)
tmp1 += const1
c0 += const2
c0 *= tmp1  // c0 = (rotate(c0) + const1)*(c0 + const2) ...squared every round

tmp2 = c1 * const1
c1 = rotate(c1)
c1 += tmp2  // c1 = rotate(c1) + c1*const1 ...doubled every round


tmp3 = c2 * const2
c2 *= c3
c2 += tmp3 // c2 = = c2*c3 + c2*const2 = c2*(c3 + const2)

c3 = c3*const1
**************/

#define DEBUG_COMPARE(C, P, M)                                                 \
  do {                                                                         \
    if (verbose) {                                                             \
      CheckCtxt(C, M);                                                         \
      debugCompare(ea, secretKey, P, C, epsilon);                              \
    }                                                                          \
  } while (0)

void testGeneralOps(const PubKey& publicKey,
                    const SecKey& secretKey,
                    const EncryptedArrayCx& ea,
                    double epsilon,
                    long nRounds)
{
  long nslots = ea.size();
  char buffer[32];

  vector<cx_double> p0, p1, p2, p3;
  ea.random(p0);
  ea.random(p1);
  ea.random(p2);
  ea.random(p3);

  Ctxt c0(publicKey), c1(publicKey), c2(publicKey), c3(publicKey);
  ea.encrypt(c0, publicKey, p0, /*size=*/1.0);
  ea.encrypt(c1, publicKey, p1, /*size=*/1.0);
  ea.encrypt(c2, publicKey, p2, /*size=*/1.0);
  ea.encrypt(c3, publicKey, p3, /*size=*/1.0);

  resetAllTimers();
  HELIB_NTIMER_START(Circuit);

  for (long i = 0; i < nRounds; i++) {

    if (verbose)
      std::cout << "*** round " << i << "..." << endl;

    if (reset) {
      resetPtxtMag(c0, p0);
      resetPtxtMag(c1, p1);
      resetPtxtMag(c2, p2);
      resetPtxtMag(c3, p3);
    }

    long rotamt = RandomBnd(2 * nslots - 1) - (nslots - 1);
    // random number in [-(nslots-1)..nslots-1]

    // two random constants
    vector<cx_double> const1, const2;
    ea.random(const1);
    ea.random(const2);

    ZZX const1_poly, const2_poly;
    ea.encode(const1_poly, const1, /*size=*/1.0);
    ea.encode(const2_poly, const2, /*size=*/1.0);

    vector<cx_double> tmp1_p(p0);
    rotate(tmp1_p, rotamt);
    Ctxt tmp1(c0);
    ea.rotate(tmp1, rotamt);
    DEBUG_COMPARE(tmp1, tmp1_p, "tmp1 = rotate(c0)");

    add(tmp1_p, const1);
    tmp1.addConstant(const1_poly);
    DEBUG_COMPARE(tmp1, tmp1_p, "tmp1 += const1");

    add(p0, const2);
    c0.addConstant(const2_poly);
    DEBUG_COMPARE(c0, p0, "c0 += const2");

    mul(p0, tmp1_p);
    c0.multiplyBy(tmp1);
    DEBUG_COMPARE(c0, p0, "c0 *= tmp1");

    vector<cx_double> tmp2_p(p1);
    mul(tmp2_p, const1);
    Ctxt tmp2(c1);
    tmp2.multByConstant(const1_poly);
    DEBUG_COMPARE(tmp2, tmp2_p, "tmp2 = c1 * const1");

    rotate(p1, rotamt);
    ea.rotate(c1, rotamt);
    DEBUG_COMPARE(c1, p1, "c1 = rotate(c1)");

    add(p1, tmp2_p);
    c1 += tmp2;
    DEBUG_COMPARE(c1, p1, "c1 += tmp2");

    vector<cx_double> tmp3_p(p2);
    mul(tmp3_p, const2);
    Ctxt tmp3(c2);
    tmp3.multByConstant(const2_poly);
    DEBUG_COMPARE(tmp3, tmp3_p, "tmp3 = c2 * const2");

    mul(p2, p3);
    c2.multiplyBy(c3);
    DEBUG_COMPARE(c2, p2, "c2 *= c3");

    add(p2, tmp3_p);
    c2.addCtxt(tmp3);
    DEBUG_COMPARE(c2, p2, "c2 += tmp3");

    mul(p3, const1);
    c3.multByConstant(const1_poly);
    DEBUG_COMPARE(c3, p3, "c3 *= const1");

    if (verbose) {
      // Check correctness after each round
      vector<cx_double> pp0, pp1, pp2, pp3;

      ea.decrypt(c0, secretKey, pp0);
      ea.decrypt(c1, secretKey, pp1);
      ea.decrypt(c2, secretKey, pp2);
      ea.decrypt(c3, secretKey, pp3);

      if (!(cx_equals(pp0, p0, epsilon) && cx_equals(pp1, p1, epsilon) &&
            cx_equals(pp2, p2, epsilon) && cx_equals(pp3, p3, epsilon))) {
        std::cout << "FAIL AT ROUND " << i << "\n";
        break;
      }
    }
  }

  c0.cleanUp();
  c1.cleanUp();
  c2.cleanUp();
  c3.cleanUp();

  HELIB_NTIMER_STOP(Circuit);

  vector<cx_double> pp0, pp1, pp2, pp3;

  ea.decrypt(c0, secretKey, pp0);
  ea.decrypt(c1, secretKey, pp1);
  ea.decrypt(c2, secretKey, pp2);
  ea.decrypt(c3, secretKey, pp3);

  std::cout << "Test " << nRounds << " rounds of mixed operations, ";
  if (cx_equals(pp0, p0, epsilon) && cx_equals(pp1, p1, epsilon) &&
      cx_equals(pp2, p2, epsilon) && cx_equals(pp3, p3, epsilon))
    std::cout << "PASS\n\n";
  else {
    std::cout << "FAIL\n\n";
    //    std::cout << "  max(p0)="<<largestCoeff(p0)
    //              << ", max(pp0)="<<largestCoeff(pp0)
    //              << ", maxDiff="<<calcMaxDiff(p0,pp0) << endl;
    //    std::cout << "  max(p1)="<<largestCoeff(p1)
    //              << ", max(pp1)="<<largestCoeff(pp1)
    //              << ", maxDiff="<<calcMaxDiff(p1,pp1) << endl;
    //    std::cout << "  max(p2)="<<largestCoeff(p2)
    //              << ", max(pp2)="<<largestCoeff(pp2)
    //              << ", maxDiff="<<calcMaxDiff(p2,pp2) << endl;
    //    std::cout << "  max(p3)="<<largestCoeff(p3)
    //              << ", max(pp3)="<<largestCoeff(pp3)
    //              << ", maxDiff="<<calcMaxDiff(p3,pp3) << endl<<endl;
  }

  if (verbose) {
    std::cout << endl;
    // printAllTimers();
    std::cout << endl;
  }
  resetAllTimers();
}

int main(int argc, char* argv[])
{

  // Commandline setup

  ArgMap amap;

  long m = 16;
  long r = 8;
  long L = 0;
  double epsilon = 0.01; // Accepted accuracy
  long R = 1;
  long seed = 0;
  bool debug = false;

  amap.arg("m", m, "Cyclotomic index");
  amap.note("e.g., m=1024, m=2047");
  amap.arg("r", r, "Bits of precision");
  amap.arg("R", R, "number of rounds");
  amap.arg("L", L, "Number of bits in modulus", "heuristic");
  amap.arg("ep", epsilon, "Accepted accuracy");
  amap.arg("seed", seed, "PRG seed");
  amap.arg("verbose", verbose, "more printouts");
  amap.arg("debug", debug, "for debugging");
  amap.arg("reset", reset, "forces calls to setPtxtMag each round");

  amap.parse(argc, argv);

  if (seed)
    NTL::SetSeed(ZZ(seed));

  if (R <= 0)
    R = 1;
  if (L == 0) {
    if (R <= 2)
      L = 100 * R;
    else
      L = 220 * (R - 1);
  }

  if (verbose) {
    cout << "** m=" << m << ", #rounds=" << R << ", |q|=" << L
         << ", epsilon=" << epsilon << endl;
  }
  try {

    // FHE setup keys, context, SKMs, etc

    Context context(m, /*p=*/-1, r);
    context.scale = 4;
    buildModChain(context, L, /*c=*/2);

    SecKey secretKey(context);
    secretKey.GenSecKey();        // A +-1/0 secret key
    addSome1DMatrices(secretKey); // compute key-switching matrices

    const PubKey publicKey = secretKey;
    const EncryptedArrayCx& ea = context.ea->getCx();

    if (verbose) {
      std::cout << "security=" << context.securityLevel() << endl;
      ea.getPAlgebra().printout();
      cout << "r = " << context.alMod.getR() << endl;
      cout << "ctxtPrimes=" << context.ctxtPrimes
           << ", specialPrimes=" << context.specialPrimes << endl
           << endl;
    }
    if (debug) {
      dbgKey = &secretKey;
      dbgEa = context.ea;
    }
#ifdef HELIB_DEBUG
    dbgKey = &secretKey;
    dbgEa = context.ea;
#endif // HELIB_DEBUG

    // Run the tests.
    testBasicArith(publicKey, secretKey, ea, epsilon);
    testComplexArith(publicKey, secretKey, ea, epsilon);
    testRotsNShifts(publicKey, secretKey, ea, epsilon);
    testGeneralOps(publicKey, secretKey, ea, epsilon, R);
  } catch (exception& e) {
    cerr << e.what() << endl;
    cerr << "***Major FAIL***" << endl;
  }

  return 0;
}

void testBasicArith(const PubKey& publicKey,
                    const SecKey& secretKey,
                    const EncryptedArrayCx& ea,
                    double epsilon)
{
  if (verbose)
    cout << "Test Arithmetic ";
  // Test objects

  Ctxt c1(publicKey), c2(publicKey), c3(publicKey);

  vector<cx_double> vd;
  vector<cx_double> vd1, vd2, vd3;
  ea.random(vd1);
  ea.random(vd2);

  // test encoding of shorter vectors
  vd1.resize(vd1.size() - 2);
  ea.encrypt(c1, publicKey, vd1, /*size=*/1.0);
  vd1.resize(vd1.size() + 2, 0.0);

  ea.encrypt(c2, publicKey, vd2, /*size=*/1.0);

  // Test - Multiplication
  c1 *= c2;
  for (long i = 0; i < lsize(vd1); i++)
    vd1[i] *= vd2[i];

  ZZX poly;
  ea.random(vd3);
  ea.encode(poly, vd3, /*size=*/1.0);
  c1.addConstant(poly); // vd1*vd2 + vd3
  for (long i = 0; i < lsize(vd1); i++)
    vd1[i] += vd3[i];

  // Test encoding, encryption of a single number
  double xx = NTL::RandomLen_long(16) / double(1L << 16); // random in [0,1]
  ea.encryptOneNum(c2, publicKey, xx);
  c1 += c2;
  for (auto& x : vd1)
    x += xx;

  // Test - Multiply by a mask
  vector<long> mask(lsize(vd1), 1);
  for (long i = 0; i * (i + 1) < lsize(mask); i++) {
    mask[i * i] = 0;
    mask[i * (i + 1)] = -1;
  }

  ea.encode(poly, mask, /*size=*/1.0);
  c1.multByConstant(poly); // mask*(vd1*vd2 + vd3)
  for (long i = 0; i < lsize(vd1); i++)
    vd1[i] *= mask[i];

  // Test - Addition
  ea.random(vd3);
  ea.encrypt(c3, publicKey, vd3, /*size=*/1.0);
  c1 += c3;
  for (long i = 0; i < lsize(vd1); i++)
    vd1[i] += vd3[i];

  c1.negate();
  c1.addConstant(to_ZZ(1));
  for (long i = 0; i < lsize(vd1); i++)
    vd1[i] = 1.0 - vd1[i];

  // Diff between approxNums HE scheme and plaintext floating
  ea.decrypt(c1, secretKey, vd);
#ifdef HELIB_DEBUG
  printVec(cout << "res=", vd, 10) << endl;
  printVec(cout << "vec=", vd1, 10) << endl;
#endif
  if (verbose)
    cout << "(max |res-vec|_{infty}=" << calcMaxDiff(vd, vd1) << "): ";

  if (cx_equals(vd, vd1, conv<double>(epsilon * c1.getPtxtMag())))
    cout << "GOOD\n";
  else {
    cout << "BAD:\n";
    std::cout << "  max(vd)=" << largestCoeff(vd)
              << ", max(vd1)=" << largestCoeff(vd1)
              << ", maxDiff=" << calcMaxDiff(vd, vd1) << endl
              << endl;
  }
}

void testComplexArith(const PubKey& publicKey,
                      const SecKey& secretKey,
                      const EncryptedArrayCx& ea,
                      double epsilon)
{

  // Test complex conjugate
  Ctxt c1(publicKey), c2(publicKey);

  vector<cx_double> vd;
  vector<cx_double> vd1, vd2;
  ea.random(vd1);
  ea.random(vd2);

  ea.encrypt(c1, publicKey, vd1, /*size=*/1.0);
  ea.encrypt(c2, publicKey, vd2, /*size=*/1.0);

  if (verbose)
    cout << "Test Conjugate: ";
  for_each(vd1.begin(), vd1.end(), [](cx_double& d) { d = std::conj(d); });
  c1.complexConj();
  ea.decrypt(c1, secretKey, vd);
#ifdef HELIB_DEBUG
  printVec(cout << "vd1=", vd1, 10) << endl;
  printVec(cout << "res=", vd, 10) << endl;
#endif
  if (cx_equals(vd, vd1, conv<double>(epsilon * c1.getPtxtMag())))
    cout << "GOOD\n";
  else {
    cout << "BAD:\n";
    std::cout << "  max(vd)=" << largestCoeff(vd)
              << ", max(vd1)=" << largestCoeff(vd1)
              << ", maxDiff=" << calcMaxDiff(vd, vd1) << endl
              << endl;
  }

  // Test that real and imaginary parts are actually extracted.
  Ctxt realCtxt(c2), imCtxt(c2);
  vector<cx_double> realParts(vd2), real_dec;
  vector<cx_double> imParts(vd2), im_dec;

  if (verbose)
    cout << "Test Real and Im parts: ";
  for_each(realParts.begin(), realParts.end(), [](cx_double& d) {
    d = std::real(d);
  });
  for_each(imParts.begin(), imParts.end(), [](cx_double& d) {
    d = std::imag(d);
  });

  ea.extractRealPart(realCtxt);
  ea.decrypt(realCtxt, secretKey, real_dec);

  ea.extractImPart(imCtxt);
  ea.decrypt(imCtxt, secretKey, im_dec);

#ifdef HELIB_DEBUG
  printVec(cout << "vd2=", vd2, 10) << endl;
  printVec(cout << "real=", realParts, 10) << endl;
  printVec(cout << "res=", real_dec, 10) << endl;
  printVec(cout << "im=", imParts, 10) << endl;
  printVec(cout << "res=", im_dec, 10) << endl;
#endif
  if (cx_equals(realParts,
                real_dec,
                conv<double>(epsilon * realCtxt.getPtxtMag())) &&
      cx_equals(imParts, im_dec, conv<double>(epsilon * imCtxt.getPtxtMag())))
    cout << "GOOD\n";
  else {
    cout << "BAD:\n";
    std::cout << "  max(re)=" << largestCoeff(realParts)
              << ", max(re1)=" << largestCoeff(real_dec)
              << ", maxDiff=" << calcMaxDiff(realParts, real_dec) << endl;
    std::cout << "  max(im)=" << largestCoeff(imParts)
              << ", max(im1)=" << largestCoeff(im_dec)
              << ", maxDiff=" << calcMaxDiff(imParts, im_dec) << endl
              << endl;
  }
}

void testRotsNShifts(const PubKey& publicKey,
                     const SecKey& secretKey,
                     const EncryptedArrayCx& ea,
                     double epsilon)
{

  long nplaces = NTL::RandomBnd(ea.size() / 2) + 1;

  if (verbose)
    cout << "Test Rotation of " << nplaces << ": ";

  Ctxt c1(publicKey);
  vector<cx_double> vd1;
  vector<cx_double> vd_dec;
  ea.random(vd1);
  ea.encrypt(c1, publicKey, vd1, /*size=*/1.0);

#ifdef HELIB_DEBUG
  printVec(cout << "vd1=", vd1, 10) << endl;
#endif
  std::rotate(vd1.begin(), vd1.end() - nplaces, vd1.end());
  ea.rotate(c1, nplaces);
  c1.reLinearize();
  ea.decrypt(c1, secretKey, vd_dec);
#ifdef HELIB_DEBUG
  printVec(cout << "vd1(rot)=", vd1, 10) << endl;
  printVec(cout << "res: ", vd_dec, 10) << endl;
#endif

  if (cx_equals(vd1, vd_dec, conv<double>(epsilon * c1.getPtxtMag())))
    cout << "GOOD\n";
  else {
    cout << "BAD:\n";
    std::cout << "  max(vd)=" << largestCoeff(vd_dec)
              << ", max(vd1)=" << largestCoeff(vd1)
              << ", maxDiff=" << calcMaxDiff(vd_dec, vd1) << endl
              << endl;
  }
}

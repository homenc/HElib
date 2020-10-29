/* Copyright (C) 2020 HElib Project
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

// This program suggests parameter-sets for CKKS, it takes as input the
// target bit-size of the modulus for fresh ciphertexts (which is loosely
// related to the bit-capacity of fresh ciphertexts), and potentially also
// target security level in the ranage 70-270 (default is 128) and target
// accuracy parameter in the range 8-32 (default is 12).
// It returns one or more suggestions for the parameters m, r, L, c to
// use to get this capacity and accuracy at the given security level.
// The parameters r, L will be roughly equal to the target accuracy and
// bitsize arguments in the input, but may not be exactly equal to them.
#include <cmath>
#include <iostream>
#include <helib/helib.h>
#include <helib/ArgMap.h>

using namespace std;


struct CKKSparams {
  int m, L, c;
  int nCtxtPrimes;
  double ctxtBits;

  // Look for the largest value of L (modulus bitlength) for the given
  // m,c, and security, Returns true if snything is found, else false.
  bool tryParameters(int security, int theM, int c2try) {

    // Set the "accruacy parameter" r to 16, it should not have any effect
    // The "plaintext space" parameter p is set to -1 for CKKS
    helib::Context theContext(theM, /*p=*/-1, /*r=*/16);
    
    // Note: teContext has an empty modulus chain, we now try to build many
    // chains with different values of L and check their security

    // Get largest modulus size (including special primes) for these value
    // of m and security. The formula from Context.h is
    //    security ~ 3.8*(phi(m)/log2(Q/sigma)) -20,
    // we overshoot |Q| by ignoring sigma and the -20 and rounding up 3.8 to 4

    int phim = theContext.zMStar.getPhiM();
    int totalBits = floor(double(4*phim)/security); // an over-estimate
    while (helib::lweEstimateSecurity(phim, double(totalBits-3), /*hwt=*/0) < security)
      totalBits--;

    if (totalBits < 50) // no point in doing this, can't suport even 50 bits
      return false;

    // count down until you find something that actually works
    this->L = 0;
    this->ctxtBits = 0;
    for (int L2try=totalBits; L2try>40; --L2try) {
      // create a new context and build a chain for it
      //helib::Context cntxt(theM, /*p=*/-1, /*r=*/16);
      helib::Context& cntxt = theContext;
      cntxt.clearModChain();
      helib::buildModChain(cntxt, L2try, c2try);
      if (cntxt.securityLevel() >= security) {
        this->ctxtBits = cntxt.logOfProduct(cntxt.ctxtPrimes)/log(2.0);
        this->L = L2try;
        this->nCtxtPrimes = cntxt.ctxtPrimes.card();
        break;
      }
    }
    if (this->L==0) // nothing was found
      return false;

    // Keep trying a few more values of L, since it is not clear that
    // reducing L always reduces the actual number of bits
    int betterL = 0;
    for (int delta=1; delta<10; delta++) {
      //helib::Context cntxt(theM, /*p=*/-1, /*r=*/16);
      helib::Context& cntxt = theContext;
      cntxt.clearModChain();
      helib::buildModChain(cntxt, this->L -delta, c2try);
      double realBits = cntxt.logOfProduct(cntxt.ctxtPrimes)/log(2.0);
      if (realBits > this->ctxtBits) {
        betterL = this->L -delta;
        this->ctxtBits = realBits;
        this->nCtxtPrimes = cntxt.ctxtPrimes.card();
      }
    }

    // Print a warning if a smaller L was found
    if (betterL > 0) {
      std::cerr << "** weird: L="<<betterL<<" gives more real bits than L="<<L
        << " (m="<<theM<<", c="<<c2try <<')'<< std::endl;
      L = betterL;
    }

    // Record the m,c values that were used here
    this->m = theM;
    this->c = c2try;
    return true;
  }
};

int main(int argc, char* argv[])
{
  // get optional security level from command-line arguments
  helib::ArgMap amap;
  int sec = 128;
  amap.arg("security", sec, "target security level in [70,270]");
  amap.parse(argc, argv);

  // ensure that security level is in the range [70,270]
  if (sec<70)       sec=70;
  else if (sec>270) sec=270;

  // print resutls in CSV format
  std::cout << "sec,m,c,L,bits" << std::endl;

  // Try power-of-two values of m from 2^11 upto 2^18
  for (int m=2048; m<524288; m*=2) {
    CKKSparams pprm; // pprm is "previous parameters"

    // Start with a ridicoulusly large value of c
    if (!pprm.tryParameters(sec, m, /*c=*/100))
      continue; // even this large c value does not work

    // get a more sensible upper-bound on c
    int maxC = std::min(pprm.nCtxtPrimes, 100);
    if (maxC == 1) // FIXME: this should not happen, but it does
      continue;
    if (maxC<100 && !pprm.tryParameters(sec, m, maxC))
        continue; // even this large c value does not work

    // Try successvely smaller values of c, down to c=2. Print the
    // previous params pprms every time the smaller c values results
    // is less modulus bits
    for (int c2try=maxC*0.86; c2try>1; c2try=0.86*c2try) {
      CKKSparams prm;
      if (prm.tryParameters(sec, m, c2try)) {
        if (pprm.L > prm.L) // the smaller c values is worse than before
          // print the previous (better) set of values
          std::cout<<sec<<','<<pprm.m<<','<<pprm.c<<','<<pprm.L<<','
            << round(pprm.ctxtBits)<<std::endl;
        pprm = prm; // record the current values for next iteration
      }
      else { // c value too small, no solution
        break;
      }
    }
    // print the last value that was recorded
    std::cout<<sec<<','<<pprm.m<<','<<pprm.c<<','<<pprm.L<<','
      << round(pprm.ctxtBits)<<std::endl;
  }
  return 0;
}

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
/* KeySwitchign.cpp - A few strategies for generating key-switching matrices
 *
 * Copyright IBM Corporation 2012 All rights reserved.
 */
#include <unordered_set>
#include "NTL/ZZ.h"
NTL_CLIENT
#include "FHE.h"
#include "permutations.h"

long KSGiantStepSize(long D)
{
  assert(D > 0);
  long g = SqrRoot(D);
  if (g*g < D) g++;  // g = ceiling(sqrt(D))
  return g;
}

// A maximalistic approach: generate matrices s(X^e)->s(X) for all e \in Zm*
void addAllMatrices(FHESecKey& sKey, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  long m = context.zMStar.getM();

  // key-switching matrices for the automorphisms
  for (long i = 0; i < m; i++) {
    if (!context.zMStar.inZmStar(i)) continue;
    sKey.GenKeySWmatrix(1, i, keyID, keyID);
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

// generate matrices s.t. you can reLinearize each s(X^e) in at most two steps
void addFewMatrices(FHESecKey& sKey, long keyID)
{
  NTL::Error("addFewMatrices Not implemented yet");
}

// This code block appears at least twice below
#define computeParams(context,m,i)\
  bool native;\
  long ord, gi, g2md;\
  NTL::mulmod_precon_t g2mdminv;\
  if (i==context.zMStar.numOfGens()) { /* Frobenius matrices */\
    ord = context.zMStar.getOrdP();\
    gi = context.zMStar.getP();\
    native = true;\
  }\
  else { /* one of the "regular" dimensions */\
    ord = context.zMStar.OrderOf(i);\
    gi = context.zMStar.ZmStarGen(i);\
    native = context.zMStar.SameOrd(i);\
    if (!native) {\
      g2md = PowerMod(gi,-ord,m); /* g^{-ord} mod m */\
      g2mdminv = PrepMulModPrecon(g2md, m);\
    }\
  }\
  NTL::mulmod_precon_t giminv = PrepMulModPrecon(gi, m);


#if 0
static void add1Dmats4dim(FHESecKey& sKey, long i, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  long m = context.zMStar.getM();
  computeParams(context,m,i); // defines vars: native, ord, gi, g2md, giminv, g2mdminv

  /* MAUTO vector<long> vals; */
  for (long j=1,val=gi; j < ord; j++) {
    // From s(X^val) to s(X)
    sKey.GenKeySWmatrix(1, val, keyID, keyID);
    if (!native) { // also from s(X^{g^{i-ord}}) to s(X)
      long val2 = MulModPrecon(val,g2md,m,g2mdminv);
      sKey.GenKeySWmatrix(1, val2, keyID, keyID);
      /* MAUTO vals.push_back(val2); */
    }
    /* MAUTO vals.push_back(val); */
    val = MulModPrecon(val, gi, m, giminv); // val *= g mod m (= g^{j+1})
  }

  if (!native) {
    sKey.GenKeySWmatrix(1, context.zMStar.genToPow(i, -ord), keyID, keyID);
  }


/* MAUTO
  sKey.resetTree(i,keyID); // remove existing tree, if any
  sKey.add2tree(i, 1, vals, keyID);
*/
}
#else
// adds all matrices for dim i.
// i == -1 => Frobenius (NOTE: in matmul1D, i ==#gens means something else,
//   so it is best to avoid that).
static void add1Dmats4dim(FHESecKey& sKey, long i, long keyID)
{
  const PAlgebra& zMStar = sKey.getContext().zMStar;
  long ord;
  bool native;

  if (i != -1) {
    ord = zMStar.OrderOf(i);
    native = zMStar.SameOrd(i);\
  }
  else {
    // Frobenius
    ord = zMStar.getOrdP();
    native = true;
  }

  for (long j = 1; j < ord; j++) 
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, j), keyID, keyID);

  if (!native)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, -ord), keyID, keyID);

  sKey.setKSStrategy(i, FHE_KSS_FULL);
}


#endif

static std::pair<long,long> computeSteps(long ord, long bound, bool native)
{
  long baby,giant;
  if (native) { // using giant+baby matrices
    if (bound*bound >= 4*ord)
      giant = ceil((bound - sqrt((double)bound*bound -4*ord))/2.0);
    else
      giant = sqrt((double)ord);
  }
  else { // using giant+2*baby matrices
    if (bound*bound >= 8*ord)
      giant = ceil((bound - sqrt((double)bound*bound -8*ord))/2.0);
    else
      giant = sqrt((double)ord);
  }
  baby = ord/giant;
  if (baby*giant<ord) baby++;

  //cerr << "*** giant steps = " << giant << "\n";
  return std::pair<long,long>(baby,giant);
}

#if 0
static void addSome1Dmats4dim(FHESecKey& sKey, long i, long bound, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  long m = context.zMStar.getM();
  computeParams(context,m,i); // defines vars: native, ord, gi, g2md, giminv, g2mdminv

  long baby, giant;
  std::tie(baby,giant) = computeSteps(ord, bound, native);

  for (long j=1,val=gi; j<=baby; j++) { // Add matrices for baby steps
    sKey.GenKeySWmatrix(1, val, keyID, keyID);
    if (!native) {
      long val2 = MulModPrecon(val,g2md,m,g2mdminv);
      sKey.GenKeySWmatrix(1, val2, keyID, keyID);
    }
    val = MulModPrecon(val, gi, m, giminv); // val *= g mod m (= g^{j+1})
   }

  long gb = PowerMod(gi,baby,m); // g^baby
  NTL::mulmod_precon_t gbminv = PrepMulModPrecon(gb, m);  
  for (long j=2,val=gb; j < giant; j++) { // Add matrices for giant steps
    val = MulModPrecon(val, gb, m, gbminv); // val = g^{(j+1)*baby}
    sKey.GenKeySWmatrix(1, val, keyID, keyID);
  }

  if (!native) {
    sKey.GenKeySWmatrix(1, context.zMStar.genToPow(i, -ord), keyID, keyID);
  }

  // VJS: experimantal feature...because the replication code
  // uses rotations by -1, -2, -4, -8, we add a few
  // of these as well...only the small ones are important,
  // and we only need them if SameOrd(i)...
  // Note: we do indeed get a nontrivial speed-up

  if (native && i<context.zMStar.numOfGens()) {
    for (long k = 1; k < giant; k = 2*k) {
      long j = ord - k;
      long val = PowerMod(gi, j, m); // val = g^j
      sKey.GenKeySWmatrix(1, val, keyID, keyID);
    }
  }

#if 0 
MAUTO

  // build the tree for this dimension, the internal nodes are 1 and
  // (subset of) gi^{giant}, gi^{2*giant}, ..., gi^{baby*giant}. We

  MAUTO sKey.resetTree(i,keyID); // remove existing tree, if any 

  // keep a list of all the elements that are covered by the tree so far,
  // initialized to only the root (=1).
  std::unordered_set<long> covered({1});

  // Make a list of the automorphisms for this dimension
  std::vector<long> autos;
  for (long j=1,val=gi; j<ord; j++) {
    // Do we have matrices for val and/or val/gi^{di}?
    if (!native) {
      long val2 = MulModPrecon(val, g2md, m, g2mdminv);
      if (sKey.haveKeySWmatrix(1,val2,keyID,keyID)) {
        autos.push_back(val2);
      }
    }
    if (sKey.haveKeySWmatrix(1,val,keyID,keyID)) {
      autos.push_back(val);
    }
    val = MulModPrecon(val, gi, m, giminv); // g^{j+1}
  }

  // Insert internal nodes and their children to tree
  for (long j=0,fromVal=1; j<giant; j++) {
    NTL::mulmod_precon_t fromminv = PrepMulModPrecon(fromVal, m);      
    vector<long> children;
    for (long k: autos) {
      long toVal = MulModPrecon(k, fromVal, m, fromminv);
      if (covered.count(toVal)==0) { // toVal not covered yet
        covered.insert(toVal);
        children.push_back(toVal);
      }
    }
    if (!children.empty()) { // insert fromVal with its children
      sKey.add2tree(i, fromVal, children, keyID);
    }
    fromVal = MulModPrecon(fromVal, gb, m, gbminv); // g^{(j+1)*baby}
  }

  // Sanity-check, did we cover everything?
  long toCover = native? ord: (2*ord-1);
  if (covered.size()<toCover)
    cerr << "**Warning: order-"<<ord<<" dimension, covered "<<covered.size()
         << " of "<<toCover<<endl;
#endif
}


#else
// same as above, but uses BS/GS strategy
static void addSome1Dmats4dim(FHESecKey& sKey, long i, long bound, long keyID)
{
  const PAlgebra& zMStar = sKey.getContext().zMStar;
  long ord;
  bool native;

  if (i != -1) {
    ord = zMStar.OrderOf(i);
    native = zMStar.SameOrd(i);\
  }
  else {
    // Frobenius
    ord = zMStar.getOrdP();
    native = true;
  }

  long g = KSGiantStepSize(ord);

  // baby steps
  for (long j = 1; j < g; j++)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, j), keyID, keyID);

  // giant steps
  for (long j = g; j < ord; j += g)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, j), keyID, keyID);

  if (!native)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, -ord), keyID, keyID);

  sKey.setKSStrategy(i, FHE_KSS_BSGS);

  // NOTE: the old code also added matrices for ord-2^k for small k,
  // in the case of (native && i<context.zMStar.numOfGens()).
  // This supposedly speeds up the replication code, but for now
  // we are leaving this out, for simplicity.   Also, it is a waste
  // of space for applications that don't use replication.

}


#endif

// generate only matrices of the form s(X^{g^i})->s(X), but not all of them.
// For a generator g whose order is larger than bound, generate only enough
// matrices for the giant-step/baby-step procedures (2*sqrt(ord(g))of them).
void addSome1DMatrices(FHESecKey& sKey, long bound, long keyID)
{
  const FHEcontext &context = sKey.getContext();

  // key-switching matrices for the automorphisms
  for (long i: range(context.zMStar.numOfGens())) {
          // For generators of small order, add all the powers
    if (bound >= context.zMStar.OrderOf(i))
      add1Dmats4dim(sKey, i, keyID);
    else  // For generators of large order, add only some of the powers
      addSome1Dmats4dim(sKey, i, bound, keyID);
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

// Generate all Frobenius matrices of the form s(X^{p^i})->s(X)
void addSomeFrbMatrices(FHESecKey& sKey, long bound, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  if (bound >= LONG(context.zMStar.getOrdP()))
    add1Dmats4dim(sKey, -1, keyID);
  else  // For generators of large order, add only some of the powers
    addSome1Dmats4dim(sKey, -1, bound, keyID);

  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

static void addMinimal1Dmats4dim(FHESecKey& sKey, long i, long keyID)
{
  const PAlgebra& zMStar = sKey.getContext().zMStar;
  long ord;
  bool native;

  if (i != -1) {
    ord = zMStar.OrderOf(i);
    native = zMStar.SameOrd(i);\
  }
  else {
    // Frobenius
    ord = zMStar.getOrdP();
    native = true;
  }

  sKey.GenKeySWmatrix(1, zMStar.genToPow(i, 1), keyID, keyID);

  if (!native)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, -ord), keyID, keyID);

  if (ord > FHE_KEYSWITCH_MIN_THRESH) {
    long g = KSGiantStepSize(ord);
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, g), keyID, keyID);
  }

  sKey.setKSStrategy(i, FHE_KSS_MIN);

}

void addMinimal1DMatrices(FHESecKey& sKey, long keyID)
{
  const FHEcontext &context = sKey.getContext();

  // key-switching matrices for the automorphisms
  for (long i: range(context.zMStar.numOfGens())) {
    addMinimal1Dmats4dim(sKey, i, keyID);
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

// Generate all Frobenius matrices of the form s(X^{p^i})->s(X)
void addMinimalFrbMatrices(FHESecKey& sKey, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  addMinimal1Dmats4dim(sKey, -1, keyID);
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

// Generate all key-switching matrices for a given permutation network
void addMatrices4Network(FHESecKey& sKey, const PermNetwork& net, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  long m = context.zMStar.getM();

  for (long i=0; i<net.depth(); i++) {
    long e = net.getLayer(i).getE();
    long gIdx = net.getLayer(i).getGenIdx();
    long g = context.zMStar.ZmStarGen(gIdx);
    long g2e = PowerMod(g, e, m); // g^e mod m
    const Vec<long>&shamts = net.getLayer(i).getShifts();
    for (long j=0; j<shamts.length(); j++) {
      if (shamts[j]==0) continue;
      long val = PowerMod(g2e, shamts[j], m);
      sKey.GenKeySWmatrix(1, val, keyID, keyID);
    }
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

void addTheseMatrices(FHESecKey& sKey,
		      const std::set<long>& automVals, long keyID)
{
  std::set<long>::iterator it;
  for (it=automVals.begin(); it!=automVals.end(); ++it) {
    long k = *it;
    sKey.GenKeySWmatrix(1, k, keyID, keyID);
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

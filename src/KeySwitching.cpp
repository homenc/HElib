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
#include "NTL/ZZ.h"
NTL_CLIENT
#include "FHE.h"
#include "permutations.h"

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
  long ord, gi, ginv;\
  NTL::mulmod_precon_t gminv;\
  if (i==context.zMStar.numOfGens()) { /* Frobenius matrices */\
    ord = context.zMStar.getOrdP();\
    gi = context.zMStar.getP();\
    native = true;\
  }\
  else { /* one of the "regular" dimensions */\
    ord = context.zMStar.OrderOf(i);\
    gi = context.zMStar.ZmStarGen(i), ginv;\
    native = context.zMStar.SameOrd(i);\
    if (!native) {\
      ginv = NTL::InvMod(PowerMod(gi,ord,m), m); /* g^{-ord} mod m */\
      gminv = PrepMulModPrecon(ginv, m);\
    }\
  }

static void add1Dmats4dim(FHESecKey& sKey, long i, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  long m = context.zMStar.getM();
  computeParams(context,m,i); // defines vars: native, ord, gi, ginv, gminv

  for (long j = 1; j < ord; j++) {
    long val = PowerMod(gi, j, m); // val = g^j
    // From s(X^val) to s(X)
    sKey.GenKeySWmatrix(1, val, keyID, keyID);
    if (!native) { // also from s(X^{g^{i-ord}}) to s(X)
      long val2 = MulModPrecon(val,ginv,m,gminv);
      sKey.GenKeySWmatrix(1, val2, keyID, keyID);
    }
  }
}

static std::pair<long,long>
computeSteps(long m, long ord, long bound, bool native)
{
  long baby,giant;
  if (native) { // using giant+baby matrices
    if (bound*bound >= 4*m)
      giant = ceil(bound - sqrt((double)bound*bound -4*m)/2.0);
    else
      giant = sqrt(double(m));
  }
  else { // using giant+2*baby matrices
    if (bound*bound >= 8*m)
      giant = ceil(bound - sqrt((double)bound*bound -8*m)/2.0);
    else
      giant = sqrt(double(m));
  }
  baby = m/giant;
  if (baby*giant<m) baby++;
  return std::pair<long,long>(baby,giant);
}

static void addSome1Dmats4dim(FHESecKey& sKey, long i, long bound, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  long m = context.zMStar.getM();
  computeParams(context,m,i); // defines vars: native, ord, gi, ginv, gminv

  long baby, giant;
  std::tie(baby,giant) = computeSteps(m, ord, bound, native);

  for (long j=1; j <= baby; j++) { // Add matrices for baby steps
    long val = PowerMod(gi, j, m);  // g^j
    sKey.GenKeySWmatrix(1, val, keyID, keyID);
    if (!native) {
      long val2 = MulModPrecon(val,ginv,m,gminv);
      sKey.GenKeySWmatrix(1, val2, keyID, keyID);
    }
  }
  for (long j=1; j <= giant; j++) { // Add matrices for giant steps
    long val = PowerMod(gi, j*baby, m);  // g^{j*baby}
    sKey.GenKeySWmatrix(1, val, keyID, keyID);
  }

  // VJS: experimantal feature...because the replication code
  // uses rotations by -1, -2, -4, -8, we add a few
  // of these as well...only the small ones are important,
  // and we only need them if SameOrd(i)...
  // Note: we do indeed get a nontrivial speed-up

  if (native && i<context.zMStar.numOfGens()) {
    for (long k = 1; k <= giant; k = 2*k) {
      long j = ord - k;
      long val = PowerMod(gi, j, m); // val = g^j
      sKey.GenKeySWmatrix(1, val, keyID, keyID);
    }
  }
}

// generate only matrices of the form s(X^{g^i})->s(X), but not all of them.
// For a generator g whose order is larger than bound, generate only enough
// matrices for the giant-step/baby-step procedures (2*sqrt(ord(g))of them).
void addSome1DMatrices(FHESecKey& sKey, long bound, long keyID)
{
  const FHEcontext &context = sKey.getContext();

  // key-switching matrices for the automorphisms
  for (long i = 0; i < (long)context.zMStar.numOfGens(); i++) {
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
  if (bound >= (long)context.zMStar.getOrdP())
    add1Dmats4dim(sKey, context.zMStar.numOfGens(), keyID);
  else  // For generators of large order, add only some of the powers
    addSome1Dmats4dim(sKey, context.zMStar.numOfGens(), bound, keyID);

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

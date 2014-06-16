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

// generate all matrices of the form s(X^{g^i})->s(X) for generators g of
// Zm* /<2> and i<ord(g). If g has different orders in Zm* and Zm* /<2>
// then generate also matrices of the form s(X^{g^{-i}})->s(X)
void add1DMatrices(FHESecKey& sKey, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  long m = context.zMStar.getM();

  // key-switching matrices for the automorphisms
  for (long i = 0; i < (long)context.zMStar.numOfGens(); i++) {
    for (long j = 1; j < (long)context.zMStar.OrderOf(i); j++) {
      long val = PowerMod(context.zMStar.ZmStarGen(i), j, m); // val = g^j
      // From s(X^val) to s(X)
      sKey.GenKeySWmatrix(1, val, keyID, keyID);
      if (!context.zMStar.SameOrd(i))
	// also from s(X^{1/val}) to s(X)
	sKey.GenKeySWmatrix(1, InvMod(val,m), keyID, keyID);
    }
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

// generate only matrices of the form s(X^{g^i})->s(X), but not all of them.
// For a generator g whose order is larger than bound, generate only enough
// matrices for the giant-step/baby-step procedures (2*sqrt(ord(g))of them).
void addSome1DMatrices(FHESecKey& sKey, long bound, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  long m = context.zMStar.getM();

  // key-switching matrices for the automorphisms
  for (long i = 0; i < (long)context.zMStar.numOfGens(); i++) {
    // For generators of small order, add all the powers
    if (bound >= (long)context.zMStar.OrderOf(i))
      for (long j = 1; j < (long)context.zMStar.OrderOf(i); j++) {
	long val = PowerMod(context.zMStar.ZmStarGen(i), j, m); // val = g^j
	// From s(X^val) to s(X)
	sKey.GenKeySWmatrix(1, val, keyID, keyID);
	if (!context.zMStar.SameOrd(i))
	  // also from s(X^{1/val}) to s(X)
	  sKey.GenKeySWmatrix(1, InvMod(val,m), keyID, keyID);
      }
    else { // For generators of large order, add only some of the powers
      long num = SqrRoot(context.zMStar.OrderOf(i)); // floor(ord^{1/2})
      if (num*num < (long) context.zMStar.OrderOf(i)) num++; // ceil(ord^{1/2})

      // VJS: the above two lines replaces the following inexact calculation
      // with an exact calculation
      // long num = ceil(sqrt((double)context.zMStar.OrderOf(i)));

      for (long j=1; j <= num; j++) { // Add matrices for g^j and g^{j*num}
	long val1 = PowerMod(context.zMStar.ZmStarGen(i), j, m);  // g^j
	long val2 = PowerMod(context.zMStar.ZmStarGen(i),num*j,m);// g^{j*num}
	if (j < num) {
	  sKey.GenKeySWmatrix(1, val1, keyID, keyID);
	  sKey.GenKeySWmatrix(1, val2, keyID, keyID);
	}
	if (!context.zMStar.SameOrd(i)) {
	  //	  sKey.GenKeySWmatrix(1, InvMod(val1,m), keyID, keyID);
	  sKey.GenKeySWmatrix(1, InvMod(val2,m), keyID, keyID);
	}
      }

      // VJS: experimantal feature...because the replication code
      // uses rotations by -1, -2, -4, -8, we add a few
      // of these as well...only the small ones are important,
      // and we only need them if SameOrd(i)...
      // Note: we do indeed get a nontrivial speed-up

      if (context.zMStar.SameOrd(i)) {
        for (long k = 1; k <= num; k = 2*k) {
          long j = context.zMStar.OrderOf(i) - k;
          long val = PowerMod(context.zMStar.ZmStarGen(i), j, m); // val = g^j
          sKey.GenKeySWmatrix(1, val, keyID, keyID);
        }
      }
    }
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

// Generate all Frobenius matrices of the form s(X^{2^i})->s(X)
void addFrbMatrices(FHESecKey& sKey, long keyID)
{
  const FHEcontext &context = sKey.getContext();
  long m = context.zMStar.getM();

  for (long j = 1; j < (long)context.zMStar.getOrdP(); j++) {
    long val = PowerMod(context.zMStar.getP(), j, m); // val = p^j mod m
    sKey.GenKeySWmatrix(1, val, keyID, keyID);
  }
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

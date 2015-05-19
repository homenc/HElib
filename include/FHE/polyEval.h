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
/**
 * @file polyEval.h
 * @brief Homomorphic Polynomial Evaluation
 */


#include "Ctxt.h"

//! @brief Evaluate a cleartext polynomial on an encrypted input
//! @param[out] res  to hold the return value
//! @param[in]  poly the degree-d polynomial to evaluate
//! @param[in]  x    the point on which to evaluate
//! @param[in]  k    optional optimization parameter, defaults to sqrt(d/2) rounded up or down to a power of two
void polyEval(Ctxt& ret, ZZX poly, const Ctxt& x, long k=0);
     // Note: poly is passed by value, so caller keeps the original

//! @brief Evaluate an encrypted polynomial on an encrypted input
//! @param[out] res  to hold the return value
//! @param[in]  poly the degree-d polynomial to evaluate
//! @param[in]  x    the point on which to evaluate
void polyEval(Ctxt& ret, const Vec<Ctxt>& poly, const Ctxt& x);


// A useful helper class

//! @brief Store powers of X, compute them synamically as needed.
// This implementation assumes that the size (# of powers) is determine
// at initialization time, it is not hard to grow the vector as needed,
// but not clear if there is any application that needs it.
class DynamicCtxtPowers {
private:
  vector<Ctxt> v;   // A vector storing the powers themselves

public:
  DynamicCtxtPowers(const Ctxt& c, long nPowers)
  {
    assert (!c.isEmpty() && nPowers>0); // Sanity-check

    Ctxt tmp(c.getPubKey(), c.getPtxtSpace());
    v.resize(nPowers, tmp); // Initializes nPowers empty cipehrtexts
    v[0] = c;               // store X itself in v[0]
  }

  //! @brief Returns the e'th power, computing it as needed
  Ctxt& getPower(long e); // must use e >= 1, else throws an exception

  //! dp.at(i) and dp[i] both return the i+1st power
  Ctxt& at(long i) { return getPower(i+1); }
  Ctxt& operator[](long i) { return getPower(i+1); }

  const vector<Ctxt>& getVector() const { return v; }
  long size() const { return v.size(); }
  bool isPowerComputed(long i)
  { return (i>0 && i<=(long)v.size() && !v[i-1].isEmpty()); }
};

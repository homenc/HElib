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
#ifndef _RECRYPTION_H_
#define _RECRYPTION_H_
/** @file recryption.h
 *  @brief Define some data structures to hold recryption data
 */

#include "NumbTh.h"

class  PAlgebraMod;
class  EncryptedArray;
class  EvalMap;
class  PowerfulDCRT;
class  FHEcontext;
class  FHEPubKey;


//! @class RecryptData
//! @brief A structure to hold recryption-related data inside the FHEcontext
class RecryptData {
public:
  //! default Hamming weight of recryption key
  static const long defSkHwt=100;

  //! Some data members that are only used for I/O
  Vec<long> mvec;     //! partition of m into co-prime factors
  long hwt;           //! Hamming weight of recryption secret-key
  bool conservative;  //! flag for choosing more conservatice parameters

  //! skey encrypted wrt space p^{e-e'+r}
  long e, ePrime;

  //! Hamming weight of recryption secret key
  long skHwt;

  //! an optimization parameter
  double alpha;

  //! for plaintext space p^{e-e'+r}
  PAlgebraMod *alMod;

  //! for plaintext space p^{e-e'+r}
  EncryptedArray *ea;

  bool build_cache;


  //! linear maps
  EvalMap *firstMap, *secondMap;

  //! conversion between ZZX and Powerful
  PowerfulDCRT *p2dConv;

  //! linPolys for uppacking the slots
  vector<ZZX> unpackSlotEncoding;

  RecryptData() {
    hwt=0; conservative=false; e=ePrime=0; alpha=0.0;
    alMod=NULL; ea=NULL; firstMap=NULL; secondMap=NULL; p2dConv=NULL;
    build_cache = false;
  }
  ~RecryptData();

  //! Initialize the recryption data in the context
  void init(const FHEcontext& context, const Vec<long>& mvec_,
            long t=0/*min Hwt for sk*/, 
            bool consFlag=false,
            bool build_cache=false,
            bool minimal=false);

  bool operator==(const RecryptData& other) const;
  bool operator!=(const RecryptData& other) const {
    return !(operator==(other));
  }

  //! Helper function for computing the recryption parameters
  static double
    setAlphaE(double& alpha, long& e, long& ePrime,
              const FHEcontext& context, bool conservative=false, long=0);
  /** To get the smallest value of e-e', the params need to satisfy:
   *  (p^e +1)/4 =>
   *   max { (t+1)( 1+ (alpha/2)*(p^e/p^{ceil(log_p(t+2))}) ) + noise      }
   *       { (t+1)( 1+ ((1-alpha)/2)*(p^e/p^{ceil(log_p(t+2))}) +p^r/2) +1 },
   *
   * where noise is taken to be twice the mod-switching additive term,
   * namely noise = p^r *sqrt((t+1)*phi(m)/3).
   * Denoting rho=(t+1)/p^{ceil(log_p(t+2))} (and ignoring fome +1 terms),
   * this is equivalent to:
   *
   *   p^e > max{4(t+noise)/(1-2*alpha*rho), 2(t+1)p^r/(1-2(1-alpha)rho)}.
   *
   * We first compute the optimal value for alpha (which must be in [0,1]),
   * that makes the two terms in the max{...} as close as possible, and
   * then compute the smallest value of e satisfying this constraint.
   *
   * If this value is too big then we try again with e-e' one larger,
   * which means that rho is a factor of p smaller.
   */
};

#endif /* _RECRYPTION_H_ */

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

extern long thinRecrypt_initial_level;

class  PAlgebraMod;
class  EncryptedArray;
class  EvalMap;
class  ThinEvalMap;
class  PowerfulDCRT;
class  FHEcontext;
class  FHEPubKey;


//! @class RecryptData
//! @brief A structure to hold recryption-related data inside the FHEcontext
class RecryptData {
public:
  //! default Hamming weight of recryption key
  static constexpr long defSkHwt=100;

  // NOTE: here's to hoping: for "random enough" x we expect
  //   |x|_powerful < |x|_canonical * sqrt(someConstant/phi(m))
  //   we set magicConst = sqrt(phi(m)/3) * sqrt(someConstant/phi(m))
  //   (sqrt(phi(m)/3) comes from context.noiseBoundForUniform())
  static /*constexpr*/ double magicConst; // defaults to 2

  //! Some data members that are only used for I/O
  NTL::Vec<long> mvec;     //! partition of m into co-prime factors

  //! skey encrypted wrt space p^{e-e'+r}
  long e, ePrime;

  //! Hamming weight of recryption secret key
  long skHwt;

  //! an optimization parameter
  long a;

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
  std::vector<NTL::ZZX> unpackSlotEncoding;

  RecryptData() {
    skHwt=0; e=ePrime=0; a=0;
    alMod=NULL; ea=NULL; firstMap=NULL; secondMap=NULL; p2dConv=NULL;
    build_cache = false;
  }
  ~RecryptData();

  //! Initialize the recryption data in the context
  void init(const FHEcontext& context, const NTL::Vec<long>& mvec_,
            long t=0/*min Hwt for sk*/, 
            bool build_cache=false,
            bool minimal=false);

  bool operator==(const RecryptData& other) const;
  bool operator!=(const RecryptData& other) const {
    return !(operator==(other));
  }

  //! Helper function for computing the recryption parameters
  static long setAE(long& a, long& e, long& ePrime,
                    const FHEcontext& context, long t=0);
  /** We want to get the smallest value of e-e', subject to a few
   * constraints: For the magicConst from above, an exponent e,
   * a norm-t secret key, and plaintext space mod p^r, we need to
   * find integers a and b so that:
   *
   *    (1) (4a+ 8p^r)(t+1)*magicConst <= q  = p^e +1
   *    (2) (4b + 5)(t+1) * magicConst <= q-4= p^e -3
   *
   * Then e' is the largest exponent such that p^{e'} <= 2(a+b).
   *
   * Note that if we let e,e' tend to infinity and set a=b=p^{e'}/4,
   * then the two constraints above degenerate to
   *
   *    2(t+1)*magicConst < p^{e-e'}
   *
   * so the smallest value of e-e' that we can hope for is
   *
   *    e-e' = ceiling( log_p( 2(t+1)*magicConst ) )
   *
   * The setAE procedure tries to find a setting of e,e',a (and
   * b = p^{e'}/2 -a) that satisfies the constraints (1,2) and
   * yeilds the smallest e-e', so long as e is "not too big".
   * Specifically, it minimizes e-e' while ensuring p^e < 2^{30},
   * if possible (and failing that it will just satisfy the
   * constraints (1,2)).
   *
   * Once e, e', a are set, it copmutes and returns the largest 
   * Hamming-weight for the key for which constraints (1,2) still hold.
   * NOTE: it returns the Hamming weight, *not* the size of the key.
   *       The size can be computed by calling the function
   *       sampleHWtBoundedEffectiveBound(context, weight)
   */
};


//! @class ThinRecryptData
//! @brief Same as above, but for "thin" bootstrapping, where the slots 
//! are assumed to contain constants
class ThinRecryptData {
public:
  //! default Hamming weight of recryption key
  static const long defSkHwt=100;

  //! Some data members that are only used for I/O
  NTL::Vec<long> mvec;     //! partition of m into co-prime factors

  //! skey encrypted wrt space p^{e-e'+r}
  long e, ePrime;

  //! Hamming weight of recryption secret key
  long skHwt;

  //! an optimization parameter
  long a;

  //! for plaintext space p^{e-e'+r}
  PAlgebraMod *alMod;

  //! for plaintext space p^{e-e'+r}
  EncryptedArray *ea;

  bool build_cache;


  //! linear maps
  ThinEvalMap *coeffToSlot, *slotToCoeff;

  ThinRecryptData() {
    skHwt=0; e=ePrime=0; a=0;
    alMod=NULL; ea=NULL; coeffToSlot=NULL; slotToCoeff=NULL; 
    build_cache = false;
  }
  ~ThinRecryptData();

  //! Initialize the recryption data in the context
  void init(const FHEcontext& context, const NTL::Vec<long>& mvec_,
            long t=0/*min Hwt for sk*/, 
            bool build_cache=false,
            bool minimal=false);

  bool operator==(const ThinRecryptData& other) const;
  bool operator!=(const ThinRecryptData& other) const {
    return !(operator==(other));
  }
};


#endif /* _RECRYPTION_H_ */

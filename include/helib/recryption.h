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
#ifndef HELIB_RECRYPTION_H
#define HELIB_RECRYPTION_H
/** @file recryption.h
 *  @brief Define some data structures to hold recryption data
 */

#include <helib/NumbTh.h>

namespace helib {

extern long thinRecrypt_initial_level;
extern long fhe_force_chen_han;
extern long printFlag;

class PAlgebraMod;
class EncryptedArray;
class EvalMap;
class ThinEvalMap;
class PowerfulDCRT;
class Context;
class PubKey;

//! @class RecryptData
//! @brief A structure to hold recryption-related data inside the Context
class RecryptData
{
public:
  //! Some data members that are only used for I/O
  NTL::Vec<long> mvec; //! partition of m into co-prime factors

  //! skey encrypted wrt space p^{e-e'+r}
  long e, ePrime;

  //! Hamming weight of recryption secret key
  long skHwt;

  //! for plaintext space p^{e-e'+r}
  std::shared_ptr<const PAlgebraMod> alMod = nullptr;

  //! for plaintext space p^{e-e'+r}
  std::shared_ptr<const EncryptedArray> ea = nullptr;

  bool build_cache;

  bool alsoThick;

  //! linear maps
  std::shared_ptr<const EvalMap> firstMap = nullptr, secondMap = nullptr;

  //! conversion between ZZX and Powerful
  std::shared_ptr<const PowerfulDCRT> p2dConv = nullptr;

  //! linPolys for unpacking the slots
  std::vector<NTL::ZZX> unpackSlotEncoding;

  RecryptData()
  {
    skHwt = 0;
    e = ePrime = 0;
    build_cache = false;
    alsoThick = false;
  }

  //! Initialize the recryption data in the context
  void init(const Context& context,
            const NTL::Vec<long>& mvec_,
            bool enableThick, /*init linear transforms for non-thin*/
            bool build_cache = false,
            bool minimal = false);

  bool operator==(const RecryptData& other) const;
  bool operator!=(const RecryptData& other) const
  {
    return !(operator==(other));
  }

  //! Helper function for computing the recryption parameters
  static void setAE(long& e, long& ePrime, const Context& context);
  // VJS-FIXME: this needs to be documented.
  // It is based on the most recent version of our bootstrapping
  // paper (see Section 6.2)
};

//! @class ThinRecryptData
//! @brief Same as above, but for "thin" bootstrapping, where the slots
//! are assumed to contain constants
class ThinRecryptData : public RecryptData
{
public:
  //! linear maps
  std::shared_ptr<const ThinEvalMap> coeffToSlot, slotToCoeff;

  //! Initialize the recryption data in the context
  void init(const Context& context,
            const NTL::Vec<long>& mvec_,
            bool alsoThick, /*init linear transforms also for non-thin*/
            bool build_cache = false,
            bool minimal = false);
};

#define HELIB_MIN_CAP_FRAC (2.0 / 3.0)
// Used in calculation of "min capacity".
// This could be set to 1.0, but just to be on the safe side,
// it is set to 2/3.  If we did set it to 1, the min capacity
// would increase by less than 6/10 of a bit.

} // namespace helib

#endif // HELIB_RECRYPTION_H
